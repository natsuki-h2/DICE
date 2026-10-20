"""Re-image the downloaded GHOST F1 MS and compare with its dirty FITS.

Run from the repository root with CASA 6. Products and the comparison summary
are written below ``data/ghost_f1_5min/`` and ignored by Git.
"""

from pathlib import Path
import json
import shutil

import numpy as np
from casatasks import exportfits, imhead, imstat, tclean
from casatools import image


REPO_ROOT = Path(__file__).resolve().parents[1]
WORKDIR = REPO_ROOT / "data" / "ghost_f1_5min"
MS = WORKDIR / "output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms"
GHOST_DIRTY = WORKDIR / "ALMA-GHOST-CCK6FJ4E.dirty.image.fits"
BLISS_DIRTY = WORKDIR / "bliss_f1_dirty"
BLISS_DIRTY_FITS = WORKDIR / "bliss_f1_dirty.image.fits"
SUMMARY = WORKDIR / "f1_reimage_comparison.json"

IMAGE_SIZE = 360
CELL = "0.12454065519280681arcsec"


def read_pixels(path):
    ia = image()
    try:
        ia.open(str(path))
        values = np.asarray(ia.getchunk(), dtype=np.float64)
    finally:
        ia.done()
    values = np.squeeze(values)
    if values.shape != (IMAGE_SIZE, IMAGE_SIZE):
        raise RuntimeError(f"Unexpected image shape: {values.shape}")
    return values


def beam(path):
    return {
        "major_arcsec": float(
            imhead(str(path), mode="get", hdkey="beammajor")["value"]
        ),
        "minor_arcsec": float(
            imhead(str(path), mode="get", hdkey="beamminor")["value"]
        ),
        "pa_deg": float(imhead(str(path), mode="get", hdkey="beampa")["value"]),
    }


def main():
    for path in WORKDIR.glob(BLISS_DIRTY.name + ".*"):
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()
    tclean(
        vis=str(MS),
        imagename=str(BLISS_DIRTY),
        specmode="mfs",
        imsize=[IMAGE_SIZE, IMAGE_SIZE],
        cell=CELL,
        weighting="natural",
        robust=0.5,
        niter=0,
        pblimit=0.2,
        interactive=False,
    )
    exportfits(
        imagename=str(BLISS_DIRTY) + ".image",
        fitsimage=str(BLISS_DIRTY_FITS),
        overwrite=True,
        dropdeg=True,
    )

    ghost = read_pixels(GHOST_DIRTY)
    bliss = read_pixels(BLISS_DIRTY_FITS)
    common_finite = np.isfinite(ghost) & np.isfinite(bliss)
    if not np.any(common_finite):
        raise RuntimeError("The dirty images have no common finite pixels")
    difference = bliss[common_finite] - ghost[common_finite]
    ghost_stats = imstat(str(GHOST_DIRTY), listit=False)
    bliss_stats = imstat(str(BLISS_DIRTY_FITS), listit=False)

    summary = {
        "settings": {
            "specmode": "mfs",
            "imsize": [IMAGE_SIZE, IMAGE_SIZE],
            "cell": CELL,
            "weighting": "natural",
            "robust": 0.5,
            "niter": 0,
            "pblimit": 0.2,
        },
        "ghost_rms_jy_per_beam": float(ghost_stats["rms"][0]),
        "bliss_rms_jy_per_beam": float(bliss_stats["rms"][0]),
        "rms_ratio_bliss_over_ghost": float(
            bliss_stats["rms"][0] / ghost_stats["rms"][0]
        ),
        "maximum_absolute_pixel_difference": float(np.max(np.abs(difference))),
        "difference_rms": float(np.sqrt(np.mean(difference**2))),
        "n_common_finite_pixels": int(np.count_nonzero(common_finite)),
        "ghost_beam": beam(GHOST_DIRTY),
        "bliss_beam": beam(BLISS_DIRTY_FITS),
    }
    with SUMMARY.open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Wrote {SUMMARY}")


if __name__ == "__main__":
    main()
