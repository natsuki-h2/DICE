"""Create full-primary-beam native dirty images for GHOST F1--F6.

Run with CASA 6. CASA image products are removed after compact summary metrics
and display arrays have been saved.
"""

from pathlib import Path
import csv
import json
import math
import shutil

import numpy as np
from casatasks import imhead, imstat, tclean
from casatools import image


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/ghost_matrix/native_full_pb_blank"
CELLS = {
    "F1": 0.12454065519280681,
    "F2": 0.093125,
    "F3": 0.034151,
    "F4": 0.020579,
    "F5": 0.007234,
    "F6": 0.003668,
}
# FFT-friendly sizes that cover the PB >= 0.2 diameter with margin.
IMSIZES = {"F1": 360, "F2": 480, "F3": 1280,
           "F4": 2112, "F5": 6000, "F6": 12000}
MS_PATHS = {
    "F1": ROOT / "data/ghost_f1_5min" /
          "output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms",
    **{
        f"F{i}": ROOT / "data/ghost_matrix" / f"f{i}_5min" / "extracted" /
        "output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms"
        for i in range(2, 7)
    },
}


def remove_products(prefix):
    for path in prefix.parent.glob(prefix.name + ".*"):
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def read_plane(path):
    ia = image()
    try:
        ia.open(str(path))
        values = np.squeeze(np.asarray(ia.getchunk(), dtype=np.float32))
        mask = np.squeeze(np.asarray(ia.getchunk(getmask=True), dtype=bool))
    finally:
        ia.done()
    return values, mask


def thumbnail(values, mask, maximum=720):
    """Block-average a square plane while respecting its CASA mask."""
    n = values.shape[0]
    factor = max(1, int(math.ceil(n / maximum)))
    usable = n - n % factor
    data = values[:usable, :usable]
    valid = mask[:usable, :usable] & np.isfinite(data)
    shape = (usable // factor, factor, usable // factor, factor)
    numerator = np.where(valid, data, 0).reshape(shape).sum(axis=(1, 3))
    denominator = valid.reshape(shape).sum(axis=(1, 3))
    result = np.full(numerator.shape, np.nan, dtype=np.float32)
    np.divide(numerator, denominator, out=result, where=denominator > 0)
    return result, factor


def header_value(path, key):
    return float(imhead(str(path), mode="get", hdkey=key)["value"])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for configuration in CELLS:
        cell = CELLS[configuration]
        imsize = IMSIZES[configuration]
        prefix = OUT / configuration.lower()
        remove_products(prefix)
        print("Imaging", configuration, imsize, "pixels at", cell, "arcsec", flush=True)
        tclean(vis=str(MS_PATHS[configuration]), imagename=str(prefix),
               specmode="mfs", imsize=[imsize, imsize],
               cell=f"{cell}arcsec", weighting="natural", robust=0.5,
               niter=0, pblimit=0.2, gridder="standard", interactive=False)
        image_path = Path(str(prefix) + ".image")
        psf_path = Path(str(prefix) + ".psf")
        pb_path = Path(str(prefix) + ".pb")
        dirty, dirty_mask = read_plane(image_path)
        psf, psf_mask = read_plane(psf_path)
        pb, pb_mask = read_plane(pb_path)
        major = header_value(image_path, "beammajor")
        minor = header_value(image_path, "beamminor")
        pa = header_value(image_path, "beampa")
        stats = imstat(str(image_path), listit=False)
        center = np.asarray(psf.shape) // 2
        yy, xx = np.indices(psf.shape)
        exclusion = np.hypot(xx - center[1], yy - center[0]) <= 2 * major / cell
        sidelobe_values = np.abs(psf[psf_mask & np.isfinite(psf) & ~exclusion])
        dirty_thumb, factor = thumbnail(dirty, dirty_mask)
        psf_thumb, psf_factor = thumbnail(psf, psf_mask)
        np.savez_compressed(OUT / f"{configuration.lower()}_display.npz",
                            dirty=dirty_thumb, psf=psf_thumb,
                            dirty_factor=factor, psf_factor=psf_factor)
        row = {
            "configuration": configuration,
            "cell_arcsec": cell,
            "imsize": imsize,
            "fov_arcsec": cell * imsize,
            "beam_major_arcsec": major,
            "beam_minor_arcsec": minor,
            "beam_pa_deg": pa,
            "minor_beam_pixels": minor / cell,
            "rms_jy_per_beam": float(stats["rms"][0]),
            "peak_abs_psf_outside_2beam": float(np.max(sidelobe_values)),
            "pb_min_valid": float(np.nanmin(pb[pb_mask & np.isfinite(pb)])),
            "pb_max": float(np.nanmax(pb[pb_mask & np.isfinite(pb)])),
            "gridder": "standard",
            "weighting": "natural",
            "uvtaper": "none",
            "niter": 0,
        }
        rows.append(row)
        (OUT / f"{configuration.lower()}_summary.json").write_text(
            json.dumps(row, indent=2) + "\n")
        remove_products(prefix)
    with (OUT / "native_full_pb_blank_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    print(OUT / "native_full_pb_blank_summary.csv")


if __name__ == "__main__":
    main()
