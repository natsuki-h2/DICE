"""Image one 5-minute GHOST MS per configuration on a common tapered grid.

This pilot verifies that a 0.7-arcsec uv taper makes the long-baseline
configurations sufficiently sampled on the F1 full-primary-beam grid.  It
does not inject sources or regenerate noise.
"""

from pathlib import Path
import csv
import shutil

import numpy as np
from casatasks import imstat, tclean
from casatools import image


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "ghost_matrix" / "common_resolution_pilot"
CELL = "0.12454065519280681arcsec"
IMSIZE = 360
UVTAPER = ["0.7arcsec"]
PB_LIMIT = 0.2
MS_PATHS = {
    "F1": ROOT / "data" / "ghost_f1_5min" / "output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms",
    **{
        f"F{index}": ROOT / "data" / "ghost_matrix" / f"f{index}_5min" / "extracted" /
        "output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms"
        for index in range(2, 7)
    },
}


def remove_products(prefix):
    for path in OUT.glob(prefix.name + ".*"):
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def plane_and_beam(path):
    ia = image()
    try:
        ia.open(str(path))
        data = np.squeeze(np.asarray(ia.getchunk(), dtype=float))
        mask = np.squeeze(np.asarray(ia.getchunk(getmask=True), dtype=bool))
        beam = ia.restoringbeam()
    finally:
        ia.done()
    return data, mask, beam


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for configuration, ms in MS_PATHS.items():
        prefix = OUT / configuration.lower()
        remove_products(prefix)
        print("Imaging", configuration, flush=True)
        tclean(vis=str(ms), imagename=str(prefix), specmode="mfs",
               imsize=[IMSIZE, IMSIZE], cell=CELL, weighting="natural",
               uvtaper=UVTAPER, robust=0.5, niter=0, pblimit=PB_LIMIT,
               interactive=False)
        stats = imstat(str(prefix) + ".image", listit=False)
        _, _, beam = plane_and_beam(str(prefix) + ".image")
        pb, valid, _ = plane_and_beam(str(prefix) + ".pb")
        finite_pb = pb[valid & np.isfinite(pb)]
        rows.append({
            "configuration": configuration,
            "uvtaper": UVTAPER[0],
            "cell_arcsec": float(CELL.removesuffix("arcsec")),
            "imsize": IMSIZE,
            "beam_major_arcsec": beam["major"]["value"],
            "beam_minor_arcsec": beam["minor"]["value"],
            "beam_pa_deg": beam["positionangle"]["value"],
            "rms_jy_per_beam": float(stats["rms"][0]),
            "pb_min_valid": float(np.min(finite_pb)),
            "pb_max": float(np.max(finite_pb)),
            "fraction_pixels_pb_ge_0p2": float(np.mean(valid & np.isfinite(pb) & (pb >= PB_LIMIT))),
        })
    with (OUT / "common_resolution_pilot.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader(); writer.writerows(rows)
    print(OUT / "common_resolution_pilot.csv")


if __name__ == "__main__":
    main()
