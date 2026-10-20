"""Create the small, deterministic FITS inputs for the public line-search demo."""

import csv
from pathlib import Path

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "tutorials/line_search_demo"
SEED = 2040
RMS_JY_BEAM = 1.0e-3


def add_line(cube, *, x, y, channel, peak_snr, spectral_sigma, beam_sigma=1.2):
    zz, yy, xx = np.indices(cube.shape)
    spatial = np.exp(-0.5 * (((xx - x) / beam_sigma) ** 2 + ((yy - y) / beam_sigma) ** 2))
    spectral = np.exp(-0.5 * ((zz - channel) / spectral_sigma) ** 2)
    cube += peak_snr * RMS_JY_BEAM * spatial * spectral


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    shape = (72, 64, 64)
    rng = np.random.default_rng(SEED)
    cube = rng.normal(0.0, RMS_JY_BEAM, shape)

    truth = [
        {"truth_id": "DEMO-T01", "x": 20, "y": 31, "channel": 17, "peak_snr": 7.5, "spectral_sigma": 1.7},
        {"truth_id": "DEMO-T02", "x": 43, "y": 21, "channel": 38, "peak_snr": 3.0, "spectral_sigma": 2.5},
        {"truth_id": "DEMO-T03", "x": 34, "y": 45, "channel": 56, "peak_snr": 2.4, "spectral_sigma": 3.4},
    ]
    for source in truth:
        add_line(cube, **{key: value for key, value in source.items() if key != "truth_id"})

    yy, xx = np.indices(shape[1:])
    radius = np.hypot(xx - 31.5, yy - 31.5)
    pb = np.exp(-0.5 * (radius / 25.0) ** 2)

    wcs = WCS(naxis=4)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN", "FREQ", "STOKES"]
    wcs.wcs.cunit = ["deg", "deg", "Hz", ""]
    wcs.wcs.crpix = [32.5, 32.5, 1.0, 1.0]
    wcs.wcs.crval = [334.5, 0.25, 100.0e9, 1.0]
    wcs.wcs.cdelt = [-0.2 / 3600.0, 0.2 / 3600.0, 15.625e6, 1.0]
    header = wcs.to_header()
    header["BUNIT"] = "Jy/beam"
    header["BMAJ"] = 0.6 / 3600.0
    header["BMIN"] = 0.5 / 3600.0
    header["BPA"] = 0.0
    header["OBJECT"] = "DICE synthetic line-search demo"
    fits.PrimaryHDU(cube[np.newaxis].astype("float32"), header=header).writeto(
        OUT / "demo.contsub.dirty.fits", overwrite=True
    )

    pb_wcs = wcs.celestial
    pb_header = pb_wcs.to_header()
    pb_header["BUNIT"] = ""
    pb_header["OBJECT"] = "DICE synthetic primary beam"
    fits.PrimaryHDU(pb.astype("float32"), header=pb_header).writeto(
        OUT / "demo.pb.fits", overwrite=True
    )

    with (OUT / "injected_truth.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=truth[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(truth)

    print(f"Wrote deterministic demo inputs to {OUT}")


if __name__ == "__main__":
    main()
