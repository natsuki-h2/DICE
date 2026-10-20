"""Inspect the downloaded GHOST F1 measurement set and dirty products.

Run from the repository root with CASA 6. The resulting JSON summary is kept
under ``data/ghost_f1_5min/`` with the downloaded products.
"""

from pathlib import Path
import json
import math

import numpy as np
from casatasks import imhead, imstat
from casatools import table


REPO_ROOT = Path(__file__).resolve().parents[1]
WORKDIR = REPO_ROOT / "data" / "ghost_f1_5min"
MS = WORKDIR / "output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms"
DIRTY_IMAGE = WORKDIR / "ALMA-GHOST-CCK6FJ4E.dirty.image.fits"
DIRTY_PSF = WORKDIR / "ALMA-GHOST-CCK6FJ4E.dirty.psf.fits"
SUMMARY = WORKDIR / "f1_inspection.json"
LIGHT_SPEED_M_S = 299792458.0


def read_table(path, columns):
    tb = table()
    try:
        tb.open(str(path))
        return {column: tb.getcol(column) for column in columns}
    finally:
        tb.done()


def inspect_ms():
    tb = table()
    try:
        tb.open(str(MS))
        columns = list(tb.colnames())
        n_rows = int(tb.nrows())
        uvw = np.asarray(tb.getcol("UVW"), dtype=np.float64)
        interval = np.asarray(tb.getcol("INTERVAL"), dtype=np.float64)
        exposure = np.asarray(tb.getcol("EXPOSURE"), dtype=np.float64)
        times = np.asarray(tb.getcol("TIME"), dtype=np.float64)
        flags = np.asarray(tb.getcol("FLAG"), dtype=bool)
        data_shape = list(tb.getcolshapestring("DATA", 0, 1))[0]
    finally:
        tb.done()

    antenna = read_table(MS / "ANTENNA", ["NAME", "DISH_DIAMETER"])
    spw = read_table(
        MS / "SPECTRAL_WINDOW",
        ["CHAN_FREQ", "CHAN_WIDTH", "TOTAL_BANDWIDTH", "NUM_CHAN"],
    )
    field = read_table(MS / "FIELD", ["NAME", "PHASE_DIR"])

    phase_dir = np.asarray(field["PHASE_DIR"], dtype=np.float64).reshape(2, -1)
    ra_rad = float(phase_dir[0, 0])
    dec_rad = float(phase_dir[1, 0])
    center_freq_hz = float(np.mean(spw["CHAN_FREQ"]))
    uv_distance_m = np.sqrt(uvw[0] ** 2 + uvw[1] ** 2)
    uv_distance_lambda = uv_distance_m * center_freq_hz / LIGHT_SPEED_M_S

    return {
        "path": str(MS),
        "n_rows": n_rows,
        "columns": columns,
        "data_shape_per_row": data_shape,
        "n_antennas": int(len(antenna["NAME"])),
        "antenna_names_first_last": [
            str(antenna["NAME"][0]),
            str(antenna["NAME"][-1]),
        ],
        "dish_diameters_m": sorted(
            float(value) for value in np.unique(antenna["DISH_DIAMETER"])
        ),
        "n_spw": int(np.asarray(spw["NUM_CHAN"]).size),
        "n_channels": [int(value) for value in np.ravel(spw["NUM_CHAN"])],
        "center_frequency_hz": center_freq_hz,
        "channel_width_hz": [
            float(value) for value in np.ravel(spw["CHAN_WIDTH"])
        ],
        "total_bandwidth_hz": [
            float(value) for value in np.ravel(spw["TOTAL_BANDWIDTH"])
        ],
        "field_names": [str(value) for value in field["NAME"]],
        "phase_centre_ra_deg": math.degrees(ra_rad) % 360.0,
        "phase_centre_dec_deg": math.degrees(dec_rad),
        "interval_s": sorted(float(value) for value in np.unique(interval)),
        "exposure_s": sorted(float(value) for value in np.unique(exposure)),
        "time_span_s": float(np.max(times) - np.min(times)),
        "uv_distance_m_min_max": [
            float(np.min(uv_distance_m)),
            float(np.max(uv_distance_m)),
        ],
        "uv_distance_klambda_min_max": [
            float(np.min(uv_distance_lambda) / 1000.0),
            float(np.max(uv_distance_lambda) / 1000.0),
        ],
        "flagged_fraction": float(np.count_nonzero(flags) / flags.size),
    }


def inspect_fits(path):
    stats = imstat(str(path), listit=False)
    return {
        "path": str(path),
        "shape": [int(value) for value in imhead(str(path), mode="get", hdkey="shape")],
        "unit": str(imhead(str(path), mode="get", hdkey="bunit")),
        "minimum": float(stats["min"][0]),
        "maximum": float(stats["max"][0]),
        "rms": float(stats["rms"][0]),
        "beam_major": imhead(str(path), mode="get", hdkey="beammajor"),
        "beam_minor": imhead(str(path), mode="get", hdkey="beamminor"),
        "beam_pa": imhead(str(path), mode="get", hdkey="beampa"),
    }


def main():
    summary = {
        "measurement_set": inspect_ms(),
        "ghost_dirty_image": inspect_fits(DIRTY_IMAGE),
        "ghost_dirty_psf": inspect_fits(DIRTY_PSF),
    }
    with SUMMARY.open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Wrote {SUMMARY}")


if __name__ == "__main__":
    main()
