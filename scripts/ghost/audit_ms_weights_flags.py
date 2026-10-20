"""Audit visibility weights, sigmas, and flags in the downloaded F2 MS.

Run with CASA Python, for example::

    casa --nologger --nogui -c scripts/ghost/28_audit_ghost_ms_weights_flags.py
"""

from pathlib import Path
import json

import numpy as np
from casatools import table


ROOT = Path(__file__).resolve().parents[2]
MS = (ROOT / "data" / "ghost_matrix" / "f2_5min" / "extracted" /
      "output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms")
OUT = ROOT / "data" / "ghost_matrix" / "f2_5min" / "ms_weight_flag_audit.json"
CHUNK_ROWS = 100_000


def percentiles(values):
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    return {
        "min": float(np.min(array)),
        "p16": float(np.percentile(array, 16)),
        "median": float(np.median(array)),
        "p84": float(np.percentile(array, 84)),
        "max": float(np.max(array)),
        "mean": float(np.mean(array)),
        "std": float(np.std(array)),
    }


def main():
    tb = table()
    tb.open(str(MS))
    try:
        columns = tb.colnames()
        nrows = tb.nrows()
        weights = []
        sigmas = []
        uv_distance = []
        data_real = []
        data_imag = []
        n_flagged = 0
        n_samples = 0
        n_flag_rows = 0
        for start in range(0, nrows, CHUNK_ROWS):
            count = min(CHUNK_ROWS, nrows - start)
            flag = tb.getcol("FLAG", startrow=start, nrow=count)
            flag_row = tb.getcol("FLAG_ROW", startrow=start, nrow=count)
            weight = tb.getcol("WEIGHT", startrow=start, nrow=count)
            sigma = tb.getcol("SIGMA", startrow=start, nrow=count)
            uvw = tb.getcol("UVW", startrow=start, nrow=count)
            data = tb.getcol("DATA", startrow=start, nrow=count)
            n_flagged += int(np.count_nonzero(flag))
            n_samples += int(flag.size)
            n_flag_rows += int(np.count_nonzero(flag_row))
            weights.append(weight.ravel())
            sigmas.append(sigma.ravel())
            uv_distance.append(np.hypot(uvw[0], uvw[1]))
            valid = ~flag
            data_real.append(data.real[valid])
            data_imag.append(data.imag[valid])
    finally:
        tb.close()

    weight = np.concatenate(weights)
    sigma = np.concatenate(sigmas)
    relation = weight * sigma**2
    result = {
        "measurement_set": str(MS.relative_to(ROOT)),
        "n_rows": int(nrows),
        "columns": columns,
        "flagged_visibility_fraction": n_flagged / n_samples,
        "flagged_row_fraction": n_flag_rows / nrows,
        "weight": percentiles(weight),
        "sigma": percentiles(sigma),
        "weight_times_sigma_squared": percentiles(relation),
        "uv_distance_m": percentiles(np.concatenate(uv_distance)),
        "unflagged_data_real_jy": percentiles(np.concatenate(data_real)),
        "unflagged_data_imag_jy": percentiles(np.concatenate(data_imag)),
        "interpretation": {
            "weight_sigma_consistent": bool(np.allclose(relation, 1.0, rtol=1e-5, atol=1e-8)),
            "weights_effectively_constant": bool(np.std(weight) <= 1e-6 * max(np.mean(weight), 1e-30)),
            "contains_flags": bool(n_flagged > 0 or n_flag_rows > 0),
        },
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
