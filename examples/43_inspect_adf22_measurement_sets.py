"""Read-only comparison of one or more measurement sets."""

import argparse
from pathlib import Path
import json

import numpy as np
from casatools import msmetadata as msmetadata_tool
from casatools import table as table_tool


def inspect(path):
    tb = table_tool()
    msmd = msmetadata_tool()
    tb.open(str(path))
    columns = tb.colnames()
    nrows = tb.nrows()
    time = tb.getcol("TIME")
    data_shape = None
    for column in ("CORRECTED_DATA", "DATA"):
        if column in columns:
            cell = tb.getcell(column, 0)
            data_shape = list(cell.shape)
            break
    tb.close()

    msmd.open(str(path))
    nspw = msmd.nspw()
    spws = []
    for spw in range(nspw):
        frequencies = np.asarray(msmd.chanfreqs(spw), dtype=float)
        spws.append({
            "spw": spw,
            "nchan": int(msmd.nchan(spw)),
            "frequency_min_ghz": float(np.min(frequencies) / 1e9),
            "frequency_max_ghz": float(np.max(frequencies) / 1e9),
            "channel_width_mhz": float(
                np.median(np.abs(np.diff(frequencies))) / 1e6
            ) if frequencies.size > 1 else None,
        })
    result = {
        "path": str(path),
        "size_gib": sum(
            item.stat().st_size for item in path.rglob("*") if item.is_file()
        ) / 1024**3,
        "nrows": nrows,
        "data_cell_shape_polarization_channel": data_shape,
        "columns": columns,
        "nfields": int(msmd.nfields()),
        "field_names": list(msmd.fieldnames()),
        "nspw": nspw,
        "spws": spws,
        "nantennas": int(msmd.nantennas()),
        "nscans": int(msmd.nscans()),
        "time_min_mjd_seconds": float(np.min(time)),
        "time_max_mjd_seconds": float(np.max(time)),
        "duration_hours": float((np.max(time) - np.min(time)) / 3600),
    }
    msmd.close()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "measurement_sets",
        nargs="+",
        type=Path,
        help="Measurement-set paths to inspect without modifying them.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("measurement_set_audit.json"),
        help="Output JSON path (default: measurement_set_audit.json).",
    )
    args = parser.parse_args()

    missing = [str(path) for path in args.measurement_sets if not path.exists()]
    if missing:
        parser.error("Measurement set not found: " + ", ".join(missing))

    results = [inspect(path) for path in args.measurement_sets]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))
    print("Wrote", args.output)


if __name__ == "__main__":
    main()
