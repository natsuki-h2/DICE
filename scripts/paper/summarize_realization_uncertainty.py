"""Summarize a metric across complete independent noise realizations.

This is a reproducibility utility, not part of the stable DICE API.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--realization-column", default="realization")
    parser.add_argument("--group-columns", default="")
    parser.add_argument("--value-column")
    parser.add_argument("--numerator")
    parser.add_argument("--denominator")
    args = parser.parse_args()

    if bool(args.value_column) == bool(args.numerator and args.denominator):
        parser.error("use either --value-column or --numerator with --denominator")
    groups = [item for item in args.group_columns.split(",") if item]
    with args.input_csv.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {args.realization_column, *groups}
    required.add(args.value_column if args.value_column else args.numerator)
    if args.denominator:
        required.add(args.denominator)
    missing = required - (set(rows[0]) if rows else set())
    if missing:
        raise ValueError(f"Missing input columns: {sorted(missing)}")

    cells: dict[tuple[str, ...], dict[str, list[dict[str, str]]]] = {}
    for row in rows:
        key = tuple(row[name] for name in groups)
        cells.setdefault(key, {}).setdefault(row[args.realization_column], []).append(row)

    output = []
    for key, realizations in sorted(cells.items()):
        values = []
        for subset in realizations.values():
            if args.value_column:
                samples = [float(row[args.value_column]) for row in subset]
                values.append(float(np.mean(samples)))
            else:
                numerator = sum(float(row[args.numerator]) for row in subset)
                denominator = sum(float(row[args.denominator]) for row in subset)
                if denominator > 0:
                    values.append(numerator / denominator)
        if not values:
            continue
        low, median, high = np.percentile(values, [16, 50, 84])
        output.append({
            **dict(zip(groups, key)),
            "n_realizations": len(values),
            "median": float(median),
            "range_low": float(low),
            "range_high": float(high),
        })

    fields = groups + ["n_realizations", "median", "range_low", "range_high"]
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)


if __name__ == "__main__":
    main()
