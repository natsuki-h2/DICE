"""Fit paper-level completeness/contamination curves and classify S/N ranges.

This is a reproducibility utility, not part of the stable DICE API.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf


def increasing_erf(x, midpoint, width):
    return 0.5 * (1.0 + erf((np.asarray(x) - midpoint) / width))


def decreasing_erf(x, midpoint, width):
    return 0.5 * (1.0 - erf((np.asarray(x) - midpoint) / width))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def fit_curve(x, y, function):
    valid = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(valid) < 3:
        return np.nan, np.nan
    xv, yv = x[valid], y[valid]
    midpoint = float(xv[np.argmin(np.abs(yv - 0.5))])
    span = max(float(np.ptp(xv)), 0.1)
    parameters, _ = curve_fit(
        function,
        xv,
        yv,
        p0=(midpoint, span / 4.0),
        bounds=([float(np.min(xv)) - span, 1e-6],
                [float(np.max(xv)) + span, 10.0 * span]),
        maxfev=20000,
    )
    return float(parameters[0]), float(parameters[1])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("fit_csv", type=Path)
    parser.add_argument("classification_csv", type=Path)
    parser.add_argument("--x-column", default="x")
    parser.add_argument("--completeness-column", default="completeness")
    parser.add_argument("--contamination-column", default="contamination")
    parser.add_argument("--group-columns", default="")
    args = parser.parse_args()

    groups = [item for item in args.group_columns.split(",") if item]
    rows = read_rows(args.input_csv)
    required = {args.x_column, args.completeness_column,
                args.contamination_column, *groups}
    missing = required - (set(rows[0]) if rows else set())
    if missing:
        raise ValueError(f"Missing input columns: {sorted(missing)}")

    grouped: dict[tuple[str, ...], list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(tuple(row[name] for name in groups), []).append(row)

    fit_rows, class_rows = [], []
    for key, subset in sorted(grouped.items()):
        subset.sort(key=lambda row: float(row[args.x_column]))
        x = np.array([float(row[args.x_column]) for row in subset])
        completeness = np.array([
            float(row[args.completeness_column]) for row in subset
        ])
        contamination = np.array([
            float(row[args.contamination_column]) for row in subset
        ])
        c_mid, c_width = fit_curve(x, completeness, increasing_erf)
        q_mid, q_width = fit_curve(x, contamination, decreasing_erf)
        identity = dict(zip(groups, key))
        fit_rows.append({
            **identity,
            "completeness_midpoint": c_mid,
            "completeness_width": c_width,
            "contamination_midpoint": q_mid,
            "contamination_width": q_width,
        })
        lower, upper = sorted((c_mid, q_mid))
        for value in x:
            category = "I" if value >= upper else "II" if value >= lower else "III"
            class_rows.append({**identity, args.x_column: value,
                               "category": category,
                               "lower_transition": lower,
                               "upper_transition": upper})

    fit_fields = groups + ["completeness_midpoint", "completeness_width",
                           "contamination_midpoint", "contamination_width"]
    class_fields = groups + [args.x_column, "category", "lower_transition",
                             "upper_transition"]
    write_rows(args.fit_csv, fit_rows, fit_fields)
    write_rows(args.classification_csv, class_rows, class_fields)


if __name__ == "__main__":
    main()

