"""Combine Hayatsu-2025 blank pilots and plot Field10 contamination.

The uncertainty unit is one complete visibility-noise realization.  The
resulting Field10 contamination is a representative-field validation; it is
kept distinct from the published full-mosaic contamination curve.
"""

import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import beta

from bliss_alma.provenance import write_provenance


ROOT = Path(__file__).resolve().parents[2]
INPUT = (
    ROOT / "data/hayatsu_reproduction/hayatsu2025_blank_realizations"
)
OBSERVED = (
    ROOT / "data/hayatsu_reproduction/hayatsu2025_reproduction"
    / "hayatsu2025_field10_contamination_diagnostic.csv"
)
OUTPUT = (
    ROOT / "data/hayatsu_reproduction/hayatsu2025_blank_realizations"
)
WIDTHS = tuple(range(2, 24, 2))
SELECTED_WIDTHS = (2, 6, 14, 22)

def read_rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def write_rows(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def main():
    paths = [INPUT / "blank_clump_counts.csv"]
    paths.extend(sorted(INPUT.glob("*/blank_clump_counts.csv")))
    all_rows = []
    for path in paths:
        if path.exists():
            all_rows.extend(read_rows(path))

    keyed = {}
    for row in all_rows:
        key = (
            int(row["realization"]),
            int(row["width_channels"]),
            float(row["snr_threshold"]),
            row["polarity"],
        )
        keyed[key] = row
    rows = list(keyed.values())
    realizations = sorted({int(row["realization"]) for row in rows})
    if len(realizations) < 2:
        raise RuntimeError("At least two blank realizations are required")

    observed = {}
    for row in read_rows(OBSERVED):
        observed[
            (int(row["width_channels"]),
             float(row["measured_snr_threshold"]))
        ] = int(row["observed_field10_clumps"])

    summary = []
    thresholds = sorted({float(row["snr_threshold"]) for row in rows})
    for width in WIDTHS:
        for threshold in thresholds:
            counts = {}
            for polarity in ("positive", "negative"):
                values = np.array([
                    int(row["n_clumps"]) for row in rows
                    if int(row["width_channels"]) == width
                    and float(row["snr_threshold"]) == threshold
                    and row["polarity"] == polarity
                ], dtype=float)
                counts[polarity] = values

            n_real = observed.get((width, threshold), 0)
            for polarity, values in counts.items():
                q_values = (
                    np.minimum(1.0, values / n_real)
                    if n_real > 0 else np.full(values.shape, np.nan)
                )
                summary.append({
                    "width_channels": width,
                    "snr_threshold": threshold,
                    "polarity": polarity,
                    "n_realizations": len(values),
                    "mean_blank_clumps": float(np.mean(values)),
                    "blank_clumps_p16": float(np.percentile(values, 16)),
                    "blank_clumps_p84": float(np.percentile(values, 84)),
                    "blank_occurrence_fraction": float(np.mean(values > 0)),
                    "blank_occurrence_p16": float(beta.ppf(
                        0.16, np.sum(values > 0) + 0.5,
                        np.sum(values == 0) + 0.5,
                    )),
                    "blank_occurrence_p84": float(beta.ppf(
                        0.84, np.sum(values > 0) + 0.5,
                        np.sum(values == 0) + 0.5,
                    )),
                    "observed_field10_clumps": n_real,
                    "field10_q_median": (
                        float(np.nanmedian(q_values))
                        if np.any(np.isfinite(q_values)) else ""
                    ),
                    "field10_q_p16": (
                        float(np.nanpercentile(q_values, 16))
                        if np.any(np.isfinite(q_values)) else ""
                    ),
                    "field10_q_p84": (
                        float(np.nanpercentile(q_values, 84))
                        if np.any(np.isfinite(q_values)) else ""
                    ),
                })

    write_rows(OUTPUT / "blank_clump_summary.csv", summary)

    colors = plt.cm.viridis(np.linspace(0.12, 0.88, len(SELECTED_WIDTHS)))
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.3))
    for width, color in zip(SELECTED_WIDTHS, colors):
        for polarity, linestyle, marker in (
            ("positive", "-", "o"), ("negative", "--", "v")
        ):
            selected = [
                row for row in summary
                if row["width_channels"] == width
                and row["polarity"] == polarity
            ]
            x = np.array([row["snr_threshold"] for row in selected])
            mean = np.array([row["mean_blank_clumps"] for row in selected])
            low = np.array([row["blank_clumps_p16"] for row in selected])
            high = np.array([row["blank_clumps_p84"] for row in selected])
            label = f"{width} ch" if polarity == "positive" else None
            axes[0].plot(
                x, mean, marker=marker, linestyle=linestyle,
                color=color, label=label,
            )
            axes[0].fill_between(x, low, high, color=color, alpha=0.10)

            occurrence = np.array([
                row["blank_occurrence_fraction"] for row in selected
            ])
            occurrence_low = np.array([
                row["blank_occurrence_p16"] for row in selected
            ])
            occurrence_high = np.array([
                row["blank_occurrence_p84"] for row in selected
            ])
            axes[1].plot(
                x, occurrence, marker=marker, linestyle=linestyle,
                color=color, label=label,
            )
            axes[1].fill_between(
                x, occurrence_low, occurrence_high, color=color, alpha=0.10
            )

    axes[0].set(
        xlabel="Measured S/N threshold",
        ylabel="Mean blank positive-clump count",
        yscale="log",
        title=f"Field10 blank realizations (N={len(realizations)})",
    )
    axes[1].set(
        xlabel="Measured S/N threshold",
        ylabel="Fraction with at least one blank clump",
        ylim=(-0.03, 1.03),
        title="Field10 false-positive occurrence",
    )
    axes[1].text(
        0.03, 0.04,
        "Solid/circles: positive; dashed/triangles: negative",
        transform=axes[1].transAxes, fontsize=8,
    )
    for axis in axes:
        axis.grid(alpha=0.2)
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    png = OUTPUT / "hayatsu2025_blank_field10_contamination.png"
    pdf = OUTPUT / "hayatsu2025_blank_field10_contamination.pdf"
    csv_output = OUTPUT / "blank_clump_summary.csv"
    fig.savefig(png, dpi=220)
    fig.savefig(pdf)
    write_provenance(
        OUTPUT / "hayatsu2025_blank_field10_contamination.provenance.json",
        root=ROOT, script=Path(__file__),
        command=[sys.executable, str(Path(__file__))],
        inputs=paths + [OBSERVED], outputs=[csv_output, png, pdf],
        parameters={
            "width_channels": list(WIDTHS),
            "selected_width_channels": list(SELECTED_WIDTHS),
            "uncertainty_unit": "whole visibility-noise realization",
        },
        random_seeds=realizations,
    )
    print("realizations", len(realizations))
    print("wrote", csv_output)


if __name__ == "__main__":
    main()
