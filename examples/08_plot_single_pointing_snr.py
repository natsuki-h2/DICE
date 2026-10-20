"""Plot S/N-dependent completeness and contamination from Monte Carlo CSVs.

This script does not require CASA. It reads the summary tables produced by
``07_casa_single_pointing_monte_carlo.py`` and creates differential and
cumulative curves split by primary-beam response.
"""

from pathlib import Path
import argparse
import csv

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = REPO_ROOT / "data" / "casa_single_pointing_monte_carlo"
PB_BINS = [(0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.01)]
PB_LABELS = ["0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"]
COLORS = ["#4477AA", "#66CCEE", "#228833", "#CCBB44"]
COMPLETENESS_EDGES = np.arange(2.0, 10.01, 1.0)
DETECTION_EDGES = np.arange(5.0, 10.01, 0.5)


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def as_bool(value):
    return str(value).lower() == "true"


def in_pb(row, lower, upper):
    return lower <= float(row["pb"]) < upper


def fraction(rows, success):
    if not rows:
        return np.nan
    return sum(bool(success(row)) for row in rows) / float(len(rows))


def realization_ids(source_rows):
    return sorted({int(row["realization"]) for row in source_rows})


def completeness_curves(source_rows, cumulative):
    curves = {}
    realizations = realization_ids(source_rows)
    for pb_index, (pb_lower, pb_upper) in enumerate(PB_BINS):
        by_realization = []
        numerator_counts = []
        denominator_counts = []
        for realization in realizations:
            base = [
                row
                for row in source_rows
                if int(row["realization"]) == realization
                and in_pb(row, pb_lower, pb_upper)
            ]
            values = []
            numerators = []
            denominators = []
            if cumulative:
                x_values = COMPLETENESS_EDGES[:-1]
                for threshold in x_values:
                    selected = [
                        row for row in base if float(row["expected_snr"]) >= threshold
                    ]
                    values.append(fraction(selected, lambda row: as_bool(row["recovered"])))
                    numerators.append(sum(as_bool(row["recovered"]) for row in selected))
                    denominators.append(len(selected))
            else:
                x_values = 0.5 * (COMPLETENESS_EDGES[:-1] + COMPLETENESS_EDGES[1:])
                for lower, upper in zip(
                    COMPLETENESS_EDGES[:-1], COMPLETENESS_EDGES[1:]
                ):
                    selected = [
                        row
                        for row in base
                        if lower <= float(row["expected_snr"]) < upper
                    ]
                    values.append(fraction(selected, lambda row: as_bool(row["recovered"])))
                    numerators.append(sum(as_bool(row["recovered"]) for row in selected))
                    denominators.append(len(selected))
            by_realization.append(values)
            numerator_counts.append(numerators)
            denominator_counts.append(denominators)
        curves[pb_index] = (
            np.asarray(x_values),
            np.asarray(by_realization),
            np.asarray(numerator_counts),
            np.asarray(denominator_counts),
        )
    return curves


def q_curves(detection_rows, cumulative, estimator):
    """Return injection-matched or blank-derived Q curves."""
    curves = {}
    realizations = sorted({int(row["realization"]) for row in detection_rows})
    for pb_index, (pb_lower, pb_upper) in enumerate(PB_BINS):
        by_realization = []
        numerator_counts = []
        denominator_counts = []
        for realization in realizations:
            base = [
                row
                for row in detection_rows
                if int(row["realization"]) == realization
                and in_pb(row, pb_lower, pb_upper)
            ]
            populated = [row for row in base if row["detection_kind"] == "populated"]
            if estimator in {"injection", "truth"}:
                numerator = [row for row in populated if not as_bool(row["matched"])]
            elif estimator == "blank":
                numerator = [
                    row for row in base if row["detection_kind"] == "blank_positive"
                ]
            else:
                raise ValueError("estimator must be injection or blank")

            values = []
            if cumulative:
                x_values = DETECTION_EDGES[:-1]
                selections = [
                    (
                        [row for row in populated if float(row["peak_snr"]) >= threshold],
                        [row for row in numerator if float(row["peak_snr"]) >= threshold],
                    )
                    for threshold in x_values
                ]
            else:
                x_values = 0.5 * (DETECTION_EDGES[:-1] + DETECTION_EDGES[1:])
                selections = [
                    (
                        [
                            row
                            for row in populated
                            if lower <= float(row["peak_snr"]) < upper
                        ],
                        [
                            row
                            for row in numerator
                            if lower <= float(row["peak_snr"]) < upper
                        ],
                    )
                    for lower, upper in zip(DETECTION_EDGES[:-1], DETECTION_EDGES[1:])
                ]
            for denominator, selected_numerator in selections:
                values.append(
                    len(selected_numerator) / float(len(denominator))
                    if denominator
                    else np.nan
                )
            numerator_counts.append([len(item[1]) for item in selections])
            denominator_counts.append([len(item[0]) for item in selections])
            by_realization.append(values)
        curves[pb_index] = (
            np.asarray(x_values),
            np.asarray(by_realization),
            np.asarray(numerator_counts),
            np.asarray(denominator_counts),
        )
    return curves


def summarize(values, numerators, denominators):
    valid_count = np.sum(np.isfinite(values), axis=0)
    mean = np.full(values.shape[1], np.nan)
    lower = np.full(values.shape[1], np.nan)
    upper = np.full(values.shape[1], np.nan)
    for column in range(values.shape[1]):
        finite = values[:, column][np.isfinite(values[:, column])]
        if finite.size:
            total_denominator = np.sum(denominators[:, column])
            if total_denominator > 0:
                mean[column] = (
                    np.sum(numerators[:, column]) / float(total_denominator)
                )
            lower[column], upper[column] = np.percentile(finite, [16, 84])
    return mean, lower, upper, valid_count


def plot_family(
    ax, curves, xlabel, ylabel, individual_style="-", aggregate_style="-"
):
    for pb_index, (x_values, values, numerators, denominators) in curves.items():
        for realization_values in values:
            ax.plot(
                x_values,
                realization_values,
                color=COLORS[pb_index],
                alpha=0.12,
                linewidth=0.8,
                linestyle=individual_style,
                marker="o",
                markersize=2.5,
            )
        mean, lower, upper, valid_count = summarize(values, numerators, denominators)
        usable = valid_count > 0
        ax.plot(
            x_values[usable],
            mean[usable],
            color=COLORS[pb_index],
            linewidth=2.0,
            linestyle=aggregate_style,
            marker="o",
            markersize=4.0,
        )
        ax.fill_between(
            x_values[usable],
            lower[usable],
            upper[usable],
            color=COLORS[pb_index],
            alpha=0.16,
            linewidth=0,
        )
    ax.set_ylim(-0.04, 1.04)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.2)


def write_metric_table(path, families):
    rows = []
    for metric, curves in families.items():
        for pb_index, (x_values, values, numerators, denominators) in curves.items():
            mean, lower, upper, valid_count = summarize(
                values, numerators, denominators
            )
            for x_value, y, lo, hi, count in zip(
                x_values, mean, lower, upper, valid_count
            ):
                rows.append(
                    {
                        "metric": metric,
                        "pb_bin": PB_LABELS[pb_index],
                        "snr": x_value,
                        "mean": y,
                        "p16": lo,
                        "p84": hi,
                        "n_valid_realizations": int(count),
                    }
                )
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--title",
        default="BliSS-ALMA single-pointing smoke test: 10 realizations",
    )
    args = parser.parse_args()
    input_path = args.input.resolve()
    output_path = (
        args.output.resolve() if args.output is not None else input_path / "figures"
    )

    source_rows = read_csv(input_path / "source_recovery_all.csv")
    detection_rows = read_csv(input_path / "detections_all.csv")
    differential_completeness = completeness_curves(source_rows, cumulative=False)
    cumulative_completeness = completeness_curves(source_rows, cumulative=True)
    differential_truth = q_curves(detection_rows, cumulative=False, estimator="injection")
    differential_blank = q_curves(detection_rows, cumulative=False, estimator="blank")
    cumulative_truth = q_curves(detection_rows, cumulative=True, estimator="injection")
    cumulative_blank = q_curves(detection_rows, cumulative=True, estimator="blank")

    output_path.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.6))
    plot_family(
        axes[0, 0], differential_completeness, "Expected S/N bin", "Completeness"
    )
    plot_family(
        axes[0, 1],
        cumulative_completeness,
        "Expected S/N threshold (≥)",
        "Completeness",
    )
    plot_family(
        axes[1, 0],
        differential_truth,
        "Measured peak S/N bin",
        "Contamination Q",
        aggregate_style="-",
    )
    plot_family(
        axes[1, 0],
        differential_blank,
        "Measured peak S/N bin",
        "Contamination Q",
        aggregate_style="--",
    )
    plot_family(
        axes[1, 1],
        cumulative_truth,
        "Measured peak S/N threshold (≥)",
        "Contamination Q",
        aggregate_style="-",
    )
    plot_family(
        axes[1, 1],
        cumulative_blank,
        "Measured peak S/N threshold (≥)",
        "Contamination Q",
        aggregate_style="--",
    )

    axes[0, 0].set_title("Differential completeness")
    axes[0, 1].set_title("Cumulative completeness")
    axes[1, 0].set_title("Differential contamination")
    axes[1, 1].set_title("Cumulative contamination")
    pb_handles = [
        Line2D([0], [0], color=color, linewidth=2, label="PB " + label)
        for color, label in zip(COLORS, PB_LABELS)
    ]
    q_handles = [
        Line2D([0], [0], color="0.25", linewidth=2, label="Q injection-matched"),
        Line2D([0], [0], color="0.25", linewidth=2, linestyle="--", label="Q blank"),
    ]
    fig.legend(
        handles=pb_handles + q_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.015),
        ncol=3,
        frameon=False,
    )
    fig.suptitle(args.title, y=0.975)
    fig.subplots_adjust(
        left=0.09, right=0.98, bottom=0.17, top=0.87, hspace=0.40, wspace=0.28
    )
    fig.savefig(output_path / "snr_pb_metrics.png", dpi=180, bbox_inches="tight")
    fig.savefig(output_path / "snr_pb_metrics.pdf", bbox_inches="tight")
    plt.close(fig)

    write_metric_table(
        output_path / "snr_pb_metrics.csv",
        {
            "completeness_differential": differential_completeness,
            "completeness_cumulative": cumulative_completeness,
            "q_truth_differential": differential_truth,
            "q_blank_differential": differential_blank,
            "q_truth_cumulative": cumulative_truth,
            "q_blank_cumulative": cumulative_blank,
        },
    )
    print("Wrote {0}".format(output_path / "snr_pb_metrics.png"))
    print("Wrote {0}".format(output_path / "snr_pb_metrics.pdf"))
    print("Wrote {0}".format(output_path / "snr_pb_metrics.csv"))


if __name__ == "__main__":
    main()
