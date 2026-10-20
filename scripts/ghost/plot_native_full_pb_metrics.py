"""Plot native full-PB GHOST completeness and contamination with provenance."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

from bliss_alma.provenance import write_provenance


ROOT = Path(__file__).resolve().parents[2]
COLORS = ["#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377"]
MARKERS = ["o", "s", "^", "D", "v", "P"]
LINESTYLES = ["-", "--", "-.", ":", (0, (5, 2)), (0, (3, 1, 1, 1))]


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "configs/ghost/native_full_pb.json",
    )
    args = parser.parse_args()
    config_path = args.config.resolve()
    settings = json.loads(config_path.read_text())
    input_root = ROOT / settings["input_directory"]
    output_root = ROOT / settings["figure_directory"]
    configurations = settings["configurations"]
    snrs = np.arange(
        settings["snr_min"],
        settings["snr_max"] + settings["snr_step"] / 2,
        settings["snr_step"],
    )
    output_root.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(settings["resample_seed"])
    metrics: list[dict] = []
    inputs: list[Path] = []
    figure, axes = plt.subplots(2, 2, figsize=(12.7, 9.2), constrained_layout=True)
    ax_c, ax_q, ax_n, ax_cq = axes.ravel()

    for configuration, color, marker, linestyle in zip(
        configurations, COLORS, MARKERS, LINESTYLES
    ):
        source_path = input_root / configuration.lower() / "source_recovery_all.csv"
        count_path = input_root / configuration.lower() / "threshold_counts.csv"
        inputs.extend([source_path, count_path])
        sources = read_rows(source_path)
        counts = read_rows(count_path)
        realization_ids = sorted({int(row["realization"]) for row in sources})
        n_realizations = len(realization_ids)
        draws = rng.integers(
            0,
            n_realizations,
            size=(settings["n_resamples"], n_realizations),
        )
        c_num = np.zeros((n_realizations, len(snrs)), int)
        c_den = np.zeros_like(c_num)
        q_num = np.zeros_like(c_num)
        q_den = np.zeros_like(c_num)
        for row_index, realization in enumerate(realization_ids):
            source_rows = [
                row for row in sources if int(row["realization"]) == realization
            ]
            count_rows = [
                row for row in counts if int(row["realization"]) == realization
            ]
            for index, snr in enumerate(snrs):
                cell = [
                    row for row in source_rows if float(row["expected_snr"]) == snr
                ]
                c_den[row_index, index] = len(cell)
                c_num[row_index, index] = sum(
                    row["recovered"].lower() == "true" for row in cell
                )
                count = next(
                    row for row in count_rows if float(row["threshold_snr"]) == snr
                )
                q_num[row_index, index] = int(count["n_negative"])
                q_den[row_index, index] = int(count["n_positive"])

        completeness = c_num.sum(0) / c_den.sum(0)
        completeness_draws = c_num[draws].sum(1) / c_den[draws].sum(1)
        c_low, c_high = np.percentile(completeness_draws, [16, 84], axis=0)
        contamination = np.minimum(q_num.sum(0) / np.maximum(q_den.sum(0), 1), 1)
        contamination_draws = np.minimum(
            q_num[draws].sum(1) / np.maximum(q_den[draws].sum(1), 1), 1
        )
        q_low, q_high = np.percentile(contamination_draws, [16, 84], axis=0)
        style = {
            "color": color,
            "marker": marker,
            "linestyle": linestyle,
            "markersize": 4.8,
            "markeredgecolor": "#222222",
            "markeredgewidth": 0.35,
            "linewidth": 1.8,
            "label": configuration,
        }
        ax_c.plot(snrs, completeness, **style)
        ax_c.fill_between(snrs, c_low, c_high, color=color, alpha=0.1)
        ax_q.plot(snrs, contamination, **style)
        ax_q.fill_between(snrs, q_low, q_high, color=color, alpha=0.1)
        ax_n.plot(
            snrs,
            q_num.sum(0) / n_realizations,
            **style,
        )
        ax_cq.plot(completeness, contamination, **style)
        for snr, value, low, high, numerator, denominator in zip(
            snrs, completeness, c_low, c_high, c_num.sum(0), c_den.sum(0)
        ):
            metrics.append(
                {
                    "configuration": configuration,
                    "metric": "differential_completeness",
                    "snr": snr,
                    "numerator": numerator,
                    "denominator": denominator,
                    "value": value,
                    "lower_16": low,
                    "upper_84": high,
                }
            )
        for snr, value, low, high, numerator, denominator in zip(
            snrs, contamination, q_low, q_high, q_num.sum(0), q_den.sum(0)
        ):
            metrics.append(
                {
                    "configuration": configuration,
                    "metric": "cumulative_contamination",
                    "snr": snr,
                    "numerator": numerator,
                    "denominator": denominator,
                    "value": value,
                    "lower_16": low,
                    "upper_84": high,
                }
            )

    ax_c.set(
        title="Differential completeness",
        xlabel="Input apparent S/N",
        ylabel="C",
        ylim=(-0.04, 1.04),
    )
    ax_q.set(
        title="Cumulative contamination",
        xlabel="Measured peak S/N threshold",
        ylabel="Q",
        ylim=(-0.04, 1.04),
    )
    ax_n.set(
        title="Negative peaks per full-PB field",
        xlabel="Measured |peak| S/N threshold",
        ylabel="Mean count per realization",
        yscale="log",
    )
    ax_cq.set(
        title="Completeness–contamination relation",
        xlabel="Differential completeness C",
        ylabel="Cumulative contamination Q",
        xlim=(-0.04, 1.04),
        ylim=(-0.04, 1.04),
    )
    for axis in axes.ravel():
        if axis in (ax_q, ax_n):
            axis.axvline(
                settings["detection_threshold_snr"],
                color="#666666",
                linestyle="--",
                linewidth=0.8,
            )
        axis.grid(alpha=0.22)
        axis.legend(title="Configuration", fontsize=8, ncol=2)
    figure.suptitle(
        "GHOST F1–F6 native full-PB: "
        f"{n_realizations} realizations, "
        f"{settings['pixels_per_minor_beam']} pixels/minor beam, "
        f"{'no' if settings['uv_taper'] == 'none' else settings['uv_taper']} taper"
    )
    png = output_root / "native_full_pb_matrix_metrics.png"
    pdf = output_root / "native_full_pb_matrix_metrics.pdf"
    table = output_root / "native_full_pb_matrix_metrics.csv"
    figure.savefig(png, dpi=200)
    figure.savefig(pdf)
    plt.close(figure)
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics[0]))
        writer.writeheader()
        writer.writerows(metrics)

    write_provenance(
        output_root / "native_full_pb_matrix_metrics.provenance.json",
        root=ROOT,
        script=Path(__file__),
        command=sys.argv,
        config=config_path,
        inputs=inputs,
        outputs=[png, pdf, table],
        parameters=settings,
        random_seeds=[settings["resample_seed"]],
    )
    print(png)


if __name__ == "__main__":
    main()
