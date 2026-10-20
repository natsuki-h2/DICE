"""Plot common-resolution completeness and contamination for F1--F6."""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "ghost_matrix" / "common_resolution_runs"
OUT = INPUT / "figures"
CONFIGS = [f"F{index}" for index in range(1, 7)]
COLORS = ["#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377"]
TARGETS = np.arange(2.0, 7.01, 0.5)
THRESHOLDS = np.arange(2.0, 7.01, 0.5)
N_RESAMPLES = 10_000


def load(configuration, name):
    path = INPUT / f"{configuration.lower()}_5min" / name
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(2040)
    metrics = []
    fig, (ax_c, ax_q, ax_cq) = plt.subplots(1, 3, figsize=(16.5, 4.8), constrained_layout=True)
    for configuration, color in zip(CONFIGS, COLORS):
        sources = load(configuration, "source_recovery_all.csv")
        detections = load(configuration, "detections_all.csv")
        realizations = sorted({int(row["realization"]) for row in sources})
        n_real = len(realizations)
        draws = rng.integers(0, n_real, size=(N_RESAMPLES, n_real))

        c_num = np.zeros((n_real, len(TARGETS)), dtype=int)
        c_den = np.zeros_like(c_num)
        q_num = np.zeros((n_real, len(THRESHOLDS)), dtype=int)
        q_den = np.zeros_like(q_num)
        for ri, realization in enumerate(realizations):
            for si, snr in enumerate(TARGETS):
                cell = [row for row in sources
                        if int(row["realization"]) == realization
                        and float(row["expected_snr"]) == snr]
                c_num[ri, si] = sum(row["recovered"].lower() == "true" for row in cell)
                c_den[ri, si] = len(cell)
            positive = [row for row in detections
                        if int(row["realization"]) == realization
                        and row["detection_kind"] == "populated_diagnostic"]
            negative = [row for row in detections
                        if int(row["realization"]) == realization
                        and row["detection_kind"] == "blank_negative"]
            for ti, threshold in enumerate(THRESHOLDS):
                q_den[ri, ti] = sum(float(row["peak_snr"]) >= threshold for row in positive)
                q_num[ri, ti] = sum(float(row["peak_snr"]) >= threshold for row in negative)

        c = c_num.sum(0) / c_den.sum(0)
        c_boot = c_num[draws].sum(1) / c_den[draws].sum(1)
        c_lo, c_hi = np.percentile(c_boot, [16, 84], axis=0)
        q = np.minimum(q_num.sum(0) / np.maximum(q_den.sum(0), 1), 1.0)
        q_boot = np.minimum(q_num[draws].sum(1) / np.maximum(q_den[draws].sum(1), 1), 1.0)
        q_lo, q_hi = np.percentile(q_boot, [16, 84], axis=0)

        ax_c.plot(TARGETS, c, "o-", color=color, label=configuration, linewidth=1.6, markersize=4)
        ax_c.fill_between(TARGETS, c_lo, c_hi, color=color, alpha=0.10)
        ax_q.plot(THRESHOLDS, q, "o-", color=color, label=configuration, linewidth=1.6, markersize=4)
        ax_q.fill_between(THRESHOLDS, q_lo, q_hi, color=color, alpha=0.10)
        ax_cq.plot(c, q, "o-", color=color, label=configuration, linewidth=1.6, markersize=4)
        for snr, value, lo, hi, num, den in zip(TARGETS, c, c_lo, c_hi, c_num.sum(0), c_den.sum(0)):
            metrics.append({"configuration": configuration, "metric": "differential_completeness",
                            "snr": snr, "numerator": num, "denominator": den,
                            "value": value, "lower_16": lo, "upper_84": hi})
        for snr, value, lo, hi, num, den in zip(THRESHOLDS, q, q_lo, q_hi, q_num.sum(0), q_den.sum(0)):
            metrics.append({"configuration": configuration, "metric": "cumulative_contamination",
                            "snr": snr, "numerator": num, "denominator": den,
                            "value": value, "lower_16": lo, "upper_84": hi})

    ax_c.set(title="Differential completeness", xlabel="Input apparent S/N", ylabel="C",
             ylim=(-0.04, 1.04))
    ax_q.set(title="Cumulative contamination", xlabel="Measured peak S/N threshold", ylabel="Q",
             ylim=(-0.04, 1.04))
    ax_q.axvline(4.0, color="#555555", linestyle="--", linewidth=1.0)
    ax_cq.set(title="Completeness–contamination relation",
              xlabel="Differential completeness C", ylabel="Cumulative contamination Q",
              xlim=(-0.04, 1.04), ylim=(-0.04, 1.04))
    ax_cq.text(0.04, 0.05, "S/N is the curve parameter", transform=ax_cq.transAxes, fontsize=9)
    for ax in (ax_c, ax_q, ax_cq):
        ax.grid(alpha=0.22); ax.legend(title="Configuration", fontsize=8, ncol=2)
    fig.suptitle("GHOST F1–F6: 10 realizations at common 1-arcsec resolution")
    fig.savefig(OUT / "common_resolution_matrix_metrics.png", dpi=200)
    fig.savefig(OUT / "common_resolution_matrix_metrics.pdf")
    plt.close(fig)
    with (OUT / "common_resolution_matrix_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics[0]))
        writer.writeheader(); writer.writerows(metrics)
    print(OUT / "common_resolution_matrix_metrics.png")


if __name__ == "__main__":
    main()
