"""Plot native full-PB completeness and contamination for F1--F6."""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/ghost_matrix/native_full_pb_runs"
OUT = INPUT / "figures"
CONFIGS = [f"F{i}" for i in range(1, 7)]
COLORS = ["#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377"]
SNRS = np.arange(2.0, 7.01, 0.5)
N_RESAMPLES = 10_000


def load(config, name):
    with (INPUT/config.lower()/name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(2040); metrics = []
    fig, axes = plt.subplots(2, 2, figsize=(12.7, 9.2), constrained_layout=True)
    ax_c, ax_q, ax_n, ax_cq = axes.ravel()
    for config, color in zip(CONFIGS, COLORS):
        sources = load(config, "source_recovery_all.csv")
        counts = load(config, "threshold_counts.csv")
        realizations = range(10)
        draws = rng.integers(0, 10, size=(N_RESAMPLES, 10))
        c_num = np.zeros((10, len(SNRS)), int); c_den = np.zeros_like(c_num)
        q_num = np.zeros((10, len(SNRS)), int); q_den = np.zeros_like(q_num)
        for r in realizations:
            source_r = [row for row in sources if int(row["realization"]) == r]
            count_r = [row for row in counts if int(row["realization"]) == r]
            for i, snr in enumerate(SNRS):
                cell = [row for row in source_r if float(row["expected_snr"]) == snr]
                c_den[r, i] = len(cell)
                c_num[r, i] = sum(row["recovered"].lower() == "true" for row in cell)
                count = next(row for row in count_r if float(row["threshold_snr"]) == snr)
                q_num[r, i] = int(count["n_negative"])
                q_den[r, i] = int(count["n_positive"])
        c = c_num.sum(0)/c_den.sum(0)
        cb = c_num[draws].sum(1)/c_den[draws].sum(1)
        c_lo, c_hi = np.percentile(cb, [16, 84], axis=0)
        q = np.minimum(q_num.sum(0)/np.maximum(q_den.sum(0), 1), 1)
        qb = np.minimum(q_num[draws].sum(1)/np.maximum(q_den[draws].sum(1), 1), 1)
        q_lo, q_hi = np.percentile(qb, [16, 84], axis=0)
        style = dict(color=color, marker="o", markersize=4, linewidth=1.6, label=config)
        ax_c.plot(SNRS, c, **style); ax_c.fill_between(SNRS, c_lo, c_hi, color=color, alpha=.1)
        ax_q.plot(SNRS, q, **style); ax_q.fill_between(SNRS, q_lo, q_hi, color=color, alpha=.1)
        ax_n.plot(SNRS, q_num.sum(0)/10, color=color, marker="o", linewidth=1.6, label=config)
        ax_cq.plot(c, q, **style)
        for snr, val, lo, hi, num, den in zip(SNRS, c, c_lo, c_hi, c_num.sum(0), c_den.sum(0)):
            metrics.append({"configuration":config, "metric":"differential_completeness",
                "snr":snr, "numerator":num, "denominator":den,
                "value":val, "lower_16":lo, "upper_84":hi})
        for snr, val, lo, hi, num, den in zip(SNRS, q, q_lo, q_hi, q_num.sum(0), q_den.sum(0)):
            metrics.append({"configuration":config, "metric":"cumulative_contamination",
                "snr":snr, "numerator":num, "denominator":den,
                "value":val, "lower_16":lo, "upper_84":hi})
    ax_c.set(title="Differential completeness", xlabel="Input apparent S/N", ylabel="C", ylim=(-.04,1.04))
    ax_q.set(title="Cumulative contamination", xlabel="Measured peak S/N threshold", ylabel="Q", ylim=(-.04,1.04))
    ax_n.set(title="Negative peaks per full-PB field", xlabel="Measured |peak| S/N threshold",
             ylabel="Mean count per realization", yscale="log")
    ax_cq.set(title="Completeness–contamination relation", xlabel="Differential completeness C",
              ylabel="Cumulative contamination Q", xlim=(-.04,1.04), ylim=(-.04,1.04))
    for ax in axes.ravel():
        ax.axvline(4, color="#666666", linestyle="--", linewidth=.8) if ax in (ax_q, ax_n) else None
        ax.grid(alpha=.22); ax.legend(title="Configuration", fontsize=8, ncol=2)
    fig.suptitle("GHOST F1–F6 native full-PB: 10 realizations, 4 pixels/minor beam, no taper")
    fig.savefig(OUT/"native_full_pb_matrix_metrics.png", dpi=200)
    fig.savefig(OUT/"native_full_pb_matrix_metrics.pdf")
    plt.close(fig)
    with (OUT/"native_full_pb_matrix_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics[0])); writer.writeheader(); writer.writerows(metrics)
    print(OUT/"native_full_pb_matrix_metrics.png")


if __name__ == "__main__": main()
