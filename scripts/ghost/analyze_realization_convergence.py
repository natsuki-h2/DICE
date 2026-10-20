"""Estimate how F1 metric uncertainty changes with the realization count.

The resampling unit is one complete visibility-noise realization.  This keeps
all sources and detections from a simulated observation together.
"""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np
from bliss_alma.provenance import write_script_provenance


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "ghost_f1_5min" / "f1_independent_noise_run"
OUT = INPUT / "convergence"
COUNTS = np.unique(np.r_[np.arange(2, 21), 25, 30, 40, 50, 75, 100])
N_DRAWS = 20_000
TARGET_SNR = 4.0


def load(name):
    with (INPUT / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    sources = load("source_recovery_all.csv")
    detections = load("detections_all.csv")
    realizations = sorted({int(row["realization"]) for row in sources})
    n_real = len(realizations)

    c_num = np.zeros(n_real, dtype=int)
    c_den = np.zeros(n_real, dtype=int)
    q_num = np.zeros(n_real, dtype=int)
    q_den = np.zeros(n_real, dtype=int)
    for index, realization in enumerate(realizations):
        cell = [row for row in sources
                if int(row["realization"]) == realization
                and float(row["expected_snr"]) == TARGET_SNR]
        c_num[index] = sum(row["recovered"].lower() == "true" for row in cell)
        c_den[index] = len(cell)
        positive = [row for row in detections
                    if int(row["realization"]) == realization
                    and row["detection_kind"] == "populated_diagnostic"]
        negative = [row for row in detections
                    if int(row["realization"]) == realization
                    and row["detection_kind"] == "blank_negative"]
        q_den[index] = sum(float(row["peak_snr"]) >= TARGET_SNR for row in positive)
        q_num[index] = sum(float(row["peak_snr"]) >= TARGET_SNR for row in negative)

    rng = np.random.default_rng(2040)
    rows = []
    for count in COUNTS:
        draws = rng.integers(0, n_real, size=(N_DRAWS, count))
        c = c_num[draws].sum(1) / c_den[draws].sum(1)
        q = np.minimum(q_num[draws].sum(1) / np.maximum(q_den[draws].sum(1), 1), 1.0)
        c16, c50, c84 = np.percentile(c, [16, 50, 84])
        q16, q50, q84 = np.percentile(q, [16, 50, 84])
        rows.append({
            "n_realizations": int(count),
            "completeness_p16": c16,
            "completeness_median": c50,
            "completeness_p84": c84,
            "completeness_half_width": (c84 - c16) / 2,
            "contamination_p16": q16,
            "contamination_median": q50,
            "contamination_p84": q84,
            "contamination_half_width": (q84 - q16) / 2,
        })

    OUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUT / "f1_realization_convergence.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)

    counts = np.asarray([row["n_realizations"] for row in rows])
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.4), constrained_layout=True)
    series = [
        ("completeness", "Differential completeness", "#4477AA"),
        ("contamination", "Cumulative contamination", "#CC6677"),
    ]
    for ax, (key, title, color) in zip(axes, series):
        median = np.asarray([row[f"{key}_median"] for row in rows])
        lower = np.asarray([row[f"{key}_p16"] for row in rows])
        upper = np.asarray([row[f"{key}_p84"] for row in rows])
        ax.fill_between(counts, lower, upper, color=color, alpha=0.2,
                        label="16–84% range")
        ax.plot(counts, median, "o-", color=color, markersize=3.5,
                label="resampled median")
        ax.axvline(10, color="#555555", linestyle="--", linewidth=1.2,
                   label="current N=10")
        ax.set(xlabel="Number of noise realizations", ylabel=key[0].upper(),
               title=f"{title} at S/N = 4", xscale="log", ylim=(-0.03, 1.03))
        ax.grid(alpha=0.22)
        ax.legend(fontsize=8)
    fig.suptitle("GHOST F1: uncertainty versus independent-noise sample size")
    fig.savefig(OUT / "f1_realization_convergence.png", dpi=200)
    fig.savefig(OUT / "f1_realization_convergence.pdf")
    plt.close(fig)
    write_script_provenance(OUT/"f1_realization_convergence.provenance.json",root=ROOT,script=__file__,
        inputs=[INPUT/"source_recovery_all.csv",INPUT/"detections_all.csv"],
        outputs=[csv_path,OUT/"f1_realization_convergence.png",OUT/"f1_realization_convergence.pdf"],
        parameters={"draws":N_DRAWS,"target_snr":TARGET_SNR,"uncertainty_unit":"whole realization"},random_seeds=[2040])
    print(csv_path)


if __name__ == "__main__":
    main()
