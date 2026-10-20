"""Make the email-ready differential completeness/contamination figure."""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "ghost_f1_5min"
RECOVERY = BASE / "f1_snr_pb_grid"
DIAGNOSTIC = BASE / "f1_snr_pb_grid_qfloor2"
OUT = RECOVERY / "figures"
PB_LABELS = ["0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"]
COLORS = ["#4477AA", "#66CCEE", "#228833", "#CCBB44"]
TARGETS = np.arange(2.0, 7.01, 0.5)
EDGES = np.arange(2.0, 10.01, 0.5)


def load(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def in_pb(row, index):
    value = float(row["pb"])
    limits = [(0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.01)]
    return limits[index][0] <= value < limits[index][1]


def main():
    sources = load(RECOVERY / "source_recovery_all.csv")
    detections = load(DIAGNOSTIC / "detections_all.csv")
    fig, (ax_c, ax_q) = plt.subplots(1, 2, figsize=(11.6, 4.8), constrained_layout=True)

    for index, (label, color) in enumerate(zip(PB_LABELS, COLORS)):
        completeness = []
        for snr in TARGETS:
            cell = [r for r in sources if int(r["pb_bin"]) == index and float(r["expected_snr"]) == snr]
            completeness.append(sum(r["recovered"].lower() == "true" for r in cell) / len(cell))
        ax_c.plot(TARGETS, completeness, "o-", color=color, label=label)

    # There is only one independent GHOST blank realization. Pool the PB bins
    # rather than drawing four apparently independent negative-clump curves.
    positive = [r for r in detections if r["detection_kind"] == "populated"]
    negative_one_blank = [r for r in detections if r["detection_kind"] == "blank_negative"
                          and int(r["realization"]) == 0]
    q = []
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        n_pos = sum(lo <= float(r["peak_snr"]) < hi for r in positive)
        n_neg = sum(lo <= float(r["peak_snr"]) < hi for r in negative_one_blank) * 10
        q.append(min(n_neg / n_pos, 1.0) if n_pos else np.nan)
    ax_q.plot(0.5 * (EDGES[:-1] + EDGES[1:]), q, "o-", color="#7A5195",
              label="all PB pooled")

    ax_c.set(title="Differential completeness", xlabel="Input apparent S/N", ylabel="C")
    ax_c.text(0.04, 0.94, r"$C_i = N_{\rm recovered,i}/N_{\rm injected,i}$",
              transform=ax_c.transAxes, va="top", fontsize=11,
              bbox={"facecolor": "white", "edgecolor": "#BBBBBB", "alpha": 0.9})
    ax_c.set_ylim(-0.04, 1.04)
    ax_c.legend(title="PB response", ncol=2, fontsize=9)

    ax_q.set(title="Differential contamination", xlabel="Measured peak S/N", ylabel="Q")
    ax_q.text(0.04, 0.94, r"$Q_i = N_{\rm negative,i}/N_{\rm positive,i}$",
              transform=ax_q.transAxes, va="top", fontsize=11,
              bbox={"facecolor": "white", "edgecolor": "#BBBBBB", "alpha": 0.9})
    ax_q.axvline(4.0, color="#555555", linestyle="--", linewidth=1.2)
    ax_q.text(4.05, 0.82, "adopted 4σ threshold", transform=ax_q.get_xaxis_transform(), fontsize=9)
    ax_q.set_ylim(bottom=-0.04)
    ax_q.text(0.98, 0.06, "one blank-noise realization reused ×10\nPB bins pooled",
              transform=ax_q.transAxes, ha="right", va="bottom", fontsize=9,
              bbox={"facecolor": "white", "edgecolor": "#BBBBBB", "alpha": 0.9})
    ax_q.legend(fontsize=9)

    for ax in (ax_c, ax_q):
        ax.grid(alpha=0.22)
    fig.suptitle("GHOST F1 continuum injection–recovery", fontsize=15)
    fig.savefig(OUT / "f1_differential_completeness_contamination.png", dpi=200)
    fig.savefig(OUT / "f1_differential_completeness_contamination.pdf")
    plt.close(fig)
    print(OUT / "f1_differential_completeness_contamination.png")


if __name__ == "__main__":
    main()
