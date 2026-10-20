"""Plot the exact target-S/N by PB recovery grid from the GHOST F1 run."""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np
from bliss_alma.provenance import write_script_provenance


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "ghost_f1_5min" / "f1_snr_pb_grid"
Q_INPUT = ROOT / "data" / "ghost_f1_5min" / "f1_snr_pb_grid_qfloor2"
OUTPUT = INPUT / "figures"
PB_LABELS = ["0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"]
COLORS = ["#4477AA", "#66CCEE", "#228833", "#CCBB44"]
TARGETS = np.arange(2.0, 7.01, 0.5)
DETECTION_EDGES = np.arange(2.0, 10.01, 0.5)


def read_csv(name, directory=INPUT):
    with (directory / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def recovered(row):
    return row["recovered"].lower() == "true"


def main():
    sources = read_csv("source_recovery_all.csv")
    detections = read_csv("detections_all.csv", Q_INPUT)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    metric_rows = []

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.2), constrained_layout=True)
    ax_diff, ax_cum, ax_qdiff, ax_qcum = axes.ravel()

    for pb_index, (label, color) in enumerate(zip(PB_LABELS, COLORS)):
        differential = []
        cumulative = []
        for snr in TARGETS:
            cell = [r for r in sources if int(r["pb_bin"]) == pb_index and float(r["expected_snr"]) == snr]
            above = [r for r in sources if int(r["pb_bin"]) == pb_index and float(r["expected_snr"]) >= snr]
            differential.append(sum(map(recovered, cell)) / len(cell))
            cumulative.append(sum(map(recovered, above)) / len(above))
            metric_rows.extend([
                {"metric": "differential_completeness", "pb_bin": label, "snr": snr,
                 "numerator": sum(map(recovered, cell)), "denominator": len(cell), "value": differential[-1]},
                {"metric": "cumulative_completeness", "pb_bin": label, "snr": snr,
                 "numerator": sum(map(recovered, above)), "denominator": len(above), "value": cumulative[-1]},
            ])
        ax_diff.plot(TARGETS, differential, "o-", color=color, label=label)
        ax_cum.plot(TARGETS, cumulative, "o-", color=color, label=label)

        populated = [r for r in detections if r["detection_kind"] == "populated" and
                     [0.2 <= float(r["pb"]) < 0.4, 0.4 <= float(r["pb"]) < 0.6,
                      0.6 <= float(r["pb"]) < 0.8, 0.8 <= float(r["pb"]) < 1.01][pb_index]]
        blank = [r for r in detections if r["detection_kind"] == "blank_negative" and
                 [0.2 <= float(r["pb"]) < 0.4, 0.4 <= float(r["pb"]) < 0.6,
                  0.6 <= float(r["pb"]) < 0.8, 0.8 <= float(r["pb"]) < 1.01][pb_index]]
        centers = 0.5 * (DETECTION_EDGES[:-1] + DETECTION_EDGES[1:])
        qdiff, qcum = [], []
        for lo, hi in zip(DETECTION_EDGES[:-1], DETECTION_EDGES[1:]):
            den = [r for r in populated if lo <= float(r["peak_snr"]) < hi]
            blank_num = [r for r in blank if lo <= float(r["peak_snr"]) < hi]
            qdiff.append(min(len(blank_num) / len(den), 1.0) if den else np.nan)
        for threshold in DETECTION_EDGES[:-1]:
            den = [r for r in populated if float(r["peak_snr"]) >= threshold]
            blank_num = [r for r in blank if float(r["peak_snr"]) >= threshold]
            qcum.append(min(len(blank_num) / len(den), 1.0) if den else np.nan)
        ax_qdiff.plot(centers, qdiff, "o-", color=color)
        ax_qcum.plot(DETECTION_EDGES[:-1], qcum, "o-", color=color)

        for x, value in zip(centers, qdiff):
            if np.isfinite(value):
                metric_rows.append({"metric": "differential_contamination", "pb_bin": label,
                                    "snr": x, "numerator": "", "denominator": "", "value": value})
        for x, value in zip(DETECTION_EDGES[:-1], qcum):
            if np.isfinite(value):
                metric_rows.append({"metric": "cumulative_contamination", "pb_bin": label,
                                    "snr": x, "numerator": "", "denominator": "", "value": value})

    for ax in (ax_diff, ax_cum, ax_qdiff, ax_qcum):
        ax.grid(alpha=0.22)
        ax.set_xlabel("S/N")
    for ax in (ax_diff, ax_cum):
        ax.set_ylim(-0.04, 1.04)
    for ax in (ax_qdiff, ax_qcum):
        ax.axvline(4.0, color="#555555", linestyle="--", linewidth=1.2)
        ax.text(4.04, 0.96, "adopted 4σ", transform=ax.get_xaxis_transform(), va="top", fontsize=8)
        ax.set_ylim(bottom=-0.04)
    ax_diff.set(title="Differential completeness (10 sources/cell)", ylabel="Completeness")
    ax_cum.set(title="Cumulative completeness", ylabel="Completeness")
    ax_qdiff.set(title="Differential contamination", ylabel="Q")
    ax_qcum.set(title="Cumulative contamination", ylabel="Q")
    ax_diff.legend(title="PB response", ncol=2, fontsize=9)
    fig.suptitle("GHOST F1: production recovery (4σ/3σ) and negative-clump Q diagnostic (floor 2σ)")
    fig.savefig(OUTPUT / "f1_snr_pb_metrics.png", dpi=180)
    fig.savefig(OUTPUT / "f1_snr_pb_metrics.pdf")
    plt.close(fig)

    with (OUTPUT / "f1_snr_pb_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "pb_bin", "snr", "numerator", "denominator", "value"])
        writer.writeheader()
        writer.writerows(metric_rows)
    write_script_provenance(OUTPUT/"f1_snr_pb_metrics.provenance.json",root=ROOT,script=__file__,
        inputs=[INPUT/"source_recovery_all.csv",Q_INPUT/"detections_all.csv"],
        outputs=[OUTPUT/"f1_snr_pb_metrics.png",OUTPUT/"f1_snr_pb_metrics.pdf",OUTPUT/"f1_snr_pb_metrics.csv"],
        parameters={"pb_bins":PB_LABELS,"target_snrs":list(map(float,TARGETS))})
    print(OUTPUT / "f1_snr_pb_metrics.png")


if __name__ == "__main__":
    main()
