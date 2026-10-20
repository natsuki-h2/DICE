"""Plot differential F1 metrics from independent visibility noise."""

from pathlib import Path
import argparse
import csv
import matplotlib.pyplot as plt
import numpy as np
from bliss_alma.provenance import write_script_provenance

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "ghost_f1_5min" / "f1_independent_noise_run"
OUT = INPUT / "figures"
PB_LABELS = ["0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"]
COLORS = ["#4477AA", "#66CCEE", "#228833", "#CCBB44"]
MARKERS = ["o", "s", "^", "D"]
LINESTYLES = ["-", "--", "-.", ":"]
LIMITS = [(0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.01)]
TARGETS = np.arange(2.0, 7.01, 0.5)
EDGES = np.arange(2.0, 10.01, 0.5)
N_BOOTSTRAP = 10000

def load(name):
    with (INPUT / name).open(newline="") as handle:
        return list(csv.DictReader(handle))

def main():
    global INPUT,OUT
    parser=argparse.ArgumentParser()
    parser.add_argument("--input-root",type=Path,default=INPUT)
    args=parser.parse_args()
    INPUT=args.input_root.resolve(); OUT=INPUT/"figures"
    sources, detections = load("source_recovery_all.csv"), load("detections_all.csv")
    OUT.mkdir(parents=True, exist_ok=True)
    fig, (ax_c, ax_q, ax_cq) = plt.subplots(1, 3, figsize=(17.2, 4.8), constrained_layout=True)
    metrics = []
    for index, (label, color, marker, linestyle) in enumerate(
        zip(PB_LABELS, COLORS, MARKERS, LINESTYLES)
    ):
        cvals = []
        for snr in TARGETS:
            cell = [r for r in sources if int(r["pb_bin"]) == index and float(r["expected_snr"]) == snr]
            nrec = sum(r["recovered"].lower() == "true" for r in cell)
            value = nrec / len(cell); cvals.append(value)
            metrics.append({"metric":"completeness", "pb_bin":label, "snr":snr,
                            "numerator":nrec, "denominator":len(cell), "value":value})
        style = dict(color=color, marker=marker, linestyle=linestyle,
                     alpha=.82, linewidth=1.5, markersize=4.5,
                     markeredgecolor="#222222", markeredgewidth=.35,
                     label=label)
        ax_c.plot(TARGETS, cvals, **style)
        lo_pb, hi_pb = LIMITS[index]
        positive = [r for r in detections if r["detection_kind"] == "populated_diagnostic"
                    and lo_pb <= float(r["pb"]) < hi_pb]
        negative = [r for r in detections if r["detection_kind"] == "blank_negative"
                    and lo_pb <= float(r["pb"]) < hi_pb]
        qvals = []
        for threshold in EDGES[:-1]:
            npos = sum(float(r["peak_snr"]) >= threshold for r in positive)
            nneg = sum(float(r["peak_snr"]) >= threshold for r in negative)
            raw = nneg / npos if npos else np.nan
            value = min(raw, 1.0) if np.isfinite(raw) else np.nan
            qvals.append(value)
            metrics.append({"metric":"cumulative_contamination", "pb_bin":label, "snr":threshold,
                            "numerator":nneg, "denominator":npos, "value":value})
        ax_q.plot(EDGES[:-1], qvals, **style)
        qtraj = np.asarray(qvals[:len(TARGETS)])
        usable_pb = np.isfinite(qtraj)
        ax_cq.plot(np.asarray(cvals)[usable_pb], qtraj[usable_pb], **style)

    # Aggregate all PB bins, resampling whole noise realizations so points from
    # one observation are never treated as independent draws.
    rng = np.random.default_rng(2040)
    realizations = np.array(sorted({int(r["realization"]) for r in sources}))
    n_realizations = len(realizations)
    draws = rng.choice(realizations, size=(N_BOOTSTRAP, len(realizations)), replace=True)
    c_num = np.zeros((n_realizations, len(TARGETS)), dtype=int)
    c_den = np.zeros_like(c_num)
    for row_index, realization in enumerate(realizations):
        for column, snr in enumerate(TARGETS):
            cell = [r for r in sources if int(r["realization"]) == realization and float(r["expected_snr"]) == snr]
            c_num[row_index, column] = sum(r["recovered"].lower() == "true" for r in cell)
            c_den[row_index, column] = len(cell)
    c_point = c_num.sum(axis=0) / c_den.sum(axis=0)
    c_boot = c_num[draws].sum(axis=1) / c_den[draws].sum(axis=1)
    c_lo, c_hi = np.percentile(c_boot, [16, 84], axis=0)
    ax_c.errorbar(TARGETS, c_point, yerr=[c_point-c_lo, c_hi-c_point], fmt="o-",
                  color="#111111", linewidth=2.2, markersize=5.5, capsize=3,
                  label="all PB (16–84% range)", zorder=10)

    q_num = np.zeros((n_realizations, len(EDGES)-1), dtype=int)
    q_den = np.zeros_like(q_num)
    for row_index, realization in enumerate(realizations):
        positive = [r for r in detections if r["detection_kind"] == "populated_diagnostic"
                    and int(r["realization"]) == realization]
        negative = [r for r in detections if r["detection_kind"] == "blank_negative"
                    and int(r["realization"]) == realization]
        for column, threshold in enumerate(EDGES[:-1]):
            q_den[row_index, column] = sum(float(r["peak_snr"]) >= threshold for r in positive)
            q_num[row_index, column] = sum(float(r["peak_snr"]) >= threshold for r in negative)
    q_point = np.minimum(q_num.sum(axis=0) / np.maximum(q_den.sum(axis=0), 1), 1.0)
    boot_num, boot_den = q_num[draws].sum(axis=1), q_den[draws].sum(axis=1)
    q_boot = np.minimum(boot_num / np.maximum(boot_den, 1), 1.0)
    q_lo, q_hi = np.percentile(q_boot, [16, 84], axis=0)
    thresholds = EDGES[:-1]
    usable = q_den.sum(axis=0) > 0
    ax_q.errorbar(thresholds[usable], q_point[usable],
                  yerr=[q_point[usable]-q_lo[usable], q_hi[usable]-q_point[usable]],
                  fmt="o-", color="#111111", linewidth=2.2, markersize=5.5,
                  capsize=3, label="all PB (16–84% range)", zorder=10)

    # Hayatsu (2025)-style completeness-contamination relation. S/N is the shared
    # parameter: differential completeness is evaluated in the injected-S/N
    # bin, while contamination is cumulative above the measured-S/N threshold.
    qtraj_point = q_point[:len(TARGETS)]
    qtraj_lo = q_lo[:len(TARGETS)]
    qtraj_hi = q_hi[:len(TARGETS)]
    ax_cq.errorbar(c_point, qtraj_point,
                   xerr=[c_point-c_lo, c_hi-c_point],
                   yerr=[qtraj_point-qtraj_lo, qtraj_hi-qtraj_point],
                   fmt="o-", color="#111111", linewidth=2.2, markersize=5.5,
                   capsize=3, label="all PB (16–84% range)", zorder=10)
    labels = {
        0: ("2σ", (6, 6)),
        2: ("3σ", (6, 6)),
        4: ("4σ", (6, 6)),
        6: ("5σ", (6, 6)),
        8: ("6–7σ", (-8, 7)),
    }
    for index, (cvalue, qvalue) in enumerate(zip(c_point, qtraj_point)):
        if index in labels:
            text, offset = labels[index]
            ax_cq.annotate(text, (cvalue, qvalue), xytext=offset,
                           textcoords="offset points", fontsize=8)
    for x, value, lo, hi, num, den in zip(TARGETS, c_point, c_lo, c_hi, c_num.sum(0), c_den.sum(0)):
        metrics.append({"metric":"completeness", "pb_bin":"all PB", "snr":x,
                        "numerator":num, "denominator":den, "value":value,
                        "lower_16":lo, "upper_84":hi})
    for x, value, lo, hi, num, den in zip(thresholds, q_point, q_lo, q_hi, q_num.sum(0), q_den.sum(0)):
        metrics.append({"metric":"cumulative_contamination", "pb_bin":"all PB", "snr":x,
                        "numerator":num, "denominator":den, "value":value,
                        "lower_16":lo, "upper_84":hi})
    ax_c.set(title="Differential completeness", xlabel="Input apparent S/N", ylabel="C", ylim=(-.04,1.04))
    ax_c.text(.04,.94,r"$C_i=N_{\rm recovered,i}/N_{\rm injected,i}$",transform=ax_c.transAxes,va="top",fontsize=11,
              bbox={"facecolor":"white","edgecolor":"#BBBBBB","alpha":.9})
    ax_c.legend(title="PB response",ncol=2,fontsize=8)
    ax_q.set(title="Cumulative contamination",xlabel="Measured peak S/N threshold",ylabel="Q",ylim=(-.04,1.04))
    ax_q.text(.98,.94,r"$Q_i=\min(N_{\rm negative,i}/N_{\rm positive,i},1)$",transform=ax_q.transAxes,ha="right",va="top",fontsize=11,
              bbox={"facecolor":"white","edgecolor":"#BBBBBB","alpha":.9})
    ax_q.axvline(4,color="#555555",linestyle="--",linewidth=1.2)
    ax_q.text(4.05,.82,"adopted 4σ threshold",transform=ax_q.get_xaxis_transform(),fontsize=9)
    ax_cq.set(title="Completeness–contamination relation",
              xlabel="Differential completeness C",
              ylabel="Cumulative contamination Q",
              xlim=(-.04, 1.04), ylim=(-.04, 1.04))
    ax_cq.legend(title="PB response", ncol=2, fontsize=8)
    ax_cq.text(.04, .05, "S/N is the curve parameter", transform=ax_cq.transAxes, fontsize=9)
    for ax in (ax_c,ax_q,ax_cq): ax.grid(alpha=.22)
    fig.suptitle(f"GHOST F1: {n_realizations} independent visibility-noise realizations",fontsize=15)
    fig.savefig(OUT/"f1_independent_differential_metrics.png",dpi=200)
    fig.savefig(OUT/"f1_independent_differential_metrics.pdf")
    plt.close(fig)
    with (OUT/"f1_independent_differential_metrics.csv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=["metric","pb_bin","snr","numerator","denominator","value","lower_16","upper_84"])
        writer.writeheader(); writer.writerows(metrics)
    write_script_provenance(OUT/"f1_independent_differential_metrics.provenance.json",
        root=ROOT,script=__file__,inputs=[INPUT/"source_recovery_all.csv",INPUT/"detections_all.csv"],
        outputs=[OUT/"f1_independent_differential_metrics.png",OUT/"f1_independent_differential_metrics.pdf",OUT/"f1_independent_differential_metrics.csv"],
        parameters={"uncertainty_unit":"whole visibility-noise realization"})
    print(OUT/"f1_independent_differential_metrics.png")

if __name__ == "__main__": main()
