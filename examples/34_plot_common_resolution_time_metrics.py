"""Plot common-resolution time dependence for GHOST F1 and F6."""

from pathlib import Path
import csv
import json

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "ghost_matrix" / "common_resolution_runs"
OUT = INPUT / "figures"
SCENARIOS = [(configuration, time) for configuration in ("F1", "F6")
             for time in ("5m", "1h", "5h")]
COLORS = {"F1": "#4477AA", "F6": "#CC6677"}
LINESTYLES = {"5m": ":", "1h": "--", "5h": "-"}
TARGETS = np.arange(2.0, 7.01, 0.5)
THRESHOLDS = np.arange(2.0, 7.01, 0.5)
TIME_MINUTES = {"5m": 5.0, "1h": 60.0, "5h": 300.0}
RUN_DIRECTORY_TIME = {"5m": "5min", "1h": "1h", "5h": "5h"}
N_RESAMPLES = 10_000


def load(configuration, time, name):
    run_time = RUN_DIRECTORY_TIME[time]
    with (INPUT / f"{configuration.lower()}_{run_time}" / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(2040)
    metrics = []
    fig, axes = plt.subplots(2, 2, figsize=(12.6, 9.0), constrained_layout=True)
    ax_c, ax_q, ax_cq, ax_rms = axes.ravel()
    for configuration, time in SCENARIOS:
        sources = load(configuration, time, "source_recovery_all.csv")
        detections = load(configuration, time, "detections_all.csv")
        realizations = sorted({int(row["realization"]) for row in sources})
        draws = rng.integers(0, len(realizations), size=(N_RESAMPLES, len(realizations)))
        c_num = np.zeros((len(realizations), len(TARGETS)), dtype=int)
        c_den = np.zeros_like(c_num)
        q_num = np.zeros((len(realizations), len(THRESHOLDS)), dtype=int)
        q_den = np.zeros_like(q_num)
        for ri, realization in enumerate(realizations):
            for si, snr in enumerate(TARGETS):
                cell = [row for row in sources if int(row["realization"]) == realization
                        and float(row["expected_snr"]) == snr]
                c_num[ri, si] = sum(row["recovered"].lower() == "true" for row in cell)
                c_den[ri, si] = len(cell)
            positive = [row for row in detections if int(row["realization"]) == realization
                        and row["detection_kind"] == "populated_diagnostic"]
            negative = [row for row in detections if int(row["realization"]) == realization
                        and row["detection_kind"] == "blank_negative"]
            for ti, threshold in enumerate(THRESHOLDS):
                q_den[ri, ti] = sum(float(row["peak_snr"]) >= threshold for row in positive)
                q_num[ri, ti] = sum(float(row["peak_snr"]) >= threshold for row in negative)
        c = c_num.sum(0) / c_den.sum(0)
        cb = c_num[draws].sum(1) / c_den[draws].sum(1)
        c_lo, c_hi = np.percentile(cb, [16, 84], axis=0)
        q = np.minimum(q_num.sum(0) / np.maximum(q_den.sum(0), 1), 1.0)
        qb = np.minimum(q_num[draws].sum(1) / np.maximum(q_den[draws].sum(1), 1), 1.0)
        q_lo, q_hi = np.percentile(qb, [16, 84], axis=0)
        label = f"{configuration} {time}"
        style = dict(color=COLORS[configuration], linestyle=LINESTYLES[time],
                     marker="o", markersize=3.5, linewidth=1.7, label=label)
        ax_c.plot(TARGETS, c, **style); ax_c.fill_between(TARGETS, c_lo, c_hi, color=COLORS[configuration], alpha=.07)
        ax_q.plot(THRESHOLDS, q, **style); ax_q.fill_between(THRESHOLDS, q_lo, q_hi, color=COLORS[configuration], alpha=.07)
        ax_cq.plot(c, q, **style)
        for snr, value, lo, hi in zip(TARGETS, c, c_lo, c_hi):
            metrics.append({"configuration": configuration, "integration_time": time,
                            "metric": "differential_completeness", "snr": snr,
                            "value": value, "lower_16": lo, "upper_84": hi})
        for snr, value, lo, hi in zip(THRESHOLDS, q, q_lo, q_hi):
            metrics.append({"configuration": configuration, "integration_time": time,
                            "metric": "cumulative_contamination", "snr": snr,
                            "value": value, "lower_16": lo, "upper_84": hi})

    for configuration in ("F1", "F6"):
        rms = []
        for time in ("5m", "1h", "5h"):
            run_time = RUN_DIRECTORY_TIME[time]
            summary = json.loads((INPUT / f"{configuration.lower()}_{run_time}" / "summary.json").read_text())
            rms.append(summary["mean_blank_rms_jy_per_beam"] * 1e6)
        minutes = np.asarray([5.0, 60.0, 300.0])
        rms = np.asarray(rms)
        ax_rms.plot(minutes, rms, "o-", color=COLORS[configuration], label=configuration)
        ax_rms.plot(minutes, rms[0] * np.sqrt(minutes[0] / minutes), "--",
                    color=COLORS[configuration], alpha=.45)

    ax_c.set(title="Differential completeness", xlabel="Input apparent S/N", ylabel="C", ylim=(-.04, 1.04))
    ax_q.set(title="Cumulative contamination", xlabel="Measured peak S/N threshold", ylabel="Q", ylim=(-.04, 1.04))
    ax_q.axvline(4, color="#555555", linestyle="--", linewidth=1)
    ax_cq.set(title="Completeness–contamination relation", xlabel="C", ylabel="Q",
              xlim=(-.04, 1.04), ylim=(-.04, 1.04))
    ax_cq.text(.04, .05, "S/N is the curve parameter", transform=ax_cq.transAxes, fontsize=9)
    ax_rms.set(title="Common-resolution RMS", xlabel="Integration time (min)",
               ylabel="RMS (microJy/beam)", xscale="log", yscale="log")
    for ax in (ax_c, ax_q, ax_cq): ax.legend(fontsize=8, ncol=2)
    ax_rms.legend()
    for ax in axes.ravel(): ax.grid(alpha=.22)
    fig.suptitle("GHOST F1/F6 time grid: 10 realizations at common 1-arcsec resolution")
    fig.savefig(OUT / "common_resolution_time_metrics.png", dpi=200)
    fig.savefig(OUT / "common_resolution_time_metrics.pdf")
    plt.close(fig)
    with (OUT / "common_resolution_time_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics[0]))
        writer.writeheader(); writer.writerows(metrics)
    print(OUT / "common_resolution_time_metrics.png")


if __name__ == "__main__":
    main()
