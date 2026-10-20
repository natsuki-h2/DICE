"""Compare three catalogue-level detection thresholds for the F1 run.

The saved 2-sigma diagnostic catalogues are filtered by measured peak S/N.
This isolates the seed-threshold effect without regenerating visibilities.  It
does not re-segment islands at a different flood threshold, so the result is a
fast threshold pilot rather than a replacement for a full extraction rerun.
"""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "ghost_f1_5min" / "f1_independent_noise_run"
OUT = INPUT / "thresholds"
THRESHOLDS = np.asarray([3.5, 4.0, 5.0])
TARGET_SNRS = np.arange(2.0, 7.01, 0.5)
MAX_DISTANCE_PIXELS = 4.0
N_RESAMPLES = 10_000


def load(name):
    with (INPUT / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    sources = load("source_recovery_all.csv")
    detections = load("detections_all.csv")
    realizations = sorted({int(row["realization"]) for row in sources})
    n_real = len(realizations)
    c_num = np.zeros((n_real, len(THRESHOLDS), len(TARGET_SNRS)), dtype=int)
    c_den = np.zeros_like(c_num)
    q_num = np.zeros((n_real, len(THRESHOLDS)), dtype=int)
    q_den = np.zeros_like(q_num)

    for ri, realization in enumerate(realizations):
        src = [row for row in sources if int(row["realization"]) == realization]
        positive = [row for row in detections
                    if int(row["realization"]) == realization
                    and row["detection_kind"] == "populated_diagnostic"]
        negative = [row for row in detections
                    if int(row["realization"]) == realization
                    and row["detection_kind"] == "blank_negative"]
        for ti, threshold in enumerate(THRESHOLDS):
            selected = [row for row in positive if float(row["peak_snr"]) >= threshold]
            for si, snr in enumerate(TARGET_SNRS):
                cell = [row for row in src if float(row["expected_snr"]) == snr]
                c_den[ri, ti, si] = len(cell)
                c_num[ri, ti, si] = sum(
                    any(np.hypot(float(det["x_pixel"]) - float(row["x_pixel"]),
                                 float(det["y_pixel"]) - float(row["y_pixel"]))
                        <= MAX_DISTANCE_PIXELS for det in selected)
                    for row in cell
                )
            q_den[ri, ti] = len(selected)
            q_num[ri, ti] = sum(float(row["peak_snr"]) >= threshold for row in negative)

    rng = np.random.default_rng(2040)
    draws = rng.integers(0, n_real, size=(N_RESAMPLES, n_real))
    c_point = c_num.sum(0) / c_den.sum(0)
    c_boot = c_num[draws].sum(1) / c_den[draws].sum(1)
    c_lo, c_hi = np.percentile(c_boot, [16, 84], axis=0)
    q_point = np.minimum(q_num.sum(0) / np.maximum(q_den.sum(0), 1), 1.0)
    q_boot = np.minimum(q_num[draws].sum(1) / np.maximum(q_den[draws].sum(1), 1), 1.0)
    q_lo, q_hi = np.percentile(q_boot, [16, 84], axis=0)

    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for ti, threshold in enumerate(THRESHOLDS):
        for si, snr in enumerate(TARGET_SNRS):
            rows.append({"metric": "differential_completeness",
                         "detection_threshold_sigma": threshold,
                         "input_snr": snr,
                         "numerator": int(c_num[:, ti, si].sum()),
                         "denominator": int(c_den[:, ti, si].sum()),
                         "value": c_point[ti, si],
                         "lower_16": c_lo[ti, si], "upper_84": c_hi[ti, si]})
        rows.append({"metric": "cumulative_contamination",
                     "detection_threshold_sigma": threshold,
                     "input_snr": "",
                     "numerator": int(q_num[:, ti].sum()),
                     "denominator": int(q_den[:, ti].sum()),
                     "value": q_point[ti], "lower_16": q_lo[ti], "upper_84": q_hi[ti]})
    with (OUT / "f1_threshold_comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader(); writer.writerows(rows)

    colors = ["#4477AA", "#228833", "#CC6677"]
    fig, (ax_c, ax_q) = plt.subplots(1, 2, figsize=(10.8, 4.5), constrained_layout=True)
    for ti, (threshold, color) in enumerate(zip(THRESHOLDS, colors)):
        ax_c.errorbar(TARGET_SNRS, c_point[ti],
                      yerr=[c_point[ti] - c_lo[ti], c_hi[ti] - c_point[ti]],
                      fmt="o-", color=color, capsize=2,
                      label=f"{threshold:g} sigma threshold")
    ax_q.errorbar(THRESHOLDS, q_point, yerr=[q_point-q_lo, q_hi-q_point],
                  fmt="o-", color="#111111", capsize=3)
    ax_c.set(title="Differential completeness", xlabel="Input apparent S/N",
             ylabel="C", ylim=(-0.04, 1.04))
    ax_q.set(title="Cumulative contamination", xlabel="Measured peak S/N threshold",
             ylabel="Q", ylim=(-0.01, max(0.12, float(q_hi.max()) * 1.15)))
    ax_c.legend(fontsize=9)
    for ax in (ax_c, ax_q): ax.grid(alpha=0.22)
    fig.suptitle("GHOST F1: catalogue-level threshold pilot (10 realizations)")
    fig.savefig(OUT / "f1_threshold_comparison.png", dpi=200)
    fig.savefig(OUT / "f1_threshold_comparison.pdf")
    plt.close(fig)
    print(OUT / "f1_threshold_comparison.csv")


if __name__ == "__main__":
    main()
