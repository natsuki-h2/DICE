"""Plot the F1/F6 non-Gaussian visibility-noise diagnostic."""

from pathlib import Path
import argparse
import csv

import matplotlib.pyplot as plt
import numpy as np
from bliss_alma.provenance import write_script_provenance


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "data/ghost_matrix/non_gaussian_visibility_pilot"
FIG = BASE / "figures"
SCENARIOS = ("gaussian", "student_t3", "sparse_outliers", "time_coherent")
LABELS = {"gaussian":"Gaussian", "student_t3":"Student-t (df=3)",
          "sparse_outliers":"0.1% ×10 outliers", "time_coherent":"10% time-coherent"}
COLORS = {"gaussian":"black", "student_t3":"#0072B2",
          "sparse_outliers":"#E69F00", "time_coherent":"#CC79A7"}
MARKERS = {"gaussian":"o", "student_t3":"s",
           "sparse_outliers":"^", "time_coherent":"D"}
LINESTYLES = {"gaussian":"-", "student_t3":"--",
              "sparse_outliers":"-.", "time_coherent":":"}


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def band(values):
    return np.percentile(np.asarray(values, float), [16, 50, 84])


def main():
    global BASE,FIG
    parser=argparse.ArgumentParser()
    parser.add_argument("--input-root",type=Path,default=BASE)
    args=parser.parse_args()
    BASE=args.input_root.resolve(); FIG=BASE/"figures"
    FIG.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 4, figsize=(15.5, 7.8), constrained_layout=True)
    summary_rows = []
    for row_index, config in enumerate(("F1", "F6")):
        stats = read(BASE/config.lower()/"image_statistics.csv")
        recovery = read(BASE/config.lower()/"source_recovery.csv")
        gaussian_rms = np.median([float(r["rms_jy_per_beam"]) for r in stats if r["scenario"] == "gaussian"])
        for scenario in SCENARIOS:
            selected = [r for r in stats if r["scenario"] == scenario]
            rms = band([float(r["rms_jy_per_beam"])/gaussian_rms for r in selected])
            kurt = band([float(r["pixel_excess_kurtosis"]) for r in selected])
            x = np.arange(len(SCENARIOS)); xi = SCENARIOS.index(scenario)
            axes[row_index, 0].errorbar(xi, rms[1], yerr=[[rms[1]-rms[0]], [rms[2]-rms[1]]],
                                        fmt=MARKERS[scenario], color=COLORS[scenario],
                                        markeredgecolor="#222222", markeredgewidth=.35, capsize=3)
            axes[row_index, 1].errorbar(xi, kurt[1], yerr=[[kurt[1]-kurt[0]], [kurt[2]-kurt[1]]],
                                        fmt=MARKERS[scenario], color=COLORS[scenario],
                                        markeredgecolor="#222222", markeredgewidth=.35, capsize=3)
            thresholds = np.arange(3.0, 7.01, .5)
            peak_med = []
            for threshold in thresholds:
                key = f"n_positive_peaks_ge_{threshold:.1f}"
                peak_med.append(band([float(r[key]) for r in selected]))
            peak_med = np.asarray(peak_med)
            style = dict(color=COLORS[scenario], marker=MARKERS[scenario],
                         linestyle=LINESTYLES[scenario], linewidth=1.7,
                         markersize=5, markeredgecolor="#222222",
                         markeredgewidth=.35)
            axes[row_index, 2].plot(thresholds, peak_med[:, 1],
                                    label=LABELS[scenario], **style)
            axes[row_index, 2].fill_between(thresholds, peak_med[:, 0], peak_med[:, 2], color=COLORS[scenario], alpha=.15)
            snrs = np.array((3., 4., 5., 6.)); comp = []
            for snr in snrs:
                picked = [float(r["completeness"]) for r in recovery
                          if r["scenario"] == scenario and float(r["input_snr"]) == snr]
                comp.append(band(picked))
            comp = np.asarray(comp)
            axes[row_index, 3].plot(snrs, comp[:, 1], **style)
            axes[row_index, 3].fill_between(snrs, comp[:, 0], comp[:, 2], color=COLORS[scenario], alpha=.15)
            summary_rows.append({"configuration":config, "scenario":scenario,
                "median_rms_ratio_to_gaussian":rms[1], "median_excess_kurtosis":kurt[1],
                "median_n_positive_peaks_ge_4sigma":peak_med[2, 1],
                "median_completeness_input_snr_4":comp[1, 1]})

        axes[row_index, 0].axhline(1, color="0.7", lw=1)
        axes[row_index, 1].axhline(0, color="0.7", lw=1)
        for col in (0, 1):
            axes[row_index, col].set_xticks(range(len(SCENARIOS)), ["G", "t3", "out", "corr"])
        axes[row_index, 0].set_ylabel(f"{config}\nRMS / Gaussian RMS")
        axes[row_index, 1].set_ylabel("Pixel excess kurtosis")
        axes[row_index, 2].set_yscale("symlog", linthresh=1)
        axes[row_index, 2].set_ylabel("Positive local peaks per image")
        axes[row_index, 2].set_xlabel("Measured peak S/N threshold")
        axes[row_index, 3].set_ylim(-.02, 1.02)
        axes[row_index, 3].set_ylabel("Completeness at 4σ detection")
        axes[row_index, 3].set_xlabel("Input S/N")
        for ax in axes[row_index]: ax.grid(alpha=.2)
    axes[0, 0].set_title("Absolute sensitivity")
    axes[0, 1].set_title("Image-pixel non-Gaussianity")
    axes[0, 2].set_title("False-positive tail (blank image)")
    axes[0, 3].set_title("Isolated-source recovery")
    axes[0, 2].legend(fontsize=8)
    n_realizations=len({int(r["realization"]) for r in read(BASE/"f1"/"image_statistics.csv")})
    fig.suptitle("Non-Gaussian visibility-noise diagnostic\n"
                 f"Native beam, natural weighting, {n_realizations} realizations; bands are realization 16–84% ranges")
    png = FIG/"non_gaussian_visibility_pilot.png"
    fig.savefig(png, dpi=200); fig.savefig(png.with_suffix(".pdf")); plt.close(fig)
    with (FIG/"non_gaussian_visibility_pilot_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0])); writer.writeheader(); writer.writerows(summary_rows)
    write_script_provenance(FIG/"non_gaussian_visibility_pilot.provenance.json",root=ROOT,script=__file__,
        inputs=[BASE],outputs=[png,png.with_suffix(".pdf"),FIG/"non_gaussian_visibility_pilot_summary.csv"],
        parameters={"scenarios":list(SCENARIOS),"scope":"diagnostic comparison with fixed marginal variance"})
    print(png)


if __name__ == "__main__": main()
