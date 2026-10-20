"""Summarize the idealized GHOST dirty-image baseline across the matrix."""

from pathlib import Path
import csv
import json
import sys

import matplotlib.pyplot as plt
import numpy as np

from bliss_alma.provenance import write_provenance


ROOT = Path(__file__).resolve().parents[2]
MATRIX = ROOT / "data" / "ghost_matrix"
OUT = MATRIX / "summary"
RUNS = {
    ("F1", "5m"): ROOT / "data" / "ghost_f1_5min" / "pipeline_log.json",
    ("F2", "5m"): MATRIX / "f2_5min" / "pipeline_log.json",
    ("F3", "5m"): MATRIX / "f3_5min" / "pipeline_log.json",
    ("F4", "5m"): MATRIX / "f4_5min" / "pipeline_log.json",
    ("F5", "5m"): MATRIX / "f5_5min" / "pipeline_log.json",
    ("F6", "5m"): MATRIX / "f6_5min" / "pipeline_log.json",
    ("F1", "1h"): MATRIX / "f1_1h" / "pipeline_log.json",
    ("F1", "5h"): MATRIX / "f1_5h" / "pipeline_log.json",
    ("F6", "1h"): MATRIX / "f6_1h" / "pipeline_log.json",
    ("F6", "5h"): MATRIX / "f6_5h" / "pipeline_log.json",
}
TIME_MINUTES = {"5m": 5.0, "1h": 60.0, "5h": 300.0}


def main():
    rows = []
    for (configuration, integration), path in RUNS.items():
        payload = json.loads(path.read_text())
        result = payload["results"][0]
        geometry = payload["geometry"]
        noise = payload["noise_validation"]["results"]["all"]
        rows.append({
            "configuration": configuration,
            "integration_time": integration,
            "integration_minutes": TIME_MINUTES[integration],
            "rms_jy_per_beam": result["rms"],
            "beam_major_arcsec": result["beam_major"],
            "beam_minor_arcsec": result["beam_minor"],
            "beam_pa_deg": result["beam_pa"],
            "max_baseline_klambda": geometry["max_baseline_uv_klam"],
            "n_good_vis": noise["N_good_vis"],
            "noise_validation_ratio": noise["ratio_meas_over_naive"],
        })
    rows.sort(key=lambda row: (TIME_MINUTES[row["integration_time"]], row["configuration"]))
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "ghost_matrix_baseline_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader(); writer.writerows(rows)

    five = [row for row in rows if row["integration_time"] == "5m"]
    labels = [row["configuration"] for row in five]
    x = np.arange(len(labels))
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.6), constrained_layout=True)
    axes[0].plot(x, [row["beam_major_arcsec"] for row in five], "o-", label="major")
    axes[0].plot(x, [row["beam_minor_arcsec"] for row in five], "o-", label="minor")
    axes[0].set(xticks=x, xticklabels=labels, yscale="log", ylabel="Beam FWHM (arcsec)",
                title="Synthesized beam at 5 min")
    axes[0].legend()
    axes[1].plot(x, np.asarray([row["rms_jy_per_beam"] for row in five]) * 1e6, "o-")
    axes[1].set(xticks=x, xticklabels=labels, ylabel="RMS (microJy/beam)",
                title="Dirty-image RMS at 5 min")

    for configuration, color in (("F1", "#4477AA"), ("F6", "#CC6677")):
        series = sorted((row for row in rows if row["configuration"] == configuration),
                        key=lambda row: row["integration_minutes"])
        time = np.asarray([row["integration_minutes"] for row in series])
        rms = np.asarray([row["rms_jy_per_beam"] for row in series]) * 1e6
        axes[2].plot(time, rms, "o-", color=color, label=configuration)
        axes[2].plot(time, rms[0] * np.sqrt(time[0] / time), "--", color=color, alpha=0.55)
    axes[2].set(xscale="log", yscale="log", xlabel="Integration time (min)",
                ylabel="RMS (microJy/beam)", title="Time scaling (dashed: t^-1/2)")
    axes[2].legend()
    for ax in axes: ax.grid(alpha=0.22)
    fig.suptitle("GHOST idealized thermal-noise baseline")
    fig.savefig(OUT / "ghost_matrix_baseline_summary.png", dpi=200)
    fig.savefig(OUT / "ghost_matrix_baseline_summary.pdf")
    plt.close(fig)
    write_provenance(
        OUT / "ghost_matrix_baseline_summary.provenance.json",
        root=ROOT,
        script=Path(__file__),
        command=sys.argv,
        inputs=RUNS.values(),
        outputs=[
            OUT / "ghost_matrix_baseline_summary.csv",
            OUT / "ghost_matrix_baseline_summary.png",
            OUT / "ghost_matrix_baseline_summary.pdf",
        ],
        parameters={"integration_minutes": TIME_MINUTES},
    )
    print(OUT / "ghost_matrix_baseline_summary.csv")


if __name__ == "__main__":
    main()
