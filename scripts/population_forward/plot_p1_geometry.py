"""Plot P1a/P1b completeness, RMS, and injection-matched contamination."""

from __future__ import annotations

import csv
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

from bliss_alma.provenance import write_provenance


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "population_forward" / "p1_geometry"
FIGURES = INPUT / "figures"
TRACKS = ["P1a_geometry_only", "P1b_array_total"]
CONFIGS = ["current_ALMA_C7", "ALMA2040_F4"]
LABELS = {"current_ALMA_C7": "Current ALMA C-7 (43)", "ALMA2040_F4": "ALMA2040 F4 (255)"}
COLORS = {"current_ALMA_C7": "#d62728", "ALMA2040_F4": "#1f77b4"}


def read_rows(name):
    with (INPUT / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def completeness(rows, track, config):
    selected = [
        row for row in rows
        if row["track"] == track and row["configuration"] == config
    ]
    snrs = sorted({float(row["target_current_alma_snr"]) for row in selected})
    numerator = np.zeros((10, len(snrs)))
    denominator = np.zeros_like(numerator)
    for realization in range(10):
        for index, snr in enumerate(snrs):
            cell = [
                row for row in selected
                if int(row["realization"]) == realization
                and float(row["target_current_alma_snr"]) == snr
            ]
            denominator[realization, index] = len(cell)
            numerator[realization, index] = sum(
                row["recovered"].lower() == "true" for row in cell
            )
    value = numerator.sum(0) / denominator.sum(0)
    rng = np.random.default_rng(2040)
    draws = rng.integers(0, 10, size=(10000, 10))
    samples = numerator[draws].sum(1) / denominator[draws].sum(1)
    low, high = np.percentile(samples, [16, 84], axis=0)
    return np.asarray(snrs), value, low, high


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    sources = read_rows("source_recovery_all.csv")
    realizations = read_rows("realization_summary.csv")
    figure, axes = plt.subplots(2, 2, figsize=(11.2, 8.5), constrained_layout=True)

    for axis, track, title in zip(
        axes[0],
        TRACKS,
        ["P1a: common image RMS (geometry only)", "P1b: common per-visibility noise"],
    ):
        for config in CONFIGS:
            snr, value, low, high = completeness(sources, track, config)
            axis.plot(
                snr, value, marker="o", color=COLORS[config], label=LABELS[config]
            )
            axis.fill_between(snr, low, high, color=COLORS[config], alpha=0.15)
        axis.axvline(4.0, color="#666666", linestyle="--", linewidth=1)
        axis.set(
            title=title,
            xlabel="Fixed source flux expressed as current-ALMA fiducial S/N",
            ylabel="Differential completeness",
            ylim=(-0.04, 1.04),
        )
        axis.grid(alpha=0.2)
        axis.legend(frameon=False)

    axis = axes[1, 0]
    positions = np.arange(len(CONFIGS))
    for offset, track in zip([-0.15, 0.15], TRACKS):
        values = []
        errors = []
        for config in CONFIGS:
            cell = [
                float(row["analysis_rms_jy_beam"]) * 1e6
                for row in realizations
                if row["track"] == track and row["configuration"] == config
            ]
            values.append(np.mean(cell))
            errors.append(np.std(cell, ddof=1))
        axis.bar(
            positions + offset,
            values,
            width=0.28,
            yerr=errors,
            label=track.replace("_", " "),
            alpha=0.8,
        )
    axis.set(
        title="Dirty-image RMS",
        ylabel=r"RMS [$\mu$Jy beam$^{-1}$]",
        xticks=positions,
        xticklabels=["Current ALMA", "ALMA2040"],
    )
    axis.legend(frameon=False, fontsize=8)
    axis.grid(axis="y", alpha=0.2)

    axis = axes[1, 1]
    for track, linestyle in zip(TRACKS, ["-", "--"]):
        means = []
        stds = []
        for config in CONFIGS:
            values = [
                float(row["truth_contamination_4sigma"])
                for row in realizations
                if row["track"] == track and row["configuration"] == config
            ]
            means.append(np.mean(values))
            stds.append(np.std(values, ddof=1))
        axis.errorbar(
            positions,
            means,
            yerr=stds,
            marker="o",
            linestyle=linestyle,
            label=track.replace("_", " "),
        )
    axis.set(
        title="Injection-matched 4σ false-peak fraction",
        ylabel="False peaks / all positive peaks",
        xticks=positions,
        xticklabels=["Current ALMA", "ALMA2040"],
        ylim=(-0.04, 1.04),
    )
    axis.legend(frameon=False, fontsize=8)
    axis.grid(alpha=0.2)

    figure.suptitle(
        "P1 fixed-7.5-GHz comparison: native dirty beams, 10 noise realizations"
    )
    figure.savefig(FIGURES / "p1_geometry_and_array_comparison.png", dpi=200)
    figure.savefig(FIGURES / "p1_geometry_and_array_comparison.pdf")
    plt.close(figure)
    write_provenance(
        FIGURES / "p1_geometry_and_array_comparison.provenance.json",
        root=ROOT,
        script=Path(__file__),
        command=sys.argv,
        inputs=[INPUT / "source_recovery_all.csv", INPUT / "realization_summary.csv"],
        outputs=[
            FIGURES / "p1_geometry_and_array_comparison.png",
            FIGURES / "p1_geometry_and_array_comparison.pdf",
        ],
        random_seeds=[2040],
    )
    print(FIGURES / "p1_geometry_and_array_comparison.png")


if __name__ == "__main__":
    main()
