"""P2 continuum bandwidth-only transfer at fixed current-ALMA geometry."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bliss_alma.population_forward import modified_blackbody_flux_jy, saunders_phi


ROOT = Path(__file__).resolve().parents[1]
P1 = ROOT / "data" / "population_forward" / "p1_geometry"
PILOT = ROOT / "data" / "population_forward" / "alma2_wideband_pilot"
OUT = ROOT / "data" / "population_forward" / "p2_continuum"
CASES = {
    "7.5 GHz": 1.0,
    "15 GHz": 0.7446189415828508,
    "30 GHz": 0.5181432139492808,
}
COLORS = {"7.5 GHz": "#d62728", "15 GHz": "#ff7f0e", "30 GHz": "#1f77b4"}
LOG_EDGES = np.arange(8.0, 13.51, 0.25)


def load_rows():
    with (P1 / "source_recovery_all.csv").open() as handle:
        return [
            row for row in csv.DictReader(handle)
            if row["track"] == "P1b_array_total"
            and row["configuration"] == "current_ALMA_C7"
        ]


def realization_curve(rows, realization):
    selected = [row for row in rows if int(row["realization"]) == realization]
    x = sorted({float(row["target_current_alma_snr"]) for row in selected})
    y = []
    for value in x:
        cell = [row for row in selected if float(row["target_current_alma_snr"]) == value]
        y.append(np.mean([row["recovered"].lower() == "true" for row in cell]))
    return np.asarray(x), np.asarray(y)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    fiducial_rms = json.loads((P1 / "summary.json").read_text())[
        "fiducial_current_alma_rms_jy_beam"
    ]
    rng = np.random.default_rng(20402040)
    count = 250_000
    log_l = rng.uniform(LOG_EDGES[0], LOG_EDGES[-1], count)
    z = rng.uniform(4.5, 6.0, count)
    pb = np.exp(rng.uniform(np.log(0.2), 0.0, count))
    flux = modified_blackbody_flux_jy(
        np.ones(1201), np.linspace(4.5, 6.0, 1201), 230.0
    )
    apparent = 10**log_l * np.interp(z, np.linspace(4.5, 6.0, 1201), flux) * pb
    weights = saunders_phi(log_l)
    centres = (LOG_EDGES[:-1] + LOG_EDGES[1:]) / 2
    width = np.diff(LOG_EDGES)
    proposal = 1.0 / (LOG_EDGES[-1] - LOG_EDGES[0])
    input_phi = (
        np.histogram(log_l, LOG_EDGES, weights=weights)[0]
        / count / width / proposal
    )
    products = {}
    table_rows = []
    for label, rms_factor in CASES.items():
        curves, outputs = [], []
        for realization in range(10):
            x, y = realization_curve(rows, realization)
            snr = apparent / (fiducial_rms * rms_factor)
            probability = np.interp(snr, x, y, left=0.0, right=1.0)
            curves.append(np.interp(np.linspace(0, 8, 161) / rms_factor, x, y, left=0, right=1))
            outputs.append(
                np.histogram(log_l, LOG_EDGES, weights=weights * probability)[0]
                / count / width / proposal
            )
        outputs = np.asarray(outputs)
        products[label] = (np.asarray(curves), outputs)
        lo, med, hi = np.percentile(outputs, [16, 50, 84], axis=0)
        for i, centre in enumerate(centres):
            table_rows.append({
                "bandwidth": label,
                "log10_lir_lsun": centre,
                "input_phi": input_phi[i],
                "output_phi_median": med[i],
                "output_phi_lower16": lo[i],
                "output_phi_upper84": hi[i],
                "output_input_median": med[i] / input_phi[i],
            })
    with (OUT / "p2_continuum_lf.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table_rows[0]))
        writer.writeheader(); writer.writerows(table_rows)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    snr_axis = np.linspace(0, 8, 161)
    for label in CASES:
        curves, outputs = products[label]
        lo, med, hi = np.percentile(curves, [16, 50, 84], axis=0)
        axes[0].plot(snr_axis, med, color=COLORS[label], label=label)
        axes[0].fill_between(snr_axis, lo, hi, color=COLORS[label], alpha=.15)
        lo, med, hi = np.percentile(outputs, [16, 50, 84], axis=0)
        axes[1].plot(centres, med, color=COLORS[label], marker="o", ms=3, label=label)
        axes[1].fill_between(centres, lo, hi, color=COLORS[label], alpha=.15)
    axes[0].axvline(4, color="#666", ls="--")
    axes[0].set(xlabel="Fixed flux in 7.5-GHz fiducial S/N", ylabel="Completeness",
                title="P2 bandwidth-only completeness", ylim=(-.03, 1.03))
    axes[1].plot(centres, input_phi, color="black", lw=2, label="Input ALPINE LF")
    axes[1].set(xlabel=r"$\log_{10}(L_{\rm IR}/L_\odot)$",
                ylabel=r"$\phi$ [Mpc$^{-3}$ dex$^{-1}$]",
                title="Directly recovered LF", yscale="log", ylim=(1e-7, 8e-3))
    for ax in axes:
        ax.grid(alpha=.2); ax.legend(frameon=False)
    fig.savefig(OUT / "p2_continuum_bandwidth.png", dpi=200)
    fig.savefig(OUT / "p2_continuum_bandwidth.pdf")
    print(OUT / "p2_continuum_bandwidth.png")


if __name__ == "__main__":
    main()
