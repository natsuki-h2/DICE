"""P2 [CII] LF: common-volume sensitivity and full-bandwidth yield."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

from bliss_alma.population_forward import (
    double_power_law_phi,
    gaussian_line_peak_jy,
    line_luminosity_to_flux_jykms,
    survey_volume_mpc3,
)
from bliss_alma.provenance import write_provenance


ROOT = Path(__file__).resolve().parents[2]
P1 = ROOT / "data" / "population_forward" / "p1_geometry"
OUT = ROOT / "data" / "population_forward" / "p2_cii"
REST_GHZ = 1900.5369
CHANNEL_MHZ = 15.625
CENTRE_GHZ = 260.0
DV = 299792.458 * CHANNEL_MHZ / 1000.0 / CENTRE_GHZ
CURRENT_CONT_RMS = json.loads((P1 / "summary.json").read_text())[
    "fiducial_current_alma_rms_jy_beam"
]
CHANNEL_RMS = CURRENT_CONT_RMS * np.sqrt(7500.0 / CHANNEL_MHZ)
AREA = 10.0
CASES = {
    "7.5 GHz": (256.25, 263.75),
    "15 GHz": (252.5, 267.5),
    "30 GHz": (245.0, 275.0),
}
COLORS = {"7.5 GHz": "#d62728", "15 GHz": "#ff7f0e", "30 GHz": "#1f77b4"}
LOG_EDGES = np.arange(6.0, 10.51, 0.2)


def frequency_to_z(bounds):
    low, high = bounds
    return REST_GHZ / high - 1.0, REST_GHZ / low - 1.0


def matched_snr(luminosity, redshift, fwhm, pb):
    integrated = line_luminosity_to_flux_jykms(luminosity, redshift)
    peak = gaussian_line_peak_jy(integrated, fwhm) * pb
    samples = 0.75269185 * fwhm / DV
    return peak / CHANNEL_RMS * np.sqrt(samples)


def completeness_curve(fwhm):
    log_l = np.linspace(6, 10.5, 500)
    pb = np.exp(np.linspace(np.log(.2), 0, 400))
    snr = matched_snr(10**log_l[:, None], 6.31, fwhm, pb[None, :])
    # Gaussian matched-filter statistic with a 4-sigma threshold.
    from scipy.special import ndtr
    return log_l, np.mean(ndtr(snr - 4.0), axis=1)


def run_case(label, bounds, realization, master):
    log_l, fwhm, pb, uniform = master
    zmin, zmax = frequency_to_z(bounds)
    rng = np.random.default_rng(20600000 + realization + 100 * list(CASES).index(label))
    z = rng.uniform(zmin, zmax, len(log_l))
    intrinsic_snr = matched_snr(10**log_l, z, fwhm, pb)
    measured_snr = intrinsic_snr + rng.normal(size=len(log_l))
    detected = measured_snr >= 4.0
    output_log_l = log_l[detected] + np.log10(
        np.maximum(measured_snr[detected] / intrinsic_snr[detected], 1e-6)
    )
    phi = double_power_law_phi(log_l)
    volume = survey_volume_mpc3(AREA, zmin, zmax)
    proposal = 1.0 / (LOG_EDGES[-1] - LOG_EDGES[0])
    width = np.diff(LOG_EDGES)
    input_density = (
        np.histogram(log_l, LOG_EDGES, weights=phi)[0]
        / len(log_l) / width / proposal
    )
    output_density = (
        np.histogram(output_log_l, LOG_EDGES, weights=phi[detected])[0]
        / len(log_l) / width / proposal
    )
    # Importance-sampled integral over logL; fwhm, PB, and z are drawn from
    # their target distributions.
    recovered_count = (
        volume
        * (LOG_EDGES[-1] - LOG_EDGES[0])
        * np.mean(phi * detected)
    )
    input_count = volume * (LOG_EDGES[-1] - LOG_EDGES[0]) * np.mean(phi)
    return input_density, output_density, input_count, recovered_count, volume


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20500000)
    n = 250_000
    log_l = rng.uniform(LOG_EDGES[0], LOG_EDGES[-1], n)
    fwhm = 10 ** rng.normal(np.log10(200.0), 0.25, n)
    pb = np.exp(rng.uniform(np.log(.2), 0, n))
    master = (log_l, fwhm, pb, rng.random(n))
    products = {}
    rows = []
    for label, bounds in CASES.items():
        outputs, recovered, inputs, volumes = [], [], [], []
        input_density = None
        for realization in range(10):
            input_density, output, n_input, n_recovered, volume = run_case(
                label, bounds, realization, master
            )
            outputs.append(output); recovered.append(n_recovered)
            inputs.append(n_input); volumes.append(volume)
        products[label] = {
            "input": input_density,
            "outputs": np.asarray(outputs),
            "recovered": np.asarray(recovered),
            "input_counts": np.asarray(inputs),
            "volume": volumes[0],
        }
        lo, med, hi = np.percentile(outputs, [16, 50, 84], axis=0)
        centres = (LOG_EDGES[:-1] + LOG_EDGES[1:]) / 2
        for i, centre in enumerate(centres):
            rows.append({
                "bandwidth": label, "log10_lcii_lsun": centre,
                "input_phi": input_density[i], "output_phi_median": med[i],
                "output_phi_lower16": lo[i], "output_phi_upper84": hi[i],
            })
    with (OUT / "p2_cii_lf.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    centres = (LOG_EDGES[:-1] + LOG_EDGES[1:]) / 2
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    for fwhm, style in zip([100, 200, 400], ["-", "--", ":"]):
        x, y = completeness_curve(fwhm)
        axes[0].plot(x, y, color="#333333", ls=style, label=f"{fwhm} km/s")
    axes[0].set(title="Common 7.5-GHz volume", xlabel=r"$\log_{10}(L_{\rm [CII]}/L_\odot)$",
                   ylabel="Completeness", ylim=(-.03,1.03))
    for label in CASES:
        p = products[label]; lo, med, hi = np.percentile(p["outputs"], [16,50,84], axis=0)
        axes[1].plot(centres, med, color=COLORS[label], marker="o", ms=3, label=label)
        axes[1].fill_between(centres, lo, hi, color=COLORS[label], alpha=.15)
    axes[1].plot(centres, products["7.5 GHz"]["input"], color="black", lw=2, label="Input benchmark LF")
    axes[1].set(title="Directly recovered [CII] LF", xlabel=r"$\log_{10}(L_{\rm [CII]}/L_\odot)$",
                   ylabel=r"$\phi$ [Mpc$^{-3}$ dex$^{-1}$]", yscale="log", ylim=(1e-8,5e-2))
    for ax in axes:
        ax.grid(alpha=.2)
        ax.legend(frameon=False)
    fig.suptitle(f"P2 [CII]: {CHANNEL_MHZ:.3f} MHz ({DV:.1f} km/s), channel RMS {CHANNEL_RMS*1e3:.2f} mJy")
    fig.savefig(OUT / "p2_cii_lf_and_volume.png", dpi=200)
    fig.savefig(OUT / "p2_cii_lf_and_volume.pdf")
    summary = {
        "channel_width_mhz": CHANNEL_MHZ, "channel_width_kms": DV,
        "channel_rms_jy": CHANNEL_RMS,
        "normalization_warning": "Hayatsu high-z amplitude is not fixed; counts use the local benchmark normalization",
        "cases": {k: {"z_range": frequency_to_z(CASES[k]),
                      "volume_mpc3": products[k]["volume"],
                      "recovered_count_median": float(np.median(products[k]["recovered"]))}
                  for k in CASES},
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    write_provenance(
        OUT / "p2_cii_lf_and_volume.provenance.json",
        root=ROOT,
        script=Path(__file__),
        command=sys.argv,
        inputs=[P1 / "source_recovery_all.csv", P1 / "summary.json"],
        outputs=[
            OUT / "p2_cii_lf.csv",
            OUT / "p2_cii_lf_and_volume.png",
            OUT / "p2_cii_lf_and_volume.pdf",
            OUT / "summary.json",
        ],
        parameters={
            "channel_width_mhz": CHANNEL_MHZ,
            "line_widths_kms": [100, 200, 400],
            "bandwidth_cases": CASES,
        },
        random_seeds=[20500000, 20600000],
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
