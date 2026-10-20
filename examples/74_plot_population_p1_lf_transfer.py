"""Propagate the ALPINE input LF through the empirical P1b response.

This is a Hayatsu-et-al.-style input/output demonstration for the five-minute
P1 control experiment.  It is not yet the 100-hour survey forecast.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u

from bliss_alma.population_forward import (
    FIDUCIAL_COSMOLOGY,
    modified_blackbody_flux_jy,
    saunders_phi,
)


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "population_forward" / "p1_geometry"
FIGURES = INPUT / "figures"
CONFIGS = ["current_ALMA_C7", "ALMA2040_F4"]
LABELS = {
    "current_ALMA_C7": "Current ALMA C-7 (43 antennas)",
    "ALMA2040_F4": "ALMA2040 F4 (255 antennas)",
}
COLORS = {"current_ALMA_C7": "#d62728", "ALMA2040_F4": "#1f77b4"}
LOG_L_EDGES = np.arange(8.0, 13.51, 0.25)
N_MONTE_CARLO = 200_000


def read_rows(name):
    with (INPUT / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def volume_redshifts(rng, count):
    grid = np.linspace(4.5, 6.0, 4097)
    weights = FIDUCIAL_COSMOLOGY.differential_comoving_volume(grid).to_value(
        u.Mpc**3 / u.sr
    )
    cdf = np.cumsum(weights)
    cdf /= cdf[-1]
    return np.interp(rng.random(count), cdf, grid)


def make_master_samples():
    rng = np.random.default_rng(20400000)
    log_l = rng.uniform(LOG_L_EDGES[0], LOG_L_EDGES[-1], N_MONTE_CARLO)
    redshift = volume_redshifts(rng, N_MONTE_CARLO)
    # A circular Gaussian PB and sources uniform in area within PB >= 0.2
    # imply that log(PB) is uniform.
    pb = np.exp(rng.uniform(np.log(0.2), 0.0, N_MONTE_CARLO))
    z_grid = np.linspace(4.5, 6.0, 801)
    flux_per_lsun = modified_blackbody_flux_jy(
        np.ones_like(z_grid), z_grid, 230.0
    )
    flux = 10**log_l * np.interp(redshift, z_grid, flux_per_lsun)
    weights = saunders_phi(log_l)
    return log_l, flux * pb, weights


def response_classes(rows, config, realization):
    selected = [
        row
        for row in rows
        if row["track"] == "P1b_array_total"
        and row["configuration"] == config
        and int(row["realization"]) == realization
    ]
    classes = sorted({float(row["target_current_alma_snr"]) for row in selected})
    grouped = {
        value: [
            row for row in selected
            if float(row["target_current_alma_snr"]) == value
        ]
        for value in classes
    }
    return np.asarray(classes), grouped


def propagate_one(log_l, apparent_flux, weights, rows, config, realization, fiducial_rms):
    rng = np.random.default_rng(20460000 + realization + 100 * CONFIGS.index(config))
    target_snr = apparent_flux / fiducial_rms
    classes, grouped = response_classes(rows, config, realization)
    midpoints = (classes[:-1] + classes[1:]) / 2.0
    class_index = np.searchsorted(midpoints, target_snr)
    detected = np.zeros(len(log_l), dtype=bool)
    output_log_l = np.full(len(log_l), np.nan)
    analysis_rms = float(
        next(
            row["rms_jy_beam"]
            for row in rows
            if row["track"] == "P1b_array_total"
            and row["configuration"] == config
            and int(row["realization"]) == realization
        )
    )
    for index, class_snr in enumerate(classes):
        members = np.flatnonzero(class_index == index)
        if members.size == 0:
            continue
        if class_snr == classes[0]:
            members = members[target_snr[members] >= classes[0]]
        pool = grouped[class_snr]
        choices = rng.integers(0, len(pool), len(members))
        for member, choice in zip(members, choices):
            response = pool[choice]
            if response["recovered"].lower() != "true":
                continue
            measured_flux = float(response["recovered_peak_snr"]) * analysis_rms
            reference_input_flux = class_snr * fiducial_rms
            ratio = measured_flux / reference_input_flux
            detected[member] = True
            output_log_l[member] = log_l[member] + np.log10(max(ratio, 1e-6))
    # Above the sampled response grid the high-S/N limit is complete with
    # negligible fractional flux scatter for this purpose.
    high = target_snr > classes[-1]
    detected[high] = True
    output_log_l[high] = log_l[high]

    width = np.diff(LOG_L_EDGES)
    proposal_density = 1.0 / (LOG_L_EDGES[-1] - LOG_L_EDGES[0])
    input_density = (
        np.histogram(log_l, bins=LOG_L_EDGES, weights=weights)[0]
        / len(log_l)
        / width
        / proposal_density
    )
    output_density = (
        np.histogram(
            output_log_l[detected],
            bins=LOG_L_EDGES,
            weights=weights[detected],
        )[0]
        / len(log_l)
        / width
        / proposal_density
    )
    input_bin_weight = np.histogram(log_l, bins=LOG_L_EDGES, weights=weights)[0]
    detected_bin_weight = np.histogram(
        log_l[detected], bins=LOG_L_EDGES, weights=weights[detected]
    )[0]
    completeness = np.divide(
        detected_bin_weight,
        input_bin_weight,
        out=np.zeros_like(detected_bin_weight),
        where=input_bin_weight > 0,
    )
    return input_density, output_density, completeness


def empirical_completeness(rows, config):
    classes = sorted(
        {
            float(row["target_current_alma_snr"])
            for row in rows
            if row["track"] == "P1b_array_total"
            and row["configuration"] == config
        }
    )
    values = []
    for value in classes:
        selected = [
            row
            for row in rows
            if row["track"] == "P1b_array_total"
            and row["configuration"] == config
            and float(row["target_current_alma_snr"]) == value
        ]
        values.append(
            np.mean([row["recovered"].lower() == "true" for row in selected])
        )
    return np.asarray(classes), np.asarray(values)


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    rows = read_rows("source_recovery_all.csv")
    summary = __import__("json").loads((INPUT / "summary.json").read_text())
    fiducial_rms = summary["fiducial_current_alma_rms_jy_beam"]
    log_l, apparent_flux, weights = make_master_samples()
    centres = (LOG_L_EDGES[:-1] + LOG_L_EDGES[1:]) / 2.0

    products = {}
    for config in CONFIGS:
        outputs = []
        completeness = []
        input_density = None
        for realization in range(10):
            input_density, output, recovery = propagate_one(
                log_l,
                apparent_flux,
                weights,
                rows,
                config,
                realization,
                fiducial_rms,
            )
            outputs.append(output)
            completeness.append(recovery)
        products[config] = {
            "input": input_density,
            "output": np.asarray(outputs),
            "completeness": np.asarray(completeness),
        }

    lf_rows = []
    response_rows = []
    for config in CONFIGS:
        output_low, output_median, output_high = np.percentile(
            products[config]["output"], [16, 50, 84], axis=0
        )
        ratio = np.divide(
            products[config]["output"],
            products[config]["input"],
            out=np.zeros_like(products[config]["output"]),
            where=products[config]["input"] > 0,
        )
        ratio_low, ratio_median, ratio_high = np.percentile(
            ratio, [16, 50, 84], axis=0
        )
        for index, centre in enumerate(centres):
            lf_rows.append(
                {
                    "configuration": config,
                    "log10_lir_lsun": float(centre),
                    "input_phi_mpc3_dex": float(products[config]["input"][index]),
                    "output_phi_median": float(output_median[index]),
                    "output_phi_lower16": float(output_low[index]),
                    "output_phi_upper84": float(output_high[index]),
                    "output_input_median": float(ratio_median[index]),
                    "output_input_lower16": float(ratio_low[index]),
                    "output_input_upper84": float(ratio_high[index]),
                }
            )
        classes = sorted(
            {
                float(row["target_current_alma_snr"])
                for row in rows
                if row["track"] == "P1b_array_total"
                and row["configuration"] == config
            }
        )
        for value in classes:
            per_realization = []
            for realization in range(10):
                selected = [
                    row
                    for row in rows
                    if row["track"] == "P1b_array_total"
                    and row["configuration"] == config
                    and int(row["realization"]) == realization
                    and float(row["target_current_alma_snr"]) == value
                ]
                per_realization.append(
                    np.mean(
                        [row["recovered"].lower() == "true" for row in selected]
                    )
                )
            low, median, high = np.percentile(
                per_realization, [16, 50, 84]
            )
            response_rows.append(
                {
                    "configuration": config,
                    "target_current_alma_snr": value,
                    "completeness_median": float(median),
                    "completeness_lower16": float(low),
                    "completeness_upper84": float(high),
                    "classification_axes": "input_SNR,PB_position,noise_realization",
                }
            )
    for path, output_rows in (
        (INPUT / "p1_alpine_lf_transfer.csv", lf_rows),
        (INPUT / "p1_response_classes.csv", response_rows),
    ):
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
            writer.writeheader()
            writer.writerows(output_rows)

    figure, axes = plt.subplots(2, 2, figsize=(11.2, 8.5), constrained_layout=True)
    for axis, config in zip(axes[0], CONFIGS):
        product = products[config]
        low, median, high = np.percentile(product["output"], [16, 50, 84], axis=0)
        axis.plot(centres, product["input"], color="black", linewidth=2, label="Input ALPINE LF")
        axis.plot(centres, median, color=COLORS[config], marker="o", markersize=3, label="Directly recovered LF")
        axis.fill_between(centres, low, high, color=COLORS[config], alpha=0.2)
        axis.set(
            title=LABELS[config],
            xlabel=r"$\log_{10}(L_{\rm IR}/L_\odot)$",
            ylabel=r"$\phi$ [Mpc$^{-3}$ dex$^{-1}$]",
            yscale="log",
            ylim=(1e-7, 8e-3),
        )
        axis.grid(alpha=0.2)
        axis.legend(frameon=False)

    axis = axes[1, 0]
    for config in CONFIGS:
        ratio = np.divide(
            products[config]["output"],
            products[config]["input"],
            out=np.zeros_like(products[config]["output"]),
            where=products[config]["input"] > 0,
        )
        low, median, high = np.percentile(ratio, [16, 50, 84], axis=0)
        axis.plot(centres, median, color=COLORS[config], marker="o", markersize=3, label=LABELS[config])
        axis.fill_between(centres, low, high, color=COLORS[config], alpha=0.18)
    axis.axhline(1.0, color="black", linewidth=1)
    axis.set(
        title="Recovered / input LF",
        xlabel=r"$\log_{10}(L_{\rm IR}/L_\odot)$",
        ylabel="Output / input",
        ylim=(-0.04, 1.3),
    )
    axis.grid(alpha=0.2)
    axis.legend(frameon=False)

    axis = axes[1, 1]
    for config in CONFIGS:
        snr, value = empirical_completeness(rows, config)
        axis.plot(snr, value, color=COLORS[config], marker="o", label=LABELS[config])
    axis.axvline(4.0, color="#666666", linestyle="--")
    axis.set(
        title="Empirical P1b transfer function",
        xlabel="Fixed flux in current-ALMA fiducial S/N",
        ylabel="Completeness",
        xlim=(0, 7.1),
        ylim=(-0.04, 1.04),
    )
    axis.grid(alpha=0.2)
    axis.legend(frameon=False)

    figure.suptitle(
        "ALPINE input/output LF through the five-minute P1 transfer function"
    )
    figure.savefig(FIGURES / "p1_alpine_input_output_lf.png", dpi=200)
    figure.savefig(FIGURES / "p1_alpine_input_output_lf.pdf")
    plt.close(figure)
    print(FIGURES / "p1_alpine_input_output_lf.png")


if __name__ == "__main__":
    main()
