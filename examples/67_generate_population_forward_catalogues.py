"""Generate the fiducial physical and weighted selection catalogues.

Run from the repository root with::

    PYTHONPATH=src python examples/67_generate_population_forward_catalogues.py

This stage does not invoke CASA.  It freezes and validates the common sky
population that will later be passed through each facility operator.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bliss_alma.population_forward import (
    double_power_law_phi,
    draw_poisson_catalogue,
    expected_count,
    gaussian_line_peak_jy,
    line_luminosity_to_flux_jykms,
    modified_blackbody_flux_jy,
    saunders_phi,
    survey_volume_mpc3,
)


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs" / "population_forward" / "fiducial_survey.json"
OUTPUT = ROOT / "data" / "population_forward" / "catalogues"


def write_rows(path, rows):
    rows = list(rows)
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def physical_rows(config, population):
    survey = config["survey"]
    catalogue_config = config["catalogues"]
    area = survey["area_arcmin2"]
    field_count = survey["independent_fields"]
    base_seed = catalogue_config["master_seed"]
    realizations = catalogue_config["physical_realizations"]
    rows = []
    counts = []
    for realization in range(realizations):
        if population == "continuum":
            setup = config["continuum"]
            catalogue = draw_poisson_catalogue(
                saunders_phi,
                setup["log10_lir_lsun_min"],
                setup["log10_lir_lsun_max"],
                area,
                setup["redshift_min"],
                setup["redshift_max"],
                seed=base_seed + realization,
                field_count=field_count,
            )
        else:
            setup = config["line"]
            catalogue = draw_poisson_catalogue(
                double_power_law_phi,
                setup["log10_l_lsun_min"],
                setup["log10_l_lsun_max"],
                area,
                setup["redshift_min"],
                setup["redshift_max"],
                seed=base_seed + 10000 + realization,
                field_count=field_count,
            )
        counts.append(len(catalogue["redshift"]))
        rng = np.random.default_rng(base_seed + 20000 + realization)
        for source_index in range(len(catalogue["redshift"])):
            redshift = float(catalogue["redshift"][source_index])
            log_l = float(catalogue["log10_l_lsun"][source_index])
            row = {
                "population": population,
                "realization": realization,
                "source_id": f"{population[0].upper()}{realization:03d}-{source_index:04d}",
                "field_id": int(catalogue["field_id"][source_index]),
                "offset_ra_arcsec": float(catalogue["offset_ra_arcsec"][source_index]),
                "offset_dec_arcsec": float(catalogue["offset_dec_arcsec"][source_index]),
                "redshift": redshift,
                "log10_l_lsun": log_l,
                "lf_importance_weight": 1.0,
            }
            if population == "continuum":
                flux = modified_blackbody_flux_jy(
                    10**log_l,
                    redshift,
                    setup["observed_frequency_ghz"],
                    dust_temperature_k=setup["dust_temperature_k"],
                    beta=setup["emissivity_beta"],
                    include_cmb=setup["include_cmb"],
                )
                row["flux_jy"] = float(flux)
            else:
                fwhm = float(
                    10
                    ** rng.normal(
                        np.log10(setup["fwhm_kms_lognormal_median"]),
                        setup["fwhm_kms_lognormal_sigma_dex"],
                    )
                )
                integrated = float(line_luminosity_to_flux_jykms(10**log_l, redshift))
                row["observed_frequency_ghz"] = setup["rest_frequency_ghz"] / (
                    1.0 + redshift
                )
                row["fwhm_kms"] = fwhm
                row["integrated_flux_jykms"] = integrated
                row["peak_flux_jy"] = float(gaussian_line_peak_jy(integrated, fwhm))
            rows.append(row)
    return rows, np.asarray(counts), catalogue["expected_count"]


def selection_rows(config, population):
    """Create an oversampled catalogue and retain LF importance weights."""

    survey = config["survey"]
    catalogue_config = config["catalogues"]
    setup = config[population]
    field_count = survey["independent_fields"]
    per_field = catalogue_config["selection_catalogue_per_field"]
    count = field_count * per_field
    rng = np.random.default_rng(
        catalogue_config["master_seed"] + (30000 if population == "continuum" else 40000)
    )
    if population == "continuum":
        low, high = setup["log10_lir_lsun_min"], setup["log10_lir_lsun_max"]
        phi = saunders_phi
    else:
        low, high = setup["log10_l_lsun_min"], setup["log10_l_lsun_max"]
        phi = double_power_law_phi
    log_l = rng.uniform(low, high, count)
    redshift = rng.uniform(setup["redshift_min"], setup["redshift_max"], count)
    # Uniform proposal in logL and z.  The common multiplicative normalization
    # is irrelevant for completeness; retain normalized weights for recovery.
    weights = phi(log_l)
    weights /= np.mean(weights)
    side_arcsec = np.sqrt(survey["area_arcmin2"] / field_count) * 60.0
    field_id = np.repeat(np.arange(field_count), per_field)
    rows = []
    for index in range(count):
        row = {
            "population": population,
            "source_id": f"S{population[0].upper()}-{index:05d}",
            "field_id": int(field_id[index]),
            "offset_ra_arcsec": float(rng.uniform(-side_arcsec / 2, side_arcsec / 2)),
            "offset_dec_arcsec": float(rng.uniform(-side_arcsec / 2, side_arcsec / 2)),
            "redshift": float(redshift[index]),
            "log10_l_lsun": float(log_l[index]),
            "lf_importance_weight": float(weights[index]),
        }
        if population == "continuum":
            row["flux_jy"] = float(
                modified_blackbody_flux_jy(
                    10 ** log_l[index],
                    redshift[index],
                    setup["observed_frequency_ghz"],
                    dust_temperature_k=setup["dust_temperature_k"],
                    beta=setup["emissivity_beta"],
                    include_cmb=setup["include_cmb"],
                )
            )
        else:
            fwhm = float(
                10
                ** rng.normal(
                    np.log10(setup["fwhm_kms_lognormal_median"]),
                    setup["fwhm_kms_lognormal_sigma_dex"],
                )
            )
            integrated = float(
                line_luminosity_to_flux_jykms(10 ** log_l[index], redshift[index])
            )
            row.update(
                {
                    "observed_frequency_ghz": setup["rest_frequency_ghz"]
                    / (1.0 + redshift[index]),
                    "fwhm_kms": fwhm,
                    "integrated_flux_jykms": integrated,
                    "peak_flux_jy": float(
                        gaussian_line_peak_jy(integrated, fwhm)
                    ),
                }
            )
        rows.append(row)
    return rows


def make_validation_figure(config, continuum_rows, line_rows, counts):
    figure, axes = plt.subplots(2, 2, figsize=(10.5, 8.2))
    ax = axes[0, 0]
    grid = np.linspace(8.0, 13.5, 500)
    ax.semilogy(grid, saunders_phi(grid), color="#d62728")
    ax.set(xlabel=r"$\log_{10}(L_{\rm IR}/L_\odot)$", ylabel=r"$\phi$ [Mpc$^{-3}$ dex$^{-1}$]")
    ax.set_title("ALPINE total-IR LF (4.5 < z < 6.0)")

    ax = axes[0, 1]
    grid = np.linspace(6.0, 10.5, 500)
    ax.semilogy(grid, double_power_law_phi(grid), color="#1f77b4")
    ax.set(xlabel=r"$\log_{10}(L_{\rm [CII]}/L_\odot)$", ylabel=r"$\phi$ [Mpc$^{-3}$ dex$^{-1}$]")
    ax.set_title("[CII] benchmark LF (normalization uncertain)")

    ax = axes[1, 0]
    continuum_flux = np.array([float(row["flux_jy"]) for row in continuum_rows]) * 1e3
    line_flux = np.array([float(row["integrated_flux_jykms"]) for row in line_rows])
    ax.hist(continuum_flux, bins=np.geomspace(max(continuum_flux.min(), 1e-7), continuum_flux.max(), 35), histtype="step", color="#d62728", label="continuum (mJy)")
    ax.hist(line_flux, bins=np.geomspace(max(line_flux.min(), 1e-7), line_flux.max(), 35), histtype="step", color="#1f77b4", label=r"[CII] (Jy km s$^{-1}$)")
    ax.set_xscale("log")
    ax.set(xlabel="Observed flux (population-specific units)", ylabel="Sources / bin")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    bins = np.arange(-0.5, max(counts["continuum"].max(), counts["line"].max()) + 1.5)
    ax.hist(counts["continuum"], bins=bins, histtype="step", color="#d62728", label="continuum")
    ax.hist(counts["line"], bins=bins, histtype="step", color="#1f77b4", label="[CII] benchmark")
    ax.axvline(counts["continuum"].mean(), color="#d62728", linestyle=":")
    ax.axvline(counts["line"].mean(), color="#1f77b4", linestyle=":")
    ax.set(xlabel="Physical sources per 10 arcmin$^2$ realization", ylabel="Realizations")
    ax.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(OUTPUT / "population_input_validation.png", dpi=200)
    figure.savefig(OUTPUT / "population_input_validation.pdf")
    plt.close(figure)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    config = json.loads(CONFIG_PATH.read_text())
    continuum_rows, continuum_counts, continuum_expected = physical_rows(
        config, "continuum"
    )
    line_rows, line_counts, line_expected = physical_rows(config, "line")
    write_rows(OUTPUT / "physical_continuum_catalogues.csv", continuum_rows)
    write_rows(OUTPUT / "physical_line_catalogues.csv", line_rows)
    write_rows(
        OUTPUT / "selection_continuum_catalogue.csv",
        selection_rows(config, "continuum"),
    )
    write_rows(
        OUTPUT / "selection_line_catalogue.csv", selection_rows(config, "line")
    )
    summary = {
        "continuum": {
            "survey_volume_mpc3": survey_volume_mpc3(10.0, 4.5, 6.0),
            "expected_sources_per_realization": continuum_expected,
            "realized_mean": float(continuum_counts.mean()),
            "realized_std": float(continuum_counts.std(ddof=1)),
        },
        "line": {
            "survey_volume_mpc3": survey_volume_mpc3(
                10.0,
                config["line"]["redshift_min"],
                config["line"]["redshift_max"],
            ),
            "expected_sources_per_realization": line_expected,
            "realized_mean": float(line_counts.mean()),
            "realized_std": float(line_counts.std(ddof=1)),
            "normalization_warning": config["line"]["normalization_status"],
        },
        "selection_sources_per_population": (
            config["survey"]["independent_fields"]
            * config["catalogues"]["selection_catalogue_per_field"]
        ),
    }
    (OUTPUT / "population_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    make_validation_figure(
        config,
        continuum_rows,
        line_rows,
        {"continuum": continuum_counts, "line": line_counts},
    )
    print(json.dumps(summary, indent=2))
    print(f"Wrote catalogues and validation plots to {OUTPUT}")


if __name__ == "__main__":
    main()
