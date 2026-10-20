"""Numerical checks for source-population forward modelling."""

import numpy as np

from bliss_alma.population_forward import (
    ALPINE_IR_Z45_60,
    HAYATSU_CII_BENCHMARK,
    double_power_law_phi,
    gaussian_line_peak_jy,
    line_luminosity_to_flux_jykms,
    modified_blackbody_flux_jy,
    saunders_phi,
    survey_volume_mpc3,
)


def test_lf_normalizations_at_lstar():
    alpine = saunders_phi(ALPINE_IR_Z45_60.log10_lstar_lsun)
    expected_alpine = 10**ALPINE_IR_Z45_60.log10_phistar_mpc3_dex * np.exp(
        -(np.log10(2.0) ** 2) / (2.0 * ALPINE_IR_Z45_60.sigma**2)
    )
    assert np.isclose(alpine, expected_alpine)
    assert np.isclose(
        double_power_law_phi(HAYATSU_CII_BENCHMARK.log10_lstar_lsun),
        0.5 * 10**HAYATSU_CII_BENCHMARK.log10_phistar_mpc3_dex,
    )


def test_survey_volume_scales_with_area():
    one = survey_volume_mpc3(1.0, 6.1, 6.5)
    ten = survey_volume_mpc3(10.0, 6.1, 6.5)
    assert np.isclose(ten, 10.0 * one)


def test_line_flux_round_trip_reference_value():
    flux = line_luminosity_to_flux_jykms(1e9, 6.3)
    assert 0.1 < flux < 10.0
    peak = gaussian_line_peak_jy(flux, 200.0)
    assert np.isclose(flux, peak * 200.0 * np.sqrt(np.pi / (4.0 * np.log(2.0))))


def test_modified_blackbody_is_positive_and_linear():
    flux = modified_blackbody_flux_jy(
        np.array([1e11, 2e11]), np.array([5.0, 5.0]), 230.0
    )
    assert np.all(flux > 0)
    assert np.isclose(flux[1] / flux[0], 2.0)
