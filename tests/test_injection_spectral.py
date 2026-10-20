"""Tests for CASA-independent image injection and absorption helpers."""

import numpy as np
import pytest

from bliss_alma.injection import (
    gaussian_line_integrated_flux,
    gaussian_line_flux_density,
    inject_source_image,
    plan_line_injection_positions,
)
from bliss_alma.recovery import line_injection_crosstalk
from bliss_alma.spectral import absorption_cube, extract_absorption_lines


def test_image_injection_places_model_without_modifying_input():
    data = np.zeros((7, 7))
    model = np.ones((3, 3))
    result = inject_source_image(data, model, position=(3, 4))
    assert np.all(data == 0)
    assert np.sum(result) == 9
    assert result[3, 4] == 1


def test_cube_absorption_injection_and_edge_clipping():
    cube = np.ones((5, 6, 6)) * 10
    model = np.ones((3, 3, 3)) * 2
    result = inject_source_image(cube, model, position=(0, 0, 0), operation="subtract")
    assert result.min() == 8
    assert np.count_nonzero(result == 8) == 8


def test_gaussian_line_profile_has_requested_peak_and_integral():
    frequencies = np.linspace(229.5e9, 230.5e9, 20001)
    profile = gaussian_line_flux_density(
        frequencies,
        center_frequency_hz=230.0e9,
        fwhm_kms=300.0,
        integrated_flux_jy_kms=0.1,
    )
    velocity_step = 299792.458 * np.diff(frequencies).mean() / 230.0e9
    expected_peak = 0.1 / (300.0 * np.sqrt(np.pi / (4 * np.log(2))))
    assert np.isclose(profile.max(), expected_peak)
    assert np.isclose(np.sum(profile) * velocity_step, 0.1, rtol=2e-3)


def test_gaussian_line_profile_rejects_centre_outside_channels():
    with pytest.raises(ValueError):
        gaussian_line_flux_density(
            [229.0e9, 229.1e9, 229.2e9],
            center_frequency_hz=230.0e9,
            fwhm_kms=300.0,
            integrated_flux_jy_kms=0.1,
        )


def test_peak_snr_converts_to_gaussian_integrated_flux():
    integrated = gaussian_line_integrated_flux(
        peak_snr=5.0,
        rms_jy_per_beam=20e-6,
        fwhm_kms=300.0,
    )
    expected = 5.0 * 20e-6 * 300.0 * np.sqrt(np.pi / (4 * np.log(2)))
    assert integrated == pytest.approx(expected)


def test_line_injection_crosstalk_measures_off_target_response():
    cubes = np.zeros((3, 12, 8, 8))
    positions = [(2, 2, 2), (6, 4, 4), (9, 6, 6)]
    for index, position in enumerate(positions):
        cubes[(index, *position)] = 5.0
    cubes[(1, *positions[0])] = 0.12
    cubes[(2, *positions[0])] = -0.03
    cubes[(0, *positions[1])] = 0.25
    result = line_injection_crosstalk(
        cubes, positions, rms=1.0, max_total_crosstalk_snr=0.2,
    )
    assert result.matrix_snr[0, 1] == pytest.approx(0.12)
    assert result.off_diagonal_sum_snr[0] == pytest.approx(0.15)
    assert result.safe.tolist() == [True, False, True]


def test_line_placement_is_reproducible_and_avoids_overlap():
    widths = [3.0, 15.0, 45.0, 7.0, 25.0]
    kwargs = dict(
        cube_shape=(256, 40, 40),
        min_spatial_separation_pixels=8.0,
        spectral_extent_fwhm=3.0,
        spatial_edge_margin_pixels=3,
        seed=2040,
    )
    first = plan_line_injection_positions(widths, **kwargs)
    second = plan_line_injection_positions(widths, **kwargs)
    assert first == second
    assert [row.source_index for row in first] == list(range(len(widths)))
    for index, left in enumerate(first):
        half_left = 1.5 * left.fwhm_channels
        assert half_left <= left.channel <= 255 - half_left
        assert 3 <= left.x < 37 and 3 <= left.y < 37
        for right in first[index + 1:]:
            sky = np.hypot(left.x - right.x, left.y - right.y)
            overlap = abs(left.channel - right.channel) < 1.5 * (
                left.fwhm_channels + right.fwhm_channels
            )
            assert not (sky < 8.0 and overlap)


def test_line_placement_rejects_line_broader_than_cube():
    with pytest.raises(ValueError, match="too broad"):
        plan_line_injection_positions([40], cube_shape=(50, 10, 10))


def test_line_placement_creates_interleaved_helices_when_needed():
    rows = plan_line_injection_positions(
        [10.0] * 8,
        cube_shape=(80, 60, 60),
        min_spatial_separation_pixels=6.0,
        spectral_extent_fwhm=3.0,
        seed=7,
    )
    helices = sorted({row.helix_index for row in rows})
    assert len(helices) > 1
    for helix in helices:
        channels = [row.channel for row in rows if row.helix_index == helix]
        assert channels == sorted(channels)


def test_absorption_representations():
    continuum = np.ones((5, 5)) * 10
    cube = np.ones((3, 5, 5)) * 10
    cube[1, 2, 2] = 5
    fractional = absorption_cube(cube, continuum, "fractional_depth")
    tau = absorption_cube(cube, continuum, "optical_depth")
    assert fractional[1, 2, 2] == pytest.approx(0.5)
    assert tau[1, 2, 2] == pytest.approx(np.log(2))


def test_absorption_extraction_uses_positive_depth():
    continuum = np.ones((5, 7, 7)) * 10
    cube = continuum.copy()
    cube[2, 3, 3] = 4
    result = extract_absorption_lines(
        cube, continuum, representation="residual_depth", rms=1,
        seed_sigma=5, flood_sigma=5, min_pixels=1,
    )
    assert len(result.catalogue) == 1
    assert result.catalogue[0].peak_index == (2, 3, 3)
