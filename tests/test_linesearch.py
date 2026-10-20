import numpy as np

from bliss_alma.linesearch import (
    local_robust_maps,
    merge_candidates,
    robust_sigma,
    search_cube,
)


def test_signed_cube_search_and_merge():
    rng = np.random.default_rng(2040)
    cube = rng.normal(0.0, 1.0, (16, 24, 24))
    cube[7:10, 8, 9] += 12.0
    cube[3:6, 16, 14] -= 12.0
    rows = search_cube(cube, np.ones((24, 24), bool), [1, 3], 5.0)
    merged = merge_candidates(rows, spatial_radius_pixels=2.0)
    assert any(row["snr"] > 0 and row["x"] == 9 for row in merged)
    assert any(row["snr"] < 0 and row["x"] == 14 for row in merged)


def test_robust_sigma_rejects_single_extreme_tail():
    values = np.r_[np.linspace(-1, 1, 101), 1e6]
    assert robust_sigma(values) < 2.0


def test_local_robust_map_tracks_spatial_noise_change():
    rng = np.random.default_rng(2041)
    image = rng.normal(size=(128, 128))
    image[:, 64:] *= 3.0
    mask = np.ones(image.shape, dtype=bool)
    median, sigma = local_robust_maps(image, mask, 31)
    normalized = (image - median) / sigma
    left = robust_sigma(normalized[16:-16, 16:48])
    right = robust_sigma(normalized[16:-16, 80:112])
    assert 0.75 < left < 1.25
    assert 0.75 < right < 1.25
    assert np.median(sigma[16:-16, 80:112]) > 2.0 * np.median(
        sigma[16:-16, 16:48]
    )


def test_local_search_requires_window_size():
    cube = np.zeros((3, 16, 16), dtype=float)
    mask = np.ones((16, 16), dtype=bool)
    with np.testing.assert_raises(ValueError):
        search_cube(cube, mask, [1], 4.0, noise_mode="local")
