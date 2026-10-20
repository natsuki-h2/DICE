"""Target-neutral helpers for signed, multi-width spectral-cube searches."""

from __future__ import annotations

from collections import deque
from typing import Iterable, Mapping

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def robust_sigma(values) -> float:
    """Return the Gaussian-equivalent median absolute deviation."""

    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if not array.size:
        return float("nan")
    median = np.median(array)
    return float(1.4826 * np.median(np.abs(array - median)))


def merge_candidates(rows: Iterable[Mapping], spatial_radius_pixels=5.0):
    """Merge same-polarity detections overlapping in space and frequency."""

    kept = []
    for source in sorted(rows, key=lambda item: abs(item["snr"]), reverse=True):
        row = dict(source)
        duplicate = any(
            row["spw"] == old["spw"]
            and np.sign(row["snr"]) == np.sign(old["snr"])
            and np.hypot(row["x"] - old["x"], row["y"] - old["y"])
            <= spatial_radius_pixels
            and abs(row["channel"] - old["channel"])
            <= 0.5 * (row["width_channels"] + old["width_channels"])
            for old in kept
        )
        if not duplicate:
            kept.append(row)
    return sorted(kept, key=lambda item: abs(item["snr"]), reverse=True)


def local_robust_maps(image, spatial_mask, box_size_pixels):
    """Estimate spatially varying median and Gaussian-equivalent MAD maps.

    The maps are evaluated on the already spectrally summed image, so the
    noise normalization follows the selected line width. Robust estimates are
    measured on an overlapping grid using only unmasked pixels and then
    interpolated. This avoids the artificially small MAD values produced when
    a PB-mask boundary is filled with a constant before median filtering.
    """

    array = np.asarray(image, dtype=float)
    mask = np.asarray(spatial_mask, dtype=bool) & np.isfinite(array)
    size = int(box_size_pixels)
    if size < 3:
        raise ValueError("noise_box_pixels must be at least 3")
    if size % 2 == 0:
        size += 1
    values = array[mask]
    if not values.size:
        shape = array.shape
        return np.full(shape, np.nan), np.full(shape, np.nan)
    global_median = float(np.median(values))
    global_sigma = robust_sigma(values)
    ny, nx = array.shape
    stride = max(1, size // 2)
    y_grid = np.unique(np.r_[np.arange(0, ny, stride), ny - 1]).astype(int)
    x_grid = np.unique(np.r_[np.arange(0, nx, stride), nx - 1]).astype(int)
    median_grid = np.full((len(y_grid), len(x_grid)), global_median)
    sigma_grid = np.full((len(y_grid), len(x_grid)), global_sigma)
    half = size // 2
    minimum_samples = max(20, int(0.1 * size * size))
    for iy, y in enumerate(y_grid):
        y0, y1 = max(0, y - half), min(ny, y + half + 1)
        for ix, x in enumerate(x_grid):
            x0, x1 = max(0, x - half), min(nx, x + half + 1)
            local_values = array[y0:y1, x0:x1][mask[y0:y1, x0:x1]]
            if local_values.size < minimum_samples:
                continue
            local_median = float(np.median(local_values))
            local_sigma = robust_sigma(local_values)
            if np.isfinite(local_sigma) and local_sigma > 0:
                median_grid[iy, ix] = local_median
                sigma_grid[iy, ix] = local_sigma
    yy, xx = np.indices(array.shape)
    points = np.column_stack((yy.ravel(), xx.ravel()))
    local_median = RegularGridInterpolator(
        (y_grid, x_grid), median_grid, bounds_error=False, fill_value=None
    )(points).reshape(array.shape)
    local_sigma = RegularGridInterpolator(
        (y_grid, x_grid), sigma_grid, bounds_error=False, fill_value=None
    )(points).reshape(array.shape)
    return local_median, local_sigma


def search_cube(
    cube,
    spatial_mask,
    widths,
    threshold,
    *,
    spw=0,
    noise_mode="global",
    noise_box_pixels=None,
):
    """Search a ``(channel, y, x)`` cube for positive and negative peaks.

    Returned rows contain pixel/channel measurements. Coordinate conversion is
    deliberately left to the caller because WCS conventions vary among cubes.
    """

    data = np.asarray(cube)
    mask = np.asarray(spatial_mask, dtype=bool)
    if data.ndim != 3 or mask.shape != data.shape[1:]:
        raise ValueError("cube must be (channel, y, x) and mask must be (y, x)")
    widths = tuple(sorted({int(width) for width in widths}))
    if not widths or widths[0] < 1:
        raise ValueError("widths must contain positive integers")
    if noise_mode not in {"global", "local"}:
        raise ValueError("noise_mode must be 'global' or 'local'")
    if noise_mode == "local" and noise_box_pixels is None:
        raise ValueError("noise_box_pixels is required for local noise mode")
    channel_rms = np.array([
        robust_sigma(np.asarray(plane)[mask & np.isfinite(plane) & (plane != 0)])
        for plane in data
    ])
    valid = np.isfinite(channel_rms) & (channel_rms > 0)
    if not np.any(valid):
        return []
    first_valid, last_valid = np.where(valid)[0][[0, -1]]
    rolling = {width: np.zeros(data.shape[1:], dtype=float) for width in widths}
    ring = deque(maxlen=max(widths))
    rows = []
    for channel, input_plane in enumerate(data):
        plane = np.where(np.isfinite(input_plane), input_plane, 0.0)
        for width in widths:
            rolling[width] += plane
            if len(ring) >= width:
                rolling[width] -= ring[-width]
            start = channel - width + 1
            if start < first_valid or channel > last_valid:
                continue
            values = rolling[width][mask]
            if noise_mode == "global":
                median, rms = float(np.median(values)), robust_sigma(values)
                if not np.isfinite(rms) or rms <= 0:
                    continue
                snr = (rolling[width] - median) / rms
            else:
                median_map, rms_map = local_robust_maps(
                    rolling[width], mask, noise_box_pixels
                )
                snr = (rolling[width] - median_map) / rms_map
            centre = channel - (width - 1) / 2.0
            for polarity in (1, -1):
                yy, xx = np.where(mask & (polarity * snr >= threshold))
                rows.extend({
                    "spw": spw, "x": int(x), "y": int(y),
                    "channel": float(centre), "width_channels": width,
                    "snr": float(snr[y, x]),
                } for y, x in zip(yy, xx))
        ring.append(plane)
    return rows
