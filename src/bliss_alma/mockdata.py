"""
Mock-data utilities for DICE.

This module will provide helper functions for describing and organizing
mock ALMA/GHOST datasets used for source-injection and recovery tests.
"""

import math


def describe_mock_setup(
    configuration,
    output_type,
    frequency=None,
    bandwidth=None,
    observing_time=None,
    noise_level=None,
    field_type="blank",
    notes=None,
):
    """
    Store a simple description of a mock ALMA/GHOST setup.

    Parameters
    ----------
    configuration : str
        Name of the array configuration, e.g. 'F4'.
    output_type : str
        Type of downloaded or generated output, e.g. 'FITS image' or 'MS'.
    frequency : str, optional
        Frequency setup.
    bandwidth : str, optional
        Bandwidth setup.
    observing_time : str, optional
        Observing time.
    noise_level : str, optional
        Noise level if specified.
    field_type : str, optional
        Type of field, e.g. 'blank' or 'noisy'.
    notes : str, optional
        Additional notes.

    Returns
    -------
    dict
        Dictionary describing the mock setup.
    """

    return {
        "configuration": configuration,
        "output_type": output_type,
        "frequency": frequency,
        "bandwidth": bandwidth,
        "observing_time": observing_time,
        "noise_level": noise_level,
        "field_type": field_type,
        "notes": notes,
    }


def ghost_output_placeholder():
    """
    Placeholder for future GHOST output handling.

    Future versions may include utilities to inspect downloaded GHOST
    FITS images, measurement sets, or source catalogs.
    """

    raise NotImplementedError(
        "GHOST output handling will be implemented after the first test dataset is available."
    )


def uniform_grid_positions(center, spacing_pixels, radius_pixels):
    """Return deterministic, uniformly spaced 2D positions within a circle."""
    if spacing_pixels <= 0:
        raise ValueError("spacing_pixels must be positive.")
    if radius_pixels < 0:
        raise ValueError("radius_pixels must be non-negative.")

    extent = int(math.floor(radius_pixels / spacing_pixels))
    positions = []
    for ix in range(-extent, extent + 1):
        for iy in range(-extent, extent + 1):
            dx = ix * spacing_pixels
            dy = iy * spacing_pixels
            if math.hypot(dx, dy) <= radius_pixels:
                positions.append((int(center[0] + dx), int(center[1] + dy)))
    return positions
