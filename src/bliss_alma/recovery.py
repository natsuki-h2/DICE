"""
Source-recovery metrics for DICE.

This module provides utilities for measuring recovery rate, recovered-to-
injected flux ratio, and positional offsets for injected sources.
"""

import math
from collections import namedtuple

import numpy as np


SourceMatch = namedtuple(
    "SourceMatch",
    ["input_index", "detection_index", "distance_pixels"],
)

LineSourceMatch = namedtuple(
    "LineSourceMatch",
    [
        "input_index", "detection_index", "spatial_distance_pixels",
        "channel_distance", "normalized_distance",
    ],
)

LineCompletenessResult = namedtuple(
    "LineCompletenessResult",
    ["completeness", "n_injected", "n_recovered", "matches", "recovered"],
)

LineInjectionCrosstalk = namedtuple(
    "LineInjectionCrosstalk",
    ["matrix_snr", "off_diagonal_sum_snr", "max_off_diagonal_snr", "safe"],
)


def line_injection_crosstalk(
    isolated_response_cubes,
    target_positions,
    *,
    rms,
    max_total_crosstalk_snr=0.2,
):
    """Measure source-to-source cross-talk from isolated response cubes.

    ``isolated_response_cubes[j]`` must contain only artificial source ``j``
    after the same imaging operation used for the batched experiment.
    ``target_positions[i]`` is the ``(channel, y, x)`` peak position of source
    ``i``. Matrix element ``[i, j]`` is therefore the response of source ``j``
    at target ``i``, divided by the supplied image RMS. The diagonal is the
    intended source response; off-diagonal values quantify dirty-beam and
    spectral leakage. A target is safe when the sum of absolute off-diagonal
    contributions is no larger than ``max_total_crosstalk_snr``.
    """
    cubes = np.asarray(isolated_response_cubes, dtype=float)
    positions = [tuple(int(value) for value in position) for position in target_positions]
    if cubes.ndim != 4:
        raise ValueError(
            "isolated_response_cubes must have shape (n_source, n_channel, n_y, n_x)."
        )
    if cubes.shape[0] != len(positions) or not positions:
        raise ValueError("Provide one target position for each isolated response cube.")
    if not np.isfinite(rms) or rms <= 0:
        raise ValueError("rms must be finite and positive.")
    if not np.isfinite(max_total_crosstalk_snr) or max_total_crosstalk_snr < 0:
        raise ValueError("max_total_crosstalk_snr must be finite and non-negative.")
    matrix = np.empty((len(positions), len(positions)), dtype=float)
    for target_index, position in enumerate(positions):
        if len(position) != 3:
            raise ValueError("Each target position must be (channel, y, x).")
        if any(index < 0 or index >= size for index, size in zip(position, cubes.shape[1:])):
            raise ValueError("A target position lies outside the response cube.")
        matrix[target_index] = cubes[(slice(None), *position)] / float(rms)
    off_diagonal = matrix.copy()
    np.fill_diagonal(off_diagonal, 0.0)
    totals = np.sum(np.abs(off_diagonal), axis=1)
    maxima = np.max(np.abs(off_diagonal), axis=1)
    return LineInjectionCrosstalk(
        matrix_snr=matrix,
        off_diagonal_sum_snr=totals,
        max_off_diagonal_snr=maxima,
        safe=totals <= float(max_total_crosstalk_snr),
    )


def line_completeness(
    injected_catalogue,
    detection_catalogue,
    *,
    max_spatial_distance_pixels,
    max_channel_distance,
):
    """Match artificial lines to detections and return catalogue completeness.

    Rows must provide ``x``, ``y``, and ``channel`` values and may provide an
    ``spw`` value.  When both rows provide ``spw``, only equal SPWs can match.
    Candidate pairs are assigned in increasing normalized 3-D distance, and
    each injected source and detection can be used at most once.
    """
    inputs = list(injected_catalogue)
    detections = list(detection_catalogue)
    if not inputs:
        raise ValueError("injected_catalogue must contain at least one source.")
    if max_spatial_distance_pixels <= 0:
        raise ValueError("max_spatial_distance_pixels must be positive.")
    if max_channel_distance <= 0:
        raise ValueError("max_channel_distance must be positive.")

    def value(row, key):
        if isinstance(row, dict):
            return row[key]
        return getattr(row, key)

    def optional_value(row, key):
        if isinstance(row, dict):
            return row.get(key)
        return getattr(row, key, None)

    pairs = []
    for input_index, injected in enumerate(inputs):
        ix, iy = float(value(injected, "x")), float(value(injected, "y"))
        ichannel = float(value(injected, "channel"))
        ispw = optional_value(injected, "spw")
        for detection_index, detected in enumerate(detections):
            dspw = optional_value(detected, "spw")
            if ispw is not None and dspw is not None and ispw != dspw:
                continue
            spatial = math.hypot(
                float(value(detected, "x")) - ix,
                float(value(detected, "y")) - iy,
            )
            channel = abs(float(value(detected, "channel")) - ichannel)
            if (
                spatial <= max_spatial_distance_pixels
                and channel <= max_channel_distance
            ):
                normalized = math.hypot(
                    spatial / max_spatial_distance_pixels,
                    channel / max_channel_distance,
                )
                pairs.append(
                    (normalized, spatial, channel, input_index, detection_index)
                )

    matched_inputs = set()
    matched_detections = set()
    matches = []
    for normalized, spatial, channel, input_index, detection_index in sorted(pairs):
        if input_index in matched_inputs or detection_index in matched_detections:
            continue
        matched_inputs.add(input_index)
        matched_detections.add(detection_index)
        matches.append(
            LineSourceMatch(
                input_index=input_index,
                detection_index=detection_index,
                spatial_distance_pixels=float(spatial),
                channel_distance=float(channel),
                normalized_distance=float(normalized),
            )
        )
    recovered = tuple(index in matched_inputs for index in range(len(inputs)))
    return LineCompletenessResult(
        completeness=len(matches) / len(inputs),
        n_injected=len(inputs),
        n_recovered=len(matches),
        matches=tuple(sorted(matches, key=lambda row: row.input_index)),
        recovered=recovered,
    )


def flux_recovery_ratio(input_flux, recovered_flux):
    """
    Return the recovered/input flux ratio.

    Parameters
    ----------
    input_flux : float
        Input flux density.
    recovered_flux : float
        Recovered flux density.

    Returns
    -------
    float
        recovered_flux / input_flux
    """

    if input_flux == 0:
        raise ValueError("input_flux must be non-zero.")

    return recovered_flux / input_flux


def flux_fidelity(input_flux, recovered_flux):
    """Return the recovered/input flux ratio (legacy function name).

    ``flux_fidelity`` is retained for compatibility with existing DICE
    scripts. New code should use :func:`flux_recovery_ratio`, whose name
    distinguishes this measurement from catalogue fidelity/reliability.
    """
    return flux_recovery_ratio(input_flux, recovered_flux)


def positional_offset_arcsec(input_position, recovered_position):
    """
    Return the positional offset between input and recovered positions.

    Parameters
    ----------
    input_position : tuple of float
        Input position as (x_arcsec, y_arcsec).
    recovered_position : tuple of float
        Recovered position as (x_arcsec, y_arcsec).

    Returns
    -------
    float
        Positional offset in arcsec.
    """

    dx = recovered_position[0] - input_position[0]
    dy = recovered_position[1] - input_position[1]

    return math.sqrt(dx**2 + dy**2)


def is_recovered(input_position, recovered_position, max_offset_arcsec):
    """
    Check whether a recovered source is close enough to an injected source.

    Parameters
    ----------
    input_position : tuple of float
        Input position as (x_arcsec, y_arcsec).
    recovered_position : tuple of float
        Recovered position as (x_arcsec, y_arcsec).
    max_offset_arcsec : float
        Maximum allowed positional offset in arcsec.

    Returns
    -------
    bool
        True if the source is counted as recovered.
    """

    offset = positional_offset_arcsec(input_position, recovered_position)

    return offset <= max_offset_arcsec


def completeness(n_injected, n_recovered):
    """
    Return the recovery fraction.

    Parameters
    ----------
    n_injected : int
        Number of injected sources.
    n_recovered : int
        Number of recovered sources.

    Returns
    -------
    float
        n_recovered / n_injected
    """

    if n_injected == 0:
        raise ValueError("n_injected must be non-zero.")

    return n_recovered / n_injected


def match_catalogue_to_positions(input_positions, catalogue, max_distance_pixels):
    """Greedily match input positions to unique catalogue peak positions.

    Candidate pairs are considered in order of increasing distance, ensuring
    that each input and each detection is used at most once.
    """
    if max_distance_pixels < 0:
        raise ValueError("max_distance_pixels must be non-negative.")

    candidates = []
    for input_index, input_position in enumerate(input_positions):
        for detection_index, source in enumerate(catalogue):
            distance = math.hypot(
                float(source.peak_index[0]) - float(input_position[0]),
                float(source.peak_index[1]) - float(input_position[1]),
            )
            if distance <= max_distance_pixels:
                candidates.append((distance, input_index, detection_index))

    matched_inputs = set()
    matched_detections = set()
    matches = []
    for distance, input_index, detection_index in sorted(candidates):
        if input_index in matched_inputs or detection_index in matched_detections:
            continue
        matched_inputs.add(input_index)
        matched_detections.add(detection_index)
        matches.append(
            SourceMatch(
                input_index=input_index,
                detection_index=detection_index,
                distance_pixels=float(distance),
            )
        )
    return sorted(matches, key=lambda match: match.input_index)
