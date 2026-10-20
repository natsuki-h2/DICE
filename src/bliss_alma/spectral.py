"""Spectral-cube helpers, including emission and absorption searches."""

import numpy as np

from .extraction import extract_sources


def absorption_cube(cube, continuum, representation="fractional_depth"):
    """Convert a spectral cube to a positive absorption-strength cube.

    Parameters
    ----------
    cube : array-like, 3-D
        Observed intensity cube.
    continuum : array-like
        Positive continuum model, either 2-D (broadcast along axis 0) or the
        same 3-D shape as ``cube``.
    representation : {"residual_depth", "fractional_depth", "optical_depth"}
        ``residual_depth`` returns ``continuum - cube``;
        ``fractional_depth`` returns ``1 - cube/continuum``; and
        ``optical_depth`` returns ``-ln(cube/continuum)``.
    """
    data = np.asarray(cube, dtype=np.float64)
    cont = np.asarray(continuum, dtype=np.float64)
    if data.ndim != 3:
        raise ValueError("cube must be three-dimensional.")
    if cont.ndim == 2:
        if cont.shape != data.shape[1:]:
            raise ValueError("2-D continuum must match the cube spatial shape.")
        cont = cont[None, :, :]
    elif cont.shape != data.shape:
        raise ValueError("3-D continuum must have the same shape as cube.")
    if not np.all(np.isfinite(cont)) or np.any(cont <= 0):
        raise ValueError("continuum must be finite and strictly positive.")
    if representation == "residual_depth":
        return cont - data
    ratio = data / cont
    if representation == "fractional_depth":
        return 1.0 - ratio
    if representation == "optical_depth":
        if np.any(ratio <= 0):
            raise ValueError("optical depth requires cube/continuum > 0.")
        return -np.log(ratio)
    raise ValueError("unknown absorption representation.")


def extract_absorption_lines(cube, continuum, representation="fractional_depth", **kwargs):
    """Detect absorption features after explicit continuum normalization.

    The converted absorption cube is positive, so the standard positive
    extractor is used. Threshold/RMS values must be expressed in the selected
    representation.
    """
    depth = absorption_cube(cube, continuum, representation=representation)
    if "polarity" in kwargs and kwargs["polarity"] != "positive":
        raise ValueError("absorption depth is positive; polarity must be 'positive'.")
    kwargs["polarity"] = "positive"
    return extract_sources(depth, **kwargs)
