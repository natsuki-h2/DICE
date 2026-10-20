"""Configurable source extraction for DICE.

The default finder uses SciPy connected-component labelling with separate
seed and flood thresholds.  It supports both 2D continuum images and 3D
spectral cubes, and can be run identically on positive and negative data.
"""

from collections import namedtuple

import numpy as np
from scipy import ndimage


Clump = namedtuple(
    "Clump",
    [
        "label",
        "n_pixels",
        "peak_value",
        "peak_index",
        "centroid_index",
        "summed_value",
    ],
)

ExtractionResult = namedtuple(
    "ExtractionResult",
    [
        "labels",
        "catalogue",
        "rms",
        "seed_threshold",
        "flood_threshold",
        "polarity",
    ],
)


def estimate_rms(data, method="mad", mask=None):
    """Estimate image/cube noise using a selectable method.

    Parameters
    ----------
    data : array-like
        Image or cube values.
    method : {"mad", "std", "rms"}, optional
        ``mad`` uses 1.4826 times the median absolute deviation, ``std`` uses
        the standard deviation about the mean, and ``rms`` uses the square
        root of the mean squared value.
    mask : array-like of bool, optional
        Pixels marked ``True`` are excluded from the estimate.
    """
    array = np.asarray(data, dtype=np.float64)
    valid = np.isfinite(array)
    if mask is not None:
        mask_array = np.asarray(mask, dtype=bool)
        if mask_array.shape != array.shape:
            raise ValueError("mask must have the same shape as data.")
        valid &= ~mask_array
    values = array[valid]
    if values.size == 0:
        raise ValueError("No finite, unmasked values are available.")

    method = method.lower()
    if method == "mad":
        median = np.median(values)
        noise = 1.4826 * np.median(np.abs(values - median))
    elif method == "std":
        noise = np.std(values)
    elif method == "rms":
        noise = np.sqrt(np.mean(values**2))
    else:
        raise ValueError("method must be one of: mad, std, rms.")

    noise = float(noise)
    if not np.isfinite(noise) or noise <= 0:
        raise ValueError("The estimated RMS must be finite and positive.")
    return noise


def catalog_from_labels(data, labels, background_label=0):
    """Convert an integer assignment array into a small source catalogue."""
    array = np.asarray(data, dtype=np.float64)
    label_array = np.asarray(labels)
    if array.shape != label_array.shape:
        raise ValueError("data and labels must have the same shape.")
    if not np.issubdtype(label_array.dtype, np.integer):
        raise TypeError("labels must have an integer dtype.")

    catalogue = []
    if background_label != 0:
        # ``ndimage.find_objects`` assumes zero is the background label.
        working_labels = label_array.copy()
        working_labels[label_array == background_label] = 0
    else:
        working_labels = label_array
    for label_value, label_slice in enumerate(ndimage.find_objects(working_labels), start=1):
        if label_slice is None:
            continue
        local_labels = working_labels[label_slice]
        source_mask = local_labels == label_value
        local_indices = np.argwhere(source_mask)
        offsets = np.asarray([axis.start for axis in label_slice], dtype=int)
        indices = local_indices + offsets
        values = array[label_slice][source_mask]
        peak_row = int(np.argmax(values))
        peak_index = tuple(int(value) for value in indices[peak_row])

        positive_weights = np.clip(values, 0.0, None)
        if float(positive_weights.sum()) > 0.0:
            centroid = np.average(indices, axis=0, weights=positive_weights)
        else:
            centroid = indices.mean(axis=0)

        catalogue.append(
            Clump(
                label=label_value,
                n_pixels=int(source_mask.sum()),
                peak_value=float(values[peak_row]),
                peak_index=peak_index,
                centroid_index=tuple(float(value) for value in centroid),
                summed_value=float(values.sum()),
            )
        )
    return catalogue


def extract_sources(
    data,
    rms=None,
    rms_method="mad",
    rms_mask=None,
    seed_sigma=5.0,
    flood_sigma=3.0,
    seed_threshold=None,
    flood_threshold=None,
    min_pixels=5,
    connectivity=1,
    polarity="positive",
    detection_mask=None,
    exclude_border=False,
):
    """Detect sources using seed-and-flood connected components.

    A flood island is retained when it contains at least one seed pixel and
    contains at least ``min_pixels`` pixels. Thresholds may be specified as
    absolute values or as multiples of the selected RMS estimate.

    ``polarity="negative"`` applies the same algorithm to ``-data``. This is
    intended for symmetric blank-field contamination tests. Pixels marked
    ``True`` in ``detection_mask`` are excluded. If ``exclude_border`` is
    enabled, islands touching an array edge are rejected.
    """
    array = np.asarray(data, dtype=np.float64)
    if array.ndim not in (2, 3):
        raise ValueError("Source extraction requires a 2D image or 3D cube.")
    if not np.all(np.isfinite(array)):
        raise ValueError("Source extraction input contains NaN or infinity.")
    if min_pixels < 1:
        raise ValueError("min_pixels must be at least one.")
    if connectivity < 1 or connectivity > array.ndim:
        raise ValueError("connectivity must be between one and data.ndim.")
    if polarity not in ("positive", "negative"):
        raise ValueError("polarity must be 'positive' or 'negative'.")

    usable = np.ones(array.shape, dtype=bool)
    if detection_mask is not None:
        detection_mask = np.asarray(detection_mask, dtype=bool)
        if detection_mask.shape != array.shape:
            raise ValueError("detection_mask must have the same shape as data.")
        usable &= ~detection_mask

    if rms is None:
        rms = estimate_rms(array, method=rms_method, mask=rms_mask)
    rms = float(rms)
    if not np.isfinite(rms) or rms <= 0:
        raise ValueError("rms must be finite and positive.")

    if seed_threshold is None:
        seed_threshold = float(seed_sigma) * rms
    if flood_threshold is None:
        flood_threshold = float(flood_sigma) * rms
    seed_threshold = float(seed_threshold)
    flood_threshold = float(flood_threshold)
    if seed_threshold < flood_threshold:
        raise ValueError("seed_threshold must be >= flood_threshold.")

    detection_data = array if polarity == "positive" else -array
    structure = ndimage.generate_binary_structure(array.ndim, connectivity)
    flood_labels, n_flood = ndimage.label(
        (detection_data >= flood_threshold) & usable,
        structure=structure,
    )
    seed_mask = (detection_data >= seed_threshold) & usable

    labels = np.zeros(array.shape, dtype=np.int32)
    next_label = 1
    for flood_label, flood_slice in enumerate(ndimage.find_objects(flood_labels), start=1):
        if flood_slice is None:
            continue
        local_flood = flood_labels[flood_slice]
        island = local_flood == flood_label
        if int(island.sum()) < min_pixels:
            continue
        if not np.any(seed_mask[flood_slice] & island):
            continue
        if exclude_border:
            touches_lower = any(axis.start == 0 for axis in flood_slice)
            touches_upper = any(
                axis.stop == array.shape[index]
                for index, axis in enumerate(flood_slice)
            )
            if touches_lower or touches_upper:
                continue
        local_output = labels[flood_slice]
        local_output[island] = next_label
        next_label += 1

    catalogue_data = detection_data
    catalogue = catalog_from_labels(catalogue_data, labels)
    return ExtractionResult(
        labels=labels,
        catalogue=catalogue,
        rms=rms,
        seed_threshold=seed_threshold,
        flood_threshold=flood_threshold,
        polarity=polarity,
    )
