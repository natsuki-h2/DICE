"""Generate independent blank ADF22 dirty cubes and count noise clumps.

This CASA 6 script is a storage-bounded pilot for the Hayatsu et al. (2025)
contamination experiment.  It reuses one working copy of the archived
Gaussian-replaced Field10 measurement set, replaces its visibilities with an
independent Gaussian realization, makes a natural-weighted dirty cube, and
records positive and negative clumps after matched boxcar smoothing.

Only catalogues and per-realization summaries are retained.  CASA image
products and the working measurement set are removed after each successful
measurement, so the run can be extended to 100 realizations without keeping
100 measurement sets or cubes.

Run a timing pilot:

    PYTHONPATH=src casa --nogui --nologger --log2term -c \
      scripts/adf22/generate_hayatsu2025_blank_realizations.py \
      --realizations 5
"""

import argparse
import csv
import json
from pathlib import Path
import shutil
import sys
import time

import numpy as np
from scipy import ndimage
from scipy.ndimage import uniform_filter1d

from casatasks import tclean
from casatools import image as image_tool
from casatools import table as table_tool


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

DEFAULT_TEMPLATE = (
    ROOT / "data/hayatsu_reproduction/adf22_mockdata_comparison"
    / "Field10_Spw0_rep.ms"
)
DEFAULT_OUTPUT = (
    ROOT / "data/hayatsu_reproduction/hayatsu2025_blank_realizations"
)
WIDTHS = tuple(range(2, 24, 2))
THRESHOLDS = tuple(np.arange(2.0, 7.01, 0.5))


def robust_sigma(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    median = np.median(values)
    return float(1.4826 * np.median(np.abs(values - median)))


def remove_path(path):
    path = Path(path)
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def remove_image_products(prefix):
    for suffix in (
        ".image", ".pb", ".psf", ".model", ".residual", ".sumwt",
        ".mask", ".weight", ".image.pbcor",
    ):
        remove_path(str(prefix) + suffix)


def template_statistics(template):
    tb = table_tool()
    tb.open(str(template))
    column = "CORRECTED_DATA" if "CORRECTED_DATA" in tb.colnames() else "DATA"
    data = tb.getcol(column)
    flags = tb.getcol("FLAG")
    tb.close()

    scales = np.zeros(data.shape[:2] + (2,), dtype=float)
    for pol in range(data.shape[0]):
        for channel in range(data.shape[1]):
            valid = ~flags[pol, channel]
            scales[pol, channel, 0] = np.std(data[pol, channel, valid].real)
            scales[pol, channel, 1] = np.std(data[pol, channel, valid].imag)
    return column, flags, scales


def independent_gaussian_data(shape, flags, scales, rng):
    data = np.zeros(shape, dtype=np.complex128)
    for pol in range(shape[0]):
        real = rng.normal(size=(shape[1], shape[2]))
        imag = rng.normal(size=(shape[1], shape[2]))

        # Approximate the archived Hanning-smoothed channel covariance, then
        # renormalize each channel to its measured visibility dispersion.
        real = uniform_filter1d(real, size=3, axis=0, mode="nearest")
        imag = uniform_filter1d(imag, size=3, axis=0, mode="nearest")
        for channel in range(shape[1]):
            valid = ~flags[pol, channel]
            real_std = np.std(real[channel, valid])
            imag_std = np.std(imag[channel, valid])
            if real_std > 0:
                real[channel] *= scales[pol, channel, 0] / real_std
            if imag_std > 0:
                imag[channel] *= scales[pol, channel, 1] / imag_std
        data[pol] = real + 1j * imag
    return data


def read_cube(path):
    ia = image_tool()
    ia.open(str(path))
    raw = np.squeeze(ia.getchunk()).astype(np.float32)
    ia.close()
    if raw.ndim != 3:
        raise ValueError("Expected a 3-D CASA cube, received {}".format(raw.shape))
    # CASA image-tool order is (x, y, channel); DICE uses (channel, y, x).
    return np.moveaxis(raw, -1, 0)


def normalized_boxcar_cube(cube, width):
    finite = np.isfinite(cube) & (cube != 0)
    filled = np.where(finite, cube, 0.0)
    smoothed = uniform_filter1d(
        filled, size=width, axis=0, mode="constant", cval=0.0
    )
    valid_fraction = uniform_filter1d(
        finite.astype(np.float32), size=width, axis=0,
        mode="constant", cval=0.0,
    )
    valid = valid_fraction >= 1.0 - 1e-6
    normalized = np.zeros(cube.shape, dtype=np.float32)
    channel_rms = np.full(cube.shape[0], np.nan)
    for channel in range(cube.shape[0]):
        use = valid[channel]
        if not np.any(use):
            continue
        channel_rms[channel] = robust_sigma(smoothed[channel][use])
        if channel_rms[channel] > 0:
            normalized[channel][use] = (
                smoothed[channel][use] / channel_rms[channel]
            )
    return normalized, ~valid, channel_rms


def count_clumps(cube, realization, seed):
    rows = []
    for width in WIDTHS:
        normalized, mask, channel_rms = normalized_boxcar_cube(cube, width)
        for threshold in THRESHOLDS:
            for polarity in ("positive", "negative"):
                detection = normalized if polarity == "positive" else -normalized
                labels, n_labels = ndimage.label(
                    (detection >= threshold) & ~mask,
                    structure=ndimage.generate_binary_structure(3, 1),
                )
                border_labels = np.unique(np.concatenate([
                    labels[0].ravel(), labels[-1].ravel(),
                    labels[:, 0].ravel(), labels[:, -1].ravel(),
                    labels[:, :, 0].ravel(), labels[:, :, -1].ravel(),
                ]))
                border_count = int(np.sum(border_labels != 0))
                rows.append({
                    "realization": realization,
                    "noise_seed": seed,
                    "width_channels": width,
                    "snr_threshold": float(threshold),
                    "polarity": polarity,
                    "n_clumps": int(n_labels) - border_count,
                    "median_smoothed_rms_jy_beam":
                        float(np.nanmedian(channel_rms)),
                })
    return rows


def write_rows(path, rows):
    exists = path.exists()
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        if not exists:
            writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--template-ms", type=Path, default=DEFAULT_TEMPLATE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--realizations", type=int, default=5)
    parser.add_argument("--start-index", type=int, default=0)
    parser.add_argument(
        "--indices",
        help="Comma-separated realization IDs; overrides start/count.",
    )
    parser.add_argument("--seed-base", type=int, default=202510000)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    work_ms = args.output_dir / "working_blank.ms"
    image_prefix = args.output_dir / "working_blank_dirtycube"
    catalogue_path = args.output_dir / "blank_clump_counts.csv"
    timing_path = args.output_dir / "blank_realization_timings.csv"

    column, flags, scales = template_statistics(args.template_ms)
    tb = table_tool()
    timing_rows = []
    try:
        remove_path(work_ms)
        shutil.copytree(args.template_ms, work_ms)
        tb.open(str(work_ms))
        shape = tb.getcol(column).shape
        tb.close()

        realization_ids = (
            [int(value) for value in args.indices.split(",")]
            if args.indices else
            list(range(args.start_index, args.start_index + args.realizations))
        )
        for realization in realization_ids:
            started = time.time()
            seed = args.seed_base + realization
            rng = np.random.default_rng(seed)
            noise = independent_gaussian_data(shape, flags, scales, rng)

            tb.open(str(work_ms), nomodify=False)
            tb.putcol(column, noise)
            if column != "DATA" and "DATA" in tb.colnames():
                tb.putcol("DATA", noise)
            tb.flush()
            tb.close()
            del noise

            remove_image_products(image_prefix)
            image_started = time.time()
            tclean(
                vis=str(work_ms),
                imagename=str(image_prefix),
                spw="0",
                specmode="cube",
                weighting="natural",
                datacolumn="corrected" if column == "CORRECTED_DATA" else "data",
                width="15.625MHz",
                niter=0,
                cell="0.1arcsec",
                imsize=256,
                pblimit=0.8,
            )
            image_seconds = time.time() - image_started
            cube = read_cube(str(image_prefix) + ".image")
            rows = count_clumps(cube, realization, seed)
            write_rows(catalogue_path, rows)
            elapsed = time.time() - started
            timing_rows.append({
                "realization": realization,
                "noise_seed": seed,
                "image_seconds": image_seconds,
                "total_seconds": elapsed,
                "cube_shape": "x".join(str(value) for value in cube.shape),
            })
            write_rows(timing_path, [timing_rows[-1]])
            remove_image_products(image_prefix)
            print(
                "realization={} seed={} image_s={:.1f} total_s={:.1f}".format(
                    realization, seed, image_seconds, elapsed
                )
            )
    finally:
        try:
            tb.close()
        except Exception:
            pass
        remove_image_products(image_prefix)
        remove_path(work_ms)

    metadata = {
        "template_ms": str(args.template_ms),
        "visibility_column": column,
        "realizations_requested": args.realizations,
        "realization_indices": realization_ids,
        "start_index": args.start_index,
        "seed_base": args.seed_base,
        "width_channels": list(WIDTHS),
        "snr_thresholds": list(THRESHOLDS),
        "imaging": {
            "specmode": "cube",
            "weighting": "natural",
            "niter": 0,
            "cell": "0.1arcsec",
            "imsize": 256,
            "channel_width": "15.625MHz",
            "pblimit": 0.8,
        },
        "scope": (
            "Independent Field10 blank realizations used as a timing and "
            "method pilot; not yet the complete 100-field ADF22 mosaic."
        ),
    }
    with (args.output_dir / "blank_realization_metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2)


if __name__ == "__main__":
    main()
