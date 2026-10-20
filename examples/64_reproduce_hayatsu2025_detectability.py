"""Reproduce the Hayatsu et al. (2025) line-completeness transition.

This CASA-independent script uses a Gaussian-replaced ADF22 dirty cube as the
noise template, injects Gaussian line profiles with the public DICE image/cube
API, applies the published matched boxcar widths, and compares the recovered
completeness with Tables 2 and 3 of Hayatsu et al. (2025).

The paper's contamination curve used 100 blank mock mosaics. Those mosaics are
not distributed with DICE, so this script reconstructs the published Q fits
but does not claim to regenerate their underlying blank-mosaic clump counts.

The default FITS file can be exported from the locally archived CASA image:

    casa --nogui --nologger --log2term -c \
      "from casatasks import exportfits; exportfits(\
imagename='data/hayatsu_reproduction/adf22_mockdata_comparison/\
Field10_Spw0_rep.ms.dirtycube.image', \
fitsimage='data/hayatsu_reproduction/adf22_mockdata_comparison/\
Field10_Spw0_rep.ms.dirtycube.image.fits', overwrite=True, dropdeg=True)"
"""

import argparse
import csv
import json
from pathlib import Path

from astropy.io import fits
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.optimize import curve_fit
from scipy.special import erf

from bliss_alma.extraction import extract_sources
from bliss_alma.injection import inject_source_image


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CUBE = (
    ROOT / "data/hayatsu_reproduction/adf22_mockdata_comparison"
    / "Field10_Spw0_rep.ms.dirtycube.image.fits"
)
DEFAULT_OBSERVED_CUBE = (
    ROOT / "data/hayatsu_reproduction/adf22_mockdata_comparison"
    / "Field10_Spw0.ms.dirtycube.image.fits"
)
DEFAULT_OUTPUT = ROOT / "data/hayatsu_reproduction/hayatsu2025_reproduction"

WIDTHS = np.arange(2, 24, 2)
INPUT_SNRS = np.arange(2.0, 7.01, 0.5)

# Hayatsu et al. (2025), Table 2. The printed error-function signs in Eqs. 4
# and 5 are opposite to the directions of the tabulated and plotted curves.
PUBLISHED_FIT = {
    "a_q": np.array([6.50, 6.42, 6.50, 6.50, 6.50, 6.50,
                     6.18, 5.87, 5.67, 5.54, 4.94]),
    "b_q": np.array([0.05, 0.03, 0.07, 0.08, 0.07, 0.08,
                     0.59, 0.12, 0.63, 0.53, 2.43]),
    "a_c": np.array([4.71, 3.61, 3.03, 2.96, 2.40, 2.15,
                     2.16, 1.98, 1.00, 2.35, 1.71]),
    "b_c": np.array([1.83, 2.05, 1.06, 1.19, 0.75, 0.74,
                     0.61, 0.22, 0.22, 0.67, 0.77]),
}

# Published differential completeness from Table 3. NaN denotes a cell shown
# as a dash because it was assigned to the unreliable category.
PUBLISHED_COMPLETENESS = np.array([
    [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
     0.65, 0.78, 0.88, 0.94, 0.97],
    [np.nan, np.nan, np.nan, np.nan,
     0.68, 0.80, 0.88, 0.93, 0.97, 0.99, 0.99],
    [np.nan, np.nan, np.nan,
     0.83, 0.95, 0.99, 1.00, 1.00, 1.00, 1.00, 1.00],
    [np.nan, np.nan,
     0.63, 0.82, 0.94, 0.98, 1.00, 1.00, 1.00, 1.00, 1.00],
    [np.nan, 0.75, 0.95, 0.99, 1.00, 1.00,
     1.00, 1.00, 1.00, 1.00, 1.00],
    [np.nan, 0.87, 0.98, 1.00, 1.00, 1.00,
     1.00, 1.00, 1.00, 1.00, 1.00],
    [np.nan, 0.92, 0.99, 1.00, 1.00, 1.00,
     1.00, 1.00, 1.00, 1.00, 1.00],
    [0.96, 1.00, 1.00, 1.00, 1.00, 1.00,
     1.00, 1.00, 1.00, 1.00, 1.00],
    [0.97, 1.00, 1.00, 1.00, 1.00, 1.00,
     1.00, 1.00, 1.00, 1.00, 1.00],
    [0.83, 0.98, 1.00, 1.00, 1.00, 1.00,
     1.00, 1.00, 1.00, 1.00, 1.00],
    [0.89, 0.98, 1.00, 1.00, 1.00, 1.00,
     1.00, 1.00, 1.00, 1.00, 1.00],
])


def robust_sigma(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    median = np.median(values)
    return float(1.4826 * np.median(np.abs(values - median)))


def increasing_erf(x, midpoint, width):
    return 0.5 * (1.0 + erf((np.asarray(x) - midpoint) / width))


def decreasing_erf(x, midpoint, width):
    return 0.5 * (1.0 - erf((np.asarray(x) - midpoint) / width))


def load_normalized_cube(path):
    with fits.open(path, memmap=True) as hdul:
        cube = np.asarray(hdul[0].data, dtype=float)
    cube = np.squeeze(cube)
    if cube.ndim != 3:
        raise ValueError(f"Expected a 3-D cube, received shape {cube.shape}")
    channel_rms = np.array([
        robust_sigma(channel) if np.any(np.isfinite(channel)) else np.nan
        for channel in cube
    ])
    normalized = np.full(cube.shape, np.nan, dtype=float)
    valid_rms = np.isfinite(channel_rms) & (channel_rms > 0)
    normalized[valid_rms] = (
        cube[valid_rms] / channel_rms[valid_rms, None, None]
    )
    return normalized, channel_rms


def draw_noise_spectra(cube, sample_count, half_window, rng):
    finite_channels = np.where(
        np.any(np.isfinite(cube), axis=(1, 2))
    )[0]
    first_channel = int(finite_channels[0])
    last_channel = int(finite_channels[-1])
    usable_channels = np.arange(
        first_channel + half_window,
        last_channel - half_window + 1,
    )
    if usable_channels.size == 0:
        raise ValueError("The cube has too few valid channels for this window")
    spatially_valid = np.all(
        np.isfinite(cube[first_channel:last_channel + 1]), axis=0
    )
    positions = np.argwhere(spatially_valid)
    if positions.size == 0:
        raise ValueError("No spatial pixel is valid across the selected channels")

    selected = positions[rng.integers(0, len(positions), sample_count)]
    centres = rng.choice(usable_channels, size=sample_count, replace=True)
    offsets = np.arange(-half_window, half_window + 1)
    spectra = cube[
        centres[:, None] + offsets[None, :],
        selected[:, 0, None],
        selected[:, 1, None],
    ]
    return spectra.T[:, :, None], selected, centres


def measure_completeness(
    blank_spectra,
    width,
    input_snr,
    threshold,
    group_count,
):
    channel_count, sample_count, _ = blank_spectra.shape
    centre = channel_count // 2
    sigma_channels = width / 2.354820045
    channel_offsets = np.arange(channel_count) - centre
    profile = np.exp(
        -0.5 * (channel_offsets / sigma_channels) ** 2
    )[:, None, None]
    model = input_snr * np.broadcast_to(profile, blank_spectra.shape)
    injected = inject_source_image(blank_spectra, model)

    blank_smoothed = uniform_filter1d(
        blank_spectra, size=width, axis=0, mode="nearest"
    )
    injected_smoothed = uniform_filter1d(
        injected, size=width, axis=0, mode="nearest"
    )
    exclusion = max(width, 3)
    noise_values = np.concatenate([
        blank_smoothed[:centre - exclusion].ravel(),
        blank_smoothed[centre + exclusion + 1:].ravel(),
    ])
    smoothed_rms = robust_sigma(noise_values)
    search_slice = slice(
        centre - width // 2,
        centre + width // 2 + 1,
    )
    detected = (
        np.max(injected_smoothed[search_slice, :, 0], axis=0)
        / smoothed_rms >= threshold
    )

    groups = np.array_split(detected, group_count)
    group_completeness = np.array([
        np.mean(group) for group in groups if len(group)
    ])
    return {
        "completeness": float(np.mean(detected)),
        "range_low": float(np.percentile(group_completeness, 16)),
        "range_high": float(np.percentile(group_completeness, 84)),
        "smoothed_rms": smoothed_rms,
    }


def fit_completeness(values):
    initial_midpoint = float(np.interp(
        0.5, values, INPUT_SNRS, left=INPUT_SNRS[0], right=INPUT_SNRS[-1]
    ))
    try:
        parameters, _ = curve_fit(
            increasing_erf,
            INPUT_SNRS,
            values,
            p0=(initial_midpoint, 1.0),
            bounds=([0.0, 0.05], [10.0, 10.0]),
            maxfev=20000,
        )
    except (RuntimeError, ValueError):
        return np.nan, np.nan
    return float(parameters[0]), float(parameters[1])


def normalized_boxcar_cube(cube, width):
    finite = np.isfinite(cube)
    filled = np.where(finite, cube, 0.0)
    smoothed = uniform_filter1d(
        filled, size=width, axis=0, mode="constant", cval=0.0
    )
    valid_fraction = uniform_filter1d(
        finite.astype(float), size=width, axis=0,
        mode="constant", cval=0.0,
    )
    valid = valid_fraction >= 1.0 - 1e-12
    normalized = np.zeros(cube.shape, dtype=np.float32)
    for channel in range(cube.shape[0]):
        channel_valid = valid[channel]
        if not np.any(channel_valid):
            continue
        noise = robust_sigma(smoothed[channel][channel_valid])
        normalized[channel][channel_valid] = (
            smoothed[channel][channel_valid] / noise
        )
    return normalized, ~valid


def single_field_contamination_diagnostic(
    observed_cube,
    blank_cube,
    output_path,
):
    rows = []
    for width in WIDTHS:
        observed_smoothed, observed_mask = normalized_boxcar_cube(
            observed_cube, int(width)
        )
        blank_smoothed, blank_mask = normalized_boxcar_cube(
            blank_cube, int(width)
        )
        for threshold in INPUT_SNRS:
            counts = {}
            for label, data, mask in (
                ("observed", observed_smoothed, observed_mask),
                ("blank", blank_smoothed, blank_mask),
            ):
                result = extract_sources(
                    data,
                    rms=1.0,
                    seed_threshold=float(threshold),
                    flood_threshold=float(threshold),
                    min_pixels=1,
                    connectivity=1,
                    polarity="positive",
                    detection_mask=mask,
                    exclude_border=True,
                )
                counts[label] = len(result.catalogue)
            q_value = (
                min(1.0, counts["blank"] / counts["observed"])
                if counts["observed"] else np.nan
            )
            rows.append({
                "width_channels": int(width),
                "measured_snr_threshold": float(threshold),
                "observed_field10_clumps": counts["observed"],
                "blank_replaced_field10_clumps": counts["blank"],
                "single_field_q": q_value,
            })
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cube-fits", type=Path, default=DEFAULT_CUBE)
    parser.add_argument(
        "--observed-cube-fits", type=Path, default=DEFAULT_OBSERVED_CUBE,
        help=(
            "Optional observed Field10 FITS cube for a noisy single-field "
            "contamination diagnostic. It does not replace the published "
            "100-mosaic calculation."
        ),
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--samples", type=int, default=10000)
    parser.add_argument("--groups", type=int, default=100)
    parser.add_argument(
        "--threshold", type=float, default=4.0,
        help=(
            "Smoothed-S/N recovery threshold. The 4.0 default best matches "
            "the historical clumpfind recovery transition; use 4.5 to test "
            "the Hayatsu et al. (2017) catalogue threshold."
        ),
    )
    parser.add_argument("--seed", type=int, default=20251002)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    cube, channel_rms = load_normalized_cube(args.cube_fits)
    rng = np.random.default_rng(args.seed)
    half_window = 36
    blank_spectra, positions, centres = draw_noise_spectra(
        cube, args.samples, half_window, rng
    )

    rows = []
    fitted = {}
    for width_index, width in enumerate(WIDTHS):
        width_values = []
        for snr_index, input_snr in enumerate(INPUT_SNRS):
            measured = measure_completeness(
                blank_spectra,
                int(width),
                float(input_snr),
                args.threshold,
                args.groups,
            )
            width_values.append(measured["completeness"])
            published = PUBLISHED_COMPLETENESS[width_index, snr_index]
            rows.append({
                "width_channels": int(width),
                "input_peak_snr": float(input_snr),
                "dice_completeness": measured["completeness"],
                "dice_range_low": measured["range_low"],
                "dice_range_high": measured["range_high"],
                "published_completeness": (
                    "" if not np.isfinite(published) else float(published)
                ),
                "dice_minus_published": (
                    "" if not np.isfinite(published)
                    else measured["completeness"] - float(published)
                ),
                "smoothed_rms": measured["smoothed_rms"],
            })
        fitted[int(width)] = fit_completeness(np.asarray(width_values))

    csv_path = args.output_dir / "hayatsu2025_dice_completeness.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    fit_rows = []
    for index, width in enumerate(WIDTHS):
        dice_midpoint, dice_width = fitted[int(width)]
        fit_rows.append({
            "width_channels": int(width),
            "published_completeness_midpoint_a2": PUBLISHED_FIT["a_c"][index],
            "dice_completeness_midpoint": dice_midpoint,
            "midpoint_difference": (
                dice_midpoint - PUBLISHED_FIT["a_c"][index]
            ),
            "published_completeness_width_b2": PUBLISHED_FIT["b_c"][index],
            "dice_completeness_width": dice_width,
            "published_contamination_midpoint_a1":
                PUBLISHED_FIT["a_q"][index],
            "published_contamination_width_b1":
                PUBLISHED_FIT["b_q"][index],
        })
    fit_path = args.output_dir / "hayatsu2025_fit_comparison.csv"
    with fit_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fit_rows[0].keys())
        writer.writeheader()
        writer.writerows(fit_rows)

    metadata = {
        "input_cube": str(args.cube_fits),
        "samples": args.samples,
        "groups_for_16_84_range": args.groups,
        "threshold_smoothed_snr": args.threshold,
        "seed": args.seed,
        "channel_rms_median_jy_per_beam": float(np.nanmedian(channel_rms)),
        "published_reference": "Hayatsu et al. 2025, Tables 2 and 3",
        "important_limitation": (
            "The 100 blank mock mosaics used for the published contamination "
            "calculation are unavailable locally; Q fits are reconstructed "
            "from the paper and not independently regenerated."
        ),
    }
    with (args.output_dir / "hayatsu2025_reproduction_metadata.json").open(
        "w"
    ) as handle:
        json.dump(metadata, handle, indent=2)

    if args.observed_cube_fits.exists():
        observed_cube, _ = load_normalized_cube(args.observed_cube_fits)
        single_field_contamination_diagnostic(
            observed_cube,
            cube,
            args.output_dir / "hayatsu2025_field10_contamination_diagnostic.csv",
        )

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.4))
    selected_widths = (2, 6, 14, 22)
    colors = plt.cm.viridis(np.linspace(0.12, 0.88, len(selected_widths)))
    for width, color in zip(selected_widths, colors):
        index = int(np.where(WIDTHS == width)[0][0])
        selected_rows = [
            row for row in rows if row["width_channels"] == width
        ]
        dice = np.array([row["dice_completeness"] for row in selected_rows])
        low = np.array([row["dice_range_low"] for row in selected_rows])
        high = np.array([row["dice_range_high"] for row in selected_rows])
        axes[0, 0].plot(
            INPUT_SNRS, dice, "o-", color=color, label=f"{width} ch"
        )
        axes[0, 0].fill_between(
            INPUT_SNRS, low, high, color=color, alpha=0.14
        )
        published = PUBLISHED_COMPLETENESS[index]
        finite = np.isfinite(published)
        axes[0, 0].plot(
            INPUT_SNRS[finite], published[finite], "s",
            markerfacecolor="none", markeredgecolor=color,
        )

    axes[0, 0].set(
        xlabel="Input peak S/N",
        ylabel="Differential completeness",
        title="DICE reproduction (filled) vs published Table 3 (open)",
        xlim=(1.9, 7.1), ylim=(-0.03, 1.03),
    )
    axes[0, 0].legend(frameon=False, ncol=2)

    dice_midpoints = np.array([
        fitted[int(width)][0] for width in WIDTHS
    ])
    axes[0, 1].plot(
        WIDTHS, PUBLISHED_FIT["a_c"], "s-", color="#c83e4d",
        label="Published $a_2$",
    )
    axes[0, 1].plot(
        WIDTHS, dice_midpoints, "o-", color="#377eb8",
        label="DICE midpoint",
    )
    axes[0, 1].axhline(args.threshold, color="0.5", linestyle=":")
    axes[0, 1].set(
        xlabel="Gaussian FWHM / boxcar width (channels)",
        ylabel="Input S/N at 50% completeness",
        title="Completeness transition",
    )
    axes[0, 1].legend(frameon=False)

    axes[1, 0].axhline(0.0, color="0.5", linewidth=0.8)
    axes[1, 0].plot(
        WIDTHS, dice_midpoints - PUBLISHED_FIT["a_c"],
        "o-", color="#4daf4a",
    )
    axes[1, 0].set(
        xlabel="Gaussian FWHM / boxcar width (channels)",
        ylabel="DICE - published midpoint (S/N)",
        title="Difference in fitted transition",
    )

    x_dense = np.linspace(2.0, 7.0, 400)
    for width, color in zip(selected_widths, colors):
        index = int(np.where(WIDTHS == width)[0][0])
        axes[1, 1].plot(
            decreasing_erf(
                x_dense,
                PUBLISHED_FIT["a_q"][index],
                PUBLISHED_FIT["b_q"][index],
            ),
            increasing_erf(
                x_dense,
                PUBLISHED_FIT["a_c"][index],
                PUBLISHED_FIT["b_c"][index],
            ),
            color=color,
            label=f"{width} ch",
        )
    axes[1, 1].set(
        xlabel="Published cumulative contamination $Q$",
        ylabel="Published completeness $C$",
        title="Reconstructed published detectability trajectories",
        xlim=(-0.03, 1.03), ylim=(-0.03, 1.03),
    )
    axes[1, 1].invert_xaxis()
    axes[1, 1].legend(frameon=False)

    for axis in axes.ravel():
        axis.grid(alpha=0.2)
    fig.suptitle(
        "Hayatsu et al. (2025) detectability reproduction\n"
        "Gaussian-replaced Field10 template; matched boxcars; "
        f"{args.samples:,} samples"
    )
    fig.tight_layout()
    figure_path = args.output_dir / "hayatsu2025_detectability_reproduction"
    fig.savefig(figure_path.with_suffix(".png"), dpi=220)
    fig.savefig(figure_path.with_suffix(".pdf"))

    print("Wrote", csv_path)
    print("Wrote", fit_path)
    print("Wrote", figure_path.with_suffix(".png"))


if __name__ == "__main__":
    main()
