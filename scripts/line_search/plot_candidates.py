"""Create target-neutral candidate-review and line-search diagnostic plots."""

from argparse import ArgumentParser
import csv
import json
from pathlib import Path
import sys

from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from bliss_alma.linesearch import robust_sigma
from bliss_alma.provenance import write_provenance


ROOT = Path(__file__).resolve().parents[2]
POSITIVE_COLOR = "#c83e4d"
NEGATIVE_COLOR = "#377eb8"


def resolve(value):
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def spectral_axis_ghz(wcs, channel_count):
    """Return spectral coordinates without assuming a 3-D or 4-D header."""

    world = wcs.spectral.pixel_to_world_values(np.arange(channel_count))
    return np.asarray(world, dtype=float) / 1e9


def load_pb(path):
    with fits.open(path, memmap=True) as hdul:
        pb = np.squeeze(hdul[0].data)
    if pb.ndim == 3:
        pb = pb[pb.shape[0] // 2]
    if pb.ndim != 2:
        raise ValueError(f"Expected a 2-D PB image after squeeze; got {pb.shape}")
    return np.asarray(pb, dtype=float)


def sampled_pixel_snr(cube, mask, limit=500_000):
    """Sample per-channel standardized pixels for a memory-safe QA histogram."""

    per_channel = max(1, limit // cube.shape[0])
    samples = []
    for plane in cube:
        values = np.asarray(plane, dtype=float)[mask]
        values = values[np.isfinite(values) & (values != 0)]
        sigma = robust_sigma(values)
        if not values.size or not np.isfinite(sigma) or sigma <= 0:
            continue
        normalized = (values - np.median(values)) / sigma
        if normalized.size > per_channel:
            step = normalized.size / per_channel
            normalized = normalized[(np.arange(per_channel) * step).astype(int)]
        samples.append(normalized)
    return np.concatenate(samples) if samples else np.array([], dtype=float)


def candidate_review_pdf(cube, wcs, rows, output, radius):
    frequencies = spectral_axis_ghz(wcs, cube.shape[0])
    with PdfPages(output) as pdf:
        for row in rows:
            x, y = int(float(row["x"])), int(float(row["y"]))
            centre = float(row["channel"])
            width = int(float(row["width_channels"]))
            snr = float(row["snr"])
            first = max(0, int(np.floor(centre - (width - 1) / 2)))
            last = min(cube.shape[0], first + width)

            moment = 1e3 * np.nansum(np.asarray(cube[first:last], dtype=float), axis=0)
            spectrum = 1e3 * np.asarray(cube[:, y, x], dtype=float)
            x0, x1 = max(0, x - radius), min(cube.shape[2], x + radius + 1)
            y0, y1 = max(0, y - radius), min(cube.shape[1], y + radius + 1)
            cutout = moment[y0:y1, x0:x1]
            scale = robust_sigma(cutout)
            if not np.isfinite(scale) or scale <= 0:
                scale = np.nanstd(cutout)
            if not np.isfinite(scale) or scale <= 0:
                scale = 1.0

            fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3))
            image = axes[0].imshow(
                cutout, origin="lower", cmap="RdBu_r", vmin=-6 * scale, vmax=6 * scale
            )
            positive_levels = scale * np.array([2.0, 3.0, 4.0])
            axes[0].contour(cutout, levels=positive_levels, colors=POSITIVE_COLOR, linewidths=0.8)
            axes[0].contour(
                cutout,
                levels=-positive_levels[::-1],
                colors=NEGATIVE_COLOR,
                linestyles="dashed",
                linewidths=0.8,
            )
            axes[0].set(title=f"Moment-0 ({width} channels)", xticks=[], yticks=[])
            fig.colorbar(image, ax=axes[0], fraction=0.046, pad=0.04, label="Channel sum (mJy beam$^{-1}$)")

            axes[1].step(
                frequencies,
                spectrum,
                where="mid",
                color=POSITIVE_COLOR if snr > 0 else NEGATIVE_COLOR,
                linestyle="-" if snr > 0 else "--",
                linewidth=1.0,
            )
            if first < last:
                axes[1].axvspan(
                    frequencies[first],
                    frequencies[last - 1],
                    color=POSITIVE_COLOR if snr > 0 else NEGATIVE_COLOR,
                    alpha=0.14,
                )
            axes[1].axhline(0, color="0.5", linewidth=0.7)
            axes[1].set(
                xlabel="Observed frequency (GHz)",
                ylabel="Flux density (mJy beam$^{-1}$)",
                title="Spectrum at candidate pixel",
            )
            axes[1].grid(alpha=0.2)
            fig.suptitle(
                f"{row['candidate_id']}  S/N={snr:.2f}\n"
                f"RA={float(row['ra_deg']):.5f} deg, "
                f"Dec={float(row['dec_deg']):.5f} deg, "
                f"ν={float(row['frequency_ghz']):.5f} GHz",
                fontsize=11,
            )
            fig.tight_layout(rect=(0, 0, 1, 0.9))
            pdf.savefig(fig)
            plt.close(fig)


def sky_distribution(pb, wcs, rows, pb_min, pdf_path, png_path):
    pixel_scale = float(np.mean(proj_plane_pixel_scales(wcs.celestial))) * 3600.0
    ny, nx = pb.shape
    x_mid, y_mid = (nx - 1) / 2.0, (ny - 1) / 2.0
    extent = ((nx - 0.5 - x_mid) * pixel_scale, (-0.5 - x_mid) * pixel_scale,
              (-0.5 - y_mid) * pixel_scale, (ny - 0.5 - y_mid) * pixel_scale)

    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    image = ax.imshow(pb, origin="lower", extent=extent, cmap="Greys", vmin=0, vmax=1)
    ax.contour(pb, levels=[pb_min], colors="#e69f00", linestyles="--", linewidths=1.2,
               extent=extent, origin="lower")
    for row in rows:
        x, y, snr = float(row["x"]), float(row["y"]), float(row["snr"])
        dx, dy = -(x - x_mid) * pixel_scale, (y - y_mid) * pixel_scale
        marker = "o" if snr > 0 else "v"
        color = POSITIVE_COLOR if snr > 0 else NEGATIVE_COLOR
        ax.scatter(dx, dy, marker=marker, s=55, color=color, edgecolor="k", linewidth=0.5)
        ax.annotate(row["candidate_id"], (dx, dy), xytext=(4, 4), textcoords="offset points", fontsize=7)
    ax.scatter([], [], marker="o", color=POSITIVE_COLOR, edgecolor="k", label="Positive candidate")
    ax.scatter([], [], marker="v", color=NEGATIVE_COLOR, edgecolor="k", label="Negative candidate")
    ax.plot([], [], color="#e69f00", linestyle="--", label=f"PB = {pb_min:g}")
    fig.colorbar(image, ax=ax, label="Primary-beam response")
    ax.set(xlabel=r"$\Delta$RA (arcsec)", ylabel=r"$\Delta$Dec (arcsec)",
           title="Line-candidate sky distribution")
    ax.legend(frameon=True, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=200)
    plt.close(fig)


def diagnostic_summary(cube, mask, rows, threshold, pdf_path, png_path):
    pixel_snr = sampled_pixel_snr(cube, mask)
    channel_rms = np.array([
        robust_sigma(np.asarray(plane, dtype=float)[mask]) for plane in cube
    ]) * 1e3
    signed_snr = np.array([float(row["snr"]) for row in rows])
    channels = np.array([float(row["channel"]) for row in rows])
    positive = signed_snr[signed_snr > 0]
    negative = -signed_snr[signed_snr < 0]

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    bins = np.linspace(-6, 6, 81)
    axes[0].hist(pixel_snr, bins=bins, density=True, histtype="step", color="0.2", label="Cube pixels")
    grid = np.linspace(-6, 6, 400)
    axes[0].plot(grid, np.exp(-0.5 * grid**2) / np.sqrt(2 * np.pi), color="#e69f00", label="Gaussian")
    axes[0].set(xlabel="Per-channel pixel S/N", ylabel="Probability density", title="Pixel distribution")
    axes[0].set_yscale("log")
    axes[0].set_ylim(1e-5, 1)
    axes[0].legend(frameon=False, fontsize=8)

    upper = max(np.max(np.abs(signed_snr)), threshold) + 0.25
    levels = np.linspace(threshold, upper, 200)
    axes[1].step(levels, [np.sum(positive >= value) for value in levels], where="post", color=POSITIVE_COLOR, label="Positive")
    axes[1].step(levels, [np.sum(negative >= value) for value in levels], where="post", color=NEGATIVE_COLOR, label="Negative")
    axes[1].set(xlabel="Peak |S/N| threshold", ylabel="Cumulative candidate count", title="Candidate S/N distribution")
    axes[1].legend(frameon=False, fontsize=8)

    colors = np.where(signed_snr > 0, POSITIVE_COLOR, NEGATIVE_COLOR)
    markers = np.where(signed_snr > 0, "o", "v")
    for channel, snr, color, marker in zip(channels, signed_snr, colors, markers):
        axes[2].scatter(channel, snr, color=color, marker=marker, s=38)
    axes[2].axhline(threshold, color="0.4", linestyle=":")
    axes[2].axhline(-threshold, color="0.4", linestyle=":")
    axes[2].set(xlabel="Channel", ylabel="Signed peak S/N", title="Channel versus candidate S/N")
    rms_axis = axes[2].twinx()
    rms_axis.plot(
        np.arange(cube.shape[0]), channel_rms, color="0.45", linewidth=1.0,
        alpha=0.75, label="Channel RMS",
    )
    rms_axis.set_ylabel("Channel RMS (mJy beam$^{-1}$)", color="0.35")
    rms_axis.tick_params(axis="y", colors="0.35")
    rms_axis.legend(frameon=False, fontsize=8, loc="upper right")
    for ax in axes:
        ax.grid(alpha=0.2)
    fig.suptitle("DICE signed line-search diagnostics")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=200)
    plt.close(fig)


def main():
    parser = ArgumentParser(description="Plot DICE candidate review and diagnostic products.")
    parser.add_argument("config", type=Path, help="JSON used for search_fits_cube.py")
    parser.add_argument("--radius-pixels", type=int, default=18)
    parser.add_argument("--max-candidates", type=int)
    args = parser.parse_args()

    config_path = resolve(args.config)
    config = json.loads(config_path.read_text())
    image_path = resolve(config["image_fits"])
    pb_path = resolve(config["pb_fits"])
    output = resolve(config["output_directory"])
    catalogue = output / "signed_line_candidates.csv"
    paths = {
        "review": output / "candidate_review.pdf",
        "sky_pdf": output / "candidate_sky_distribution.pdf",
        "sky_png": output / "candidate_sky_distribution.png",
        "diagnostics_pdf": output / "line_search_diagnostics.pdf",
        "diagnostics_png": output / "line_search_diagnostics.png",
    }

    with catalogue.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit(f"No candidates to plot in {catalogue}")
    if args.max_candidates is not None:
        rows = rows[: args.max_candidates]

    with fits.open(image_path, memmap=True) as hdul:
        cube = np.squeeze(hdul[0].data)
        if cube.ndim != 3:
            raise ValueError(f"Expected a 3-D spectral cube after squeeze; got {cube.shape}")
        cube = np.asarray(cube, dtype=float)
        wcs = WCS(hdul[0].header)
    pb = load_pb(pb_path)
    if pb.shape != cube.shape[1:]:
        raise ValueError(f"PB shape {pb.shape} does not match cube spatial shape {cube.shape[1:]}")
    mask = np.isfinite(pb) & (pb >= float(config["pb_min"]))

    candidate_review_pdf(cube, wcs, rows, paths["review"], args.radius_pixels)
    sky_distribution(pb, wcs, rows, float(config["pb_min"]), paths["sky_pdf"], paths["sky_png"])
    diagnostic_summary(cube, mask, rows, float(config["threshold_snr"]), paths["diagnostics_pdf"], paths["diagnostics_png"])

    provenance = output / "candidate_plots.provenance.json"
    write_provenance(
        provenance,
        root=ROOT,
        script=Path(__file__).resolve(),
        command=[
            sys.executable,
            str(Path(__file__).resolve().relative_to(ROOT)),
            str(args.config),
        ],
        config=config_path,
        inputs=[image_path, pb_path, catalogue],
        outputs=list(paths.values()),
        parameters={"radius_pixels": args.radius_pixels, "max_candidates": args.max_candidates},
    )
    print(f"Wrote {paths['review']} ({len(rows)} candidates)")
    print(f"Wrote {paths['sky_pdf']}")
    print(f"Wrote {paths['diagnostics_pdf']}")


if __name__ == "__main__":
    main()
