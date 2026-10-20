"""Plot the fixed-width, sparse line-injection tutorial design."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bliss_alma.injection import plan_line_injection_positions
from bliss_alma.provenance import write_provenance


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "tutorials" / "line_injection_demo"
PEAK_SNRS = np.array([2.5, 3.5, 4.5, 5.5, 7.0])
FWHM_CHANNELS = 10.0
CUBE_SHAPE = (240, 100, 100)
SEPARATION = 8.0
PROTECTED_FWHM = 4.0
MARGIN = 8
SEED = 7


def _primary_beam(shape):
    _, ny, nx = shape
    yy, xx = np.mgrid[:ny, :nx]
    radius = np.hypot(xx - (nx - 1) / 2, yy - (ny - 1) / 2)
    half_power_radius = 0.42 * min(nx, ny)
    return xx, yy, np.exp(-np.log(2.0) * (radius / half_power_radius) ** 2)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = plan_line_injection_positions(
        [FWHM_CHANNELS] * len(PEAK_SNRS), CUBE_SHAPE,
        min_spatial_separation_pixels=SEPARATION,
        spectral_extent_fwhm=PROTECTED_FWHM,
        spatial_edge_margin_pixels=MARGIN,
        seed=SEED,
    )
    csv_path = OUT / "expected_helix_placements.csv"
    with csv_path.open("w", newline="") as handle:
        fieldnames = [*rows[0]._fields, "peak_snr"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row, peak_snr in zip(rows, PEAK_SNRS):
            writer.writerow({**row._asdict(), "peak_snr": float(peak_snr)})

    norm = plt.Normalize(PEAK_SNRS.min(), PEAK_SNRS.max())
    cmap = plt.get_cmap("viridis")
    fig, (sky, spectra) = plt.subplots(
        1, 2, figsize=(11.2, 4.7), constrained_layout=True,
        gridspec_kw={"width_ratios": [1.0, 1.16]},
    )

    xx, yy, pb = _primary_beam(CUBE_SHAPE)
    background = sky.imshow(pb, origin="lower", cmap="Greys", vmin=0, vmax=1,
                            alpha=0.45)
    sky.contour(xx, yy, pb, levels=[0.5, 0.8], colors=["0.35", "0.15"],
                linestyles=["--", ":"], linewidths=1.1)
    points = sky.scatter(
        [row.x for row in rows], [row.y for row in rows], c=PEAK_SNRS,
        cmap=cmap, norm=norm, s=92, edgecolor="black", linewidth=0.8, zorder=3,
    )
    for index, row in enumerate(rows):
        sky.annotate(f"S{index + 1}", (row.x, row.y), xytext=(5, 5),
                     textcoords="offset points", fontsize=9, weight="semibold")
    sky.set(title="(a) Artificial-source positions", xlabel="x pixel", ylabel="y pixel",
            xlim=(-0.5, CUBE_SHAPE[2] - 0.5), ylim=(-0.5, CUBE_SHAPE[1] - 0.5),
            aspect="equal")
    sky.invert_xaxis()
    pb_bar = fig.colorbar(background, ax=sky, location="left", pad=0.11, shrink=0.84)
    pb_bar.set_label("Primary-beam response")
    snr_bar = fig.colorbar(points, ax=sky, pad=0.02, shrink=0.84)
    snr_bar.set_label("Input peak S/N")

    channel_grid = np.linspace(0, CUBE_SHAPE[0] - 1, 1800)
    sigma = FWHM_CHANNELS / (2 * np.sqrt(2 * np.log(2)))
    half_protected = 0.5 * PROTECTED_FWHM * FWHM_CHANNELS
    ordered = sorted(zip(rows, PEAK_SNRS), key=lambda item: item[0].channel)
    for level, (row, peak_snr) in enumerate(ordered, start=1):
        profile = peak_snr * np.exp(-0.5 * ((channel_grid - row.channel) / sigma) ** 2)
        color = cmap(norm(peak_snr))
        local = np.abs(channel_grid - row.channel) <= 3.0 * sigma
        spectra.hlines(level, 0, CUBE_SHAPE[0] - 1, color="0.82", lw=0.8)
        spectra.fill_between(
            [row.channel - half_protected, row.channel + half_protected],
            level - 0.31, level + 0.31, color=color, alpha=0.16, linewidth=0,
        )
        spectra.plot(channel_grid[local], level + 0.085 * profile[local],
                     color=color, lw=2.2)
        spectra.plot([row.channel - half_protected, row.channel + half_protected],
                     [level, level], color=color, lw=3)
        spectra.scatter(row.channel, level, color=color, edgecolor="black",
                        linewidth=0.7, s=48, zorder=3)
    spectra.set(
        title=f"(b) Fixed FWHM = {FWHM_CHANNELS:g} channels",
        xlabel="Channel", ylabel="Artificial source",
        xlim=(0, CUBE_SHAPE[0] - 1), ylim=(0.45, len(ordered) + 0.95),
        yticks=range(1, len(ordered) + 1),
        yticklabels=[f"S{row.source_index + 1}  (peak {snr:g}σ)"
                     for row, snr in ordered],
    )
    spectra.grid(axis="x", alpha=0.22)
    spectra.text(
        0.99, 0.02,
        "Shading/thick line: protected interval (4 × FWHM)\n"
        "A separate cube is used for each FWHM.",
        transform=spectra.transAxes, ha="right", va="bottom", fontsize=8.2,
        bbox={"facecolor": "white", "edgecolor": "0.75", "alpha": 0.92},
    )

    png = OUT / "line_injection_helix_example.png"
    pdf = OUT / "line_injection_helix_example.pdf"
    fig.savefig(png, dpi=220)
    fig.savefig(pdf)
    plt.close(fig)
    write_provenance(
        OUT / "line_injection_helix_example.provenance.json",
        root=ROOT, script=Path(__file__), command=["python", str(Path(__file__))],
        inputs=[], outputs=[png, pdf, csv_path],
        parameters={
            "peak_snr": PEAK_SNRS.tolist(), "fwhm_channels": FWHM_CHANNELS,
            "cube_shape": CUBE_SHAPE, "min_spatial_separation_pixels": SEPARATION,
            "spectral_extent_fwhm": PROTECTED_FWHM,
            "spatial_edge_margin_pixels": MARGIN, "seed": SEED,
        }, random_seeds=[SEED],
    )
    print(png)


if __name__ == "__main__":
    main()
