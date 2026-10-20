"""Blind matched-filter line search in the Hayatsu et al. (2017) ADF22 cubes.

The input FITS files are the historical continuum-subtracted dirty cubes.
This implementation deliberately keeps detection and comparison with the
published catalogue separate: candidates are found without using known source
positions, and are matched to the literature catalogue only after extraction.
"""

from pathlib import Path
import csv

from astropy.io import fits
from astropy.wcs import WCS
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/hayatsu_reproduction/adf22_image_fits"
OUT = ROOT / "data/hayatsu_reproduction/results"
OUT.mkdir(parents=True, exist_ok=True)

# The published catalogue was selected at approximately 6 sigma.  A 5.5-sigma
# extraction floor retains margin for implementation-level RMS differences.
THRESHOLD = 5.5
PB_LIMIT = 0.8
WIDTHS = (6, 9, 12, 15, 18, 21)

PUBLISHED = [
    ("ADF22-LineA", "22:17:37.43", "+00:17:10.7", 253.79, 6.5, 12),
    ("ADF22-LineB", "22:17:31.95", "+00:18:20.3", 269.92, 6.2, 15),
    ("ADF22-LineC", "22:17:36.97", "+00:18:20.8", 253.49, 7.8, 12),
    ("ADF22-LineD", "22:17:33.07", "+00:17:18.8", 269.70, 6.5, 21),
]


def image_sigma(values):
    values = np.asarray(values)
    values = values[np.isfinite(values) & (values != 0)]
    if values.size == 0:
        return np.nan
    return np.std(values)


def merge_candidates(rows, spatial_radius=7.0, channel_radius=12):
    """Greedily retain the strongest peak across widths and adjacent voxels."""
    kept = []
    for row in sorted(rows, key=lambda item: abs(item["snr"]), reverse=True):
        duplicate = False
        for old in kept:
            distance = np.hypot(row["x"] - old["x"], row["y"] - old["y"])
            if (row["spw"] == old["spw"] and distance <= spatial_radius
                    and abs(row["channel"] - old["channel"]) <= channel_radius
                    and np.sign(row["snr"]) == np.sign(old["snr"])):
                duplicate = True
                break
        if not duplicate:
            kept.append(row)
    return sorted(kept, key=lambda item: item["snr"], reverse=True)


pb_path = DATA / "spw0.Core_region.ms.dirtycube.pb.fits"
with fits.open(pb_path, memmap=True) as hdul:
    pb_mid = np.asarray(hdul[0].data[0, hdul[0].data.shape[1] // 2])
search_mask = np.isfinite(pb_mid) & (pb_mid >= PB_LIMIT)

raw_candidates = []
diagnostics = []
for spw in range(4):
    path = DATA / f"spw{spw}.Core_region.ms.dirtycube.image.fits"
    with fits.open(path, memmap=True) as hdul:
        cube = hdul[0].data[0]
        header = hdul[0].header.copy()
        wcs = WCS(header)
        nchan = cube.shape[0]
        plane_sigma = np.array([
            image_sigma(np.asarray(cube[channel])[search_mask])
            for channel in range(nchan)
        ])
        good = np.isfinite(plane_sigma) & (plane_sigma > 0)
        good &= plane_sigma <= 0.9 * np.nanmax(plane_sigma)
        valid_channels = np.where(good)[0]
        start_channel = int(valid_channels.min())
        end_channel = int(valid_channels.max())
        # The archived FITS cubes retain continuum emission, whereas the
        # published search used uvcontsub with fitorder=0.  Subtract the
        # corresponding constant spectrum per image pixel before filtering.
        continuum = np.zeros(cube.shape[1:], dtype=np.float64)
        continuum_count = np.zeros(cube.shape[1:], dtype=np.int16)
        for continuum_channel in valid_channels:
            continuum_plane = np.asarray(
                cube[continuum_channel], dtype=np.float64
            )
            finite_plane = np.isfinite(continuum_plane)
            continuum[finite_plane] += continuum_plane[finite_plane]
            continuum_count[finite_plane] += 1
        valid_continuum = continuum_count > 0
        continuum[valid_continuum] /= continuum_count[valid_continuum]

        # Cumulative sum permits each boxcar plane to be evaluated without
        # materialising another multi-gigabyte cube.
        ring = []
        rolling = {
            width: np.zeros(cube.shape[1:], dtype=np.float64) for width in WIDTHS
        }
        for channel in range(nchan):
            plane = np.asarray(cube[channel], dtype=np.float64)
            plane = np.where(np.isfinite(plane), plane - continuum, 0.0)
            ring.append(plane)
            if len(ring) > max(WIDTHS):
                ring.pop(0)

            for width in WIDTHS:
                rolling[width] += plane
                if channel >= width:
                    old_plane = np.asarray(cube[channel - width],
                                           dtype=np.float64)
                    old_plane = np.where(
                        np.isfinite(old_plane), old_plane - continuum, 0.0
                    )
                    rolling[width] -= old_plane
                if len(ring) < width:
                    continue
                if channel - width + 1 < start_channel or channel > end_channel:
                    continue
                summed = rolling[width]
                sigma = image_sigma(summed[search_mask])
                if not np.isfinite(sigma) or sigma <= 0:
                    continue
                snr = summed / sigma
                centre = channel - (width - 1) / 2.0
                diagnostics.append((spw, width, centre, sigma))

                for polarity in (1, -1):
                    work = polarity * snr
                    yy, xx = np.where(search_mask & (work >= THRESHOLD))
                    for y, x in zip(yy, xx):
                        world = wcs.pixel_to_world_values(
                            float(x), float(y), float(centre), 0.0
                        )
                        raw_candidates.append({
                            "spw": spw,
                            "x": int(x),
                            "y": int(y),
                            "channel": float(centre),
                            "width_channels": width,
                            "snr": float(snr[y, x]),
                            "ra_deg": float(world[0]),
                            "dec_deg": float(world[1]),
                            "frequency_ghz": float(world[2]) / 1e9,
                        })

candidates = merge_candidates(raw_candidates)
fieldnames = list(candidates[0].keys()) if candidates else [
    "spw", "x", "y", "channel", "width_channels", "snr",
    "ra_deg", "dec_deg", "frequency_ghz",
]
with (OUT / "adf22_blind_candidates.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(candidates)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
positive = sorted([row["snr"] for row in candidates if row["snr"] > 0], reverse=True)
negative = sorted([-row["snr"] for row in candidates if row["snr"] < 0], reverse=True)
for values, label, color in [
    (positive, "Positive", "#c83e4d"), (negative, "Inverted", "#377eb8")
]:
    if values:
        axes[0].step(values, np.arange(1, len(values) + 1),
                     where="post", label=label, color=color)
axes[0].set(xlabel="Peak |S/N|", ylabel="Cumulative number",
            xlim=(THRESHOLD, max([8] + positive + negative)))
axes[0].legend(frameon=False)

for row in candidates:
    axes[1].scatter(row["frequency_ghz"], row["snr"],
                    color="#c83e4d" if row["snr"] > 0 else "#377eb8", s=18)
axes[1].axhline(THRESHOLD, color="0.4", linestyle=":")
axes[1].axhline(-THRESHOLD, color="0.4", linestyle=":")
axes[1].set(xlabel="Observed frequency (GHz)", ylabel="Peak S/N")
for axis in axes:
    axis.grid(alpha=0.2)
fig.tight_layout()
fig.savefig(OUT / "adf22_blind_search_diagnostic.png", dpi=220)
print("Wrote", OUT / "adf22_blind_candidates.csv")
print("Wrote", OUT / "adf22_blind_search_diagnostic.png")
