"""Match and plot the exact uv-continuum-subtracted ADF22 recovery."""

from pathlib import Path
import csv

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/hayatsu_reproduction/adf22_contsub_images"
OUT = ROOT / "data/hayatsu_reproduction/results"
PUBLISHED = [
    ("ADF22-LineA", "22:17:37.43", "+00:17:10.7", 253.79, 6.5, 12),
    ("ADF22-LineB", "22:17:31.95", "+00:18:20.3", 269.92, 6.2, 15),
    ("ADF22-LineC", "22:17:36.97", "+00:18:20.8", 253.49, 7.8, 12),
    ("ADF22-LineD", "22:17:33.07", "+00:17:18.8", 269.70, 6.5, 21),
]

candidates = []
for spw in (0, 2):
    with (OUT / f"adf22_exact_spw{spw}_blind_candidates.csv").open() as handle:
        candidates.extend(csv.DictReader(handle))
for row in candidates:
    for key in ("snr", "ra_deg", "dec_deg", "frequency_ghz", "channel"):
        row[key] = float(row[key])
    for key in ("spw", "width_channels", "x", "y"):
        row[key] = int(row[key])

matches = []
for name, ra, dec, frequency, published_snr, published_width in PUBLISHED:
    target = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    nearby = []
    for row in candidates:
        coordinate = SkyCoord(row["ra_deg"] * u.deg, row["dec_deg"] * u.deg)
        separation = target.separation(coordinate).arcsec
        frequency_offset = abs(row["frequency_ghz"] - frequency)
        if row["snr"] > 0 and separation <= 1.0 and frequency_offset <= 0.15:
            nearby.append((frequency_offset, separation, row))
    if not nearby:
        raise RuntimeError(f"No untargeted-search candidate matches {name}")
    _, separation, recovered = min(nearby, key=lambda item: item[0])
    matches.append({
        "name": name,
        "published_ra": ra,
        "published_dec": dec,
        "published_frequency_ghz": frequency,
        "published_snr": published_snr,
        "published_width_channels": published_width,
        "recovered_ra_deg": recovered["ra_deg"],
        "recovered_dec_deg": recovered["dec_deg"],
        "separation_arcsec": separation,
        "recovered_frequency_ghz": recovered["frequency_ghz"],
        "frequency_offset_mhz": 1e3 * (recovered["frequency_ghz"] - frequency),
        "recovered_snr": recovered["snr"],
        "recovered_width_channels": recovered["width_channels"],
        "spw": recovered["spw"],
        "x": recovered["x"],
        "y": recovered["y"],
        "channel": recovered["channel"],
    })

table = OUT / "adf22_exact_published_candidate_recovery.csv"
with table.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=matches[0].keys())
    writer.writeheader()
    writer.writerows(matches)

fig, axes = plt.subplots(4, 2, figsize=(9, 11))
for row_index, match in enumerate(matches):
    path = (
        DATA
        / f"spw{match['spw']}.Core_region.contsub.dirtycube.image.fits"
    )
    with fits.open(path, memmap=True) as hdul:
        cube = hdul[0].data[0]
        wcs = WCS(hdul[0].header)
        x, y = match["x"], match["y"]
        channel = int(round(match["channel"]))
        width = match["recovered_width_channels"]
        first = channel - width // 2
        last = first + width
        moment = np.nansum(np.asarray(cube[first:last], dtype=np.float64), axis=0)
        spectrum = np.asarray(cube[:, y, x], dtype=np.float64)
        frequencies = wcs.pixel_to_world_values(
            np.full(cube.shape[0], x), np.full(cube.shape[0], y),
            np.arange(cube.shape[0]), np.zeros(cube.shape[0])
        )[2] / 1e9

    radius = 35
    cutout = moment[y-radius:y+radius+1, x-radius:x+radius+1]
    finite = cutout[np.isfinite(cutout) & (cutout != 0)]
    rms = np.std(finite)
    axes[row_index, 0].imshow(
        cutout, origin="lower", cmap="RdBu_r", vmin=-6*rms, vmax=6*rms
    )
    axes[row_index, 0].axhline(radius, color="k", alpha=0.4, linewidth=0.7)
    axes[row_index, 0].axvline(radius, color="k", alpha=0.4, linewidth=0.7)
    axes[row_index, 0].set_title(
        f"{match['name']}: {match['recovered_snr']:.2f}σ "
        f"(published {match['published_snr']:.1f}σ)"
    )
    axes[row_index, 0].set(xticks=[], yticks=[])

    axes[row_index, 1].step(
        frequencies, spectrum * 1e3, where="mid", color="0.2"
    )
    axes[row_index, 1].axvspan(
        frequencies[max(0, first)],
        frequencies[min(len(frequencies) - 1, last - 1)],
        color="#377eb8", alpha=0.18,
    )
    axes[row_index, 1].axvline(
        match["published_frequency_ghz"], color="#c83e4d", linestyle=":",
        label="Published frequency" if row_index == 0 else None,
    )
    axes[row_index, 1].set_ylabel("Flux density (mJy beam$^{-1}$)")
    axes[row_index, 1].grid(alpha=0.2)

axes[-1, 1].set_xlabel("Observed frequency (GHz)")
axes[0, 1].legend(frameon=False, fontsize=8)
fig.suptitle(
    "Hayatsu et al. (2017) ADF22 untargeted-search reproduction\n"
    "visibility-domain continuum subtraction; native dirty cubes; PB ≥ 0.8",
    y=0.995,
)
fig.tight_layout()
png = OUT / "adf22_exact_four_candidate_recovery.png"
pdf = OUT / "adf22_exact_four_candidate_recovery.pdf"
fig.savefig(png, dpi=220)
fig.savefig(pdf)
print("Wrote", table)
print("Wrote", png)
