"""Match the blind ADF22 catalogue to Hayatsu et al. (2017) and plot it."""

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
from bliss_alma.provenance import write_script_provenance


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/hayatsu_reproduction/adf22_image_fits"
OUT = ROOT / "data/hayatsu_reproduction/results"

PUBLISHED = [
    ("ADF22-LineA", "22:17:37.43", "+00:17:10.7", 253.79, 6.5, 12),
    ("ADF22-LineB", "22:17:31.95", "+00:18:20.3", 269.92, 6.2, 15),
    ("ADF22-LineC", "22:17:36.97", "+00:18:20.8", 253.49, 7.8, 12),
    ("ADF22-LineD", "22:17:33.07", "+00:17:18.8", 269.70, 6.5, 21),
]

with (OUT / "adf22_blind_candidates.csv").open() as handle:
    candidates = list(csv.DictReader(handle))
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
        if row["snr"] > 0 and separation <= 1.0:
            nearby.append((abs(row["frequency_ghz"] - frequency), separation, row))
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

with (OUT / "adf22_published_candidate_recovery.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=matches[0].keys())
    writer.writeheader()
    writer.writerows(matches)

fig, axes = plt.subplots(4, 2, figsize=(9, 11))
for row_index, match in enumerate(matches):
    spw = match["spw"]
    path = DATA / f"spw{spw}.Core_region.ms.dirtycube.image.fits"
    with fits.open(path, memmap=True) as hdul:
        cube = hdul[0].data[0]
        wcs = WCS(hdul[0].header)
        valid = np.arange(9 if spw < 2 else 10, 119 if spw < 2 else 120)
        continuum = np.zeros(cube.shape[1:], dtype=np.float64)
        continuum_count = np.zeros(cube.shape[1:], dtype=np.int16)
        for continuum_channel in valid:
            plane = np.asarray(cube[continuum_channel], dtype=np.float64)
            finite = np.isfinite(plane)
            continuum[finite] += plane[finite]
            continuum_count[finite] += 1
        finite_continuum = continuum_count > 0
        continuum[finite_continuum] /= continuum_count[finite_continuum]
        x, y = match["x"], match["y"]
        channel = int(round(match["channel"]))
        width = match["recovered_width_channels"]
        first = channel - width // 2
        last = first + width
        moment = np.nansum(
            np.asarray(cube[first:last], dtype=np.float64) - continuum, axis=0
        )
        spectrum = np.asarray(cube[:, y, x], dtype=np.float64)
        spectrum -= np.nanmean(spectrum[valid])
        frequencies = wcs.pixel_to_world_values(
            np.full(cube.shape[0], x), np.full(cube.shape[0], y),
            np.arange(cube.shape[0]), np.zeros(cube.shape[0])
        )[2] / 1e9

    radius = 35
    cutout = moment[y-radius:y+radius+1, x-radius:x+radius+1]
    rms = np.nanstd(moment)
    axes[row_index, 0].imshow(
        cutout, origin="lower", cmap="RdBu_r", vmin=-6*rms, vmax=6*rms
    )
    axes[row_index, 0].axhline(radius, color="k", alpha=0.4, linewidth=0.7)
    axes[row_index, 0].axvline(radius, color="k", alpha=0.4, linewidth=0.7)
    axes[row_index, 0].set_title(
        f"{match['name']}: recovered {match['recovered_snr']:.2f}σ"
    )
    axes[row_index, 0].set(xticks=[], yticks=[])

    axes[row_index, 1].step(frequencies, spectrum * 1e3,
                            where="mid", color="0.2")
    axes[row_index, 1].axvspan(
        frequencies[max(0, first)], frequencies[min(len(frequencies)-1, last-1)],
        color="#377eb8", alpha=0.18,
    )
    axes[row_index, 1].axvline(
        match["published_frequency_ghz"], color="#c83e4d", linestyle=":",
        label="Published peak" if row_index == 0 else None,
    )
    axes[row_index, 1].set_ylabel("Flux density (mJy beam$^{-1}$)")
    axes[row_index, 1].grid(alpha=0.2)

axes[-1, 1].set_xlabel("Observed frequency (GHz)")
axes[0, 1].legend(frameon=False, fontsize=8)
fig.suptitle(
    "ADF22 untargeted-search recovery of the four Hayatsu et al. (2017) candidates\n"
    "dirty cubes; constant continuum subtraction; PB ≥ 0.8",
    y=0.995,
)
fig.tight_layout()
fig.savefig(OUT / "adf22_four_candidate_recovery.png", dpi=220)
fig.savefig(OUT / "adf22_four_candidate_recovery.pdf")
write_script_provenance(OUT/"adf22_four_candidate_recovery.provenance.json",root=ROOT,script=__file__,
    inputs=[DATA,OUT/"adf22_blind_candidates.csv"],outputs=[OUT/"adf22_published_candidate_recovery.csv",OUT/"adf22_four_candidate_recovery.png",OUT/"adf22_four_candidate_recovery.pdf"],
    parameters={"match_radius_arcsec":1.0,"pb_min":0.8})
print("Wrote", OUT / "adf22_published_candidate_recovery.csv")
print("Wrote", OUT / "adf22_four_candidate_recovery.png")
