"""Search one exactly uv-continuum-subtracted ADF22 dirty cube."""

from collections import deque
from pathlib import Path
import csv
import sys

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/hayatsu_reproduction/adf22_contsub_images"
OUT = ROOT / "data/hayatsu_reproduction/results"
OUT.mkdir(parents=True, exist_ok=True)

spw = int(sys.argv[-1])
basename = DATA / f"spw{spw}.Core_region.contsub.dirtycube"
WIDTHS = (6, 9, 12, 15, 18, 21)
THRESHOLD = 5.5
PB_LIMIT = 0.8
PUBLISHED = [
    ("ADF22-LineA", "22:17:37.43", "+00:17:10.7", 253.79, 6.5, 12),
    ("ADF22-LineB", "22:17:31.95", "+00:18:20.3", 269.92, 6.2, 15),
    ("ADF22-LineC", "22:17:36.97", "+00:18:20.8", 253.49, 7.8, 12),
    ("ADF22-LineD", "22:17:33.07", "+00:17:18.8", 269.70, 6.5, 21),
]


def sigma(values):
    values = np.asarray(values)
    values = values[np.isfinite(values) & (values != 0)]
    return np.std(values) if values.size else np.nan


def merge(rows, spatial_radius=7.0, channel_radius=12):
    kept = []
    for row in sorted(rows, key=lambda item: abs(item["snr"]), reverse=True):
        if any(
            np.hypot(row["x"] - old["x"], row["y"] - old["y"])
            <= spatial_radius
            and abs(row["channel"] - old["channel"]) <= channel_radius
            and np.sign(row["snr"]) == np.sign(old["snr"])
            for old in kept
        ):
            continue
        kept.append(row)
    return sorted(kept, key=lambda item: item["snr"], reverse=True)


pb_path = Path(str(basename) + ".pb.fits")
with fits.open(pb_path, memmap=True) as hdul:
    pb_cube = hdul[0].data[0]
    pb_mid = np.asarray(pb_cube[pb_cube.shape[0] // 2])
if not np.any(np.isfinite(pb_mid) & (pb_mid > 0)):
    # The PB is independent of uv continuum subtraction.  A calcpsf=False
    # restart can leave CASA's new PB image blank, so use the archived PB
    # made with the same grid and measurement set in that case.
    pb_path = (
        ROOT / "data/hayatsu_reproduction/adf22_image_fits"
        / f"spw{spw}.Core_region.ms.dirtycube.pb.fits"
    )
    with fits.open(pb_path, memmap=True) as hdul:
        pb_cube = hdul[0].data[0]
        pb_mid = np.asarray(pb_cube[pb_cube.shape[0] // 2])
search_mask = np.isfinite(pb_mid) & (pb_mid >= PB_LIMIT)

with fits.open(str(basename) + ".image.fits", memmap=True) as hdul:
    cube = hdul[0].data[0]
    header = hdul[0].header.copy()
    wcs = WCS(header)
    nchan = cube.shape[0]
    plane_sigma = np.array(
        [sigma(np.asarray(cube[channel])[search_mask]) for channel in range(nchan)]
    )
    good = np.isfinite(plane_sigma) & (plane_sigma > 0)
    good &= plane_sigma <= 0.9 * np.nanmax(plane_sigma)
    valid = np.where(good)[0]
    first, last = int(valid.min()), int(valid.max())

    published_pixels = {}
    for name, ra, dec, frequency, published_snr, published_width in PUBLISHED:
        coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
        x, y, spectral_pixel, _stokes_pixel = wcs.world_to_pixel_values(
            coord.ra.deg, coord.dec.deg, frequency * 1e9, 0.0
        )
        x = float(np.asarray(x))
        y = float(np.asarray(y))
        spectral_pixel = float(np.asarray(spectral_pixel))
        if -1 <= x < cube.shape[2] and -1 <= y < cube.shape[1]:
            published_pixels[name] = (
                int(round(x)), int(round(y)), spectral_pixel, frequency,
                published_snr, published_width,
            )

    rolling = {
        width: np.zeros(cube.shape[1:], dtype=np.float64) for width in WIDTHS
    }
    ring = deque(maxlen=max(WIDTHS))
    raw = []
    measurements = []
    for channel in range(nchan):
        plane = np.asarray(cube[channel], dtype=np.float64)
        plane = np.where(np.isfinite(plane), plane, 0.0)
        for width in WIDTHS:
            rolling[width] += plane
            if len(ring) >= width:
                rolling[width] -= ring[-width]
            start = channel - width + 1
            if start < first or channel > last:
                continue
            snr_map = rolling[width] / sigma(rolling[width][search_mask])
            centre = channel - (width - 1) / 2.0

            for polarity in (1, -1):
                work = polarity * snr_map
                peaks = search_mask & (work >= THRESHOLD)
                for y, x in zip(*np.where(peaks)):
                    world = wcs.pixel_to_world_values(
                        float(x), float(y), float(centre), 0.0
                    )
                    raw.append({
                        "spw": spw, "x": int(x), "y": int(y),
                        "channel": float(centre), "width_channels": width,
                        "snr": float(snr_map[y, x]),
                        "ra_deg": float(world[0]), "dec_deg": float(world[1]),
                        "frequency_ghz": float(world[2]) / 1e9,
                    })

            for name, (
                x, y, target_channel, frequency, pub_snr, pub_width
            ) in published_pixels.items():
                if width != pub_width:
                    continue
                if abs(centre - target_channel) > 0.6:
                    continue
                cut = snr_map[max(0, y - 7):y + 8, max(0, x - 7):x + 8]
                iy, ix = np.unravel_index(np.nanargmax(cut), cut.shape)
                measurements.append({
                    "name": name, "spw": spw, "width_channels": width,
                    "published_snr": pub_snr,
                    "measured_snr_at_position": float(snr_map[y, x]),
                    "measured_local_peak_snr": float(cut[iy, ix]),
                    "position_offset_pixels": float(
                        np.hypot(ix + max(0, x - 7) - x,
                                 iy + max(0, y - 7) - y)
                    ),
                })
        ring.append(plane)

candidates = merge(raw)
for filename, rows in [
    (f"adf22_exact_spw{spw}_blind_candidates.csv", candidates),
    (f"adf22_exact_spw{spw}_published_measurements.csv", measurements),
]:
    path = OUT / filename
    fields = list(rows[0]) if rows else ["spw"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print("Wrote", path)
