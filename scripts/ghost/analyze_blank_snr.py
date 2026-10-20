"""Plot the blank GHOST F1 dirty-image S/N distribution.

The analysis is restricted to the usable primary-beam region (PB >= 0.2) and
compares empirical positive/negative tails with a unit Gaussian. It also runs
the BliSS extractor with 4- and 5-sigma seeds before selecting the population
test threshold.
"""

from pathlib import Path
import json
import math
import sys

import matplotlib.pyplot as plt
import numpy as np
from casatasks import imstat
from casatools import image
from scipy import stats


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from bliss_alma.extraction import extract_sources
from bliss_alma.provenance import write_script_provenance


WORKDIR = REPO_ROOT / "data" / "ghost_f1_5min"
BLANK_IMAGE = WORKDIR / "bliss_f1_dirty.image"
PB_IMAGE = WORKDIR / "bliss_f1_dirty.pb"
OUTPUT_DIR = WORKDIR / "figures"
SUMMARY = WORKDIR / "f1_blank_snr_distribution.json"
PB_LIMIT = 0.2
FLOOD_SIGMA = 3.0
MIN_PIXELS = 5


def read_plane(path):
    ia = image()
    try:
        ia.open(str(path))
        pixels = np.squeeze(np.asarray(ia.getchunk(), dtype=np.float64))
        valid = np.squeeze(np.asarray(ia.getchunk(getmask=True), dtype=bool))
    finally:
        ia.done()
    return pixels, valid


def gaussian_tail_probability(sigma):
    return 0.5 * math.erfc(float(sigma) / math.sqrt(2.0))


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    dirty, dirty_valid = read_plane(BLANK_IMAGE)
    pb, pb_valid = read_plane(PB_IMAGE)
    usable = dirty_valid & pb_valid & np.isfinite(dirty) & np.isfinite(pb) & (pb >= PB_LIMIT)
    rms = float(imstat(str(BLANK_IMAGE), listit=False)["rms"][0])
    snr_map = np.where(usable, dirty / rms, 0.0)
    snr = snr_map[usable]
    detection_mask = ~usable

    extractor = {}
    for seed_sigma in (4.0, 5.0):
        positive = extract_sources(
            snr_map,
            rms=1.0,
            seed_sigma=seed_sigma,
            flood_sigma=FLOOD_SIGMA,
            min_pixels=MIN_PIXELS,
            detection_mask=detection_mask,
            exclude_border=True,
        )
        negative = extract_sources(
            snr_map,
            rms=1.0,
            seed_sigma=seed_sigma,
            flood_sigma=FLOOD_SIGMA,
            min_pixels=MIN_PIXELS,
            detection_mask=detection_mask,
            exclude_border=True,
            polarity="negative",
        )
        extractor[str(seed_sigma)] = {
            "n_positive_islands": len(positive.catalogue),
            "n_negative_islands": len(negative.catalogue),
            "positive_peak_snr": [float(item.peak_value) for item in positive.catalogue],
            "negative_peak_snr": [float(item.peak_value) for item in negative.catalogue],
        }

    thresholds = {}
    for threshold in (3.0, 4.0, 5.0):
        thresholds[str(threshold)] = {
            "n_positive_pixels": int(np.count_nonzero(snr >= threshold)),
            "n_negative_pixels": int(np.count_nonzero(snr <= -threshold)),
            "gaussian_expected_one_sided_pixels": float(
                snr.size * gaussian_tail_probability(threshold)
            ),
        }

    summary = {
        "n_usable_pixels": int(snr.size),
        "pb_limit": PB_LIMIT,
        "rms_jy_per_beam": rms,
        "snr_mean": float(np.mean(snr)),
        "snr_std": float(np.std(snr)),
        "snr_median": float(np.median(snr)),
        "snr_skewness": float(stats.skew(snr)),
        "snr_excess_kurtosis": float(stats.kurtosis(snr)),
        "snr_minimum": float(np.min(snr)),
        "snr_maximum": float(np.max(snr)),
        "pixel_tails": thresholds,
        "extractor": extractor,
    }
    with SUMMARY.open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    edges = np.linspace(-6.0, 6.0, 121)
    centers = 0.5 * (edges[:-1] + edges[1:])
    histogram, _ = np.histogram(snr, bins=edges, density=True)
    gaussian = np.exp(-0.5 * centers**2) / math.sqrt(2.0 * math.pi)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    axes[0].step(centers, histogram, where="mid", label="F1 blank", color="#4477AA")
    axes[0].plot(centers, gaussian, label="N(0,1)", color="#333333", linestyle="--")
    axes[0].set_xlim(-6, 6)
    axes[0].set_xlabel("Pixel S/N")
    axes[0].set_ylabel("Probability density")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.2)

    tail_centers = np.abs(centers[centers >= 0])
    positive_density = histogram[centers >= 0]
    negative_density = histogram[centers < 0][::-1]
    axes[1].step(tail_centers, positive_density, where="mid", label="positive", color="#CC6677")
    axes[1].step(tail_centers, negative_density, where="mid", label="negative", color="#4477AA")
    axes[1].plot(tail_centers, gaussian[centers >= 0], label="N(0,1)", color="#333333", linestyle="--")
    axes[1].axvline(4.0, color="#228833", linestyle=":", label="4σ seed")
    axes[1].axvline(5.0, color="#AA3377", linestyle=":", label="5σ seed")
    axes[1].set_yscale("log")
    axes[1].set_xlim(2.5, 6.0)
    axes[1].set_ylim(1e-6, 0.1)
    axes[1].set_xlabel("|Pixel S/N|")
    axes[1].set_ylabel("Probability density")
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=0.2)
    fig.suptitle("GHOST F1 blank dirty-image S/N distribution")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "f1_blank_snr_distribution.png", dpi=180)
    fig.savefig(OUTPUT_DIR / "f1_blank_snr_distribution.pdf")
    plt.close(fig)
    write_script_provenance(OUTPUT_DIR/"f1_blank_snr_distribution.provenance.json",
        root=REPO_ROOT,script=__file__,inputs=[BLANK_IMAGE,PB_IMAGE],
        outputs=[SUMMARY,OUTPUT_DIR/"f1_blank_snr_distribution.png",OUTPUT_DIR/"f1_blank_snr_distribution.pdf"],
        parameters={"pb_min":PB_LIMIT,"flood_sigma":FLOOD_SIGMA})

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"Wrote {OUTPUT_DIR / 'f1_blank_snr_distribution.png'}")


if __name__ == "__main__":
    main()
