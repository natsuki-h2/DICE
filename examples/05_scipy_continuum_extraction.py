"""Extract positive and negative sources from the CASA continuum test."""

from pathlib import Path
import sys

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from bliss_alma.extraction import extract_sources


DATA_DIR = REPO_ROOT / "data" / "casa_continuum_e2e"
SEED_SIGMA = 5.0
FLOOD_SIGMA = 3.0
MIN_PIXELS = 5


def print_catalogue(name, result):
    print(name)
    print("  RMS: {0:.6g} Jy/beam".format(result.rms))
    print("  Seed threshold: {0:.6g} Jy/beam".format(result.seed_threshold))
    print("  Flood threshold: {0:.6g} Jy/beam".format(result.flood_threshold))
    print("  Number of sources: {0}".format(len(result.catalogue)))
    for source in result.catalogue:
        print("  {0}".format(source))


def main():
    injected = np.load(str(DATA_DIR / "injected.npy"))
    reference = np.load(str(DATA_DIR / "reference.npy"))
    reference_rms = float(np.load(str(DATA_DIR / "reference_rms_jy.npy")))

    common = {
        "rms": reference_rms,
        "seed_sigma": SEED_SIGMA,
        "flood_sigma": FLOOD_SIGMA,
        "min_pixels": MIN_PIXELS,
        "exclude_border": True,
    }
    injected_positive = extract_sources(injected, polarity="positive", **common)
    blank_positive = extract_sources(reference, polarity="positive", **common)
    blank_negative = extract_sources(reference, polarity="negative", **common)

    np.save(str(DATA_DIR / "injected_source_labels.npy"), injected_positive.labels)
    np.save(str(DATA_DIR / "blank_positive_labels.npy"), blank_positive.labels)
    np.save(str(DATA_DIR / "blank_negative_labels.npy"), blank_negative.labels)

    print_catalogue("Injected image: positive detections", injected_positive)
    print_catalogue("Blank image: positive detections", blank_positive)
    print_catalogue("Blank image: negative detections", blank_negative)

    expected_pixel = np.array([64, 64])
    recovered_central = any(
        np.linalg.norm(np.array(source.peak_index) - expected_pixel) <= 1.0
        for source in injected_positive.catalogue
    )
    if not recovered_central:
        raise RuntimeError("The phase-centre continuum source was not detected.")
    print("PASS: SciPy source extraction recovered the injected source.")


if __name__ == "__main__":
    main()
