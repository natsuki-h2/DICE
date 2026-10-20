"""Export the existing CASA diagnostic images as NumPy arrays."""

from pathlib import Path
import numpy as np
from casatools import image

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "data/ghost_matrix/common_resolution_runs/f6_5min/visualization_work"


def read(path):
    ia = image()
    try:
        ia.open(str(path))
        return np.squeeze(np.asarray(ia.getchunk(), dtype=float))
    finally:
        ia.done()


for name, casa_name in {
    "populated": "populated.smooth.image",
    "source_only": "source_only.smooth.image",
    "psf": "populated.psf",
}.items():
    np.save(WORK / (name + ".npy"), read(WORK / casa_name))
