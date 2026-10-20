"""Run the first GHOST F1 visibility-injection and recovery smoke test.

The downloaded blank MS is never modified. This script copies it, injects a
small five-source population with ``ft`` and ``uvsub(reverse=True)``, makes a
dirty image, and applies the BliSS SciPy extractor to populated and blank maps.
"""

from pathlib import Path
import json
import math
import shutil
import sys

import numpy as np
from casatasks import ft, imstat, tclean, uvsub
from casatools import componentlist, image


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from bliss_alma.extraction import extract_sources
from bliss_alma.recovery import match_catalogue_to_positions


WORKDIR = REPO_ROOT / "data" / "ghost_f1_5min"
BLANK_MS = WORKDIR / "output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms"
POPULATED_MS = WORKDIR / "f1_injected.ms"
COMPONENT_LIST = WORKDIR / "f1_injected_sources.cl"
BLANK_IMAGE = WORKDIR / "bliss_f1_dirty"
POPULATED_IMAGE = WORKDIR / "f1_injected_dirty"
SUMMARY = WORKDIR / "f1_injection_recovery.json"

IMAGE_SIZE = 360
IMAGE_CENTRE = (IMAGE_SIZE // 2, IMAGE_SIZE // 2)
CELL_ARCSEC = 0.12454065519280681
CELL = f"{CELL_ARCSEC}arcsec"
PHASE_CENTRE_RA_DEG = 180.0
PHASE_CENTRE_DEC_DEG = -27.8
SOURCE_FLUX_JY = 0.0001
SOURCE_OFFSET_PIXELS = 48
SOURCE_POSITIONS = [
    IMAGE_CENTRE,
    (IMAGE_CENTRE[0] - SOURCE_OFFSET_PIXELS, IMAGE_CENTRE[1]),
    (IMAGE_CENTRE[0] + SOURCE_OFFSET_PIXELS, IMAGE_CENTRE[1]),
    (IMAGE_CENTRE[0], IMAGE_CENTRE[1] - SOURCE_OFFSET_PIXELS),
    (IMAGE_CENTRE[0], IMAGE_CENTRE[1] + SOURCE_OFFSET_PIXELS),
]
SEED_SIGMA = 5.0
FLOOD_SIGMA = 3.0
MIN_PIXELS = 5
PB_LIMIT = 0.2
MAX_MATCH_DISTANCE_PIXELS = 3.0


def remove_path(path):
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def direction_for_pixel(position):
    """Approximate the local tangent-plane pixel position as a J2000 direction."""
    x_pixel, y_pixel = position
    east_arcsec = -(x_pixel - IMAGE_CENTRE[0]) * CELL_ARCSEC
    north_arcsec = (y_pixel - IMAGE_CENTRE[1]) * CELL_ARCSEC
    ra_deg = PHASE_CENTRE_RA_DEG + east_arcsec / (
        3600.0 * math.cos(math.radians(PHASE_CENTRE_DEC_DEG))
    )
    dec_deg = PHASE_CENTRE_DEC_DEG + north_arcsec / 3600.0
    return f"J2000 {ra_deg}deg {dec_deg}deg"


def inject_population(pb):
    remove_path(POPULATED_MS)
    remove_path(COMPONENT_LIST)
    shutil.copytree(BLANK_MS, POPULATED_MS)

    cl = componentlist()
    try:
        for position in SOURCE_POSITIONS:
            apparent_flux_jy = SOURCE_FLUX_JY * float(pb[position])
            cl.addcomponent(
                dir=direction_for_pixel(position),
                flux=apparent_flux_jy,
                fluxunit="Jy",
                shape="point",
            )
        cl.rename(str(COMPONENT_LIST))
    finally:
        cl.done()

    ft(
        vis=str(POPULATED_MS),
        spw="0:0",
        complist=str(COMPONENT_LIST),
        usescratch=True,
    )
    uvsub(vis=str(POPULATED_MS), reverse=True)


def make_dirty_image():
    for path in WORKDIR.glob(POPULATED_IMAGE.name + ".*"):
        remove_path(path)
    tclean(
        vis=str(POPULATED_MS),
        imagename=str(POPULATED_IMAGE),
        specmode="mfs",
        imsize=[IMAGE_SIZE, IMAGE_SIZE],
        cell=CELL,
        weighting="natural",
        robust=0.5,
        niter=0,
        pblimit=PB_LIMIT,
        interactive=False,
    )


def read_plane(path):
    ia = image()
    try:
        ia.open(str(path))
        pixels = np.squeeze(np.asarray(ia.getchunk(), dtype=np.float64))
        valid = np.squeeze(np.asarray(ia.getchunk(getmask=True), dtype=bool))
    finally:
        ia.done()
    if pixels.shape != (IMAGE_SIZE, IMAGE_SIZE):
        raise RuntimeError(f"Unexpected image shape: {pixels.shape}")
    return np.ascontiguousarray(pixels), np.ascontiguousarray(valid)


def main():
    blank, blank_valid = read_plane(str(BLANK_IMAGE) + ".image")
    blank_pb, blank_pb_valid = read_plane(str(BLANK_IMAGE) + ".pb")
    inject_population(blank_pb)
    make_dirty_image()

    populated, populated_valid = read_plane(str(POPULATED_IMAGE) + ".image")
    pb, pb_valid = read_plane(str(POPULATED_IMAGE) + ".pb")
    valid = blank_valid & populated_valid & pb_valid & np.isfinite(pb) & (pb >= PB_LIMIT)
    detection_mask = ~valid
    blank_for_extraction = np.where(valid, blank, 0.0)
    populated_for_extraction = np.where(valid, populated, 0.0)
    blank_rms = float(imstat(str(BLANK_IMAGE) + ".image", listit=False)["rms"][0])

    common = dict(
        rms=blank_rms,
        seed_sigma=SEED_SIGMA,
        flood_sigma=FLOOD_SIGMA,
        min_pixels=MIN_PIXELS,
        detection_mask=detection_mask,
        exclude_border=True,
    )
    detected = extract_sources(populated_for_extraction, **common)
    blank_positive = extract_sources(blank_for_extraction, **common)
    blank_negative = extract_sources(
        blank_for_extraction, polarity="negative", **common
    )
    matches = match_catalogue_to_positions(
        SOURCE_POSITIONS,
        detected.catalogue,
        max_distance_pixels=MAX_MATCH_DISTANCE_PIXELS,
    )

    recovered = []
    for match in matches:
        source = detected.catalogue[match.detection_index]
        position = SOURCE_POSITIONS[match.input_index]
        apparent_flux_jy = SOURCE_FLUX_JY * float(blank_pb[position])
        response_jy = float(populated[position] - blank[position])
        recovered.append(
            {
                "input_index": int(match.input_index),
                "input_position": [int(position[0]), int(position[1])],
                "direction": direction_for_pixel(position),
                "pb": float(pb[position]),
                "injected_apparent_flux_jy": apparent_flux_jy,
                "peak_jy_per_beam": float(source.peak_value),
                "peak_snr": float(source.peak_value / blank_rms),
                "noise_cancelled_response_jy_per_beam": response_jy,
                # Legacy serialized column name retained for reproducibility.
                "response_fidelity": response_jy / apparent_flux_jy,
                "offset_pixels": float(match.distance_pixels),
            }
        )

    summary = {
        "status": "PASS" if len(matches) == len(SOURCE_POSITIONS) else "FAIL",
        "blank_ms_preserved": str(BLANK_MS),
        "injected_ms": str(POPULATED_MS),
        "source_flux_jy": SOURCE_FLUX_JY,
        "n_injected": len(SOURCE_POSITIONS),
        "n_detected_positive": len(detected.catalogue),
        "n_matched": len(matches),
        "n_unmatched_positive": len(detected.catalogue) - len(matches),
        "n_blank_positive": len(blank_positive.catalogue),
        "n_blank_negative": len(blank_negative.catalogue),
        "blank_rms_jy_per_beam": blank_rms,
        "recovered": recovered,
    }
    with SUMMARY.open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    if summary["status"] != "PASS":
        raise RuntimeError("Not all injected F1 sources were recovered")
    if not all(abs(row["response_fidelity"] - 1.0) < 0.02 for row in recovered):
        raise RuntimeError("Injected F1 source response differs from unity by more than 2%")
    print("PASS: GHOST F1 visibility injection and recovery succeeded")


if __name__ == "__main__":
    main()
