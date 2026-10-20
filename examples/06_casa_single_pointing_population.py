"""Simulate and recover a uniform source population in one ALMA pointing.

This example samples different primary-beam responses efficiently:

1. place equal-flux point sources on a uniform grid within one pointing;
2. run one populated and one blank noisy ``simobserve`` simulation;
3. image each MS once;
4. extract every source in the populated image at once;
5. match detections to the input positions;
6. report completeness as a function of the measured PB response.

Run from the repository root with CASA 6:

    casa --nogui --nologger --log2term \
        -c examples/06_casa_single_pointing_population.py
"""

from pathlib import Path
import csv
import glob
import math
import os
import shutil
import sys

import numpy as np

from casatasks import imstat, simobserve, tclean
from casatools import coordsys, image, measures


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from bliss_alma.extraction import extract_sources
from bliss_alma.mockdata import uniform_grid_positions
from bliss_alma.recovery import match_catalogue_to_positions


WORKDIR = REPO_ROOT / "data" / "casa_single_pointing_population"
PHASE_CENTRE = "J2000 03h32m39.0s -27d46m00.0s"
CENTRE_FREQUENCY = "230GHz"
CONTINUUM_BANDWIDTH = "2GHz"
ANTENNA_LIST = "alma.cycle8.1.cfg"
CELL_ARCSEC = 0.2
CELL = "{0}arcsec".format(CELL_ARCSEC)
IMAGE_SIZE = 256
IMAGE_CENTRE = (IMAGE_SIZE // 2, IMAGE_SIZE // 2)
SOURCE_SPACING_ARCSEC = 8.0
SOURCE_RADIUS_ARCSEC = 18.0
INTRINSIC_FLUX_JY = 0.001
SEED_SIGMA = 5.0
FLOOD_SIGMA = 3.0
MIN_PIXELS = 5
PB_LIMIT = 0.2
MAX_MATCH_DISTANCE_PIXELS = 3.0
NOISE_SEED = 11111


def reset_workdir():
    if WORKDIR.exists():
        shutil.rmtree(WORKDIR)
    WORKDIR.mkdir(parents=True)


def create_sky_model(path, source_positions):
    """Create a valid one-channel CASA model containing pixel point sources."""
    cs_tool = coordsys()
    cs = cs_tool.newcoordsys(direction=True, stokes=["I"], spectral=True)
    me = measures()

    direction = me.direction("J2000", "03h32m39.0s", "-27d46m00.0s")
    cs.setreferencevalue([direction["m0"], direction["m1"]], type="direction")
    cs.setreferencepixel(list(IMAGE_CENTRE), type="direction")
    cs.setincrement(["-{0}".format(CELL), CELL], type="direction")
    cs.setreferencevalue(CENTRE_FREQUENCY, type="spectral")
    cs.setincrement(CONTINUUM_BANDWIDTH, type="spectral")

    ia = image()
    try:
        ia.fromshape(
            outfile=str(path),
            shape=[IMAGE_SIZE, IMAGE_SIZE, 1, 1],
            csys=cs.torecord(),
            overwrite=True,
        )
        pixels = np.zeros((IMAGE_SIZE, IMAGE_SIZE, 1, 1), dtype=np.float32)
        for x_pixel, y_pixel in source_positions:
            pixels[x_pixel, y_pixel, 0, 0] = INTRINSIC_FLUX_JY
        ia.putchunk(pixels)
        ia.setbrightnessunit("Jy/pixel")
    finally:
        ia.done()
        cs.done()
        cs_tool.done()
        me.done()


def write_single_pointing_file(path):
    path.write_text(
        "#Epoch RA DEC TIME(optional)\n"
        "J2000 03h32m39.0s -27d46m00.0s 10.0\n"
    )


def simulate(project, sky_model, pointing_file):
    simobserve(
        project=project,
        skymodel=str(sky_model),
        inbright="",
        indirection=PHASE_CENTRE,
        incell=CELL,
        incenter=CENTRE_FREQUENCY,
        inwidth=CONTINUUM_BANDWIDTH,
        setpointings=False,
        ptgfile=str(pointing_file),
        integration="10s",
        totaltime="120s",
        antennalist=ANTENNA_LIST,
        obsmode="int",
        thermalnoise="tsys-atm",
        user_pwv=0.5,
        seed=NOISE_SEED,
        overwrite=True,
        graphics="none",
    )


def find_noisy_ms(project):
    matches = glob.glob(str(WORKDIR / project / "*.noisy.ms"))
    if len(matches) != 1:
        raise RuntimeError("Expected one noisy MS for {0}: {1}".format(project, matches))
    return Path(matches[0])


def image_continuum(vis, imagename):
    tclean(
        vis=str(vis),
        imagename=str(imagename),
        specmode="mfs",
        imsize=[IMAGE_SIZE, IMAGE_SIZE],
        cell=CELL,
        phasecenter=PHASE_CENTRE,
        weighting="natural",
        niter=0,
        pblimit=0.1,
        interactive=False,
    )


def read_casa_plane(path):
    ia = image()
    try:
        ia.open(str(path))
        pixels = np.asarray(ia.getchunk(), dtype=np.float64)
    finally:
        ia.done()
    pixels = np.squeeze(pixels)
    if pixels.ndim != 2:
        raise RuntimeError("Expected a 2D CASA plane: {0}".format(pixels.shape))
    return np.ascontiguousarray(pixels)


def summarize_pb_bins(rows):
    edges = [0.2, 0.4, 0.6, 0.8, 1.01]
    print("Completeness by primary-beam response")
    for lower, upper in zip(edges[:-1], edges[1:]):
        selected = [row for row in rows if lower <= row["pb"] < upper]
        if not selected:
            continue
        recovered = sum(row["recovered"] for row in selected)
        print(
            "  {0:.1f} <= PB < {1:.1f}: {2}/{3} = {4:.3f}".format(
                lower,
                upper,
                recovered,
                len(selected),
                recovered / float(len(selected)),
            )
        )


def write_recovery_table(rows):
    output = WORKDIR / "source_recovery.csv"
    fieldnames = [
        "source_id",
        "x_pixel",
        "y_pixel",
        "radius_arcsec",
        "pb",
        "intrinsic_flux_jy",
        "expected_attenuated_flux_jy",
        "expected_snr",
        "recovered",
        "recovered_peak_jy_per_beam",
        "recovered_peak_snr",
        "offset_pixels",
    ]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return output


def write_detection_table(filename, catalogue, pb, rms, matched_indices=None):
    """Write a compact detection catalogue for later Q calculations."""
    matched_indices = set() if matched_indices is None else set(matched_indices)
    output = WORKDIR / filename
    fieldnames = [
        "detection_id",
        "x_pixel",
        "y_pixel",
        "pb",
        "peak_jy_per_beam",
        "peak_snr",
        "n_pixels",
        "matched",
    ]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for detection_id, source in enumerate(catalogue):
            x_pixel, y_pixel = source.peak_index
            writer.writerow(
                {
                    "detection_id": detection_id,
                    "x_pixel": x_pixel,
                    "y_pixel": y_pixel,
                    "pb": float(pb[x_pixel, y_pixel]),
                    "peak_jy_per_beam": source.peak_value,
                    "peak_snr": source.peak_value / rms,
                    "n_pixels": source.n_pixels,
                    "matched": detection_id in matched_indices,
                }
            )
    return output


def main():
    reset_workdir()
    spacing_pixels = int(round(SOURCE_SPACING_ARCSEC / CELL_ARCSEC))
    radius_pixels = int(round(SOURCE_RADIUS_ARCSEC / CELL_ARCSEC))
    source_positions = uniform_grid_positions(
        IMAGE_CENTRE,
        spacing_pixels=spacing_pixels,
        radius_pixels=radius_pixels,
    )

    blank_model = WORKDIR / "blank_sky.image"
    populated_model = WORKDIR / "populated_sky.image"
    pointing_file = WORKDIR / "single_pointing.txt"
    create_sky_model(blank_model, [])
    create_sky_model(populated_model, source_positions)
    write_single_pointing_file(pointing_file)

    os.chdir(WORKDIR)
    simulate("blank", blank_model, pointing_file)
    simulate("population", populated_model, pointing_file)

    blank_ms = find_noisy_ms("blank")
    population_ms = find_noisy_ms("population")
    blank_image = WORKDIR / "blank_dirty"
    population_image = WORKDIR / "population_dirty"
    image_continuum(blank_ms, blank_image)
    image_continuum(population_ms, population_image)

    blank = read_casa_plane(str(blank_image) + ".image")
    population = read_casa_plane(str(population_image) + ".image")
    pb = read_casa_plane(str(population_image) + ".pb")
    blank_rms = float(imstat(imagename=str(blank_image) + ".image")["rms"][0])
    detection_mask = (~np.isfinite(pb)) | (pb < PB_LIMIT)

    extraction = extract_sources(
        population,
        rms=blank_rms,
        seed_sigma=SEED_SIGMA,
        flood_sigma=FLOOD_SIGMA,
        min_pixels=MIN_PIXELS,
        detection_mask=detection_mask,
        exclude_border=True,
    )
    blank_positive = extract_sources(
        blank,
        rms=blank_rms,
        seed_sigma=SEED_SIGMA,
        flood_sigma=FLOOD_SIGMA,
        min_pixels=MIN_PIXELS,
        detection_mask=detection_mask,
        exclude_border=True,
    )
    blank_negative = extract_sources(
        blank,
        rms=blank_rms,
        seed_sigma=SEED_SIGMA,
        flood_sigma=FLOOD_SIGMA,
        min_pixels=MIN_PIXELS,
        detection_mask=detection_mask,
        exclude_border=True,
        polarity="negative",
    )

    matches = match_catalogue_to_positions(
        source_positions,
        extraction.catalogue,
        max_distance_pixels=MAX_MATCH_DISTANCE_PIXELS,
    )
    match_by_input = {match.input_index: match for match in matches}
    matched_detection_indices = {match.detection_index for match in matches}

    rows = []
    for source_id, position in enumerate(source_positions):
        x_pixel, y_pixel = position
        radius_arcsec = CELL_ARCSEC * math.hypot(
            x_pixel - IMAGE_CENTRE[0],
            y_pixel - IMAGE_CENTRE[1],
        )
        pb_value = float(pb[x_pixel, y_pixel])
        match = match_by_input.get(source_id)
        recovered = match is not None
        recovered_peak = ""
        offset_pixels = ""
        if recovered:
            detection = extraction.catalogue[match.detection_index]
            recovered_peak = detection.peak_value
            offset_pixels = match.distance_pixels
        rows.append(
            {
                "source_id": source_id,
                "x_pixel": x_pixel,
                "y_pixel": y_pixel,
                "radius_arcsec": radius_arcsec,
                "pb": pb_value,
                "intrinsic_flux_jy": INTRINSIC_FLUX_JY,
                "expected_attenuated_flux_jy": INTRINSIC_FLUX_JY * pb_value,
                "expected_snr": INTRINSIC_FLUX_JY * pb_value / blank_rms,
                "recovered": recovered,
                "recovered_peak_jy_per_beam": recovered_peak,
                "recovered_peak_snr": (
                    recovered_peak / blank_rms if recovered else ""
                ),
                "offset_pixels": offset_pixels,
            }
        )

    table_path = write_recovery_table(rows)
    unmatched_detections = len(extraction.catalogue) - len(matched_detection_indices)
    populated_catalogue_path = write_detection_table(
        "populated_detections.csv",
        extraction.catalogue,
        pb,
        blank_rms,
        matched_detection_indices,
    )
    blank_positive_path = write_detection_table(
        "blank_positive_detections.csv", blank_positive.catalogue, pb, blank_rms
    )
    blank_negative_path = write_detection_table(
        "blank_negative_detections.csv", blank_negative.catalogue, pb, blank_rms
    )
    print("BliSS-ALMA single-pointing population result")
    print("  Input sources: {0}".format(len(source_positions)))
    print("  Extracted positive islands: {0}".format(len(extraction.catalogue)))
    print("  Matched sources: {0}".format(len(matches)))
    print("  Unmatched positive islands: {0}".format(unmatched_detections))
    print("  Blank positive islands: {0}".format(len(blank_positive.catalogue)))
    print("  Blank negative islands: {0}".format(len(blank_negative.catalogue)))
    print("  Blank RMS: {0:.6g} Jy/beam".format(blank_rms))
    print("  Recovery table: {0}".format(table_path))
    summarize_pb_bins(rows)

    np.save(str(WORKDIR / "population.npy"), population)
    np.save(str(WORKDIR / "blank.npy"), blank)
    np.save(str(WORKDIR / "primary_beam.npy"), pb)
    np.save(str(WORKDIR / "source_labels.npy"), extraction.labels)

    if not matches:
        raise RuntimeError("No injected sources were recovered.")
    print("PASS: the single-pointing source population was extracted and matched.")
    return {
        "noise_seed": NOISE_SEED,
        "n_injected": len(source_positions),
        "n_detected": len(extraction.catalogue),
        "n_matched": len(matches),
        "n_unmatched": unmatched_detections,
        "n_blank_positive": len(blank_positive.catalogue),
        "n_blank_negative": len(blank_negative.catalogue),
        "blank_rms_jy_per_beam": blank_rms,
        "recovery_table": str(table_path),
        "populated_catalogue": str(populated_catalogue_path),
        "blank_positive_catalogue": str(blank_positive_path),
        "blank_negative_catalogue": str(blank_negative_path),
    }


if __name__ == "__main__":
    main()
