"""Run the Hayatsu-style 2D continuum S/N-PB grid on GHOST F1.

Ten deterministic batches each contain one source for every combination of
11 target S/N values and four PB bins (44 simultaneous sources per batch,
440 injections total). Intermediate MS and CASA image products are streamed
and removed after each successful batch.
"""

from pathlib import Path
import argparse
import csv
import json
import math
import shutil
import sys

import numpy as np
from casatasks import ft, imstat, tclean, uvsub
from casatools import componentlist, image


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from bliss_alma.extraction import extract_sources
from bliss_alma.recovery import match_catalogue_to_positions


WORKDIR = REPO_ROOT / "data" / "ghost_f1_5min"
OUTPUT_DIR = WORKDIR / "f1_snr_pb_grid"
BLANK_MS = WORKDIR / "output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms"
BLANK_IMAGE = WORKDIR / "bliss_f1_dirty.image"
BLANK_PB = WORKDIR / "bliss_f1_dirty.pb"
WORKING_MS = OUTPUT_DIR / "working_injected.ms"
WORKING_COMPONENTS = OUTPUT_DIR / "working_sources.cl"
WORKING_IMAGE = OUTPUT_DIR / "working_dirty"

IMAGE_SIZE = 360
IMAGE_CENTRE = (IMAGE_SIZE // 2, IMAGE_SIZE // 2)
CELL_ARCSEC = 0.12454065519280681
CELL = f"{CELL_ARCSEC}arcsec"
PHASE_CENTRE_RA_DEG = 180.0
PHASE_CENTRE_DEC_DEG = -27.8
TARGET_SNRS = np.arange(2.0, 7.01, 0.5)
PB_BINS = [(0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.01)]
N_BATCHES = 10
MIN_SEPARATION_PIXELS = 20.0
POSITION_MARGIN_PIXELS = 10
SEED_SIGMA = 4.0
FLOOD_SIGMA = 3.0
MIN_PIXELS = 5
PB_LIMIT = 0.2
MAX_MATCH_DISTANCE_PIXELS = 4.0


def remove_path(path):
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def remove_working_products():
    remove_path(WORKING_MS)
    remove_path(WORKING_COMPONENTS)
    for path in OUTPUT_DIR.glob(WORKING_IMAGE.name + ".*"):
        remove_path(path)


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


def direction_for_pixel(position):
    x_pixel, y_pixel = position
    east_arcsec = -(x_pixel - IMAGE_CENTRE[0]) * CELL_ARCSEC
    north_arcsec = (y_pixel - IMAGE_CENTRE[1]) * CELL_ARCSEC
    ra_deg = PHASE_CENTRE_RA_DEG + east_arcsec / (
        3600.0 * math.cos(math.radians(PHASE_CENTRE_DEC_DEG))
    )
    dec_deg = PHASE_CENTRE_DEC_DEG + north_arcsec / 3600.0
    return f"J2000 {ra_deg}deg {dec_deg}deg"


def choose_positions(pb, valid, batch):
    selected = []
    by_bin = {}
    x_grid, y_grid = np.indices(pb.shape)
    interior = (
        valid
        & (x_grid >= POSITION_MARGIN_PIXELS)
        & (x_grid < IMAGE_SIZE - POSITION_MARGIN_PIXELS)
        & (y_grid >= POSITION_MARGIN_PIXELS)
        & (y_grid < IMAGE_SIZE - POSITION_MARGIN_PIXELS)
    )
    # Walk from the central PB bin outward along a reproducible spiral. Rotate
    # the whole pattern between batches so sky direction is still sampled.
    spiral_index = 0
    angle_step = math.radians(40.0)
    batch_rotation = math.radians(13.0 * batch)
    for pb_index in reversed(range(len(PB_BINS))):
        lower, upper = PB_BINS[pb_index]
        candidates = np.argwhere(interior & (pb >= lower) & (pb < upper))
        radii = np.hypot(candidates[:, 0] - IMAGE_CENTRE[0], candidates[:, 1] - IMAGE_CENTRE[1])
        radius_targets = np.quantile(radii, np.linspace(0.12, 0.88, len(TARGET_SNRS)))
        chosen = []
        for radius_target in radius_targets:
            angle = batch_rotation + spiral_index * angle_step
            target = np.array([
                IMAGE_CENTRE[0] + radius_target * math.cos(angle),
                IMAGE_CENTRE[1] + radius_target * math.sin(angle),
            ])
            order = np.argsort(np.sum((candidates - target) ** 2, axis=1))
            position = None
            for candidate_index in order:
                candidate = candidates[candidate_index]
                trial = (int(candidate[0]), int(candidate[1]))
                if all(math.hypot(trial[0] - other[0], trial[1] - other[1])
                       >= MIN_SEPARATION_PIXELS for other in selected):
                    position = trial
                    break
            if position is None:
                raise RuntimeError(f"Batch {batch}: spiral separation failed in PB {lower}-{upper}")
            chosen.append(position)
            selected.append(position)
            spiral_index += 1
        if len(chosen) != len(TARGET_SNRS):
            raise RuntimeError(
                f"Batch {batch}: found only {len(chosen)} positions in PB {lower}-{upper}"
            )
        by_bin[pb_index] = (lower, upper, chosen)
    return [by_bin[index] for index in range(len(PB_BINS))]


def create_components(source_rows):
    remove_path(WORKING_COMPONENTS)
    cl = componentlist()
    try:
        for row in source_rows:
            cl.addcomponent(
                dir=row["direction"],
                flux=row["apparent_flux_jy"],
                fluxunit="Jy",
                shape="point",
            )
        cl.rename(str(WORKING_COMPONENTS))
    finally:
        cl.done()


def inject_and_image(source_rows):
    remove_working_products()
    shutil.copytree(BLANK_MS, WORKING_MS)
    create_components(source_rows)
    ft(
        vis=str(WORKING_MS),
        spw="0:0",
        complist=str(WORKING_COMPONENTS),
        usescratch=True,
    )
    uvsub(vis=str(WORKING_MS), reverse=True)
    tclean(
        vis=str(WORKING_MS),
        imagename=str(WORKING_IMAGE),
        specmode="mfs",
        imsize=[IMAGE_SIZE, IMAGE_SIZE],
        cell=CELL,
        weighting="natural",
        robust=0.5,
        niter=0,
        pblimit=PB_LIMIT,
        interactive=False,
    )


def write_csv(path, rows, fieldnames):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    blank, blank_valid = read_plane(BLANK_IMAGE)
    pb, pb_valid = read_plane(BLANK_PB)
    usable = blank_valid & pb_valid & np.isfinite(blank) & np.isfinite(pb) & (pb >= PB_LIMIT)
    blank_rms = float(imstat(str(BLANK_IMAGE), listit=False)["rms"][0])
    blank_for_extraction = np.where(usable, blank, 0.0)
    detection_mask = ~usable
    common = dict(
        rms=blank_rms,
        seed_sigma=SEED_SIGMA,
        flood_sigma=FLOOD_SIGMA,
        min_pixels=MIN_PIXELS,
        detection_mask=detection_mask,
        exclude_border=True,
    )
    blank_positive = extract_sources(blank_for_extraction, **common)
    blank_negative = extract_sources(blank_for_extraction, polarity="negative", **common)

    all_source_rows = []
    all_detection_rows = []
    batch_summaries = []
    for batch in range(N_BATCHES):
        positions_by_bin = choose_positions(pb, usable, batch)
        source_rows = []
        source_id = 0
        for pb_index, (lower, upper, positions) in enumerate(positions_by_bin):
            # Rotate the S/N-to-position assignment between batches.
            snrs = np.roll(TARGET_SNRS, batch % len(TARGET_SNRS))
            for target_snr, position in zip(snrs, positions):
                pb_value = float(pb[position])
                apparent_flux_jy = float(target_snr * blank_rms)
                source_rows.append(
                    {
                        "realization": batch,
                        "source_id": source_id,
                        "pb_bin": pb_index,
                        "pb_lower": lower,
                        "pb_upper": upper,
                        "x_pixel": position[0],
                        "y_pixel": position[1],
                        "direction": direction_for_pixel(position),
                        "pb": pb_value,
                        "expected_snr": float(target_snr),
                        "apparent_flux_jy": apparent_flux_jy,
                        "intrinsic_flux_jy": apparent_flux_jy / pb_value,
                    }
                )
                source_id += 1

        print(f"Batch {batch + 1}/{N_BATCHES}: injecting {len(source_rows)} sources", flush=True)
        inject_and_image(source_rows)
        populated, populated_valid = read_plane(str(WORKING_IMAGE) + ".image")
        populated_pb, populated_pb_valid = read_plane(str(WORKING_IMAGE) + ".pb")
        batch_usable = usable & populated_valid & populated_pb_valid
        batch_mask = ~batch_usable
        populated_for_extraction = np.where(batch_usable, populated, 0.0)
        detected = extract_sources(
            populated_for_extraction,
            rms=blank_rms,
            seed_sigma=SEED_SIGMA,
            flood_sigma=FLOOD_SIGMA,
            min_pixels=MIN_PIXELS,
            detection_mask=batch_mask,
            exclude_border=True,
        )
        input_positions = [(row["x_pixel"], row["y_pixel"]) for row in source_rows]
        matches = match_catalogue_to_positions(
            input_positions,
            detected.catalogue,
            max_distance_pixels=MAX_MATCH_DISTANCE_PIXELS,
        )
        by_input = {match.input_index: match for match in matches}
        matched_detection_ids = {match.detection_index for match in matches}

        for input_index, row in enumerate(source_rows):
            match = by_input.get(input_index)
            recovered = match is not None
            detection = detected.catalogue[match.detection_index] if recovered else None
            position = (row["x_pixel"], row["y_pixel"])
            output_row = dict(row)
            output_row.update(
                {
                    "recovered": recovered,
                    "recovered_peak_jy_per_beam": (
                        float(detection.peak_value) if recovered else ""
                    ),
                    "recovered_peak_snr": (
                        float(detection.peak_value / blank_rms) if recovered else ""
                    ),
                    "offset_pixels": float(match.distance_pixels) if recovered else "",
                    "noise_cancelled_response_jy_per_beam": float(
                        populated[position] - blank[position]
                    ),
                }
            )
            all_source_rows.append(output_row)

        for detection_id, detection in enumerate(detected.catalogue):
            x_pixel, y_pixel = detection.peak_index
            all_detection_rows.append(
                {
                    "realization": batch,
                    "detection_kind": "populated",
                    "detection_id": detection_id,
                    "x_pixel": x_pixel,
                    "y_pixel": y_pixel,
                    "pb": float(populated_pb[x_pixel, y_pixel]),
                    "peak_jy_per_beam": float(detection.peak_value),
                    "peak_snr": float(detection.peak_value / blank_rms),
                    "n_pixels": detection.n_pixels,
                    "matched": detection_id in matched_detection_ids,
                }
            )
        for kind, catalogue in (
            ("blank_positive", blank_positive.catalogue),
            ("blank_negative", blank_negative.catalogue),
        ):
            for detection_id, detection in enumerate(catalogue):
                x_pixel, y_pixel = detection.peak_index
                all_detection_rows.append(
                    {
                        "realization": batch,
                        "detection_kind": kind,
                        "detection_id": detection_id,
                        "x_pixel": x_pixel,
                        "y_pixel": y_pixel,
                        "pb": float(pb[x_pixel, y_pixel]),
                        "peak_jy_per_beam": float(detection.peak_value),
                        "peak_snr": float(detection.peak_value / blank_rms),
                        "n_pixels": detection.n_pixels,
                        "matched": False,
                    }
                )

        batch_summaries.append(
            {
                "realization": batch,
                "n_injected": len(source_rows),
                "n_detected_positive": len(detected.catalogue),
                "n_matched": len(matches),
                "n_unmatched_positive": len(detected.catalogue) - len(matches),
                "n_blank_positive": len(blank_positive.catalogue),
                "n_blank_negative": len(blank_negative.catalogue),
            }
        )
        print(
            f"Batch {batch + 1}: matched {len(matches)}/{len(source_rows)}, "
            f"unmatched {len(detected.catalogue) - len(matches)}",
            flush=True,
        )
        remove_working_products()

    source_fields = list(all_source_rows[0])
    detection_fields = list(all_detection_rows[0])
    write_csv(OUTPUT_DIR / "source_recovery_all.csv", all_source_rows, source_fields)
    write_csv(OUTPUT_DIR / "detections_all.csv", all_detection_rows, detection_fields)
    write_csv(
        OUTPUT_DIR / "realization_summary.csv",
        batch_summaries,
        list(batch_summaries[0]),
    )
    summary = {
        "status": "PASS",
        "n_batches": N_BATCHES,
        "sources_per_batch": len(TARGET_SNRS) * len(PB_BINS),
        "n_injected": len(all_source_rows),
        "n_recovered": sum(str(row["recovered"]).lower() == "true" for row in all_source_rows),
        "n_populated_detections": sum(
            row["detection_kind"] == "populated" for row in all_detection_rows
        ),
        "n_unmatched_populated": sum(
            row["detection_kind"] == "populated" and not row["matched"]
            for row in all_detection_rows
        ),
        "blank_positive_per_batch": len(blank_positive.catalogue),
        "blank_negative_per_batch": len(blank_negative.catalogue),
        "blank_rms_jy_per_beam": blank_rms,
        "seed_sigma": SEED_SIGMA,
        "flood_sigma": FLOOD_SIGMA,
        "pb_bins": PB_BINS,
        "target_snrs": [float(value) for value in TARGET_SNRS],
        "note": "Blank catalogues repeat the same F1 noise realization in each batch.",
    }
    with (OUTPUT_DIR / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed-sigma", type=float, default=SEED_SIGMA)
    parser.add_argument("--flood-sigma", type=float, default=FLOOD_SIGMA)
    parser.add_argument("--output-name", default=OUTPUT_DIR.name)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.flood_sigma > args.seed_sigma:
        raise ValueError("flood sigma must not exceed seed sigma")
    SEED_SIGMA = args.seed_sigma
    FLOOD_SIGMA = args.flood_sigma
    OUTPUT_DIR = WORKDIR / args.output_name
    WORKING_MS = OUTPUT_DIR / "working_injected.ms"
    WORKING_COMPONENTS = OUTPUT_DIR / "working_sources.cl"
    WORKING_IMAGE = OUTPUT_DIR / "working_dirty"
    main()
