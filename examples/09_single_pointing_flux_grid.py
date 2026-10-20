"""Create a continuum flux-grid smoke test from existing dirty images.

Dirty imaging is linear. For each existing noise realization this script
separates the source response from the blank noise image and rescales that
response over a grid of intrinsic source fluxes. No new measurement sets or
CASA images are created.
"""

from pathlib import Path
import csv
import sys

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from bliss_alma.extraction import extract_sources
from bliss_alma.recovery import match_catalogue_to_positions


INPUT = REPO_ROOT / "data" / "casa_single_pointing_monte_carlo"
OUTPUT = REPO_ROOT / "data" / "casa_single_pointing_flux_grid"
BASE_FLUX_JY = 0.001
FLUX_GRID_MJY = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
PB_LIMIT = 0.2
SEED_SIGMA = 5.0
FLOOD_SIGMA = 3.0
MIN_PIXELS = 5
MAX_MATCH_DISTANCE_PIXELS = 3.0


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    if not rows:
        raise RuntimeError("No rows available for {0}".format(path))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def detection_rows(
    realization,
    noise_seed,
    flux_mjy,
    kind,
    catalogue,
    pb,
    rms,
    matched_indices=None,
):
    matched_indices = set() if matched_indices is None else set(matched_indices)
    rows = []
    for detection_id, source in enumerate(catalogue):
        x_pixel, y_pixel = source.peak_index
        rows.append(
            {
                "realization": realization,
                "noise_seed": noise_seed,
                "flux_mjy": flux_mjy,
                "detection_kind": kind,
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
    return rows


def main():
    summary_rows = read_csv(INPUT / "realization_summary.csv")
    all_recovery = []
    all_detections = []
    grid_summary = []

    for summary in summary_rows:
        realization = int(summary["realization"])
        noise_seed = int(summary["noise_seed"])
        rms = float(summary["blank_rms_jy_per_beam"])
        root = INPUT / "realization_{0:03d}".format(realization)
        blank = np.load(root / "blank.npy")
        base_population = np.load(root / "population.npy")
        source_response = base_population - blank
        pb = np.load(root / "primary_beam.npy")
        input_rows = read_csv(root / "source_recovery.csv")
        source_positions = [
            (int(row["x_pixel"]), int(row["y_pixel"])) for row in input_rows
        ]
        detection_mask = (~np.isfinite(pb)) | (pb < PB_LIMIT)

        blank_positive = extract_sources(
            blank,
            rms=rms,
            seed_sigma=SEED_SIGMA,
            flood_sigma=FLOOD_SIGMA,
            min_pixels=MIN_PIXELS,
            detection_mask=detection_mask,
            exclude_border=True,
        )

        for flux_mjy in FLUX_GRID_MJY:
            flux_jy = flux_mjy * 1.0e-3
            populated = blank + (flux_jy / BASE_FLUX_JY) * source_response
            extracted = extract_sources(
                populated,
                rms=rms,
                seed_sigma=SEED_SIGMA,
                flood_sigma=FLOOD_SIGMA,
                min_pixels=MIN_PIXELS,
                detection_mask=detection_mask,
                exclude_border=True,
            )
            matches = match_catalogue_to_positions(
                source_positions,
                extracted.catalogue,
                max_distance_pixels=MAX_MATCH_DISTANCE_PIXELS,
            )
            by_input = {match.input_index: match for match in matches}
            matched_detection_indices = {
                match.detection_index for match in matches
            }

            for source_id, row in enumerate(input_rows):
                match = by_input.get(source_id)
                recovered_peak = ""
                recovered_peak_snr = ""
                offset_pixels = ""
                if match is not None:
                    detection = extracted.catalogue[match.detection_index]
                    recovered_peak = detection.peak_value
                    recovered_peak_snr = detection.peak_value / rms
                    offset_pixels = match.distance_pixels
                pb_value = float(row["pb"])
                all_recovery.append(
                    {
                        "realization": realization,
                        "noise_seed": noise_seed,
                        "flux_mjy": flux_mjy,
                        "source_id": source_id,
                        "x_pixel": row["x_pixel"],
                        "y_pixel": row["y_pixel"],
                        "pb": pb_value,
                        "intrinsic_flux_jy": flux_jy,
                        "expected_attenuated_flux_jy": flux_jy * pb_value,
                        "expected_snr": flux_jy * pb_value / rms,
                        "recovered": match is not None,
                        "recovered_peak_jy_per_beam": recovered_peak,
                        "recovered_peak_snr": recovered_peak_snr,
                        "offset_pixels": offset_pixels,
                    }
                )

            all_detections.extend(
                detection_rows(
                    realization,
                    noise_seed,
                    flux_mjy,
                    "populated",
                    extracted.catalogue,
                    pb,
                    rms,
                    matched_detection_indices,
                )
            )
            all_detections.extend(
                detection_rows(
                    realization,
                    noise_seed,
                    flux_mjy,
                    "blank_positive",
                    blank_positive.catalogue,
                    pb,
                    rms,
                )
            )
            grid_summary.append(
                {
                    "realization": realization,
                    "noise_seed": noise_seed,
                    "flux_mjy": flux_mjy,
                    "n_injected": len(source_positions),
                    "n_detected": len(extracted.catalogue),
                    "n_matched": len(matches),
                    "n_unmatched": len(extracted.catalogue) - len(matches),
                    "n_blank_positive": len(blank_positive.catalogue),
                    "blank_rms_jy_per_beam": rms,
                }
            )

    OUTPUT.mkdir(parents=True, exist_ok=True)
    write_csv(OUTPUT / "source_recovery_all.csv", all_recovery)
    write_csv(OUTPUT / "detections_all.csv", all_detections)
    write_csv(OUTPUT / "realization_summary.csv", grid_summary)
    print("Flux grid [mJy]: {0}".format(FLUX_GRID_MJY))
    print("Source trials: {0}".format(len(all_recovery)))
    print("Wrote summary tables to {0}".format(OUTPUT))


if __name__ == "__main__":
    main()
