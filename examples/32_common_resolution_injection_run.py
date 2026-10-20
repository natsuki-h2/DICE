"""Run common-resolution injection/recovery for one GHOST configuration.

The visibilities are imaged with a 0.7-arcsec uv taper on the full-PB F1
grid, then convolved to a circular 1.0-arcsec restoring beam.  Use one
realization as a pilot before requesting the ten-realization production run.
"""

from pathlib import Path
import argparse
import atexit
import csv
import importlib.util
import json
import shutil
import sys

import numpy as np
from casatasks import ft, imsmooth, imstat, tclean, uvsub
from casatools import componentlist, table


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
from bliss_alma.extraction import extract_sources
from bliss_alma.recovery import match_catalogue_to_positions

spec = importlib.util.spec_from_file_location(
    "f1grid", ROOT / "examples" / "15_ghost_f1_snr_pb_grid.py")
grid = importlib.util.module_from_spec(spec)
spec.loader.exec_module(grid)

BASE_MS = {
    ("F1", "5m"): ROOT / "data" / "ghost_f1_5min" /
          "output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms",
    **{
        (f"F{index}", "5m"): ROOT / "data" / "ghost_matrix" / f"f{index}_5min" /
        "extracted" / "output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms"
        for index in range(2, 7)
    },
    **{
        (configuration, integration): ROOT / "data" / "ghost_matrix" /
        f"{configuration.lower()}_{integration}" / "extracted" /
        f"output_custom_skymodel_SIM_MIXED_{integration}_-27_480000_ALMA.ms"
        for configuration in ("F1", "F6") for integration in ("1h", "5h")
    },
}
COMMON_ROOT = ROOT / "data" / "ghost_matrix" / "common_resolution_runs"
UVTAPER = ["0.7arcsec"]
TARGET_BEAM = "1.0arcsec"
CHUNK_ROWS = 100_000


def remove(path):
    if path.is_dir(): shutil.rmtree(path)
    elif path.exists(): path.unlink()


def remove_images(out, prefix):
    for path in out.glob(prefix.name + ".*"):
        remove(path)


def csv_write(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", choices=[f"F{index}" for index in range(1, 7)], required=True)
    parser.add_argument("--integration", choices=["5m", "1h", "5h"], default="5m")
    parser.add_argument("--realizations", type=int, default=1)
    parser.add_argument("--seed-base", type=int, default=20410000)
    args = parser.parse_args()
    if args.realizations < 1:
        raise ValueError("--realizations must be positive")

    scenario = (args.config, args.integration)
    if scenario not in BASE_MS:
        parser.error(f"No downloaded template for {args.config} at {args.integration}")
    base_ms = BASE_MS[scenario]
    out = COMMON_ROOT / f"{args.config.lower()}_{args.integration}"
    out.mkdir(parents=True, exist_ok=True)
    for stale in out.glob("working*"):
        remove(stale)
    ms = out / "working_00.ms"
    cl_path = out / "working_00.cl"
    blank_prefix = out / "working_blank_00"
    populated_prefix = out / "working_populated_00"

    def cleanup():
        remove(ms); remove(cl_path)
        remove_images(out, blank_prefix); remove_images(out, populated_prefix)
    atexit.register(cleanup)

    # Stream the noise moments so the 5-hour templates do not require loading
    # several gigabytes of complex visibilities into memory at once.
    tb = table(); tb.open(str(base_ms))
    try:
        nrows = tb.nrows()
        sample = tb.getcol("DATA", startrow=0, nrow=1)
        npol = sample.shape[0]
        count_total = 0
        real_sum = np.zeros(npol); real_sumsq = np.zeros(npol)
        imag_sum = np.zeros(npol); imag_sumsq = np.zeros(npol)
        for start in range(0, nrows, CHUNK_ROWS):
            count = min(CHUNK_ROWS, nrows - start)
            data = tb.getcol("DATA", startrow=start, nrow=count)
            real = data.real.reshape(npol, -1)
            imag = data.imag.reshape(npol, -1)
            count_total += real.shape[1]
            real_sum += real.sum(axis=1); real_sumsq += np.square(real).sum(axis=1)
            imag_sum += imag.sum(axis=1); imag_sumsq += np.square(imag).sum(axis=1)
    finally:
        tb.close()
    real_sigma = np.sqrt(real_sumsq / count_total - np.square(real_sum / count_total))
    imag_sigma = np.sqrt(imag_sumsq / count_total - np.square(imag_sum / count_total))

    def regenerate_noise(seed):
        remove(ms); shutil.copytree(base_ms, ms)
        rng = np.random.default_rng(seed)
        tb = table(); tb.open(str(ms), nomodify=False)
        try:
            nrows = tb.nrows()
            for start in range(0, nrows, CHUNK_ROWS):
                count = min(CHUNK_ROWS, nrows - start)
                shape = (len(real_sigma), 1, count)
                values = (rng.normal(size=shape) * real_sigma[:, None, None] +
                          1j * rng.normal(size=shape) * imag_sigma[:, None, None])
                tb.putcol("DATA", values, startrow=start, nrow=count)
                tb.putcol("CORRECTED_DATA", values, startrow=start, nrow=count)
        finally:
            tb.close()

    def image_ms(prefix):
        remove_images(out, prefix)
        tclean(vis=str(ms), imagename=str(prefix), specmode="mfs",
               imsize=[grid.IMAGE_SIZE, grid.IMAGE_SIZE], cell=grid.CELL,
               weighting="natural", uvtaper=UVTAPER, robust=0.5, niter=0,
               pblimit=grid.PB_LIMIT, interactive=False)
        imsmooth(imagename=str(prefix) + ".image",
                 outfile=str(prefix) + ".smooth.image", kernel="gauss",
                 major=TARGET_BEAM, minor=TARGET_BEAM, pa="0deg", targetres=True,
                 overwrite=True)

    def make_components(rows):
        remove(cl_path); cl = componentlist()
        try:
            for row in rows:
                cl.addcomponent(dir=row["direction"], flux=row["apparent_flux_jy"],
                                fluxunit="Jy", shape="point")
            cl.rename(str(cl_path))
        finally:
            cl.done()

    all_sources, all_detections, logs = [], [], []
    config_index = int(args.config[1:])
    time_offset = {"5m": 0, "1h": 10_000, "5h": 20_000}[args.integration]
    for realization in range(args.realizations):
        # CASA can retain table descriptors after a Measurement Set is closed.
        # A unique path per realization prevents a cached descriptor from
        # referring to a directory that was removed and recreated in place.
        ms = out / f"working_{realization:02d}.ms"
        cl_path = out / f"working_{realization:02d}.cl"
        blank_prefix = out / f"working_blank_{realization:02d}"
        populated_prefix = out / f"working_populated_{realization:02d}"
        seed = args.seed_base + time_offset + config_index * 100 + realization
        print(args.config, args.integration, realization + 1, "/", args.realizations,
              "seed", seed, flush=True)
        regenerate_noise(seed)
        image_ms(blank_prefix)
        blank, blank_valid = grid.read_plane(str(blank_prefix) + ".smooth.image")
        pb, pb_valid = grid.read_plane(str(blank_prefix) + ".pb")
        usable = blank_valid & pb_valid & np.isfinite(blank) & np.isfinite(pb) & (pb >= grid.PB_LIMIT)
        rms = float(imstat(str(blank_prefix) + ".smooth.image", listit=False)["rms"][0])
        diagnostic = dict(rms=rms, seed_sigma=2.0, flood_sigma=1.5,
                          min_pixels=grid.MIN_PIXELS, detection_mask=~usable,
                          exclude_border=True)
        blank_input = np.where(usable, blank, 0.0)
        negative = extract_sources(blank_input, polarity="negative", **diagnostic)
        positive_blank = extract_sources(blank_input, **diagnostic)

        source_rows = []
        for pb_index, (lower, upper, positions) in enumerate(grid.choose_positions(pb, usable, realization)):
            snrs = np.roll(grid.TARGET_SNRS, realization % len(grid.TARGET_SNRS))
            for local_id, (snr, position) in enumerate(zip(snrs, positions)):
                p = float(pb[position]); apparent = float(snr * rms)
                source_rows.append({"realization": realization,
                    "source_id": pb_index * len(grid.TARGET_SNRS) + local_id,
                    "pb_bin": pb_index, "pb_lower": lower, "pb_upper": upper,
                    "x_pixel": position[0], "y_pixel": position[1],
                    "direction": grid.direction_for_pixel(position), "pb": p,
                    "expected_snr": float(snr), "apparent_flux_jy": apparent,
                    "intrinsic_flux_jy": apparent / p})

        make_components(source_rows)
        ft(vis=str(ms), spw="0:0", complist=str(cl_path), usescratch=True)
        uvsub(vis=str(ms), reverse=True)
        image_ms(populated_prefix)
        populated, pop_valid = grid.read_plane(str(populated_prefix) + ".smooth.image")
        pop_usable = usable & pop_valid
        pop_input = np.where(pop_usable, populated, 0.0)
        production = extract_sources(pop_input, rms=rms, seed_sigma=4.0, flood_sigma=3.0,
                                     min_pixels=grid.MIN_PIXELS, detection_mask=~pop_usable,
                                     exclude_border=True)
        diagnostic_positive = extract_sources(pop_input, **diagnostic)
        positions = [(row["x_pixel"], row["y_pixel"]) for row in source_rows]
        matches = match_catalogue_to_positions(positions, production.catalogue,
                                               max_distance_pixels=grid.MAX_MATCH_DISTANCE_PIXELS)
        by_input = {match.input_index: match for match in matches}
        for index, row in enumerate(source_rows):
            match = by_input.get(index)
            item = production.catalogue[match.detection_index] if match else None
            output = dict(row)
            output.update({"recovered": bool(match),
                "recovered_peak_jy_per_beam": float(item.peak_value) if item else "",
                "recovered_peak_snr": float(item.peak_value / rms) if item else "",
                "offset_pixels": float(match.distance_pixels) if match else ""})
            all_sources.append(output)

        for kind, result in (("populated_diagnostic", diagnostic_positive),
                             ("blank_negative", negative),
                             ("blank_positive", positive_blank)):
            for detection_id, item in enumerate(result.catalogue):
                x, y = item.peak_index
                all_detections.append({"realization": realization,
                    "detection_kind": kind, "detection_id": detection_id,
                    "x_pixel": x, "y_pixel": y, "pb": float(pb[x, y]),
                    "peak_jy_per_beam": float(item.peak_value),
                    "peak_snr": abs(float(item.peak_value)) / rms,
                    "n_pixels": item.n_pixels, "matched": ""})
        logs.append({"configuration": args.config, "integration_time": args.integration,
                     "realization": realization,
                     "noise_seed": seed, "blank_rms_jy_per_beam": rms,
                     "n_blank_negative": len(negative.catalogue),
                     "n_blank_positive": len(positive_blank.catalogue),
                     "n_recovered": len(matches), "n_injected": len(source_rows)})
        cleanup()

    csv_write(out / "source_recovery_all.csv", all_sources)
    csv_write(out / "detections_all.csv", all_detections)
    csv_write(out / "realization_summary.csv", logs)
    summary = {"status": "PASS", "configuration": args.config,
               "integration_time": args.integration,
               "track": "common_resolution", "uvtaper": UVTAPER[0],
               "target_beam": TARGET_BEAM, "n_realizations": args.realizations,
               "n_injected": len(all_sources),
               "n_recovered": sum(row["recovered"] for row in all_sources),
               "mean_blank_rms_jy_per_beam": float(np.mean([row["blank_rms_jy_per_beam"] for row in logs]))}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
