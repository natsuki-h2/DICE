"""Production F1 run with ten independent visibility-noise realizations."""

from pathlib import Path
import atexit
import csv
import importlib.util
import json
import shutil
import sys

import numpy as np
from casatasks import ft, imstat, tclean, uvsub
from casatools import componentlist, table


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
from bliss_alma.extraction import extract_sources
from bliss_alma.recovery import match_catalogue_to_positions

spec = importlib.util.spec_from_file_location("f1grid", ROOT / "examples" / "15_ghost_f1_snr_pb_grid.py")
grid = importlib.util.module_from_spec(spec)
spec.loader.exec_module(grid)

WORK = ROOT / "data" / "ghost_f1_5min"
OUT = WORK / "f1_independent_noise_run"
MS = OUT / "working.ms"
CL = OUT / "working.cl"
BLANK_PREFIX = OUT / "working_blank"
POP_PREFIX = OUT / "working_populated"
N_REALIZATIONS = 10
NOISE_SEEDS = [20400000 + value for value in range(N_REALIZATIONS)]
CHUNK_ROWS = 100000


def remove(path):
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def remove_images(prefix):
    for path in OUT.glob(prefix.name + ".*"):
        remove(path)


def cleanup_working_products():
    remove(MS); remove(CL); remove_images(BLANK_PREFIX); remove_images(POP_PREFIX)


def csv_write(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def noise_scales():
    tb = table(); tb.open(str(grid.BLANK_MS))
    try:
        data = tb.getcol("DATA")
    finally:
        tb.close()
    return np.std(data.real, axis=(1, 2)), np.std(data.imag, axis=(1, 2))


def regenerate_noise(seed, real_sigma, imag_sigma):
    remove(MS); shutil.copytree(grid.BLANK_MS, MS)
    rng = np.random.default_rng(seed)
    tb = table(); tb.open(str(MS), nomodify=False)
    try:
        nrows = tb.nrows()
        for start in range(0, nrows, CHUNK_ROWS):
            count = min(CHUNK_ROWS, nrows - start)
            shape = (len(real_sigma), 1, count)
            real = rng.normal(size=shape) * real_sigma[:, None, None]
            imag = rng.normal(size=shape) * imag_sigma[:, None, None]
            data = real + 1j * imag
            tb.putcol("DATA", data, startrow=start, nrow=count)
            tb.putcol("CORRECTED_DATA", data, startrow=start, nrow=count)
    finally:
        tb.close()


def image_ms(prefix):
    remove_images(prefix)
    tclean(vis=str(MS), imagename=str(prefix), specmode="mfs",
           imsize=[grid.IMAGE_SIZE, grid.IMAGE_SIZE], cell=grid.CELL,
           weighting="natural", robust=0.5, niter=0, pblimit=grid.PB_LIMIT,
           interactive=False)


def make_components(rows):
    remove(CL); cl = componentlist()
    try:
        for row in rows:
            cl.addcomponent(dir=row["direction"], flux=row["apparent_flux_jy"],
                            fluxunit="Jy", shape="point")
        cl.rename(str(CL))
    finally:
        cl.done()


def detection_rows(catalogue, realization, kind, pb):
    rows = []
    for index, item in enumerate(catalogue):
        x, y = item.peak_index
        rows.append({"realization": realization, "detection_kind": kind,
                     "detection_id": index, "x_pixel": x, "y_pixel": y,
                     "pb": float(pb[x, y]), "peak_jy_per_beam": float(item.peak_value),
                     "peak_snr": "", "n_pixels": item.n_pixels, "matched": ""})
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    real_sigma, imag_sigma = noise_scales()
    all_sources, all_detections, logs = [], [], []

    for realization, seed in enumerate(NOISE_SEEDS):
        print(f"Realization {realization + 1}/{N_REALIZATIONS}, seed={seed}", flush=True)
        regenerate_noise(seed, real_sigma, imag_sigma)
        image_ms(BLANK_PREFIX)
        blank, blank_valid = grid.read_plane(str(BLANK_PREFIX) + ".image")
        pb, pb_valid = grid.read_plane(str(BLANK_PREFIX) + ".pb")
        usable = blank_valid & pb_valid & np.isfinite(blank) & np.isfinite(pb) & (pb >= grid.PB_LIMIT)
        mask = ~usable
        rms = float(imstat(str(BLANK_PREFIX) + ".image", listit=False)["rms"][0])
        blank_input = np.where(usable, blank, 0.0)
        diagnostic = dict(rms=rms, seed_sigma=2.0, flood_sigma=1.5,
                          min_pixels=grid.MIN_PIXELS, detection_mask=mask, exclude_border=True)
        neg = extract_sources(blank_input, polarity="negative", **diagnostic)
        pos_blank = extract_sources(blank_input, **diagnostic)

        source_rows = []
        for pb_index, (lower, upper, positions) in enumerate(grid.choose_positions(pb, usable, realization)):
            snrs = np.roll(grid.TARGET_SNRS, realization % len(grid.TARGET_SNRS))
            for source_id_local, (snr, position) in enumerate(zip(snrs, positions)):
                p = float(pb[position]); apparent = float(snr * rms)
                source_rows.append({"realization": realization,
                    "source_id": pb_index * len(grid.TARGET_SNRS) + source_id_local,
                    "pb_bin": pb_index, "pb_lower": lower, "pb_upper": upper,
                    "x_pixel": position[0], "y_pixel": position[1],
                    "direction": grid.direction_for_pixel(position), "pb": p,
                    "expected_snr": float(snr), "apparent_flux_jy": apparent,
                    "intrinsic_flux_jy": apparent / p})

        make_components(source_rows)
        ft(vis=str(MS), spw="0:0", complist=str(CL), usescratch=True)
        uvsub(vis=str(MS), reverse=True)
        image_ms(POP_PREFIX)
        populated, pop_valid = grid.read_plane(str(POP_PREFIX) + ".image")
        pop_input = np.where(usable & pop_valid, populated, 0.0)
        prod = extract_sources(pop_input, rms=rms, seed_sigma=4.0, flood_sigma=3.0,
                               min_pixels=grid.MIN_PIXELS, detection_mask=~(usable & pop_valid),
                               exclude_border=True)
        diag_pos = extract_sources(pop_input, **diagnostic)
        positions = [(r["x_pixel"], r["y_pixel"]) for r in source_rows]
        matches = match_catalogue_to_positions(positions, prod.catalogue,
                                               max_distance_pixels=grid.MAX_MATCH_DISTANCE_PIXELS)
        by_input = {m.input_index: m for m in matches}
        for index, row in enumerate(source_rows):
            match = by_input.get(index); output = dict(row)
            item = prod.catalogue[match.detection_index] if match else None
            output.update({"recovered": bool(match),
                "recovered_peak_jy_per_beam": float(item.peak_value) if item else "",
                "recovered_peak_snr": float(item.peak_value / rms) if item else "",
                "offset_pixels": float(match.distance_pixels) if match else ""})
            all_sources.append(output)

        for kind, result in (("populated_diagnostic", diag_pos),
                             ("blank_negative", neg), ("blank_positive", pos_blank)):
            rows = detection_rows(result.catalogue, realization, kind, pb)
            for row in rows:
                row["peak_snr"] = abs(float(row["peak_jy_per_beam"])) / rms
            all_detections.extend(rows)
        logs.append({"realization": realization, "noise_seed": seed,
                     "real_sigma_pol0_jy": real_sigma[0], "real_sigma_pol1_jy": real_sigma[1],
                     "imag_sigma_pol0_jy": imag_sigma[0], "imag_sigma_pol1_jy": imag_sigma[1],
                     "blank_rms_jy_per_beam": rms, "n_blank_negative": len(neg.catalogue),
                     "n_blank_positive": len(pos_blank.catalogue), "n_recovered": len(matches)})
        remove(MS); remove(CL); remove_images(BLANK_PREFIX); remove_images(POP_PREFIX)

    source_fields = list(all_sources[0]); detection_fields = list(all_detections[0]); log_fields = list(logs[0])
    csv_write(OUT / "source_recovery_all.csv", all_sources, source_fields)
    csv_write(OUT / "detections_all.csv", all_detections, detection_fields)
    csv_write(OUT / "realization_summary.csv", logs, log_fields)
    summary = {"status": "PASS", "n_realizations": N_REALIZATIONS,
               "noise_seeds": NOISE_SEEDS, "n_injected": len(all_sources),
               "n_recovered": sum(r["recovered"] for r in all_sources),
               "mean_blank_rms_jy_per_beam": float(np.mean([r["blank_rms_jy_per_beam"] for r in logs]))}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    atexit.register(cleanup_working_products)
    main()
