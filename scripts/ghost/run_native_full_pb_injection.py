"""Native-resolution, full-PB injection/recovery for GHOST F1--F6.

Run one configuration at a time with CASA 6.  All intermediate Measurement
Sets and CASA images are deleted after each realization.
"""

from pathlib import Path
import argparse
import csv
import json
import shutil
import sys

import numpy as np
from casatasks import ft, imstat, tclean, uvsub
from casatools import componentlist, image, table
from scipy.ndimage import maximum_filter

from bliss_alma.provenance import write_provenance


ROOT = Path(__file__).resolve().parents[2]
OUT_ROOT = ROOT / "data/ghost_matrix/native_full_pb_runs"
REFERENCE = ROOT / "data/ghost_f1_5min/f1_independent_noise_run/source_recovery_all.csv"
F1_CELL = 0.12454065519280681
BEAM_MINOR = {"F1": 0.7166181802749634, "F2": 0.5299125909805298,
              "F3": 0.28328394889831543, "F4": 0.13086126744747162,
              "F5": 0.04867600277066231, "F6": 0.01248726062476635}
CELLS = {key: value / 4.0 for key, value in BEAM_MINOR.items()}
IMSIZES = {"F1": 256, "F2": 336, "F3": 640,
           "F4": 1344, "F5": 3600, "F6": 14400}
BEAM_MAJOR = {"F1": 0.73606938123703, "F2": 0.5423505306243896,
              "F3": 0.2971660792827606, "F4": 0.15921376645565033,
              "F5": 0.06545752286911011, "F6": 0.013991966843605042}
BASE_MS = {
    "F1": ROOT / "data/ghost_f1_5min/output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms",
    **{f"F{i}": ROOT / f"data/ghost_matrix/f{i}_5min/extracted/"
        "output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms" for i in range(2, 7)},
}
THRESHOLDS = np.arange(2.0, 7.01, 0.5)
CHUNK = 100_000
PB_LIMIT = 0.2


def remove(path):
    if path.is_dir(): shutil.rmtree(path)
    elif path.exists(): path.unlink()


def remove_images(prefix):
    for path in prefix.parent.glob(prefix.name + ".*"): remove(path)


def read_plane(path):
    ia = image()
    try:
        ia.open(str(path))
        data = np.squeeze(np.asarray(ia.getchunk(), dtype=np.float32))
        valid = np.squeeze(np.asarray(ia.getchunk(getmask=True), dtype=bool))
    finally:
        ia.done()
    return data, valid


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def peak_catalogue(data, usable, rms, polarity="positive"):
    values = data if polarity == "positive" else -data
    work = np.where(usable, values, -np.inf).astype(np.float32, copy=False)
    local = maximum_filter(work, size=5, mode="constant", cval=-np.inf)
    selected = usable & (work >= 2.0 * rms) & (work == local)
    positions = np.argwhere(selected).astype(np.int32)
    peaks = work[selected].astype(np.float32) / rms
    return positions, peaks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", choices=list(CELLS), required=True)
    parser.add_argument("--realizations", type=int, default=10)
    parser.add_argument("--output-root", type=Path, default=OUT_ROOT)
    parser.add_argument("--keep-working", action="store_true")
    args = parser.parse_args()
    config = args.config; cell = CELLS[config]; imsize = IMSIZES[config]
    output_root = args.output_root.resolve()
    out = output_root / config.lower(); out.mkdir(parents=True, exist_ok=True)
    work = out / "working"
    remove(work); work.mkdir()

    reference = list(csv.DictReader(REFERENCE.open()))
    reference_realizations = sorted({int(row["realization"]) for row in reference})
    by_realization = {
        r: [row for row in reference if int(row["realization"]) == r]
        for r in reference_realizations
    }

    tb = table(); tb.open(str(BASE_MS[config]))
    try:
        nrows = tb.nrows(); sample = tb.getcol("DATA", startrow=0, nrow=1)
        npol = sample.shape[0]; total = 0
        rs = np.zeros(npol); rss = np.zeros(npol); ims = np.zeros(npol); imss = np.zeros(npol)
        for start in range(0, nrows, CHUNK):
            count = min(CHUNK, nrows-start); data = tb.getcol("DATA", startrow=start, nrow=count)
            real = data.real.reshape(npol, -1); imag = data.imag.reshape(npol, -1)
            total += real.shape[1]; rs += real.sum(1); rss += np.square(real).sum(1)
            ims += imag.sum(1); imss += np.square(imag).sum(1)
    finally: tb.close()
    real_sigma = np.sqrt(rss/total - np.square(rs/total))
    imag_sigma = np.sqrt(imss/total - np.square(ims/total))

    source_rows_all, threshold_rows, realization_rows = [], [], []
    for realization in range(args.realizations):
        print(config, realization+1, "/", args.realizations, flush=True)
        ms = work / f"r{realization:02d}.ms"
        blank_prefix = work / f"blank_{realization:02d}"
        populated_prefix = work / f"populated_{realization:02d}"
        components = work / f"sources_{realization:02d}.cl"
        shutil.copytree(BASE_MS[config], ms)
        rng = np.random.default_rng(20420000 + int(config[1:])*100 + realization)
        tb = table(); tb.open(str(ms), nomodify=False)
        try:
            for start in range(0, tb.nrows(), CHUNK):
                count = min(CHUNK, tb.nrows()-start); shape = (npol, 1, count)
                noise = (rng.normal(size=shape)*real_sigma[:, None, None] +
                         1j*rng.normal(size=shape)*imag_sigma[:, None, None])
                tb.putcol("DATA", noise, startrow=start, nrow=count)
                tb.putcol("CORRECTED_DATA", noise, startrow=start, nrow=count)
        finally: tb.close()
        common_tclean = dict(vis=str(ms), specmode="mfs", imsize=[imsize, imsize],
                             cell=f"{cell}arcsec", weighting="natural", robust=0.5,
                             niter=0, pblimit=PB_LIMIT, gridder="standard", interactive=False)
        tclean(imagename=str(blank_prefix), **common_tclean)
        blank, blank_valid = read_plane(Path(str(blank_prefix)+".image"))
        pb, pb_valid = read_plane(Path(str(blank_prefix)+".pb"))
        usable = blank_valid & pb_valid & np.isfinite(blank) & np.isfinite(pb) & (pb >= PB_LIMIT)
        rms = float(imstat(str(blank_prefix)+".image", listit=False)["rms"][0])

        rows = []
        reference_realization = reference_realizations[
            realization % len(reference_realizations)
        ]
        old_center = 180; new_center = imsize//2
        cl = componentlist()
        try:
            for ref in by_realization[reference_realization]:
                x = new_center + int(round((int(ref["x_pixel"])-old_center)*F1_CELL/cell))
                y = new_center + int(round((int(ref["y_pixel"])-old_center)*F1_CELL/cell))
                target = float(ref["expected_snr"]); apparent = target*rms
                pb_value = float(pb[x, y])
                row = {"configuration": config, "realization": realization,
                       "reference_layout_realization": reference_realization,
                       "source_id": int(ref["source_id"]), "pb_bin": int(ref["pb_bin"]),
                       "pb_lower": float(ref["pb_lower"]), "pb_upper": float(ref["pb_upper"]),
                       "x_pixel": x, "y_pixel": y, "direction": ref["direction"],
                       "pb": pb_value, "expected_snr": target,
                       "apparent_flux_jy": apparent, "intrinsic_flux_jy": apparent/pb_value}
                rows.append(row)
                cl.addcomponent(dir=ref["direction"], flux=apparent, fluxunit="Jy", shape="point")
            cl.rename(str(components))
        finally: cl.done()
        ft(vis=str(ms), spw="0:0", complist=str(components), usescratch=True)
        uvsub(vis=str(ms), reverse=True)
        tclean(imagename=str(populated_prefix), **common_tclean)
        populated, populated_valid = read_plane(Path(str(populated_prefix)+".image"))
        pop_usable = usable & populated_valid
        positive_pos, positive_snr = peak_catalogue(populated, pop_usable, rms, "positive")
        negative_pos, negative_snr = peak_catalogue(blank, usable, rms, "negative")

        match_radius = BEAM_MAJOR[config] / cell
        for row in rows:
            distances = np.hypot(positive_pos[:, 0]-row["x_pixel"], positive_pos[:, 1]-row["y_pixel"])
            candidates = np.flatnonzero((distances <= match_radius) & (positive_snr >= 4.0))
            if candidates.size:
                best = candidates[np.argmax(positive_snr[candidates])]
                row.update({"recovered": True, "recovered_peak_snr": float(positive_snr[best]),
                            "recovered_peak_jy_per_beam": float(positive_snr[best]*rms),
                            "offset_pixels": float(distances[best]),
                            "offset_beams": float(distances[best]*cell/BEAM_MAJOR[config])})
            else:
                row.update({"recovered": False, "recovered_peak_snr": "",
                            "recovered_peak_jy_per_beam": "", "offset_pixels": "", "offset_beams": ""})
            source_rows_all.append(row)
        for threshold in THRESHOLDS:
            threshold_rows.append({"configuration": config, "realization": realization,
                "threshold_snr": float(threshold),
                "n_positive": int(np.count_nonzero(positive_snr >= threshold)),
                "n_negative": int(np.count_nonzero(negative_snr >= threshold))})
        realization_rows.append({"configuration": config, "realization": realization,
            "noise_seed": 20420000 + int(config[1:])*100 + realization,
            "rms_jy_per_beam": rms, "n_injected": len(rows),
            "n_recovered": sum(row["recovered"] for row in rows),
            "n_positive_4sigma": int(np.count_nonzero(positive_snr >= 4)),
            "n_negative_4sigma": int(np.count_nonzero(negative_snr >= 4))})
        remove(ms); remove(components); remove_images(blank_prefix); remove_images(populated_prefix)

    write_csv(out/"source_recovery_all.csv", source_rows_all)
    write_csv(out/"threshold_counts.csv", threshold_rows)
    write_csv(out/"realization_summary.csv", realization_rows)
    summary = {"status":"PASS", "configuration":config, "track":"native_full_pb",
        "weighting":"natural", "uvtaper":"none", "gridder":"standard",
        "pixels_per_minor_beam":4, "cell_arcsec":cell, "imsize":imsize,
        "n_realizations":args.realizations, "n_injected":len(source_rows_all),
        "n_recovered":sum(row["recovered"] for row in source_rows_all),
        "mean_rms_jy_per_beam":float(np.mean([r["rms_jy_per_beam"] for r in realization_rows]))}
    (out/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    write_provenance(
        out / "run.provenance.json",
        root=ROOT,
        script=Path(__file__),
        command=sys.argv,
        inputs=[BASE_MS[config], REFERENCE],
        outputs=[
            out / "source_recovery_all.csv",
            out / "threshold_counts.csv",
            out / "realization_summary.csv",
            out / "summary.json",
        ],
        parameters={
            "configuration": config,
            "realizations": args.realizations,
            "cell_arcsec": cell,
            "imsize": imsize,
            "pb_limit": PB_LIMIT,
            "thresholds_snr": THRESHOLDS.tolist(),
        },
        random_seeds=[
            20420000 + int(config[1:]) * 100 + realization
            for realization in range(args.realizations)
        ],
    )
    if not args.keep_working:
        remove(work)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__": main()
