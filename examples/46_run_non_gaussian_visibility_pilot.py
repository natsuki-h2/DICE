"""Non-Gaussian visibility-noise diagnostic at native resolution.

This is a diagnostic comparison, not an empirical ALMA residual model. It holds the
visibility marginal variance fixed and compares image-plane tails and isolated
point-source recovery after normalising each dirty image by its measured RMS.

Run with CASA 6, for example:
    casa --nologger --nogui -c examples/46_run_non_gaussian_visibility_pilot.py \
        --config F1 --realizations 10
"""

from pathlib import Path
import argparse
import csv
import json
import shutil

import numpy as np
from casatasks import imstat, tclean
from casatools import image, table
from scipy.ndimage import maximum_filter


ROOT = Path(__file__).resolve().parents[1]
OUT_ROOT = ROOT / "data/ghost_matrix/non_gaussian_visibility_pilot"
BASE_MS = {
    "F1": ROOT / "data/ghost_f1_5min/output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms",
    "F6": ROOT / "data/ghost_matrix/f6_5min/extracted/output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms",
}
BEAM_MINOR = {"F1": 0.7166181802749634, "F6": 0.01248726062476635}
BEAM_MAJOR = {"F1": 0.73606938123703, "F6": 0.013991966843605042}
CELL = {key: value / 4.0 for key, value in BEAM_MINOR.items()}
SCENARIOS = ("gaussian", "student_t3", "sparse_outliers", "time_coherent")
THRESHOLDS = np.arange(3.0, 7.01, 0.5)
INPUT_SNRS = (3.0, 4.0, 5.0, 6.0)
IMSIZE = 1024
CHUNK = 100_000


def remove(path):
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def remove_images(prefix):
    for path in prefix.parent.glob(prefix.name + ".*"):
        remove(path)


def read_plane(path):
    ia = image()
    try:
        ia.open(str(path))
        data = np.squeeze(np.asarray(ia.getchunk(), dtype=np.float32))
        mask = np.squeeze(np.asarray(ia.getchunk(getmask=True), dtype=bool))
    finally:
        ia.done()
    return data, mask


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def moments(z):
    values = z[np.isfinite(z)].astype(np.float64)
    mean = values.mean(); centered = values - mean
    variance = np.mean(centered**2)
    return mean, np.sqrt(variance), np.mean(centered**3) / variance**1.5, np.mean(centered**4) / variance**2 - 3.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", choices=list(BASE_MS), required=True)
    parser.add_argument("--realizations", type=int, default=10)
    args = parser.parse_args()
    config = args.config
    out = OUT_ROOT / config.lower(); out.mkdir(parents=True, exist_ok=True)
    work = out / "working"; remove(work); work.mkdir()

    tb = table(); tb.open(str(BASE_MS[config]))
    try:
        nrows = tb.nrows()
        sample = tb.getcol("DATA", startrow=0, nrow=1)
        npol = sample.shape[0]
        times = np.asarray(tb.getcol("TIME"), dtype=np.float64)
        _, time_index = np.unique(times, return_inverse=True)
        total = 0
        sums = np.zeros((2, npol)); sums2 = np.zeros((2, npol))
        for start in range(0, nrows, CHUNK):
            count = min(CHUNK, nrows-start)
            data = tb.getcol("DATA", startrow=start, nrow=count)
            for part, values in enumerate((data.real, data.imag)):
                flat = values.reshape(npol, -1)
                sums[part] += flat.sum(1); sums2[part] += np.square(flat).sum(1)
            total += data.shape[-1]
    finally:
        tb.close()
    sigma = np.sqrt(sums2/total - np.square(sums/total))

    image_rows = []; recovery_rows = []
    for realization in range(args.realizations):
        for scenario_index, scenario in enumerate(SCENARIOS):
            print(config, realization+1, args.realizations, scenario, flush=True)
            seed = 20460000 + int(config[1:])*10000 + realization*100 + scenario_index
            rng = np.random.default_rng(seed)
            ms = work / f"r{realization:02d}_{scenario}.ms"
            prefix = work / f"r{realization:02d}_{scenario}"
            shutil.copytree(BASE_MS[config], ms)
            coherent = None
            if scenario == "time_coherent":
                # A controlled common-mode residual shared by baselines at a
                # timestamp. rho fixes its fraction of marginal variance.
                rho = 0.10
                coherent = rng.normal(size=(2, npol, time_index.max()+1))
            tb = table(); tb.open(str(ms), nomodify=False)
            try:
                for start in range(0, nrows, CHUNK):
                    count = min(CHUNK, nrows-start); shape = (npol, 1, count)
                    if scenario == "gaussian":
                        parts = rng.normal(size=(2,) + shape)
                    elif scenario == "student_t3":
                        parts = rng.standard_t(3.0, size=(2,) + shape) / np.sqrt(3.0)
                    elif scenario == "sparse_outliers":
                        p = 1.0e-3; factor = 10.0
                        parts = rng.normal(size=(2,) + shape)
                        mask = rng.random(size=(2,) + shape) < p
                        parts[mask] *= factor
                        parts /= np.sqrt((1.0-p) + p*factor**2)
                    else:
                        rho = 0.10
                        independent = rng.normal(size=(2,) + shape)
                        shared = coherent[:, :, time_index[start:start+count]][:, :, None, :]
                        parts = np.sqrt(1.0-rho)*independent + np.sqrt(rho)*shared
                    noise = parts[0] * sigma[0, :, None, None] + 1j*parts[1] * sigma[1, :, None, None]
                    tb.putcol("DATA", noise, startrow=start, nrow=count)
                    tb.putcol("CORRECTED_DATA", noise, startrow=start, nrow=count)
            finally:
                tb.close()
            tclean(vis=str(ms), imagename=str(prefix), specmode="mfs",
                   imsize=[IMSIZE, IMSIZE], cell=f"{CELL[config]}arcsec",
                   weighting="natural", robust=0.5, niter=0, pblimit=0.2,
                   gridder="standard", interactive=False)
            dirty, valid = read_plane(Path(str(prefix)+".image"))
            psf, psf_valid = read_plane(Path(str(prefix)+".psf"))
            usable = valid & psf_valid & np.isfinite(dirty)
            # Exclude an edge guard region from all statistics.
            guard = 32; usable[:guard] = False; usable[-guard:] = False
            usable[:, :guard] = False; usable[:, -guard:] = False
            rms = float(imstat(str(prefix)+".image", listit=False)["rms"][0])
            z = dirty / rms
            mean, std, skew, excess = moments(z[usable])
            work_z = np.where(usable, z, -np.inf)
            local = maximum_filter(work_z, size=5, mode="constant", cval=-np.inf)
            peaks = work_z[usable & (work_z == local)]
            negative_work = np.where(usable, -z, -np.inf)
            negative_local = maximum_filter(negative_work, size=5, mode="constant", cval=-np.inf)
            negative_peaks = negative_work[usable & (negative_work == negative_local)]
            row = {"configuration":config, "realization":realization,
                   "scenario":scenario, "seed":seed, "rms_jy_per_beam":rms,
                   "pixel_mean_snr":mean, "pixel_std_snr":std,
                   "pixel_skewness":skew, "pixel_excess_kurtosis":excess}
            for threshold in THRESHOLDS:
                row[f"n_positive_peaks_ge_{threshold:.1f}"] = int(np.count_nonzero(peaks >= threshold))
                row[f"n_negative_peaks_ge_{threshold:.1f}"] = int(np.count_nonzero(negative_peaks >= threshold))
            image_rows.append(row)

            # Isolated-source recovery: sample positions, then add a shifted
            # measured PSF locally. This preserves the generated noise exactly.
            prng = np.random.default_rng(seed + 500000)
            radius_pix = int(np.ceil(BEAM_MAJOR[config] / CELL[config]))
            half = 24; center = np.array(psf.shape)//2
            psf_patch = psf[center[0]-half:center[0]+half+1, center[1]-half:center[1]+half+1]
            psf_patch = psf_patch / psf_patch[half, half]
            positions = prng.integers(guard+half, IMSIZE-guard-half, size=(200, 2))
            for input_snr in INPUT_SNRS:
                recovered = 0
                for x, y in positions:
                    patch = z[x-half:x+half+1, y-half:y+half+1] + input_snr*psf_patch
                    yy, xx = np.ogrid[-half:half+1, -half:half+1]
                    search = patch[(xx*xx + yy*yy) <= radius_pix**2]
                    recovered += bool(np.max(search) >= 4.0)
                recovery_rows.append({"configuration":config, "realization":realization,
                    "scenario":scenario, "input_snr":input_snr,
                    "n_trials":len(positions), "n_recovered":recovered,
                    "completeness":recovered/len(positions)})
            remove(ms); remove_images(prefix)

    write_csv(out/"image_statistics.csv", image_rows)
    write_csv(out/"source_recovery.csv", recovery_rows)
    summary = {"status":"PASS", "configuration":config,
        "description":"non-Gaussian visibility-noise diagnostic comparison",
        "scenarios":list(SCENARIOS), "n_realizations":args.realizations,
        "imsize":IMSIZE, "cell_arcsec":CELL[config], "pixels_per_minor_beam":4,
        "source_trials_per_snr_per_realization":200,
        "guardrail":"Synthetic diagnostic cases; not empirical ALMA residual models."}
    (out/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    remove(work)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
