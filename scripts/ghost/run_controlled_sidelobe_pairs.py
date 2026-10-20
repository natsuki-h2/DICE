"""Run a controlled native-resolution target + bright-neighbour sidelobe test.

Dirty imaging is linear.  For each independent visibility-noise realization,
the target-only and target-plus-neighbour images are constructed from the
measured dirty PSF and exactly the same dirty-noise image.  Run with CASA 6.
"""

from pathlib import Path
import argparse
import csv
import json
import math
import shutil

import numpy as np
from casatasks import imhead, imstat, tclean
from casatools import image, table
from scipy.ndimage import maximum_filter

ROOT = Path(__file__).resolve().parents[2]
OUT_ROOT = ROOT / "data/ghost_matrix/controlled_sidelobe_pair_test"
BEAM_MINOR = {"F1": 0.7166181802749634, "F5": 0.04867600277066231,
              "F6": 0.01248726062476635}
BEAM_MAJOR = {"F1": 0.73606938123703, "F5": 0.06545752286911011,
              "F6": 0.013991966843605042}
CELLS = {key: value / 4 for key, value in BEAM_MINOR.items()}
BASE_MS = {
    "F1": ROOT / "data/ghost_f1_5min/output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms",
    "F5": ROOT / "data/ghost_matrix/f5_5min/extracted/output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms",
    "F6": ROOT / "data/ghost_matrix/f6_5min/extracted/output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms",
}
TARGET_SNRS = (3.0, 4.0, 5.0)
BRIGHT_SNRS = (10.0, 30.0, 100.0)
SEPARATIONS_BEAM = (2.0, 5.0, 10.0, 20.0)
ANGLES_DEG = tuple(range(0, 360, 45))
IMSIZE = 512
CHUNK = 100_000


def remove(path):
    if path.is_dir(): shutil.rmtree(path)
    elif path.exists(): path.unlink()


def remove_images(prefix):
    for path in prefix.parent.glob(prefix.name + ".*"): remove(path)


def read_plane(path):
    ia = image()
    try:
        ia.open(str(path))
        return np.squeeze(np.asarray(ia.getchunk(), dtype=np.float32))
    finally: ia.done()


def shifted_no_wrap(array, shift):
    result = np.zeros_like(array)
    source = []; destination = []
    for size, delta in zip(array.shape, shift):
        if delta >= 0:
            source.append(slice(0, size-delta)); destination.append(slice(delta, size))
        else:
            source.append(slice(-delta, size)); destination.append(slice(0, size+delta))
    result[tuple(destination)] = array[tuple(source)]
    return result


def local_peak_snr(image_data, center, radius_pixels, rms):
    x0, y0 = center; radius = int(math.ceil(radius_pixels))
    xs = slice(x0-radius, x0+radius+1); ys = slice(y0-radius, y0+radius+1)
    patch = image_data[xs, ys]
    xx, yy = np.indices(patch.shape)
    mask = np.hypot(xx-radius, yy-radius) <= radius_pixels
    values = np.where(mask, patch, -np.inf)
    local = np.unravel_index(np.argmax(values), values.shape)
    peak = float(values[local]/rms)
    offset = float(np.hypot(local[0]-radius, local[1]-radius)/radius_pixels)
    return peak, offset


def injection_matched_catalogue(image_data, injected_positions, radius_pixels, rms):
    """Return total and injection-unmatched >=4-sigma local maxima."""
    peaks = (image_data >= 4*rms) & (image_data == maximum_filter(image_data, size=5))
    positions = np.argwhere(peaks)
    candidates = []
    for injection_index, injected_position in enumerate(injected_positions):
        for detection_index, position in enumerate(positions):
            distance = float(np.hypot(*(position-injected_position)))
            if distance <= radius_pixels:
                candidates.append((distance, injection_index, detection_index))
    matched_injections=set(); matched_detections=set()
    for _,injection_index,detection_index in sorted(candidates):
        if injection_index in matched_injections or detection_index in matched_detections: continue
        matched_injections.add(injection_index); matched_detections.add(detection_index)
    return len(positions), len(positions)-len(matched_detections)


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", choices=list(CELLS), required=True)
    parser.add_argument("--realizations", type=int, default=10)
    parser.add_argument("--output-root", type=Path, default=OUT_ROOT)
    args = parser.parse_args(); config = args.config
    out = args.output_root.resolve()/config.lower(); out.mkdir(parents=True, exist_ok=True)
    work = out/"working"; remove(work); work.mkdir()
    cell = CELLS[config]

    tb = table(); tb.open(str(BASE_MS[config]))
    try:
        nrows = tb.nrows(); sample = tb.getcol("DATA", startrow=0, nrow=1)
        npol = sample.shape[0]; total = 0
        rs=np.zeros(npol); rss=np.zeros(npol); ims=np.zeros(npol); imss=np.zeros(npol)
        for start in range(0, nrows, CHUNK):
            count=min(CHUNK,nrows-start); data=tb.getcol("DATA",startrow=start,nrow=count)
            real=data.real.reshape(npol,-1); imag=data.imag.reshape(npol,-1); total+=real.shape[1]
            rs+=real.sum(1); rss+=np.square(real).sum(1); ims+=imag.sum(1); imss+=np.square(imag).sum(1)
    finally: tb.close()
    real_sigma=np.sqrt(rss/total-np.square(rs/total)); imag_sigma=np.sqrt(imss/total-np.square(ims/total))

    rows=[]; realization_rows=[]; psf_saved=False
    for realization in range(args.realizations):
        print(config, realization+1, "/", args.realizations, flush=True)
        ms=work/f"r{realization:02d}.ms"; prefix=work/f"blank_{realization:02d}"
        shutil.copytree(BASE_MS[config],ms)
        rng=np.random.default_rng(20430000+int(config[1:])*100+realization)
        tb=table(); tb.open(str(ms),nomodify=False)
        try:
            for start in range(0,tb.nrows(),CHUNK):
                count=min(CHUNK,tb.nrows()-start); shape=(npol,1,count)
                noise=(rng.normal(size=shape)*real_sigma[:,None,None]+
                       1j*rng.normal(size=shape)*imag_sigma[:,None,None])
                tb.putcol("DATA",noise,startrow=start,nrow=count)
                tb.putcol("CORRECTED_DATA",noise,startrow=start,nrow=count)
        finally: tb.close()
        tclean(vis=str(ms),imagename=str(prefix),specmode="mfs",imsize=[IMSIZE,IMSIZE],
               cell=f"{cell}arcsec",weighting="natural",robust=.5,niter=0,
               pblimit=.2,gridder="standard",interactive=False)
        blank=read_plane(Path(str(prefix)+".image")); psf=read_plane(Path(str(prefix)+".psf"))
        rms=float(imstat(str(prefix)+".image",listit=False)["rms"][0])
        measured_major=float(imhead(str(prefix)+".image",mode="get",hdkey="beammajor")["value"])
        measured_minor=float(imhead(str(prefix)+".image",mode="get",hdkey="beamminor")["value"])
        center=(IMSIZE//2,IMSIZE//2); psf_peak=np.unravel_index(np.argmax(psf),psf.shape)
        psf=np.roll(psf,(center[0]-psf_peak[0],center[1]-psf_peak[1]),axis=(0,1))
        beam_scale=math.sqrt(measured_major*measured_minor); match_radius=measured_major/cell
        if not psf_saved:
            np.savez_compressed(out/"psf_display.npz",psf=psf,cell_arcsec=cell,
                                beam_major_arcsec=measured_major,beam_minor_arcsec=measured_minor)
            psf_saved=True
        for target_snr in TARGET_SNRS:
            target_image=blank+target_snr*rms*psf
            isolated_peak,isolated_offset=local_peak_snr(target_image,center,match_radius,rms)
            for bright_snr in BRIGHT_SNRS:
                for separation in SEPARATIONS_BEAM:
                    radius_pixels=separation*beam_scale/cell
                    for angle in ANGLES_DEG:
                        dx=int(round(radius_pixels*math.cos(math.radians(angle))))
                        dy=int(round(radius_pixels*math.sin(math.radians(angle))))
                        bright_psf=shifted_no_wrap(psf,(dx,dy))
                        bright_only=blank+bright_snr*rms*bright_psf
                        combined=target_image+bright_snr*rms*bright_psf
                        bright_center=(center[0]+dx,center[1]+dy)
                        bright_only_peak,bright_only_offset=local_peak_snr(
                            bright_only,center,match_radius,rms)
                        combined_peak,combined_offset=local_peak_snr(combined,center,match_radius,rms)
                        bright_total,bright_unmatched=injection_matched_catalogue(
                            bright_only,[bright_center],match_radius,rms)
                        combined_total,combined_unmatched=injection_matched_catalogue(
                            combined,[center,bright_center],match_radius,rms)
                        eta=float(bright_snr*bright_psf[center])
                        isolated_recovered=isolated_peak>=4
                        bright_false_association=bright_only_peak>=4
                        combined_recovered=combined_peak>=4
                        if bright_false_association:
                            recovery_class="sidelobe_only_false_association"
                        elif not isolated_recovered and combined_recovered:
                            recovery_class="sidelobe_boosted_recovery"
                        elif isolated_recovered and not combined_recovered:
                            recovery_class="sidelobe_suppressed_recovery"
                        elif isolated_recovered and combined_recovered:
                            recovery_class="robust_target_recovery"
                        else:
                            recovery_class="missed"
                        rows.append({"configuration":config,"realization":realization,
                            "target_snr":target_snr,"bright_snr":bright_snr,
                            "separation_beams":separation,"angle_deg":angle,
                            "dx_pixels":dx,"dy_pixels":dy,"psf_coupling":eta/bright_snr,
                            "eta_snr":eta,"isolated_peak_snr":isolated_peak,
                            "bright_only_peak_at_target_snr":bright_only_peak,
                            "combined_peak_snr":combined_peak,
                            "delta_peak_snr":combined_peak-isolated_peak,
                            "isolated_recovered":isolated_recovered,
                            "bright_only_false_association":bright_false_association,
                            "combined_recovered":combined_recovered,
                            "recovery_class":recovery_class,
                            "isolated_offset_beams":isolated_offset,
                            "bright_only_offset_beams":bright_only_offset,
                            "combined_offset_beams":combined_offset,
                            "bright_only_n_detections":bright_total,
                            "bright_only_n_unmatched":bright_unmatched,
                            "bright_only_truth_contamination":bright_unmatched/max(bright_total,1),
                            "combined_n_detections":combined_total,
                            "combined_n_unmatched":combined_unmatched,
                            "combined_truth_contamination":combined_unmatched/max(combined_total,1)})
        realization_rows.append({"configuration":config,"realization":realization,
            "noise_seed":20430000+int(config[1:])*100+realization,
            "rms_jy_per_beam":rms,"beam_major_arcsec":measured_major,
            "beam_minor_arcsec":measured_minor,"cell_arcsec":cell})
        remove(ms); remove_images(prefix)
    write_csv(out/"pair_results.csv",rows); write_csv(out/"realization_summary.csv",realization_rows)
    summary={"status":"PASS","configuration":config,"realizations":args.realizations,
        "n_paired_trials":len(rows),"target_snrs":TARGET_SNRS,"bright_snrs":BRIGHT_SNRS,
        "separations_beam":SEPARATIONS_BEAM,"angles_deg":ANGLES_DEG,
        "weighting":"natural","uvtaper":"none","pixels_per_minor_beam":4}
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    remove(work); print(json.dumps(summary,indent=2),flush=True)


if __name__=="__main__": main()
