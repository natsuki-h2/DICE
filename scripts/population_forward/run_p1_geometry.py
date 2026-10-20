"""P1 geometry/collecting-area comparison at fixed 7.5-GHz bandwidth.

P1a rescales each blank dirty image to a common RMS, isolating dirty-beam and
geometry effects.  P1b uses the same per-visibility noise for 43-antenna
current ALMA C-7 and 255-antenna ALMA2040 F4, retaining the collecting-area
gain.  Point-source dirty responses are formed by shifted native dirty PSFs;
this is the linear dirty-imaging equivalent of visibility injection.

Run with CASA 6::

    casa --nogui --nologger -c scripts/population_forward/run_p1_geometry.py
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import shutil

import numpy as np
from casatasks import simobserve, tclean
from casatools import coordsys, ctsys, image, measures, table
from scipy.ndimage import maximum_filter


ROOT = Path(__file__).resolve().parents[2]
POPULATION_ROOT = Path(
    os.environ.get("DICE_POPULATION_OUTPUT_ROOT", ROOT / "data/population_forward")
).resolve()
OUTPUT = POPULATION_ROOT / "p1_geometry"
WORK = OUTPUT / "working"
PHASE_CENTRE = "J2000 03h32m39.0s -27d48m00.0s"
F4_BASE = (
    ROOT
    / "data/ghost_matrix/f4_5min/extracted/"
    "output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms"
)
N_REALIZATIONS = 10
PER_VISIBILITY_SIGMA_JY = 0.0085
PB_LIMIT = 0.2
TARGET_SNRS = np.array(
    [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0,
     3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
)
CONFIGS = {
    "current_ALMA_C7": {
        "cell_arcsec": 0.025,
        "imsize": 1536,
        "beam_match_pixels": 6.0,
    },
    "ALMA2040_F4": {
        "cell_arcsec": 0.13086126744747162 / 4.0,
        "imsize": 1344,
        "beam_match_pixels": 6.0,
    },
}


def remove(path):
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def create_blank_model(path):
    cs_tool = coordsys()
    cs = cs_tool.newcoordsys(direction=True, stokes=["I"], spectral=True)
    me = measures()
    direction = me.direction("J2000", "03h32m39.0s", "-27d48m00.0s")
    cs.setreferencevalue([direction["m0"], direction["m1"]], type="direction")
    cs.setreferencepixel([31.5, 31.5], type="direction")
    cs.setincrement(["-0.025arcsec", "0.025arcsec"], type="direction")
    cs.setreferencevalue("230GHz", type="spectral")
    cs.setincrement("7.5GHz", type="spectral")
    ia = image()
    try:
        ia.fromshape(
            outfile=str(path),
            shape=[64, 64, 1, 1],
            csys=cs.torecord(),
            overwrite=True,
        )
        ia.setbrightnessunit("Jy/pixel")
    finally:
        ia.done()
        cs.done()
        cs_tool.done()
        me.done()


def current_base_ms():
    project = WORK / "current_structure"
    matches = list(project.glob("*.ms")) if project.exists() else []
    plain = [path for path in matches if not path.name.endswith(".noisy.ms")]
    if len(plain) == 1:
        return plain[0]
    model = WORK / "current_blank.image"
    create_blank_model(model)
    os.chdir(WORK)
    simobserve(
        project="current_structure",
        skymodel=str(model),
        inbright="0Jy/pixel",
        indirection=PHASE_CENTRE,
        incell="0.025arcsec",
        incenter="230GHz",
        inwidth="7.5GHz",
        direction=PHASE_CENTRE,
        setpointings=True,
        mapsize=["1.6arcsec", "1.6arcsec"],
        integration="10s",
        totaltime="300s",
        antennalist=str(
            Path(ctsys.resolve("alma/simmos")) / "alma.cycle9.7.cfg"
        ),
        obsmode="int",
        thermalnoise="",
        overwrite=True,
        graphics="none",
    )
    matches = [
        path
        for path in (WORK / "current_structure").glob("*.ms")
        if not path.name.endswith(".noisy.ms")
    ]
    if len(matches) != 1:
        raise RuntimeError(f"Could not resolve current-ALMA base MS: {matches}")
    return matches[0]


def put_noise(base_ms, output_ms, seed):
    shutil.copytree(base_ms, output_ms)
    rng = np.random.default_rng(seed)
    tb = table()
    try:
        tb.open(str(output_ms), nomodify=False)
        columns = set(tb.colnames())
        for start in range(0, tb.nrows(), 100_000):
            count = min(100_000, tb.nrows() - start)
            sample = tb.getcol("DATA", startrow=start, nrow=count)
            noise = (
                rng.normal(size=sample.shape)
                + 1j * rng.normal(size=sample.shape)
            ) * PER_VISIBILITY_SIGMA_JY
            tb.putcol("DATA", noise, startrow=start, nrow=count)
            if "CORRECTED_DATA" in columns:
                tb.putcol("CORRECTED_DATA", noise, startrow=start, nrow=count)
    finally:
        tb.done()


def read_plane(path, mask=False):
    ia = image()
    try:
        ia.open(str(path))
        values = ia.getchunk(getmask=True) if mask else ia.getchunk()
        values = np.squeeze(np.asarray(values))
    finally:
        ia.done()
    return values


def image_blank(config, measurement_set, realization):
    settings = CONFIGS[config]
    prefix = WORK / f"{config}_r{realization:02d}"
    tclean(
        vis=str(measurement_set),
        imagename=str(prefix),
        specmode="mfs",
        imsize=[settings["imsize"], settings["imsize"]],
        cell=f"{settings['cell_arcsec']}arcsec",
        weighting="natural",
        niter=0,
        pblimit=PB_LIMIT,
        interactive=False,
    )
    blank = read_plane(Path(f"{prefix}.image")).astype(np.float32)
    valid = read_plane(Path(f"{prefix}.image"), mask=True).astype(bool)
    pb = read_plane(Path(f"{prefix}.pb")).astype(np.float32)
    psf = read_plane(Path(f"{prefix}.psf")).astype(np.float32)
    usable = valid & np.isfinite(blank) & np.isfinite(pb) & (pb >= PB_LIMIT)
    rms = float(np.std(blank[usable]))
    for product in prefix.parent.glob(prefix.name + ".*"):
        remove(product)
    return blank, pb, psf, usable, rms


def source_offsets():
    # Four independent sky positions per input-S/N class.
    count = 4 * len(TARGET_SNRS)
    radius = np.linspace(2.0, 10.0, count)
    theta = np.arange(count) * np.pi * (3.0 - np.sqrt(5.0))
    return np.column_stack([radius * np.cos(theta), radius * np.sin(theta)])


def shifted_psf_population(blank, psf, pb, cell, intrinsic_fluxes, offsets):
    output = np.array(blank, copy=True)
    centre = np.array(blank.shape) // 2
    positions = []
    for flux, (dra, ddec) in zip(intrinsic_fluxes, offsets):
        x = int(round(centre[0] - dra / cell))
        y = int(round(centre[1] + ddec / cell))
        positions.append((x, y))
        apparent_flux = flux * float(pb[x, y])
        output += apparent_flux * np.roll(
            np.roll(psf, x - centre[0], axis=0), y - centre[1], axis=1
        )
    return output, np.asarray(positions, dtype=int)


def peak_catalogue(data, usable, rms, polarity=1):
    work = np.where(usable, polarity * data, -np.inf)
    local = maximum_filter(work, size=5, mode="constant", cval=-np.inf)
    selected = usable & (work == local) & (work >= 2.0 * rms)
    return np.argwhere(selected), work[selected] / rms


def evaluate(data, blank, usable, rms, positions, target_snrs, match_radius):
    positive_pos, positive_snr = peak_catalogue(data, usable, rms, 1)
    negative_pos, negative_snr = peak_catalogue(blank, usable, rms, -1)
    recovered = []
    matched_detection_indices = set()
    for position in positions:
        distances = np.hypot(
            positive_pos[:, 0] - position[0],
            positive_pos[:, 1] - position[1],
        )
        candidates = np.flatnonzero(
            (distances <= match_radius) & (positive_snr >= 4.0)
        )
        if candidates.size:
            best = int(candidates[np.argmax(positive_snr[candidates])])
            matched_detection_indices.add(best)
            recovered.append((True, float(positive_snr[best])))
        else:
            recovered.append((False, np.nan))
    rows = []
    for index, ((detected, peak_snr), target_snr) in enumerate(
        zip(recovered, target_snrs)
    ):
        rows.append(
            {
                "source_id": index,
                "target_current_alma_snr": float(target_snr),
                "recovered": bool(detected),
                "recovered_peak_snr": peak_snr,
            }
        )
    positive_4 = int(np.count_nonzero(positive_snr >= 4.0))
    negative_4 = int(np.count_nonzero(negative_snr >= 4.0))
    false_4 = max(positive_4 - len(matched_detection_indices), 0)
    return rows, {
        "positive_peaks_4sigma": positive_4,
        "negative_peaks_4sigma": negative_4,
        "truth_false_peaks_4sigma": false_4,
        "truth_contamination_4sigma": false_4 / max(positive_4, 1),
    }


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    if OUTPUT.exists():
        shutil.rmtree(OUTPUT)
    WORK.mkdir(parents=True)
    bases = {
        "current_ALMA_C7": current_base_ms(),
        "ALMA2040_F4": F4_BASE,
    }
    offsets = source_offsets()
    target_snr_template = np.resize(TARGET_SNRS, len(offsets))
    source_rows = []
    realization_rows = []
    fiducial_rms = None
    pb_values = None

    for realization in range(N_REALIZATIONS):
        products = {}
        for config_index, config in enumerate(CONFIGS):
            measurement_set = WORK / f"{config}_r{realization:02d}.ms"
            put_noise(
                bases[config],
                measurement_set,
                seed=20450000 + realization * 10 + config_index,
            )
            products[config] = image_blank(config, measurement_set, realization)
            remove(measurement_set)

        if fiducial_rms is None:
            fiducial_rms = products["current_ALMA_C7"][4]
            current_pb = products["current_ALMA_C7"][1]
            current_cell = CONFIGS["current_ALMA_C7"]["cell_arcsec"]
            centre = np.array(current_pb.shape) // 2
            pb_values = []
            for dra, ddec in offsets:
                x = int(round(centre[0] - dra / current_cell))
                y = int(round(centre[1] + ddec / current_cell))
                pb_values.append(float(current_pb[x, y]))
            pb_values = np.asarray(pb_values)

        # Cycle S/N values through sky positions so each differential bin
        # samples many PB/PSF locations across whole-field realizations.
        target_snrs = np.roll(target_snr_template, realization * 5)
        intrinsic_fluxes = target_snrs * fiducial_rms / pb_values

        for config in CONFIGS:
            blank_raw, pb, psf, usable, raw_rms = products[config]
            for track in ("P1a_geometry_only", "P1b_array_total"):
                if track == "P1a_geometry_only":
                    blank = blank_raw * (fiducial_rms / raw_rms)
                    rms = fiducial_rms
                else:
                    blank = blank_raw
                    rms = raw_rms
                populated, positions = shifted_psf_population(
                    blank,
                    psf,
                    pb,
                    CONFIGS[config]["cell_arcsec"],
                    intrinsic_fluxes,
                    offsets,
                )
                rows, counts = evaluate(
                    populated,
                    blank,
                    usable,
                    rms,
                    positions,
                    target_snrs,
                    CONFIGS[config]["beam_match_pixels"],
                )
                for row in rows:
                    row.update(
                        {
                            "track": track,
                            "configuration": config,
                            "realization": realization,
                            "rms_jy_beam": rms,
                        }
                    )
                    source_rows.append(row)
                realization_rows.append(
                    {
                        "track": track,
                        "configuration": config,
                        "realization": realization,
                        "raw_rms_jy_beam": raw_rms,
                        "analysis_rms_jy_beam": rms,
                        **counts,
                    }
                )
        print(f"P1 realization {realization + 1}/{N_REALIZATIONS}", flush=True)

    write_csv(OUTPUT / "source_recovery_all.csv", source_rows)
    write_csv(OUTPUT / "realization_summary.csv", realization_rows)
    means = {}
    for track in ("P1a_geometry_only", "P1b_array_total"):
        means[track] = {}
        for config in CONFIGS:
            selected = [
                row
                for row in realization_rows
                if row["track"] == track and row["configuration"] == config
            ]
            means[track][config] = {
                "mean_raw_rms_jy_beam": float(
                    np.mean([row["raw_rms_jy_beam"] for row in selected])
                ),
                "mean_analysis_rms_jy_beam": float(
                    np.mean([row["analysis_rms_jy_beam"] for row in selected])
                ),
                "mean_truth_contamination_4sigma": float(
                    np.mean([row["truth_contamination_4sigma"] for row in selected])
                ),
            }
    summary = {
        "status": "PASS",
        "bandwidth_ghz": 7.5,
        "per_visibility_sigma_jy": PER_VISIBILITY_SIGMA_JY,
        "realizations": N_REALIZATIONS,
        "sources_per_realization": len(offsets),
        "fiducial_current_alma_rms_jy_beam": fiducial_rms,
        "means": means,
        "method": (
            "native dirty images plus shifted native dirty PSFs; exact linear "
            "dirty response for pixel-centred point sources"
        ),
    }
    (OUTPUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    os.chdir(ROOT)
    remove(WORK)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
