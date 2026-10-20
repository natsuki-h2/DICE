"""Validate the P2 [CII] matched-filter response with visibility cubes.

Creates three 64-channel, 15.625-MHz CASA sky cubes containing central
Gaussian lines, simulates five-minute current-ALMA C-7 observations, forms
native dirty cubes, and compares measured and predicted matched-filter S/N.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil

import numpy as np
from casatasks import simobserve, tclean
from casatools import coordsys, ctsys, image, measures


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "population_forward" / "p2_cii_visibility_validation"
PHASE = "J2000 03h32m39.0s -27d48m00.0s"
FREQ_GHZ = 260.0
WIDTH_MHZ = 15.625
NCHAN = 64
DV = 299792.458 * WIDTH_MHZ / 1000.0 / FREQ_GHZ
FWHM_KMS = [100.0, 200.0, 400.0]
PEAK_JY = 0.002


def remove(path):
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def make_model(path, fwhm_kms):
    cs_tool = coordsys()
    cs = cs_tool.newcoordsys(direction=True, stokes=["I"], spectral=True)
    me = measures()
    direction = me.direction("J2000", "03h32m39.0s", "-27d48m00.0s")
    cs.setreferencevalue([direction["m0"], direction["m1"]], type="direction")
    cs.setreferencepixel([31.5, 31.5], type="direction")
    cs.setincrement(["-0.025arcsec", "0.025arcsec"], type="direction")
    cs.setreferencevalue(f"{FREQ_GHZ}GHz", type="spectral")
    cs.setreferencepixel([(NCHAN - 1) / 2.0], type="spectral")
    cs.setincrement(f"{WIDTH_MHZ}MHz", type="spectral")
    channels = np.arange(NCHAN) - (NCHAN - 1) / 2.0
    fwhm_channels = fwhm_kms / DV
    profile = np.exp(-4.0 * np.log(2.0) * channels**2 / fwhm_channels**2)
    ia = image()
    try:
        ia.fromshape(
            outfile=str(path),
            shape=[64, 64, 1, NCHAN],
            csys=cs.torecord(),
            overwrite=True,
        )
        values = np.zeros((64, 64, 1, NCHAN), dtype=np.float32)
        values[32, 32, 0, :] = PEAK_JY * profile
        ia.putchunk(values)
        ia.setbrightnessunit("Jy/pixel")
    finally:
        ia.done(); cs.done(); cs_tool.done(); me.done()
    return profile


def read_cube(path):
    ia = image()
    try:
        ia.open(str(path))
        values = np.squeeze(np.asarray(ia.getchunk(), dtype=np.float64))
    finally:
        ia.done()
    if values.shape[0] == NCHAN:
        values = np.moveaxis(values, 0, -1)
    return values


def find_noisy_ms(project):
    matches = list(project.glob("*.noisy.ms"))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one noisy MS in {project}: {matches}")
    return matches[0]


def run_case(fwhm_kms, case_index):
    label = f"fwhm_{int(fwhm_kms)}"
    model = OUT / f"{label}.image"
    profile = make_model(model, fwhm_kms)
    os.chdir(OUT)
    simobserve(
        project=label,
        skymodel=str(model),
        inbright="",
        indirection=PHASE,
        incell="0.025arcsec",
        incenter=f"{FREQ_GHZ}GHz",
        inwidth=f"{WIDTH_MHZ}MHz",
        direction=PHASE,
        setpointings=True,
        mapsize=["1.6arcsec", "1.6arcsec"],
        integration="10s",
        totaltime="300s",
        antennalist=str(
            Path(ctsys.resolve("alma/simmos")) / "alma.cycle9.7.cfg"
        ),
        obsmode="int",
        thermalnoise="tsys-atm",
        user_pwv=0.5,
        seed=20700000 + case_index,
        overwrite=True,
        graphics="none",
    )
    measurement_set = find_noisy_ms(OUT / label)
    prefix = OUT / f"{label}_dirty"
    tclean(
        vis=str(measurement_set),
        imagename=str(prefix),
        specmode="cube",
        nchan=NCHAN,
        start=0,
        width=1,
        imsize=[256, 256],
        cell="0.025arcsec",
        weighting="natural",
        niter=0,
        interactive=False,
    )
    cube = read_cube(Path(f"{prefix}.image"))
    centre = np.array(cube.shape[:2]) // 2
    yy, xx = np.indices(cube.shape[:2])
    noise_mask = (xx - centre[1]) ** 2 + (yy - centre[0]) ** 2 > 20**2
    channel_rms = float(np.std(cube[noise_mask, :]))
    # Use the brightest pixel within one beam of the expected source.
    patch = cube[
        centre[0] - 4 : centre[0] + 5,
        centre[1] - 4 : centre[1] + 5,
        :,
    ]
    summed = np.tensordot(patch, profile, axes=([-1], [0]))
    peak_position = np.unravel_index(np.argmax(summed), summed.shape)
    spectrum = patch[peak_position[0], peak_position[1], :]
    measured_snr = float(
        np.dot(spectrum, profile) / (channel_rms * np.sqrt(np.sum(profile**2)))
    )
    predicted_snr = float(
        PEAK_JY * np.sqrt(np.sum(profile**2)) / channel_rms
    )
    return {
        "fwhm_kms": float(fwhm_kms),
        "fwhm_channels": float(fwhm_kms / DV),
        "channel_rms_jy_beam": channel_rms,
        "predicted_matched_snr": predicted_snr,
        "measured_matched_snr": measured_snr,
        "measured_over_predicted": measured_snr / predicted_snr,
        "peak_pixel_in_9x9_patch": [int(value) for value in peak_position],
        "measurement_set": str(measurement_set.relative_to(ROOT)),
    }


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    results = [run_case(width, index) for index, width in enumerate(FWHM_KMS)]
    ratios = np.asarray([row["measured_over_predicted"] for row in results])
    summary = {
        "status": "PASS" if np.all(np.abs(ratios - 1.0) < 0.25) else "CHECK",
        "channel_width_mhz": WIDTH_MHZ,
        "channel_width_kms": DV,
        "peak_flux_jy": PEAK_JY,
        "results": results,
        "acceptance": "absolute measured/predicted deviation < 25%",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    os.chdir(ROOT)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
