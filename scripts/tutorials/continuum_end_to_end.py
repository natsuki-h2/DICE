"""Run a minimal noisy-continuum end-to-end test for BliSS-ALMA.

This example follows the structure of the original GenerateMockData workflow:

1. create a one-channel zero-valued sky model;
2. use simobserve to generate a noisy reference measurement set;
3. copy the reference MS and inject a continuum point source in the uv plane;
4. image both datasets with MFS;
5. report the reference RMS, recovered peak, flux recovery ratio, and S/N.

Run from the repository root with CASA 6:

    casa --nogui --nologger --log2term \
        -c scripts/tutorials/continuum_end_to_end.py

Generated products are written under ``data/`` and ignored by Git.
"""

from pathlib import Path
import glob
import os
import shutil
import sys

import numpy as np

from casatasks import imstat, simobserve, tclean
from casatools import coordsys, image, measures


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from bliss_alma.injection import inject_point_source_visibility


TUTORIAL_DATA_ROOT = Path(
    os.environ.get("DICE_TUTORIAL_DATA_ROOT", REPO_ROOT / "data")
).resolve()
WORKDIR = TUTORIAL_DATA_ROOT / "casa_continuum_e2e"
PROJECT = "reference_simulation"
PHASE_CENTRE = "J2000 03h32m39.0s -27d46m00.0s"
CENTRE_FREQUENCY = "230GHz"
CONTINUUM_BANDWIDTH = "2GHz"
CELL = "0.2arcsec"
IMAGE_SIZE = 128
INJECTED_FLUX_JY = 0.01


def reset_workdir():
    if WORKDIR.exists():
        shutil.rmtree(WORKDIR)
    WORKDIR.mkdir(parents=True)


def create_zero_sky_model(path):
    """Create a one-channel zero-valued CASA image with valid coordinates."""
    cs_tool = coordsys()
    cs = cs_tool.newcoordsys(direction=True, stokes=["I"], spectral=True)
    me = measures()

    direction = me.direction("J2000", "03h32m39.0s", "-27d46m00.0s")
    cs.setreferencevalue(
        [direction["m0"], direction["m1"]],
        type="direction",
    )
    cs.setreferencepixel(
        [(IMAGE_SIZE - 1) / 2.0, (IMAGE_SIZE - 1) / 2.0],
        type="direction",
    )
    cs.setincrement([f"-{CELL}", CELL], type="direction")
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
        ia.setbrightnessunit("Jy/pixel")
    finally:
        ia.done()
        cs.done()
        cs_tool.done()
        me.done()


def find_reference_ms():
    matches = glob.glob(str(WORKDIR / PROJECT / "*.ms"))
    noisy_matches = [path for path in matches if path.endswith(".noisy.ms")]
    if len(noisy_matches) == 1:
        return Path(noisy_matches[0])
    if len(matches) == 1:
        return Path(matches[0])
    raise RuntimeError(f"Expected one noisy reference MS, found: {matches}")


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
        interactive=False,
    )


def save_image_array(imagename, output_path):
    """Save a CASA continuum image as a 2D NumPy array."""
    ia = image()
    try:
        ia.open(str(imagename))
        pixels = np.asarray(ia.getchunk(), dtype=np.float64)
    finally:
        ia.done()
    pixels = np.squeeze(pixels)
    if pixels.ndim != 2:
        raise RuntimeError(
            "Expected a 2D continuum image after squeezing; got {0}.".format(
                pixels.shape
            )
        )
    np.save(str(output_path), np.ascontiguousarray(pixels))


def main():
    reset_workdir()
    zero_model = WORKDIR / "zero_sky.image"
    create_zero_sky_model(zero_model)

    os.chdir(WORKDIR)
    simobserve(
        project=PROJECT,
        skymodel=str(zero_model),
        inbright="0Jy/pixel",
        indirection=PHASE_CENTRE,
        incell=CELL,
        incenter=CENTRE_FREQUENCY,
        inwidth=CONTINUUM_BANDWIDTH,
        direction=PHASE_CENTRE,
        setpointings=True,
        mapsize=[f"{IMAGE_SIZE * 0.2}arcsec", f"{IMAGE_SIZE * 0.2}arcsec"],
        integration="10s",
        totaltime="120s",
        antennalist="alma.cycle8.1.cfg",
        obsmode="int",
        thermalnoise="tsys-atm",
        user_pwv=0.5,
        overwrite=True,
        graphics="none",
    )

    reference_ms = find_reference_ms()
    reference_image = WORKDIR / "reference"
    image_continuum(reference_ms, reference_image)
    reference_stats = imstat(imagename=f"{reference_image}.image")
    reference_rms_jy = float(reference_stats["rms"][0])
    save_image_array(f"{reference_image}.image", WORKDIR / "reference.npy")
    np.save(str(WORKDIR / "reference_rms_jy.npy"), reference_rms_jy)

    injected_ms = WORKDIR / "injected.ms"
    shutil.copytree(reference_ms, injected_ms)
    inject_point_source_visibility(
        vis=str(injected_ms),
        direction=PHASE_CENTRE,
        flux_jy=INJECTED_FLUX_JY,
        spw="",
        complist=str(WORKDIR / "injected_point.cl"),
    )

    injected_image = WORKDIR / "injected"
    image_continuum(injected_ms, injected_image)
    injected_stats = imstat(imagename=f"{injected_image}.image")
    save_image_array(f"{injected_image}.image", WORKDIR / "injected.npy")
    recovered_peak_jy = float(injected_stats["max"][0])
    peak_position = injected_stats["maxpos"]
    flux_recovery_ratio = recovered_peak_jy / INJECTED_FLUX_JY
    recovered_snr = recovered_peak_jy / reference_rms_jy

    print("BliSS-ALMA noisy continuum end-to-end result")
    print(f"  Reference MS:     {reference_ms}")
    print(f"  Injected MS:      {injected_ms}")
    print(f"  Injected flux:    {INJECTED_FLUX_JY:.6f} Jy")
    print(f"  Reference RMS:    {reference_rms_jy:.6f} Jy/beam")
    print(f"  Recovered peak:   {recovered_peak_jy:.6f} Jy/beam")
    print(f"  Flux recovery:    {flux_recovery_ratio:.3f}")
    print(f"  Recovered S/N:    {recovered_snr:.1f}")
    print(f"  Peak pixel:       {peak_position}")

    if recovered_snr <= 5.0:
        raise RuntimeError(
            "The injected source was not recovered above S/N=5; inspect products."
        )

    print("PASS: a continuum point source was recovered in the noisy simulation.")


if __name__ == "__main__":
    main()
