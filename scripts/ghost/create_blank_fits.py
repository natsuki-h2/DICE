"""Create the zero-valued FITS sky model for the first GHOST F1 test.

Run from the repository root with CASA 6:

    casa --nogui --nologger --log2term \
        -c scripts/ghost/10_create_ghost_blank_fits.py

The generated FITS file is written below ``data/`` and is ignored by Git.
"""

from pathlib import Path
import math
import shutil

from casatasks import exportfits, imhead, imstat
from casatools import coordsys, image


REPO_ROOT = Path(__file__).resolve().parents[2]
WORKDIR = REPO_ROOT / "data" / "ghost_f1_5min"
CASA_IMAGE = WORKDIR / "ghost_blank_230ghz.image"
FITS_IMAGE = WORKDIR / "ghost_blank_230ghz.fits"

IMAGE_SIZE = 128
CELL = "0.2arcsec"
PHASE_CENTRE_RA = "03h32m39.0s"
PHASE_CENTRE_DEC = "-27d48m00.0s"
PHASE_CENTRE_RA_DEG = 53.1625
PHASE_CENTRE_DEC_DEG = -27.8
CENTRE_FREQUENCY = "230GHz"
BANDWIDTH = "7.5GHz"


def create_blank_image():
    """Create a four-axis, zero-valued CASA image with valid celestial WCS."""
    WORKDIR.mkdir(parents=True, exist_ok=True)
    if CASA_IMAGE.exists():
        shutil.rmtree(CASA_IMAGE)
    if FITS_IMAGE.exists():
        FITS_IMAGE.unlink()

    cs_tool = coordsys()
    cs = cs_tool.newcoordsys(direction=True, stokes=["I"], spectral=True)
    try:
        cs.setreferencecode("J2000", type="direction")
        cs.setreferencevalue(
            [f"{PHASE_CENTRE_RA_DEG}deg", f"{PHASE_CENTRE_DEC_DEG}deg"],
            type="direction",
        )
        cs.setreferencepixel(
            [(IMAGE_SIZE - 1) / 2.0, (IMAGE_SIZE - 1) / 2.0],
            type="direction",
        )
        cs.setincrement([f"-{CELL}", CELL], type="direction")
        cs.setreferencevalue(CENTRE_FREQUENCY, type="spectral")
        cs.setincrement(BANDWIDTH, type="spectral")

        ia = image()
        try:
            ia.fromshape(
                outfile=str(CASA_IMAGE),
                shape=[IMAGE_SIZE, IMAGE_SIZE, 1, 1],
                csys=cs.torecord(),
                overwrite=True,
            )
            ia.setbrightnessunit("Jy/pixel")
            ia.set(0.0)
        finally:
            ia.done()
    finally:
        cs.done()
        cs_tool.done()

    exportfits(
        imagename=str(CASA_IMAGE),
        fitsimage=str(FITS_IMAGE),
        overwrite=True,
        # GHOST accepts a 2D continuum image; frequency and bandwidth are
        # supplied separately in its simulation form.
        dropdeg=True,
    )
    stats = imstat(str(FITS_IMAGE), listit=False)
    ra = imhead(str(FITS_IMAGE), mode="get", hdkey="crval1")["value"]
    dec = imhead(str(FITS_IMAGE), mode="get", hdkey="crval2")["value"]
    assert math.isclose(ra, math.radians(PHASE_CENTRE_RA_DEG), abs_tol=1e-12)
    assert math.isclose(dec, math.radians(PHASE_CENTRE_DEC_DEG), abs_tol=1e-12)
    assert float(stats["min"][0]) == 0.0
    assert float(stats["max"][0]) == 0.0
    assert float(stats["sum"][0]) == 0.0
    print(f"Created {FITS_IMAGE}")
    print(f"Phase centre: {PHASE_CENTRE_RA} {PHASE_CENTRE_DEC} (J2000)")
    print("Verified all FITS pixels are zero")


if __name__ == "__main__":
    create_blank_image()
