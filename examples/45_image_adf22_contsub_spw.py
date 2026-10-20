"""Create the native dirty cube used by the ADF22 untargeted search."""

from pathlib import Path
import sys

from casatasks import exportfits, tclean


ROOT = Path(__file__).resolve().parents[1]
MS_ROOT = ROOT / "data/hayatsu_reproduction/adf22_contsub_ms"
IMAGE_ROOT = ROOT / "data/hayatsu_reproduction/adf22_contsub_images"
IMAGE_ROOT.mkdir(parents=True, exist_ok=True)

spw = int(sys.argv[-1])
vis = MS_ROOT / f"Core_region.spw{spw}.ms.contsub"
imagename = IMAGE_ROOT / f"spw{spw}.Core_region.contsub.dirtycube"
if not vis.exists():
    raise FileNotFoundError(vis)
if Path(str(imagename) + ".image").exists():
    raise FileExistsError(str(imagename) + ".image")
psf_exists = Path(str(imagename) + ".psf").exists()

tclean(
    vis=str(vis),
    imagename=str(imagename),
    field="Core_region",
    specmode="cube",
    spw="0",
    interpolation="nearest",
    gridder="mosaic",
    imsize=[1600, 2400],
    cell="0.1arcsec",
    niter=0,
    weighting="natural",
    pbcor=True,
    outframe="BARY",
    calcpsf=not psf_exists,
    calcres=True,
    restart=True,
    parallel=False,
    interactive=False,
)
exportfits(
    imagename=str(imagename) + ".image",
    fitsimage=str(imagename) + ".image.fits",
    overwrite=False,
)
exportfits(
    imagename=str(imagename) + ".pb",
    fitsimage=str(imagename) + ".pb.fits",
    overwrite=False,
)
print("Wrote", str(imagename) + ".image.fits")
