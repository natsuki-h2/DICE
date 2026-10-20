"""Run the legacy Hayatsu et al. ADF22 continuum subtraction with CASA 6.

Usage:
    casa --nologger --nogui -c examples/44_uvcontsub_adf22_spw.py 0

The task writes a temporary ``Core_region.ms.contsub`` and atomically renames
it to a per-SPW result. Existing outputs are never overwritten.
"""

from pathlib import Path
import sys

from casatasks import uvcontsub_old


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/hayatsu_reproduction/adf22_original_ms/Core_region.ms"
OUTPUT_ROOT = ROOT / "data/hayatsu_reproduction/adf22_contsub_ms"
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

spw = int(sys.argv[-1])
if spw not in range(4):
    raise ValueError("SPW must be 0, 1, 2, or 3.")

temporary = Path(str(SOURCE) + ".contsub")
destination = OUTPUT_ROOT / f"Core_region.spw{spw}.ms.contsub"
if temporary.exists():
    raise FileExistsError(f"Temporary output already exists: {temporary}")
if destination.exists():
    raise FileExistsError(f"Destination already exists: {destination}")

uvcontsub_old(
    vis=str(SOURCE),
    field="Core_region",
    fitspw=str(spw),
    excludechans=True,
    combine="spw",
    solint="int",
    fitorder=0,
    spw=str(spw),
    want_cont=False,
)
temporary.rename(destination)
print("Wrote", destination)
