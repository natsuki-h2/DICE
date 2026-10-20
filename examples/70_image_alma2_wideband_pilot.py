"""Image the 7.5/15/30-GHz MS pilot and test radiometer scaling.

Run after ``69_build_alma2_wideband_ms_pilot.py`` with CASA 6.
"""

from __future__ import annotations

import json
from pathlib import Path
import shutil

from casatasks import imstat, tclean


ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT / "data" / "population_forward" / "alma2_wideband_pilot"
OUTPUT = PILOT / "imaging"
PRODUCTS = {
    "current_alma_7p5": (
        PILOT
        / "current_alma_7p5"
        / "current_alma_7p5.alma.cycle9.7.noisy.ms"
    ),
    "alma2_15": PILOT / "alma2_15ghz.ms",
    "alma2_30": PILOT / "alma2_30ghz.ms",
}


def main():
    if OUTPUT.exists():
        shutil.rmtree(OUTPUT)
    OUTPUT.mkdir(parents=True)
    summary = {}
    for label, measurement_set in PRODUCTS.items():
        imagename = OUTPUT / label
        tclean(
            vis=str(measurement_set),
            imagename=str(imagename),
            specmode="mfs",
            deconvolver="hogbom",
            gridder="standard",
            imsize=[512, 512],
            cell="0.02arcsec",
            phasecenter="J2000 03h32m39.0s -27d48m00.0s",
            weighting="natural",
            niter=0,
            interactive=False,
            pbcor=False,
        )
        statistics = imstat(
            imagename=f"{imagename}.image",
            box="32,32,479,479",
        )
        summary[label] = {
            "measurement_set": str(measurement_set),
            "rms_jy_beam": float(statistics["rms"][0]),
        }
    reference = summary["current_alma_7p5"]["rms_jy_beam"]
    for label in summary:
        summary[label]["rms_relative_to_7p5"] = (
            summary[label]["rms_jy_beam"] / reference
        )
    summary["expected_bandwidth_only"] = {
        "alma2_15_relative_to_7p5": 2.0**-0.5,
        "alma2_30_relative_to_7p5": 4.0**-0.5,
    }
    (OUTPUT / "wideband_imaging_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
