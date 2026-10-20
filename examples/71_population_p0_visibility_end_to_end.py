"""P0 visibility-domain identity test for source-population forward modelling.

Two facility labels receive byte-identical copies of one noisy MS.  The same
multi-source component list is injected, dirty-imaged, searched without a
preselected source position, and matched to the injected-source catalogue.
P0 passes only when the images and serialized catalogues are
identical.

Run with CASA 6 after example 69::

    casa --nogui --nologger -c \
        examples/71_population_p0_visibility_end_to_end.py
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from casatasks import ft, imstat, tclean, uvsub
from casatools import componentlist, image
from scipy.ndimage import maximum_filter


ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT / "data" / "population_forward" / "alma2_wideband_pilot"
BASE_MS = (
    PILOT
    / "current_alma_7p5"
    / "current_alma_7p5.alma.cycle9.7.noisy.ms"
)
OUTPUT = ROOT / "data" / "population_forward" / "p0_visibility_end_to_end"
LABELS = ["current_alma", "comparison_label"]
CELL_ARCSEC = 0.02
IMSIZE = 512
CENTRE_RA_DEG = 53.1625
CENTRE_DEC_DEG = -27.8


def catalogue():
    offsets = [
        (-3.5, -3.0), (-1.5, -3.0), (0.5, -3.0), (2.5, -3.0),
        (-3.5, -0.8), (-1.5, -0.8), (0.5, -0.8), (2.5, -0.8),
        (-3.5, 1.5), (-1.5, 1.5), (0.5, 1.5), (2.5, 1.5),
    ]
    fluxes = np.geomspace(8e-5, 8e-4, len(offsets))
    rows = []
    for index, ((dra, ddec), flux) in enumerate(zip(offsets, fluxes)):
        ra = CENTRE_RA_DEG + dra / (3600.0 * np.cos(np.deg2rad(CENTRE_DEC_DEG)))
        dec = CENTRE_DEC_DEG + ddec / 3600.0
        rows.append(
            {
                "source_id": f"P0-{index:02d}",
                "offset_ra_arcsec": dra,
                "offset_dec_arcsec": ddec,
                "flux_jy": float(flux),
                "direction": f"J2000 {ra:.10f}deg {dec:.10f}deg",
                "x_pixel": int(round((IMSIZE - 1) / 2.0 - dra / CELL_ARCSEC)),
                "y_pixel": int(round((IMSIZE - 1) / 2.0 + ddec / CELL_ARCSEC)),
            }
        )
    return rows


def read_plane(path):
    ia = image()
    try:
        ia.open(str(path))
        data = np.squeeze(np.asarray(ia.getchunk(), dtype=np.float64))
    finally:
        ia.done()
    return data


def image_digest(array):
    return hashlib.sha256(np.ascontiguousarray(array).view(np.uint8)).hexdigest()


def run_label(label, sources):
    measurement_set = OUTPUT / f"{label}.ms"
    shutil.copytree(BASE_MS, measurement_set)
    components = OUTPUT / f"{label}.cl"
    cl = componentlist()
    try:
        for source in sources:
            cl.addcomponent(
                dir=source["direction"],
                flux=source["flux_jy"],
                fluxunit="Jy",
                shape="point",
            )
        cl.rename(str(components))
    finally:
        cl.done()
    ft(vis=str(measurement_set), spw="", complist=str(components), usescratch=True)
    uvsub(vis=str(measurement_set), reverse=True)
    prefix = OUTPUT / label
    tclean(
        vis=str(measurement_set),
        imagename=str(prefix),
        specmode="mfs",
        imsize=[IMSIZE, IMSIZE],
        cell=f"{CELL_ARCSEC}arcsec",
        phasecenter="J2000 03h32m39.0s -27d48m00.0s",
        weighting="natural",
        niter=0,
        interactive=False,
    )
    data = read_plane(Path(f"{prefix}.image"))
    rms = float(imstat(imagename=f"{prefix}.image", box="32,32,479,479")["rms"][0])
    local = maximum_filter(data, size=5, mode="constant", cval=-np.inf)
    peaks = np.argwhere((data == local) & (data >= 4.0 * rms))
    detections = []
    for source in sources:
        distances = np.hypot(
            peaks[:, 0] - source["x_pixel"], peaks[:, 1] - source["y_pixel"]
        )
        candidates = np.flatnonzero(distances <= 6.0)
        if candidates.size:
            best = candidates[np.argmax(data[tuple(peaks[candidates].T)])]
            position = peaks[best]
            detections.append(
                {
                    "source_id": source["source_id"],
                    "detected": True,
                    "peak_x": int(position[0]),
                    "peak_y": int(position[1]),
                    "peak_snr": float(data[tuple(position)] / rms),
                }
            )
        else:
            detections.append(
                {
                    "source_id": source["source_id"],
                    "detected": False,
                    "peak_x": None,
                    "peak_y": None,
                    "peak_snr": None,
                }
            )
    return {
        "label": label,
        "rms_jy_beam": rms,
        "image_sha256": image_digest(data),
        "n_blind_peaks_4sigma": int(len(peaks)),
        "detections": detections,
    }


def canonical_without_label(result):
    comparable = dict(result)
    comparable.pop("label")
    return json.dumps(comparable, sort_keys=True, separators=(",", ":"))


def main():
    if OUTPUT.exists():
        shutil.rmtree(OUTPUT)
    OUTPUT.mkdir(parents=True)
    sources = catalogue()
    results = [run_label(label, sources) for label in LABELS]
    passed = canonical_without_label(results[0]) == canonical_without_label(results[1])
    summary = {
        "status": "PASS" if passed else "FAIL",
        "scope": "visibility injection, dirty imaging, untargeted search, injected-source matching",
        "source_count": len(sources),
        "results": results,
        "facility_label_only_outputs_identical": passed,
    }
    (OUTPUT / "p0_end_to_end_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    if not passed:
        raise AssertionError("P0 end-to-end outputs differ")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
