"""Build and audit a short current-ALMA/ALMA2 wideband MS pilot.

The current-ALMA reference is one 7.5-GHz block centred at 230 GHz.  The
ALMA2 approximations concatenate two or four independently simulated
7.5-GHz blocks using the same 43-antenna C-7 geometry:

* 15 GHz: 222.5--237.5 GHz
* 30 GHz: 215.0--245.0 GHz

This pilot validates MS structure only.  The production run will use GHOST
for each supported 7.5-GHz block and will apply separately documented WSU
digital-efficiency/receiver noise scalings.

Run with CASA 6 from the repository root::

    casa --nogui --nologger --log2term \
        -c scripts/population_forward/build_wideband_ms.py
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import sys

import numpy as np

from casatasks import concat, simobserve
from casatools import coordsys, ctsys, image, measures, table


ROOT = Path(__file__).resolve().parents[2]
POPULATION_ROOT = Path(
    os.environ.get("DICE_POPULATION_OUTPUT_ROOT", ROOT / "data/population_forward")
).resolve()
OUTPUT = POPULATION_ROOT / "alma2_wideband_pilot"
PHASE_CENTRE = "J2000 03h32m39.0s -27d48m00.0s"
CONFIGURATION = str(Path(ctsys.resolve("alma/simmos")) / "alma.cycle9.7.cfg")
BLOCK_WIDTH_GHZ = 7.5
BLOCK_CENTRES_GHZ = [218.75, 226.25, 233.75, 241.25]
REFERENCE_CENTRE_GHZ = 230.0
TOTAL_TIME = "60s"
INTEGRATION = "10s"


def reset_output():
    if OUTPUT.exists():
        shutil.rmtree(OUTPUT)
    OUTPUT.mkdir(parents=True)


def create_blank_model(path, centre_ghz):
    cs_tool = coordsys()
    cs = cs_tool.newcoordsys(direction=True, stokes=["I"], spectral=True)
    me = measures()
    direction = me.direction("J2000", "03h32m39.0s", "-27d48m00.0s")
    cs.setreferencevalue([direction["m0"], direction["m1"]], type="direction")
    cs.setreferencepixel([31.5, 31.5], type="direction")
    cs.setincrement(["-0.02arcsec", "0.02arcsec"], type="direction")
    cs.setreferencevalue(f"{centre_ghz}GHz", type="spectral")
    cs.setincrement(f"{BLOCK_WIDTH_GHZ}GHz", type="spectral")
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


def find_ms(project_dir):
    matches = sorted(project_dir.glob("*.ms"))
    noisy = [path for path in matches if path.name.endswith(".noisy.ms")]
    candidates = noisy or matches
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one MS in {project_dir}, found {candidates}")
    return candidates[0]


def simulate_block(label, centre_ghz, seed):
    model = OUTPUT / f"{label}_blank.image"
    create_blank_model(model, centre_ghz)
    simobserve(
        project=label,
        skymodel=str(model),
        inbright="0Jy/pixel",
        indirection=PHASE_CENTRE,
        incell="0.02arcsec",
        incenter=f"{centre_ghz}GHz",
        inwidth=f"{BLOCK_WIDTH_GHZ}GHz",
        direction=PHASE_CENTRE,
        setpointings=True,
        mapsize=["1.28arcsec", "1.28arcsec"],
        integration=INTEGRATION,
        totaltime=TOTAL_TIME,
        antennalist=CONFIGURATION,
        obsmode="int",
        thermalnoise="tsys-atm",
        user_pwv=0.5,
        seed=seed,
        overwrite=True,
        graphics="none",
    )
    return find_ms(OUTPUT / label)


def array_digest(array):
    return hashlib.sha256(np.ascontiguousarray(array).view(np.uint8)).hexdigest()


def inspect_ms(path):
    tb = table()
    try:
        tb.open(str(path))
        data = tb.getcol("DATA")
        uvw = tb.getcol("UVW")
        rows = tb.nrows()
    finally:
        tb.done()
    try:
        tb.open(str(path / "SPECTRAL_WINDOW"))
        chan_freq = tb.getcol("CHAN_FREQ")
        chan_width = tb.getcol("CHAN_WIDTH")
        total_bandwidth = tb.getcol("TOTAL_BANDWIDTH")
    finally:
        tb.done()
    return {
        "path": str(path),
        "rows": int(rows),
        "data_shape": list(data.shape),
        "data_sha256": array_digest(data),
        "uvw_sha256": array_digest(uvw),
        "spw_count": int(total_bandwidth.size),
        "frequency_min_ghz": float(
            np.min(chan_freq - np.abs(chan_width) / 2.0) / 1e9
        ),
        "frequency_max_ghz": float(
            np.max(chan_freq + np.abs(chan_width) / 2.0) / 1e9
        ),
        "summed_bandwidth_ghz": float(np.sum(total_bandwidth) / 1e9),
        "visibility_real_std_jy": float(np.std(data.real)),
        "visibility_imag_std_jy": float(np.std(data.imag)),
    }


def main():
    reset_output()
    os.chdir(OUTPUT)
    reference = simulate_block(
        "current_alma_7p5", REFERENCE_CENTRE_GHZ, seed=20400000
    )
    blocks = [
        simulate_block(f"alma2_block_{index}", centre, seed=20400100 + index)
        for index, centre in enumerate(BLOCK_CENTRES_GHZ)
    ]

    # P0: facility labels point to byte-identical copies of one physical MS.
    p0_current = OUTPUT / "p0_current_alma.ms"
    p0_comparison = OUTPUT / "p0_comparison_label.ms"
    shutil.copytree(reference, p0_current)
    shutil.copytree(reference, p0_comparison)

    alma2_15 = OUTPUT / "alma2_15ghz.ms"
    alma2_30 = OUTPUT / "alma2_30ghz.ms"
    concat(
        vis=[str(blocks[1]), str(blocks[2])],
        concatvis=str(alma2_15),
        copypointing=True,
    )
    concat(
        vis=[str(block) for block in blocks],
        concatvis=str(alma2_30),
        copypointing=True,
    )

    products = {
        "configuration": CONFIGURATION,
        "configuration_note": (
            "43-antenna Cycle-9 C-7 proxy; ALMA configurations have not "
            "changed significantly since Cycle 9"
        ),
        "scope": (
            "short local structural pilot; production uses independently "
            "generated GHOST-compatible 7.5-GHz blocks"
        ),
        "current_alma_7p5": inspect_ms(reference),
        "alma2_15": inspect_ms(alma2_15),
        "alma2_30": inspect_ms(alma2_30),
        "blocks": [inspect_ms(block) for block in blocks],
        "p0_current": inspect_ms(p0_current),
        "p0_comparison": inspect_ms(p0_comparison),
    }
    products["p0_byte_equivalent"] = (
        products["p0_current"]["data_sha256"]
        == products["p0_comparison"]["data_sha256"]
        and products["p0_current"]["uvw_sha256"]
        == products["p0_comparison"]["uvw_sha256"]
    )
    products["independent_block_noise"] = len(
        {block["data_sha256"] for block in products["blocks"]}
    ) == len(blocks)
    tb = table()
    noise_vectors = []
    for block in blocks:
        try:
            tb.open(str(block))
            values = np.asarray(tb.getcol("DATA").real, dtype=float).ravel()
        finally:
            tb.done()
        noise_vectors.append((values - values.mean()) / values.std())
    products["block_noise_correlation"] = np.corrcoef(noise_vectors).tolist()
    products["maximum_absolute_offdiagonal_noise_correlation"] = float(
        np.max(
            np.abs(
                np.asarray(products["block_noise_correlation"])
                - np.eye(len(noise_vectors))
            )
        )
    )
    if not products["p0_byte_equivalent"]:
        raise AssertionError("P0 visibility identity failed")
    if not products["independent_block_noise"]:
        raise AssertionError("Wideband blocks do not have independent DATA")
    if products["maximum_absolute_offdiagonal_noise_correlation"] > 0.05:
        raise AssertionError("Wideband block noise is unexpectedly correlated")
    (OUTPUT / "wideband_ms_audit.json").write_text(
        json.dumps(products, indent=2) + "\n"
    )
    print(json.dumps(products, indent=2))


if __name__ == "__main__":
    main()
