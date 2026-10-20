"""Run the CASA-independent DICE recovery and mock-setup tutorials."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from bliss_alma.mockdata import describe_mock_setup
from bliss_alma.provenance import write_provenance
from bliss_alma.recovery import (
    completeness,
    flux_recovery_ratio,
    is_recovered,
    positional_offset_arcsec,
)


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "configs/tutorials/basic_examples.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "data/tutorials/products/basic_metrics.json",
    )
    args = parser.parse_args()
    config_path = args.config.resolve()
    settings = json.loads(config_path.read_text())
    recovery = settings["recovery"]
    setup = describe_mock_setup(**settings["mock_setup"])
    result = {
        # Preserve the published JSON schema; use flux_recovery_ratio in new APIs.
        "flux_fidelity": flux_recovery_ratio(
            recovery["input_flux_jy"], recovery["recovered_flux_jy"]
        ),
        "offset_arcsec": positional_offset_arcsec(
            recovery["input_position_arcsec"],
            recovery["recovered_position_arcsec"],
        ),
        "recovered": is_recovered(
            recovery["input_position_arcsec"],
            recovery["recovered_position_arcsec"],
            recovery["max_offset_arcsec"],
        ),
        "completeness": completeness(
            recovery["n_injected"], recovery["n_recovered"]
        ),
        "mock_setup": setup,
    }
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    write_provenance(
        output.with_suffix(".provenance.json"),
        root=ROOT,
        script=Path(__file__),
        command=sys.argv,
        config=config_path,
        inputs=[config_path],
        outputs=[output],
        parameters=settings,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
