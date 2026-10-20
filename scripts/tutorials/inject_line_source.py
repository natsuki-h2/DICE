"""Inject one Gaussian-line point source into a working-copy ALMA MS."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from bliss_alma.injection import (
    gaussian_line_integrated_flux,
    inject_gaussian_line_point_source_visibility,
)
from bliss_alma.provenance import write_provenance


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", type=Path,
        default=ROOT / "configs/tutorials/basic_examples.json",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print the requested line injection without changing an MS.",
    )
    args = parser.parse_args()
    config_path = args.config.resolve()
    settings = json.loads(config_path.read_text())["visibility_line_injection"]
    vis = ROOT / settings["vis"]
    complist = ROOT / settings["complist"]
    request = {**settings, "vis": str(vis), "complist": str(complist)}
    if "peak_snr" in settings:
        request["derived_integrated_flux_jy_kms"] = gaussian_line_integrated_flux(
            peak_snr=settings["peak_snr"],
            rms_jy_per_beam=settings["rms_jy_per_beam"],
            fwhm_kms=settings["fwhm_kms"],
        )
    if args.dry_run:
        print(json.dumps({"status": "DRY-RUN", **request}, indent=2))
        return

    inject_gaussian_line_point_source_visibility(
        vis=str(vis),
        direction=settings["direction"],
        center_frequency_hz=settings["center_frequency_hz"],
        fwhm_kms=settings["fwhm_kms"],
        integrated_flux_jy_kms=settings.get("integrated_flux_jy_kms"),
        peak_snr=settings.get("peak_snr"),
        rms_jy_per_beam=settings.get("rms_jy_per_beam"),
        complist=str(complist),
        frequency_frame=settings.get("frequency_frame", "TOPO"),
    )
    manifest = ROOT / "data/tutorials/manifests/visibility_line_injection.json"
    write_provenance(
        manifest, root=ROOT, script=Path(__file__), command=sys.argv,
        config=config_path, inputs=[vis], outputs=[vis, complist],
        parameters=settings,
        casa_version="recorded by the CASA runtime when available",
    )
    print(manifest)


if __name__ == "__main__":
    main()
