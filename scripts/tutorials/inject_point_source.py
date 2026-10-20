"""Inject one continuum point source into a user-supplied ALMA MS."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from bliss_alma.injection import inject_point_source_visibility
from bliss_alma.provenance import write_provenance


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "configs/tutorials/basic_examples.json",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate and print the requested injection without changing an MS.",
    )
    args = parser.parse_args()
    config_path = args.config.resolve()
    settings = json.loads(config_path.read_text())["visibility_injection"]
    vis = ROOT / settings["vis"]
    complist = ROOT / settings["complist"]
    request = {**settings, "vis": str(vis), "complist": str(complist)}
    if args.dry_run:
        print(json.dumps({"status": "DRY-RUN", **request}, indent=2))
        return
    inject_point_source_visibility(
        vis=str(vis),
        direction=settings["direction"],
        flux_jy=settings["flux_jy"],
        spw=settings["spw"],
        complist=str(complist),
    )
    manifest = ROOT / "data/tutorials/manifests/visibility_injection.json"
    write_provenance(
        manifest,
        root=ROOT,
        script=Path(__file__),
        command=sys.argv,
        config=config_path,
        inputs=[vis],
        outputs=[vis],
        parameters=settings,
        casa_version="recorded by the CASA runtime when available",
    )
    print(manifest)


if __name__ == "__main__":
    main()
