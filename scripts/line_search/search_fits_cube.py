"""Run a target-neutral signed multi-width search on one FITS cube."""

from argparse import ArgumentParser
import csv
import json
from pathlib import Path
import sys

from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import numpy as np

from bliss_alma.linesearch import merge_candidates, search_cube
from bliss_alma.provenance import write_provenance


ROOT = Path(__file__).resolve().parents[2]


def resolve(value):
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def noise_box_pixels(config, header, wcs):
    """Resolve the local-noise window, preferably in synthesized-beam units."""

    if config.get("noise_mode", "global") != "local":
        return None
    if "noise_box_pixels" in config:
        return int(config["noise_box_pixels"])
    if "noise_box_beams" not in config:
        raise ValueError(
            "local noise mode requires noise_box_beams or noise_box_pixels"
        )
    if "BMIN" not in header:
        raise ValueError(
            "noise_box_beams requires a FITS BMIN beam keyword; "
            "use noise_box_pixels only when the beam is unavailable"
        )
    pixel_scale_deg = float(np.min(np.abs(proj_plane_pixel_scales(wcs.celestial))))
    if not np.isfinite(pixel_scale_deg) or pixel_scale_deg <= 0:
        raise ValueError("Could not determine a valid spatial pixel scale")
    pixels_per_minor_beam = float(header["BMIN"]) / pixel_scale_deg
    return max(3, int(round(float(config["noise_box_beams"]) * pixels_per_minor_beam)))


def main():
    parser = ArgumentParser()
    parser.add_argument("config", type=Path)
    args = parser.parse_args()
    config_path = resolve(args.config)
    config = json.loads(config_path.read_text())
    image_path, pb_path = resolve(config["image_fits"]), resolve(config["pb_fits"])
    output = resolve(config["output_directory"])
    output.mkdir(parents=True, exist_ok=True)
    with fits.open(image_path, memmap=True) as image_hdul, fits.open(pb_path, memmap=True) as pb_hdul:
        cube = np.squeeze(image_hdul[0].data)
        pb = np.squeeze(pb_hdul[0].data)
        if pb.ndim == 3:
            pb = pb[pb.shape[0] // 2]
        header = image_hdul[0].header
        wcs = WCS(header)
        mask = np.isfinite(pb) & (pb >= float(config["pb_min"]))
        box_pixels = noise_box_pixels(config, header, wcs)
        raw = search_cube(
            cube,
            mask,
            config["width_channels"],
            float(config["threshold_snr"]),
            spw=config.get("spw_label", 0),
            noise_mode=config.get("noise_mode", "global"),
            noise_box_pixels=box_pixels,
        )
        rows = merge_candidates(raw, float(config.get("merge_radius_pixels", 5.0)))
        positive = negative = 0
        for row in rows:
            if row["snr"] > 0:
                positive += 1; identifier = f"{config.get('candidate_prefix', 'DICE')}-E{positive:02d}"
            else:
                negative += 1; identifier = f"{config.get('candidate_prefix', 'DICE')}-A{negative:02d}"
            world = wcs.pixel_to_world_values(row["x"], row["y"], row["channel"], 0.0)
            row.update(candidate_id=identifier, ra_deg=float(world[0]), dec_deg=float(world[1]), frequency_ghz=float(world[2]) / 1e9, pb=float(pb[row["y"], row["x"]]))
    table = output / "signed_line_candidates.csv"
    fields = list(rows[0]) if rows else ["candidate_id", "spw", "snr"]
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
    script_path = Path(__file__).resolve()
    write_provenance(output / "signed_line_candidates.provenance.json", root=ROOT, script=script_path, command=[sys.executable, str(script_path.relative_to(ROOT)), str(args.config)], config=config_path, inputs=[image_path, pb_path], outputs=[table], parameters=config)
    print(f"Wrote {table} ({positive} positive, {negative} negative)")


if __name__ == "__main__":
    main()
