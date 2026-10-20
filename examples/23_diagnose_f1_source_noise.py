"""Recreate blank maps and measure local noise at every F1 injection position."""

from pathlib import Path
import atexit
import csv
import importlib.util


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("production", ROOT / "examples" / "21_ghost_f1_independent_noise_run.py")
run = importlib.util.module_from_spec(spec); spec.loader.exec_module(run)
SOURCE_CSV = run.OUT / "source_recovery_all.csv"
OUTPUT_CSV = run.OUT / "source_noise_diagnostics.csv"


def main():
    with SOURCE_CSV.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    real_sigma, imag_sigma = run.noise_scales()
    for realization, seed in enumerate(run.NOISE_SEEDS):
        print(f"Blank diagnostic {realization + 1}/10, seed={seed}", flush=True)
        run.regenerate_noise(seed, real_sigma, imag_sigma)
        run.image_ms(run.BLANK_PREFIX)
        blank, _ = run.grid.read_plane(str(run.BLANK_PREFIX) + ".image")
        current = [r for r in rows if int(r["realization"]) == realization]
        rms = float(current[0]["apparent_flux_jy"]) / float(current[0]["expected_snr"])
        for row in current:
            position = (int(row["x_pixel"]), int(row["y_pixel"]))
            local = float(blank[position] / rms)
            row["blank_pixel_snr"] = local
            row["input_plus_blank_pixel_snr"] = float(row["expected_snr"]) + local
        run.cleanup_working_products()
    fields = list(rows[0])
    with OUTPUT_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)
    print(OUTPUT_CSV)


if __name__ == "__main__":
    atexit.register(run.cleanup_working_products)
    main()
