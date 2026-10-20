"""Inspect the PB-high S/N=4 source missed despite slightly positive local noise."""

from pathlib import Path
import atexit
import csv
import importlib.util
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "production", ROOT / "scripts" / "ghost" / "run_f1_independent_noise.py"
)
run = importlib.util.module_from_spec(spec); spec.loader.exec_module(run)

def main():
    realization = 5; seed = run.NOISE_SEEDS[realization]
    with (run.OUT / "source_recovery_all.csv").open(newline="") as handle:
        rows = [r for r in csv.DictReader(handle) if int(r["realization"]) == realization]
    target = next(r for r in rows if int(r["pb_bin"]) == 3 and float(r["expected_snr"]) == 4.0)
    real_sigma, imag_sigma = run.noise_scales(); run.regenerate_noise(seed, real_sigma, imag_sigma)
    run.image_ms(run.BLANK_PREFIX); blank, _ = run.grid.read_plane(str(run.BLANK_PREFIX) + ".image")
    numeric = []
    for r in rows:
        numeric.append({"direction":r["direction"], "apparent_flux_jy":float(r["apparent_flux_jy"])})
    run.make_components(numeric)
    run.ft(vis=str(run.MS), spw="0:0", complist=str(run.CL), usescratch=True)
    run.uvsub(vis=str(run.MS), reverse=True); run.image_ms(run.POP_PREFIX)
    populated, _ = run.grid.read_plane(str(run.POP_PREFIX) + ".image")
    x, y = int(target["x_pixel"]), int(target["y_pixel"])
    rms = float(target["apparent_flux_jy"]) / 4.0
    sl = np.s_[x-3:x+4, y-3:y+4]
    peak_index = np.unravel_index(np.argmax(populated[sl]), populated[sl].shape)
    peak_xy = (x - 3 + peak_index[0], y - 3 + peak_index[1])
    print({"realization":realization, "seed":seed, "position":(x,y),
           "blank_at_input_snr":float(blank[x,y]/rms),
           "populated_at_input_snr":float(populated[x,y]/rms),
           "local_peak_xy":peak_xy, "local_peak_snr":float(populated[peak_xy]/rms),
           "source_response_at_input_snr":float((populated[x,y]-blank[x,y])/rms)})
    run.cleanup_working_products()

if __name__ == "__main__":
    atexit.register(run.cleanup_working_products); main()
