"""Run a 10-realization single-pointing Monte Carlo smoke test.

This runner reuses the validated CASA population example and varies only the
thermal-noise seed. It writes aggregate CSV files that can later
be used to calculate differential and cumulative completeness and the two
contamination estimators without reopening CASA products.

Run from the repository root with CASA 6:

    casa --nogui --nologger --log2term \
        -c scripts/tutorials/run_single_pointing_monte_carlo.py
"""

from pathlib import Path
import csv
import importlib.util
import os
import shutil


REPO_ROOT = Path(__file__).resolve().parents[2]
TUTORIAL_DATA_ROOT = Path(
    os.environ.get("DICE_TUTORIAL_DATA_ROOT", REPO_ROOT / "data")
).resolve()
OUTPUT_ROOT = TUTORIAL_DATA_ROOT / "casa_single_pointing_monte_carlo"
N_REALIZATIONS = 10
FIRST_NOISE_SEED = 21000


def load_population_example():
    path = REPO_ROOT / "scripts" / "tutorials" / "simulate_single_pointing_population.py"
    spec = importlib.util.spec_from_file_location("bliss_population_example", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_rows(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(path, rows):
    if not rows:
        raise RuntimeError("Cannot write an empty aggregate table: {0}".format(path))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def add_run_columns(rows, realization, noise_seed, detection_kind=None):
    output = []
    for row in rows:
        combined = {
            "realization": realization,
            "noise_seed": noise_seed,
        }
        if detection_kind is not None:
            combined["detection_kind"] = detection_kind
        combined.update(row)
        output.append(combined)
    return output


def main():
    if OUTPUT_ROOT.exists():
        shutil.rmtree(OUTPUT_ROOT)
    OUTPUT_ROOT.mkdir(parents=True)

    population = load_population_example()
    summaries = []
    recovery_rows = []
    detection_rows = []

    for realization in range(N_REALIZATIONS):
        noise_seed = FIRST_NOISE_SEED + realization
        population.WORKDIR = OUTPUT_ROOT / "realization_{0:03d}".format(realization)
        population.NOISE_SEED = noise_seed
        print(
            "\n=== Realization {0}/{1}; seed={2} ===".format(
                realization + 1, N_REALIZATIONS, noise_seed
            )
        )
        summary = population.main()

        summary_row = {
            "realization": realization,
            "noise_seed": noise_seed,
            "n_injected": summary["n_injected"],
            "n_detected": summary["n_detected"],
            "n_matched": summary["n_matched"],
            "n_unmatched": summary["n_unmatched"],
            "n_blank_positive": summary["n_blank_positive"],
            "n_blank_negative": summary["n_blank_negative"],
            "blank_rms_jy_per_beam": summary["blank_rms_jy_per_beam"],
        }
        summaries.append(summary_row)
        recovery_rows.extend(
            add_run_columns(
                read_rows(summary["recovery_table"]), realization, noise_seed
            )
        )
        for kind, key in (
            ("populated", "populated_catalogue"),
            ("blank_positive", "blank_positive_catalogue"),
            ("blank_negative", "blank_negative_catalogue"),
        ):
            rows = read_rows(summary[key])
            if rows:
                detection_rows.extend(
                    add_run_columns(rows, realization, noise_seed, kind)
                )

    write_rows(OUTPUT_ROOT / "realization_summary.csv", summaries)
    write_rows(OUTPUT_ROOT / "source_recovery_all.csv", recovery_rows)
    if detection_rows:
        write_rows(OUTPUT_ROOT / "detections_all.csv", detection_rows)

    total_injected = sum(int(row["n_injected"]) for row in summaries)
    total_matched = sum(int(row["n_matched"]) for row in summaries)
    total_unmatched = sum(int(row["n_unmatched"]) for row in summaries)
    total_blank_positive = sum(int(row["n_blank_positive"]) for row in summaries)
    print("\nBliSS-ALMA 10-realization smoke-test result")
    print("  Total injected: {0}".format(total_injected))
    print("  Total matched: {0}".format(total_matched))
    print("  Overall completeness: {0:.3f}".format(total_matched / total_injected))
    print("  Unmatched populated detections: {0}".format(total_unmatched))
    print("  Blank positive detections: {0}".format(total_blank_positive))
    print("  Aggregate output: {0}".format(OUTPUT_ROOT))
    print("PASS: all {0} realizations completed.".format(N_REALIZATIONS))


if __name__ == "__main__":
    main()
