"""P0 deterministic identity check for source-population forward modelling.

This inexpensive test verifies that changing only a facility label cannot
change catalogue ingestion, primary-beam attenuation, seeded thermal noise,
thresholding, or output serialization.  It precedes, but does not replace,
the corresponding CASA visibility-domain identity test.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "population_forward" / "catalogues"
OUTPUT = ROOT / "data" / "population_forward" / "p0_identity"


def read_catalogue(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def run_operator(rows, *, facility_label, seed, rms_jy, threshold_sigma):
    """Apply one label-independent analytic detection operator."""

    del facility_label
    rng = np.random.default_rng(seed)
    results = []
    for row in rows:
        # This is a mosaic-preflight PB draw, not the final CASA mosaic PB.
        # Hashing the source ID makes the response deterministic and shared.
        digest = hashlib.sha256(row["source_id"].encode()).digest()
        pb = 0.2 + 0.8 * int.from_bytes(digest[:8], "big") / (2**64 - 1)
        intrinsic_flux = float(row.get("flux_jy") or row.get("peak_flux_jy"))
        observed_flux = intrinsic_flux * pb + rng.normal(0.0, rms_jy)
        snr = observed_flux / rms_jy
        results.append(
            {
                "source_id": row["source_id"],
                "pb": pb,
                "observed_flux_jy": observed_flux,
                "snr": snr,
                "detected": bool(snr >= threshold_sigma),
            }
        )
    return results


def canonical_digest(results):
    payload = json.dumps(results, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    populations = {
        "continuum": (INPUT / "selection_continuum_catalogue.csv", 5e-6),
        "line": (INPUT / "selection_line_catalogue.csv", 1e-4),
    }
    summary = {}
    for population, (path, rms_jy) in populations.items():
        rows = read_catalogue(path)
        reference = run_operator(
            rows,
            facility_label="current_ALMA",
            seed=20400000,
            rms_jy=rms_jy,
            threshold_sigma=4.0,
        )
        comparison = run_operator(
            rows,
            facility_label="ALMA2040",
            seed=20400000,
            rms_jy=rms_jy,
            threshold_sigma=4.0,
        )
        digest_reference = canonical_digest(reference)
        digest_comparison = canonical_digest(comparison)
        if reference != comparison:
            raise AssertionError(f"P0 identity failed for {population}")
        summary[population] = {
            "sources": len(rows),
            "detected": sum(item["detected"] for item in reference),
            "sha256_reference": digest_reference,
            "sha256_comparison": digest_comparison,
            "byte_equivalent": digest_reference == digest_comparison,
        }
    summary["scope"] = (
        "analytic catalogue/operator preflight; CASA visibility-domain P0 remains"
    )
    (OUTPUT / "p0_identity_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
