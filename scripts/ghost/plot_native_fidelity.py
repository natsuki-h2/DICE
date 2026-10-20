"""Plot flux recovery and beam-normalized positional offsets for native F1--F6 runs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bliss_alma.provenance import write_script_provenance


ROOT = Path(__file__).resolve().parents[2]
CONFIGURATIONS = ["F1", "F2", "F3", "F4", "F5", "F6"]
BEAM_MAJOR = {
    "F1": 0.73606938123703, "F2": 0.5423505306243896,
    "F3": 0.2971660792827606, "F4": 0.15921376645565033,
    "F5": 0.06545752286911011, "F6": 0.013991966843605042,
}
BEAM_MINOR = {
    "F1": 0.7166181802749634, "F2": 0.5299125909805298,
    "F3": 0.28328394889831543, "F4": 0.13086126744747162,
    "F5": 0.04867600277066231, "F6": 0.01248726062476635,
}
COLORS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#000000"]
MARKERS = ["o", "s", "^", "D", "v", "P"]
LINESTYLES = ["-", "--", "-.", ":", (0, (5, 2)), (0, (3, 1, 1, 1))]


def realization_statistics(rows: list[dict[str, str]], quantity: str) -> tuple[np.ndarray, ...]:
    snr_values = np.array(sorted({float(row["expected_snr"]) for row in rows}))
    medians, lows, highs = [], [], []
    for snr in snr_values:
        realization_values = []
        realization_ids = sorted({int(row["realization"]) for row in rows})
        for realization in realization_ids:
            selected = [
                row for row in rows
                if int(row["realization"]) == realization
                and float(row["expected_snr"]) == snr
                and row["recovered"].lower() == "true"
            ]
            if not selected:
                continue
            if quantity == "flux":
                values = [
                    float(row["recovered_peak_jy_per_beam"]) / float(row["apparent_flux_jy"])
                    for row in selected
                ]
            else:
                config = selected[0]["configuration"].upper()
                cell = BEAM_MINOR[config] / 4.0
                effective_beam = np.sqrt(BEAM_MAJOR[config] * BEAM_MINOR[config])
                values = [float(row["offset_pixels"]) * cell / effective_beam for row in selected]
            realization_values.append(float(np.median(values)))
        if realization_values:
            medians.append(float(np.median(realization_values)))
            lows.append(float(np.percentile(realization_values, 16)))
            highs.append(float(np.percentile(realization_values, 84)))
        else:
            medians.append(np.nan); lows.append(np.nan); highs.append(np.nan)
    return snr_values, np.array(medians), np.array(lows), np.array(highs)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-directory",
        type=Path,
        default=ROOT / "data/ghost_matrix/native_full_pb_runs_30",
    )
    parser.add_argument(
        "--output",
        type=Path,
        # Legacy filename retained so existing paper scripts find the output.
        default=ROOT / "data/ghost_matrix/native_full_pb_runs_30/figures/native_fidelity.png",
    )
    args = parser.parse_args()
    input_directory = args.input_directory.resolve()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    inputs = []
    realization_counts = set()
    for config, color, marker, linestyle in zip(CONFIGURATIONS, COLORS, MARKERS, LINESTYLES):
        table = input_directory / config.lower() / "source_recovery_all.csv"
        summary = input_directory / config.lower() / "summary.json"
        inputs.extend([table, summary])
        rows = list(csv.DictReader(table.open(newline="")))
        realization_counts.add(json.loads(summary.read_text())["n_realizations"])
        for axis, quantity in zip(axes, ["flux", "position"]):
            snr, median, low, high = realization_statistics(rows, quantity)
            axis.plot(
                snr, median, color=color, marker=marker, linestyle=linestyle,
                linewidth=1.8, markersize=5, label=config,
            )
            axis.fill_between(snr, low, high, color=color, alpha=0.10)

    axes[0].axhline(1.0, color="0.35", linestyle="--", linewidth=1.0)
    axes[0].set(
        title="Flux recovery", xlabel="Injected apparent peak S/N",
        ylabel=r"$R_f=f_{\rm recovered}/f_{\rm injected}$", ylim=(0.75, 2.35),
    )
    axes[1].set(
        title="Astrometric accuracy", xlabel="Injected apparent peak S/N",
        ylabel=r"$R_\theta=\Delta\theta/\sqrt{B_{\rm maj}B_{\rm min}}$",
        ylim=(-0.02, 0.82),
    )
    for axis in axes:
        axis.axvline(4.0, color="0.4", linestyle=":", linewidth=1.0)
        axis.grid(alpha=0.22)
    axes[0].legend(title="Configuration", ncol=2, fontsize=8)
    n_realizations = realization_counts.pop() if len(realization_counts) == 1 else "mixed"
    fig.suptitle(
        f"Native-resolution F1--F6 recovery accuracy ({n_realizations} realizations per configuration)"
    )
    fig.savefig(output, dpi=220)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)
    write_script_provenance(
        output.with_suffix(".provenance.json"), root=ROOT, script=__file__,
        inputs=inputs, outputs=[output, output.with_suffix(".pdf")],
        parameters={
            "position_normalization": "sqrt(beam_major_fwhm * beam_minor_fwhm)",
            "summary": "median and 16--84 percentiles of realization medians",
        },
    )
    print(output)


if __name__ == "__main__":
    main()
