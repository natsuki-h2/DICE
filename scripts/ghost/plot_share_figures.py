"""Create email-ready figures for the final GHOST F1 independent-noise run."""

from pathlib import Path
import csv
import json

import matplotlib.pyplot as plt
import numpy as np
from casatools import image
from bliss_alma.provenance import write_script_provenance


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "data" / "ghost_f1_5min"
INPUT = WORK / "f1_independent_noise_run"
OUTPUT = INPUT / "share_figures"
PB_PATH = WORK / "bliss_f1_dirty.pb"
PB_LABELS = ["0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"]
COLORS = ["#4477AA", "#66CCEE", "#228833", "#CCBB44"]
MARKERS = ["o", "s", "^", "D"]


def load_csv(name):
    with (INPUT / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def load_pb():
    ia = image()
    try:
        ia.open(str(PB_PATH))
        return np.squeeze(np.asarray(ia.getchunk(), dtype=float))
    finally:
        ia.done()


def save(fig, stem):
    fig.savefig(OUTPUT / f"{stem}.png", dpi=240, bbox_inches="tight")
    fig.savefig(OUTPUT / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_setup(summary, realization_summary):
    mean_rms = 1e6 * float(summary["mean_blank_rms_jy_per_beam"])
    setup = [
        ("Array configuration", "alma.F1.cfg"),
        ("Antennas", "255 × 12 m"),
        ("Maximum baseline", "~0.39 km"),
        ("Declination", "−27°48′00.00″"),
        ("Integration time", "5 min"),
        ("Central frequency", "230 GHz"),
        ("Bandwidth", "7.5 GHz"),
        ("Pointing", "Single field"),
        ("Visibility noise", "10 independent Gaussian realizations"),
        ("Noise seeds", "20400000–20400009"),
        ("Injected sources", f"{summary['n_injected']} (44 per realization)"),
        ("Input grid", "11 S/N values × 4 PB bins"),
        ("Detection", "4σ seed / 3σ flood"),
        ("Mean dirty-image RMS", f"{mean_rms:.3f} µJy beam⁻¹"),
    ]
    fig, ax = plt.subplots(figsize=(9.4, 7.0))
    ax.axis("off")
    table = ax.table(
        cellText=setup,
        colLabels=["Parameter", "F1 setup"],
        cellLoc="left",
        colLoc="left",
        loc="center",
        colWidths=[0.38, 0.58],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 1.55)
    for (row, _), cell in table.get_celld().items():
        cell.set_edgecolor("#BBBBBB")
        if row == 0:
            cell.set_facecolor("#DCEAF7")
            cell.set_text_props(weight="bold")
        elif row % 2 == 0:
            cell.set_facecolor("#F3F3F3")
    ax.set_title("GHOST F1 continuum injection–recovery setup", fontsize=16, pad=18)
    ax.text(
        0.5, -0.015,
        "Idealized thermal-noise baseline; calibration residuals and other non-Gaussian systematics are not included.",
        transform=ax.transAxes, ha="center", va="top", fontsize=9, color="#555555",
    )
    save(fig, "01_f1_simulation_setup")


def plot_positions(sources, pb):
    fig, ax = plt.subplots(figsize=(7.2, 6.4), constrained_layout=True)
    ax.imshow(pb.T, origin="lower", cmap="Greys", vmin=0.2, vmax=1.0)
    batch = [row for row in sources if int(row["realization"]) == 0]
    points = None
    for index, (label, marker) in enumerate(zip(PB_LABELS, MARKERS)):
        subset = [row for row in batch if int(row["pb_bin"]) == index]
        points = ax.scatter(
            [float(row["x_pixel"]) for row in subset],
            [float(row["y_pixel"]) for row in subset],
            c=[float(row["expected_snr"]) for row in subset],
            cmap="viridis", vmin=2, vmax=7, marker=marker, s=70,
            edgecolor="white", linewidth=0.7, label=label,
        )
    ax.contour(pb.T, levels=[0.2, 0.4, 0.6, 0.8], colors="white", linewidths=0.8, alpha=0.8)
    ax.set(title="Injected source positions: one F1 realization", xlabel="x pixel", ylabel="y pixel")
    ax.legend(title="PB response", fontsize=9)
    fig.colorbar(points, ax=ax, label="Input apparent S/N", shrink=0.84)
    save(fig, "03_f1_injected_source_positions")


def plot_input_output(sources):
    fig, ax = plt.subplots(figsize=(7.2, 5.8), constrained_layout=True)
    for index, (label, color, marker) in enumerate(zip(PB_LABELS, COLORS, MARKERS)):
        subset = [
            row for row in sources
            if int(row["pb_bin"]) == index and row["recovered"].lower() == "true"
        ]
        ax.scatter(
            [float(row["expected_snr"]) for row in subset],
            [float(row["recovered_peak_snr"]) for row in subset],
            color=color, marker=marker, s=30, alpha=0.55, label=label,
        )
    ax.plot([2, 7.5], [2, 7.5], color="#444444", linestyle="--", label="Output = input")
    ax.axhline(4.0, color="#AA3377", linestyle=":", linewidth=1.5, label="4σ detection threshold")
    ax.set(
        xlim=(1.8, 7.3), ylim=(3.7, 10.0),
        xlabel="Input apparent S/N", ylabel="Recovered peak S/N",
        title="GHOST F1 recovered sources: input versus output S/N",
    )
    ax.grid(alpha=0.22)
    ax.legend(title="PB response", fontsize=9, ncol=2)
    save(fig, "04_f1_input_output_snr")


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    sources = load_csv("source_recovery_all.csv")
    realization_summary = load_csv("realization_summary.csv")
    with (INPUT / "summary.json").open() as handle:
        summary = json.load(handle)
    plot_setup(summary, realization_summary)
    plot_positions(sources, load_pb())
    plot_input_output(sources)
    figure_outputs = sorted(path for path in OUTPUT.iterdir() if path.suffix in {".png", ".pdf"})
    write_script_provenance(OUTPUT/"share_figures.provenance.json",root=ROOT,script=__file__,
        inputs=[INPUT/"source_recovery_all.csv",INPUT/"realization_summary.csv",INPUT/"summary.json",PB_PATH],
        outputs=figure_outputs,parameters={"purpose":"email-ready F1 validation figures"})
    for path in sorted(OUTPUT.iterdir()):
        print(path)


if __name__ == "__main__":
    main()
