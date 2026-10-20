"""Plot the F1 injection layout and input-to-output S/N selection bias."""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np
from casatools import image
from bliss_alma.provenance import write_script_provenance


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "data" / "ghost_f1_5min"
GRID = WORK / "f1_snr_pb_grid"
OUT = GRID / "figures"
PB_PATH = WORK / "bliss_f1_dirty.pb"
PB_LABELS = ["0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"]
MARKERS = ["o", "s", "^", "D"]


def rows():
    with (GRID / "source_recovery_all.csv").open(newline="") as handle:
        return list(csv.DictReader(handle))


def pb_plane():
    ia = image()
    try:
        ia.open(str(PB_PATH))
        return np.squeeze(np.asarray(ia.getchunk(), dtype=float))
    finally:
        ia.done()


def main():
    data = rows()
    pb = pb_plane()
    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), constrained_layout=True)

    ax = axes[0]
    im = ax.imshow(pb.T, origin="lower", cmap="Greys", vmin=0.2, vmax=1.0)
    batch = [r for r in data if int(r["realization"]) == 0]
    for pb_index, (label, marker) in enumerate(zip(PB_LABELS, MARKERS)):
        subset = [r for r in batch if int(r["pb_bin"]) == pb_index]
        points = ax.scatter(
            [float(r["x_pixel"]) for r in subset], [float(r["y_pixel"]) for r in subset],
            c=[float(r["expected_snr"]) for r in subset], cmap="viridis", vmin=2, vmax=7,
            marker=marker, s=58, edgecolor="white", linewidth=0.6, label=label,
        )
    ax.contour(pb.T, levels=[0.2, 0.4, 0.6, 0.8], colors="white", linewidths=0.7, alpha=0.7)
    ax.set(title="One simultaneous batch: 44 injected sources", xlabel="x pixel", ylabel="y pixel")
    ax.legend(title="PB", loc="upper right", fontsize=8)
    fig.colorbar(points, ax=ax, label="Input S/N", shrink=0.83)

    ax = axes[1]
    for pb_index, (label, marker) in enumerate(zip(PB_LABELS, MARKERS)):
        subset = [r for r in data if int(r["pb_bin"]) == pb_index and r["recovered"].lower() == "true"]
        x = np.array([float(r["expected_snr"]) for r in subset])
        y = np.array([float(r["recovered_peak_snr"]) for r in subset])
        ax.scatter(x, y, marker=marker, s=25, alpha=0.5, label=label)
    ax.plot([2, 7.5], [2, 7.5], color="#444444", linestyle="--", label="output = input")
    ax.axhline(4.0, color="#AA3377", linestyle=":", label="4σ selection")
    ax.set(xlim=(1.8, 7.3), ylim=(3.7, 10.0), xlabel="Input apparent S/N",
           ylabel="Recovered peak S/N", title="Detected sources: threshold-selection bias")
    ax.grid(alpha=0.2)
    ax.legend(title="PB", fontsize=8, ncol=2)

    fig.savefig(OUT / "f1_input_positions_and_snr_bias.png", dpi=180)
    fig.savefig(OUT / "f1_input_positions_and_snr_bias.pdf")
    plt.close(fig)
    write_script_provenance(OUT/"f1_input_positions_and_snr_bias.provenance.json",root=ROOT,script=__file__,
        inputs=[GRID/"source_recovery_all.csv",PB_PATH],outputs=[OUT/"f1_input_positions_and_snr_bias.png",OUT/"f1_input_positions_and_snr_bias.pdf"],
        parameters={"pb_bins":PB_LABELS,"detection_threshold_sigma":4.0})
    print(OUT / "f1_input_positions_and_snr_bias.png")


if __name__ == "__main__":
    main()
