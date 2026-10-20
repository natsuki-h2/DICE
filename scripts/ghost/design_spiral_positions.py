"""Plot the planned spiral placement used by future S/N-PB batches."""

from pathlib import Path
import importlib.util

import matplotlib.pyplot as plt
import numpy as np
from casatools import image
from bliss_alma.provenance import write_script_provenance


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "ghost_f1_5min" / "f1_snr_pb_grid" / "figures"


def main():
    spec = importlib.util.spec_from_file_location(
        "f1grid", ROOT / "scripts" / "ghost" / "run_snr_pb_grid.py"
    )
    grid = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(grid)
    ia = image()
    try:
        ia.open(str(grid.BLANK_PB))
        pb = np.squeeze(np.asarray(ia.getchunk(), dtype=float))
        valid = np.squeeze(np.asarray(ia.getchunk(getmask=True), dtype=bool))
    finally:
        ia.done()
    groups = grid.choose_positions(pb, valid, batch=0)
    markers = ["o", "s", "^", "D"]
    labels = ["0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"]
    fig, ax = plt.subplots(figsize=(6.4, 6.0), constrained_layout=True)
    ax.imshow(pb.T, origin="lower", cmap="Greys", vmin=0.2, vmax=1.0)
    for index, (_, _, positions) in enumerate(groups):
        points = ax.scatter([p[0] for p in positions], [p[1] for p in positions],
                            c=grid.TARGET_SNRS, cmap="viridis", vmin=2, vmax=7,
                            marker=markers[index], s=65, edgecolor="white", linewidth=0.7,
                            label=labels[index])
    ax.contour(pb.T, levels=[0.2, 0.4, 0.6, 0.8], colors="white", linewidths=0.7)
    ax.set(title="Planned spiral injection pattern (one batch)", xlabel="x pixel", ylabel="y pixel")
    ax.legend(title="PB response", fontsize=8)
    fig.colorbar(points, ax=ax, label="Input S/N", shrink=0.82)
    fig.savefig(OUT / "f1_planned_spiral_positions.png", dpi=200)
    fig.savefig(OUT / "f1_planned_spiral_positions.pdf")
    plt.close(fig)
    write_script_provenance(OUT/"f1_planned_spiral_positions.provenance.json",root=ROOT,script=__file__,
        inputs=[grid.BLANK_PB],outputs=[OUT/"f1_planned_spiral_positions.png",OUT/"f1_planned_spiral_positions.pdf"],
        parameters={"batch":0,"target_snrs":list(map(float,grid.TARGET_SNRS))})
    print(OUT / "f1_planned_spiral_positions.png")


if __name__ == "__main__":
    main()
