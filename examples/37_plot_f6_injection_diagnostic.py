"""Plot NumPy arrays exported by the F6 CASA injection diagnostic."""

from pathlib import Path
import csv
import shutil

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "data/ghost_matrix/common_resolution_runs/f6_5min"
WORK = RUN / "visualization_work"
OUT = RUN / "figures"
CELL_ARCSEC = 0.12454065519280681
IMSIZE = 360

populated = np.load(WORK / "populated.npy")
source_only = np.load(WORK / "source_only.npy")
noise_only = populated - source_only
psf = np.load(WORK / "psf.npy")
rms = float(np.nanstd(noise_only))
rows = [row for row in csv.DictReader((RUN / "source_recovery_all.csv").open())
        if int(row["realization"]) == 0]
extent = np.array([-IMSIZE/2, IMSIZE/2, -IMSIZE/2, IMSIZE/2]) * CELL_ARCSEC

fig, axes = plt.subplots(2, 2, figsize=(11.2, 9.2), constrained_layout=True)
panels = [(populated/rms, "Injected sources + noise", -5, 8, "S/N"),
          (source_only/rms, "Injected sources only", -3, 8, "Response / noise RMS"),
          (noise_only/rms, "Noise realization", -5, 5, "S/N"),
          (psf, "F6 dirty PSF after 0.7-arcsec taper", -0.15, 1, "Normalized response")]
for ax, (data, title, vmin, vmax, label) in zip(axes.ravel(), panels):
    im = ax.imshow(data.T, origin="lower", extent=extent, cmap="RdBu_r",
                   vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set(title=title, xlabel="Relative RA (arcsec)", ylabel="Relative Dec (arcsec)")
    fig.colorbar(im, ax=ax, shrink=.84, label=label)
for row in rows:
    x = (float(row["x_pixel"]) - IMSIZE/2) * CELL_ARCSEC
    y = (float(row["y_pixel"]) - IMSIZE/2) * CELL_ARCSEC
    for ax in axes.ravel()[:2]:
        ax.plot(x, y, marker="o", ms=3, mfc="none", mec="black", mew=.45)
fig.suptitle("F6 5-min realization 0 at common 1-arcsec resolution")
OUT.mkdir(parents=True, exist_ok=True)
path = OUT / "f6_injection_dynamic_range_diagnostic.png"
fig.savefig(path, dpi=200); plt.close(fig)
print(path)
shutil.rmtree(WORK)
