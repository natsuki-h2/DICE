"""Plot the native full-primary-beam blank dirty-image pilot."""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np
from bliss_alma.provenance import write_script_provenance


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data/ghost_matrix/native_full_pb_blank"


def main():
    rows = list(csv.DictReader((INPUT / "native_full_pb_blank_summary.csv").open()))
    fig, axes = plt.subplots(2, 3, figsize=(14.8, 9.2), constrained_layout=True)
    for ax, row in zip(axes.ravel(), rows):
        configuration = row["configuration"]
        data = np.load(INPUT / f"{configuration.lower()}_display.npz")
        dirty = data["dirty"]
        # Each configuration requires a different block-averaging factor for
        # display.  Normalize by the displayed plane's own RMS so visual
        # contrast is not confused with the scientifically measured native
        # Jy/beam RMS reported in the separate metrics panel and CSV.
        dirty = dirty / np.nanstd(dirty)
        fov = float(row["fov_arcsec"])
        im = ax.imshow(dirty.T, origin="lower", cmap="RdBu_r", vmin=-4, vmax=4,
                       extent=[-fov/2, fov/2, -fov/2, fov/2])
        ax.set(title=(f"{configuration}: {float(row['beam_major_arcsec']):.4f} × "
                      f"{float(row['beam_minor_arcsec']):.4f} arcsec"),
               xlabel="Relative RA (arcsec)", ylabel="Relative Dec (arcsec)")
        fig.colorbar(im, ax=ax, shrink=.78, label="Display RMS")
    fig.suptitle("Native full-PB dirty images: 5 min, no taper (block-averaged for display)")
    path = INPUT / "native_full_pb_blank_images.png"
    fig.savefig(path, dpi=180); plt.close(fig)

    labels = [row["configuration"] for row in rows]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.4), constrained_layout=True)
    axes[0].plot(x, [float(row["beam_major_arcsec"]) for row in rows], "o-", label="major")
    axes[0].plot(x, [float(row["beam_minor_arcsec"]) for row in rows], "o-", label="minor")
    axes[0].set(yscale="log", ylabel="Beam FWHM (arcsec)", title="Native beam")
    axes[0].legend()
    axes[1].plot(x, [float(row["rms_jy_per_beam"])*1e6 for row in rows], "o-")
    axes[1].set(ylabel="RMS (microJy/beam)", title="Full-PB dirty-image RMS")
    axes[2].plot(x, [float(row["peak_abs_psf_outside_2beam"]) for row in rows], "o-")
    axes[2].set(ylabel="Maximum |PSF|", title="Outside two beam FWHM")
    for ax in axes:
        ax.set(xticks=x, xticklabels=labels); ax.grid(alpha=.25)
    fig.suptitle("Native-resolution F1–F6 baseline: 5 min, natural weighting, no taper")
    path2 = INPUT / "native_full_pb_blank_metrics.png"
    fig.savefig(path2, dpi=200); plt.close(fig)
    write_script_provenance(INPUT/"native_full_pb_blank.provenance.json",
        root=ROOT,script=__file__,inputs=[INPUT/"native_full_pb_blank_summary.csv"],
        outputs=[path,path2],parameters={"display":"block-averaged; native science metrics unchanged"})
    print(path); print(path2)


if __name__ == "__main__":
    main()
