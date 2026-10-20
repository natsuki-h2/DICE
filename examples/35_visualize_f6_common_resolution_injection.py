"""Visualize one F6 common-resolution injection realization and its PSF."""

from pathlib import Path
import csv
import shutil

import matplotlib.pyplot as plt
import numpy as np
from casatasks import ft, imsmooth, tclean, uvsub
from casatools import componentlist, image, table


ROOT = Path(__file__).resolve().parents[1]
BASE_MS = (ROOT / "data/ghost_matrix/f6_5min/extracted" /
           "output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms")
RUN = ROOT / "data/ghost_matrix/common_resolution_runs/f6_5min"
OUT = RUN / "figures"
WORK = RUN / "visualization_work"
CELL = "0.12454065519280681arcsec"
IMSIZE = 360
UVTAPER = ["0.7arcsec"]
TARGET_BEAM = "1.0arcsec"
SEED = 20410600
CHUNK = 100_000


def remove(path):
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def read_plane(path):
    ia = image()
    try:
        ia.open(str(path))
        return np.squeeze(np.asarray(ia.getchunk(), dtype=float))
    finally:
        ia.done()


def image_ms(ms, prefix):
    tclean(vis=str(ms), imagename=str(prefix), specmode="mfs",
           imsize=[IMSIZE, IMSIZE], cell=CELL, weighting="natural",
           uvtaper=UVTAPER, robust=0.5, niter=0, pblimit=0.2,
           interactive=False)
    imsmooth(imagename=str(prefix) + ".image",
             outfile=str(prefix) + ".smooth.image", kernel="gauss",
             major=TARGET_BEAM, minor=TARGET_BEAM, pa="0deg",
             targetres=True, overwrite=True)


def add_sources(ms, rows, component_path):
    cl = componentlist()
    try:
        for row in rows:
            cl.addcomponent(dir=row["direction"], flux=float(row["apparent_flux_jy"]),
                            fluxunit="Jy", shape="point")
        cl.rename(str(component_path))
    finally:
        cl.done()
    ft(vis=str(ms), spw="0:0", complist=str(component_path), usescratch=True)
    uvsub(vis=str(ms), reverse=True)


def main():
    remove(WORK)
    WORK.mkdir(parents=True)
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [row for row in csv.DictReader((RUN / "source_recovery_all.csv").open())
            if int(row["realization"]) == 0]

    # Measure the template visibility dispersion without loading the full MS.
    tb = table(); tb.open(str(BASE_MS))
    try:
        nrows = tb.nrows(); sample = tb.getcol("DATA", startrow=0, nrow=1)
        npol = sample.shape[0]; count_total = 0
        sums = [np.zeros(npol) for _ in range(4)]
        for start in range(0, nrows, CHUNK):
            count = min(CHUNK, nrows - start)
            data = tb.getcol("DATA", startrow=start, nrow=count)
            real = data.real.reshape(npol, -1); imag = data.imag.reshape(npol, -1)
            count_total += real.shape[1]
            sums[0] += real.sum(1); sums[1] += np.square(real).sum(1)
            sums[2] += imag.sum(1); sums[3] += np.square(imag).sum(1)
    finally:
        tb.close()
    real_sigma = np.sqrt(sums[1] / count_total - np.square(sums[0] / count_total))
    imag_sigma = np.sqrt(sums[3] / count_total - np.square(sums[2] / count_total))

    noise_ms = WORK / "noise.ms"; source_ms = WORK / "source.ms"
    shutil.copytree(BASE_MS, noise_ms); shutil.copytree(BASE_MS, source_ms)
    rng = np.random.default_rng(SEED)
    tb = table(); tb.open(str(noise_ms), nomodify=False)
    try:
        for start in range(0, tb.nrows(), CHUNK):
            count = min(CHUNK, tb.nrows() - start); shape = (npol, 1, count)
            values = (rng.normal(size=shape) * real_sigma[:, None, None] +
                      1j * rng.normal(size=shape) * imag_sigma[:, None, None])
            tb.putcol("DATA", values, startrow=start, nrow=count)
            tb.putcol("CORRECTED_DATA", values, startrow=start, nrow=count)
    finally:
        tb.close()
    tb = table(); tb.open(str(source_ms), nomodify=False)
    try:
        for start in range(0, tb.nrows(), CHUNK):
            count = min(CHUNK, tb.nrows() - start)
            zeros = np.zeros((npol, 1, count), dtype=complex)
            tb.putcol("DATA", zeros, startrow=start, nrow=count)
            tb.putcol("CORRECTED_DATA", zeros, startrow=start, nrow=count)
    finally:
        tb.close()

    add_sources(noise_ms, rows, WORK / "noise_sources.cl")
    add_sources(source_ms, rows, WORK / "sources_only.cl")
    image_ms(noise_ms, WORK / "populated")
    image_ms(source_ms, WORK / "source_only")
    populated = read_plane(WORK / "populated.smooth.image")
    source_only = read_plane(WORK / "source_only.smooth.image")
    noise_only = populated - source_only
    psf = read_plane(WORK / "populated.psf")
    rms = float(np.nanstd(noise_only))

    extent = np.array([-IMSIZE/2, IMSIZE/2, -IMSIZE/2, IMSIZE/2]) * float(CELL[:-6])
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 9.2), constrained_layout=True)
    panels = [(populated/rms, "Injected sources + noise", -5, 8),
              (source_only/rms, "Injected sources only", -3, 8),
              (noise_only/rms, "Noise realization", -5, 5),
              (psf, "F6 dirty PSF after 0.7-arcsec taper", -0.15, 1)]
    for ax, (data, title, vmin, vmax) in zip(axes.ravel(), panels):
        im = ax.imshow(data.T, origin="lower", extent=extent, cmap="RdBu_r",
                       vmin=vmin, vmax=vmax, interpolation="nearest")
        ax.set(title=title, xlabel="Relative RA (arcsec)", ylabel="Relative Dec (arcsec)")
        fig.colorbar(im, ax=ax, shrink=.84, label="S/N" if title != panels[-1][1] else "Normalized response")
    for row in rows:
        x = (float(row["x_pixel"]) - IMSIZE/2) * float(CELL[:-6])
        y = (float(row["y_pixel"]) - IMSIZE/2) * float(CELL[:-6])
        for ax in axes.ravel()[:2]:
            ax.plot(x, y, marker="o", ms=3, mfc="none", mec="black", mew=.45)
    fig.suptitle("F6 5-min realization 0 at common 1-arcsec resolution")
    path = OUT / "f6_injection_dynamic_range_diagnostic.png"
    fig.savefig(path, dpi=200); plt.close(fig)
    remove(WORK)
    print(path)


if __name__ == "__main__":
    main()
