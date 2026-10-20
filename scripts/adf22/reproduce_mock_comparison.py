"""Reproduce the ADF22 real/Gaussian/simulated-noise comparison safely.

Run with CASA 6, for example:

    casa --nologger --nogui -c scripts/adf22/reproduce_mock_comparison.py

Unlike the historical script, this program is read-only with respect to the
input measurement sets and CASA images.
"""

from pathlib import Path
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from bliss_alma.provenance import write_script_provenance
from casatools import image as image_tool
from casatools import table as table_tool


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/hayatsu_reproduction/adf22_mockdata_comparison"
OUT = ROOT / "data/hayatsu_reproduction/results"
OUT.mkdir(parents=True, exist_ok=True)

DATASETS = [
    ("Observed", "Field10_Spw0.ms", "Field10_Spw0.ms.dirtycube.image", "#c83e4d"),
    ("Gaussian-replaced", "Field10_Spw0_rep.ms",
     "Field10_Spw0_rep.ms.dirtycube.image", "#377eb8"),
    ("Simulated", "sim.SB0_Field10_Spw0.ms",
     "sim.SB0_Field10_Spw0.ms.dirtycube.image", "#4daf4a"),
]


def finite(values):
    values = np.asarray(values)
    return values[np.isfinite(values) & (values != 0)]


def robust_sigma(values):
    values = finite(values)
    median = np.median(values)
    return 1.4826 * np.median(np.abs(values - median))


rows = []
products = {}
ia = image_tool()
tb = table_tool()

for label, ms_name, image_name, color in DATASETS:
    image_path = str(DATA / image_name)
    ia.open(image_path)
    cube = np.squeeze(ia.getchunk()).astype(np.float64)
    beam = ia.restoringbeam(channel=cube.shape[-1] // 2, polarization=0)
    ia.close()

    # CASA returns spatial axes first and frequency last for these images.
    channel_sigma = np.array([robust_sigma(cube[..., chan])
                              for chan in range(cube.shape[-1])])
    snr = cube / channel_sigma[np.newaxis, np.newaxis, :]
    snr_values = finite(snr)

    tb.open(str(DATA / ms_name))
    column = "CORRECTED_DATA" if "CORRECTED_DATA" in tb.colnames() else "DATA"
    vis = tb.getcol(column)
    flags = tb.getcol("FLAG")
    tb.close()
    vis = vis[~flags]

    major = beam.get("major", {}).get("value", np.nan)
    minor = beam.get("minor", {}).get("value", np.nan)
    rows.append({
        "dataset": label,
        "image_shape": "x".join(str(v) for v in cube.shape),
        "beam_major_arcsec": major,
        "beam_minor_arcsec": minor,
        "median_channel_sigma_jy_beam": np.nanmedian(channel_sigma),
        "snr_gt_4_pixels": int(np.sum(snr_values > 4)),
        "snr_lt_minus4_pixels": int(np.sum(snr_values < -4)),
        "snr_gt_5_pixels": int(np.sum(snr_values > 5)),
        "snr_lt_minus5_pixels": int(np.sum(snr_values < -5)),
        "visibility_real_std_jy": np.std(vis.real),
        "visibility_imag_std_jy": np.std(vis.imag),
    })
    products[label] = {
        "color": color,
        "vis_real": vis.real,
        "snr": snr_values,
        "channel_sigma": channel_sigma,
    }

with (OUT / "hayatsu_mock_comparison_summary.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.7))
for label, product in products.items():
    color = product["color"]
    axes[0].hist(product["vis_real"], bins=160, range=(-2, 2),
                 histtype="step", density=True, color=color, label=label)
    axes[1].hist(product["snr"], bins=160, range=(-8, 8),
                 histtype="step", density=True, color=color, label=label)
    axes[2].plot(product["channel_sigma"] * 1e3, color=color, label=label)

x = np.linspace(-8, 8, 400)
axes[1].plot(x, np.exp(-0.5 * x**2) / np.sqrt(2 * np.pi),
             color="0.25", linestyle=":", label="Unit Gaussian")
axes[0].set(xlabel="Visibility real part (Jy)", ylabel="Probability density",
            yscale="log", ylim=(1e-5, None))
axes[1].set(xlabel="Image-plane S/N", ylabel="Probability density",
            yscale="log", ylim=(1e-7, 1))
axes[2].set(xlabel="Frequency channel", ylabel="Robust RMS (mJy beam$^{-1}$)")
for axis in axes:
    axis.grid(alpha=0.2)
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "hayatsu_mock_comparison.png", dpi=220)
fig.savefig(OUT / "hayatsu_mock_comparison.pdf")
write_script_provenance(OUT/"hayatsu_mock_comparison.provenance.json",root=ROOT,script=__file__,
    inputs=[DATA],outputs=[OUT/"hayatsu_mock_comparison_summary.csv",OUT/"hayatsu_mock_comparison.png",OUT/"hayatsu_mock_comparison.pdf"],
    parameters={"comparison":"historical versus Gaussian-replaced visibility noise"})
print("Wrote", OUT / "hayatsu_mock_comparison_summary.csv")
print("Wrote", OUT / "hayatsu_mock_comparison.png")
