"""Combine and classify the exact ADF22 positive and negative candidates."""

from pathlib import Path
import csv

from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/hayatsu_reproduction/results"
PUBLISHED = [
    ("Line A", "22:17:37.43", "+00:17:10.7", 253.79),
    ("Line B", "22:17:31.95", "+00:18:20.3", 269.92),
    ("Line C", "22:17:36.97", "+00:18:20.8", 253.49),
    ("Line D", "22:17:33.07", "+00:17:18.8", 269.70),
]
REPEATED_POSITION = SkyCoord(334.385028 * u.deg, 0.295500 * u.deg)

rows = []
for spw in range(4):
    with (OUT / f"adf22_exact_spw{spw}_blind_candidates.csv").open() as handle:
        for row in csv.DictReader(handle):
            for key in ("snr", "ra_deg", "dec_deg", "frequency_ghz", "channel"):
                row[key] = float(row[key])
            for key in ("spw", "width_channels", "x", "y"):
                row[key] = int(row[key])
            coordinate = SkyCoord(row["ra_deg"] * u.deg, row["dec_deg"] * u.deg)
            classification = "negative"
            matched_name = ""
            if row["snr"] > 0:
                classification = "other_positive"
                for name, ra, dec, frequency in PUBLISHED:
                    target = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
                    if (
                        target.separation(coordinate).arcsec <= 1.0
                        and abs(row["frequency_ghz"] - frequency) <= 0.15
                    ):
                        classification = "published_candidate"
                        matched_name = name
                        break
                if (
                    classification == "other_positive"
                    and REPEATED_POSITION.separation(coordinate).arcsec <= 1.0
                ):
                    classification = "multi_spw_positive_residual_candidate"
            row["classification"] = classification
            row["matched_name"] = matched_name
            rows.append(row)

rows.sort(key=lambda row: abs(row["snr"]), reverse=True)
table = OUT / "adf22_exact_all_spw_candidates.csv"
with table.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
thresholds = np.linspace(5.5, 8.5, 301)
positive = np.array([row["snr"] for row in rows if row["snr"] > 0])
negative = np.array([-row["snr"] for row in rows if row["snr"] < 0])
axes[0].step(
    thresholds, [(positive >= value).sum() for value in thresholds],
    where="post", color="#c83e4d", label="Positive",
)
axes[0].step(
    thresholds, [(negative >= value).sum() for value in thresholds],
    where="post", color="#377eb8", label="Negative",
)
axes[0].axvline(6.0, color="0.35", linestyle=":")
axes[0].set(
    xlabel="Peak |S/N| threshold", ylabel="Cumulative candidate count",
    xlim=(5.5, 8.5),
)
axes[0].legend(frameon=False)

styles = {
    "published_candidate": ("#c83e4d", "o", "Published candidates"),
    "other_positive": ("#e6a532", "o", "Other positive"),
    "multi_spw_positive_residual_candidate":
        ("#7b4ab3", "s", "Repeated-position positive"),
    "negative": ("#377eb8", "v", "Negative"),
}
used = set()
for row in rows:
    color, marker, label = styles[row["classification"]]
    axes[1].scatter(
        row["frequency_ghz"], row["snr"], color=color, marker=marker, s=35,
        label=label if label not in used else None,
    )
    used.add(label)
    if row["matched_name"]:
        axes[1].annotate(
            row["matched_name"], (row["frequency_ghz"], row["snr"]),
            xytext=(3, 4), textcoords="offset points", fontsize=8,
        )
axes[1].axhline(6.0, color="0.35", linestyle=":")
axes[1].axhline(-6.0, color="0.35", linestyle=":")
axes[1].set(xlabel="Observed frequency (GHz)", ylabel="Peak S/N")
axes[1].legend(
    frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5)
)
for axis in axes:
    axis.grid(alpha=0.2)
fig.suptitle(
    "ADF22 exact visibility-domain untargeted search: all four spectral windows"
)
fig.tight_layout()
png = OUT / "adf22_exact_all_spw_candidate_summary.png"
fig.savefig(png, dpi=220)
fig.savefig(OUT / "adf22_exact_all_spw_candidate_summary.pdf")
print("Wrote", table)
print("Wrote", png)
