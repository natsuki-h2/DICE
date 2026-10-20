"""Fast controlled smoke tests for four planned DICE extensions.

These tests validate expected selection effects and output conventions. They
are idealized image/cube-domain calculations, not visibility-domain forecasts.
"""

from pathlib import Path
import csv
import json

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/ghost_matrix/extension_smoke_tests"
FIG = OUT / "figures"
RNG = np.random.default_rng(20490000)
N = 20000
THRESHOLD = 4.0
SNRS = np.arange(2.0, 8.01, 0.5)


def recovery_curve(signal_factor=1.0):
    return np.array([np.mean(RNG.normal(s*signal_factor, 1.0, N) >= THRESHOLD) for s in SNRS])


def line_profile(fwhm):
    if fwhm <= 1:
        return np.array([1.0])
    x = np.arange(-max(4, int(np.ceil(3*fwhm))), max(4, int(np.ceil(3*fwhm)))+1)
    profile = np.exp(-4*np.log(2)*x*x/fwhm**2)
    return profile / np.sqrt(np.sum(profile**2))


def write_rows(path, rows):
    fieldnames = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames); writer.writeheader(); writer.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    rows = []
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.2), constrained_layout=True)

    # 1. Extended continuum: intrinsic Gaussian FWHM in beam units. For fixed
    # point-source-equivalent integrated flux, beam dilution reduces peak S/N.
    ax = axes[0, 0]
    for size, color in zip((0., .5, 1., 2.), ("black", "#0072B2", "#E69F00", "#CC79A7")):
        factor = 1.0/(1.0 + size**2)
        curve = recovery_curve(factor)
        ax.plot(SNRS, curve, "o-", ms=3, color=color, label=f"source FWHM = {size:g} beam")
        for snr, value in zip(SNRS, curve):
            rows.append({"test":"extended_continuum", "parameter":size, "parameter_name":"source_fwhm_beams",
                         "input_snr":snr, "method":"peak", "completeness":value})
    ax.set_title("Extended continuum: peak dilution")
    ax.set_xlabel("Point-source-equivalent input S/N"); ax.set_ylabel("Completeness")
    ax.legend(fontsize=8)

    # 2. Linear polarization: Ricean P=sqrt(Q^2+U^2), simple one-sigma
    # debiasing, and a fixed fractional I->Q/U leakage vector.
    ax = axes[0, 1]; i_snr = 20.0; p_values = np.linspace(0, .20, 21)
    for leakage, color in zip((0., .005, .01, .03), ("black", "#0072B2", "#E69F00", "#CC79A7")):
        medians = []; detections = []
        for p in p_values:
            angle = RNG.uniform(0, np.pi, N); leak_angle = RNG.uniform(0, np.pi, N)
            q = i_snr*p*np.cos(2*angle) + i_snr*leakage*np.cos(2*leak_angle) + RNG.normal(size=N)
            u = i_snr*p*np.sin(2*angle) + i_snr*leakage*np.sin(2*leak_angle) + RNG.normal(size=N)
            observed = np.hypot(q, u); debiased = np.sqrt(np.maximum(observed**2 - 1.0, 0.0))/i_snr
            medians.append(np.median(debiased)); detections.append(np.mean(observed >= THRESHOLD))
            rows.append({"test":"linear_polarization", "parameter":leakage,
                         "parameter_name":"fractional_leakage", "input_snr":p*i_snr,
                         "method":"P_debiased", "completeness":detections[-1],
                         "median_recovered_fraction":medians[-1]})
        ax.plot(100*p_values, 100*np.asarray(medians), "o-", ms=3, color=color, label=f"leakage {100*leakage:g}%")
    ax.plot([0, 20], [0, 20], color="0.6", ls="--", label="one-to-one")
    ax.set_title("Linear polarization: P bias and leakage")
    ax.set_xlabel("Input polarization fraction (%)"); ax.set_ylabel("Median debiased fraction (%)")
    ax.legend(fontsize=8)

    # 3. Point spectral lines: compare a single-channel peak search with the
    # correctly matched spectral filter at fixed integrated matched-filter S/N.
    ax = axes[1, 0]
    for width, color in zip((1, 3, 5, 9), ("black", "#0072B2", "#E69F00", "#CC79A7")):
        profile = line_profile(width); peak_factor = float(profile.max())
        peak_curve = recovery_curve(peak_factor); matched_curve = recovery_curve(1.0)
        ax.plot(SNRS, peak_curve, "o-", ms=3, color=color, label=f"peak, FWHM={width} ch")
        ax.plot(SNRS, matched_curve, color=color, ls="--", alpha=.65)
        for snr, peak, matched in zip(SNRS, peak_curve, matched_curve):
            rows.append({"test":"point_spectral_line", "parameter":width,
                         "parameter_name":"line_fwhm_channels", "input_snr":snr,
                         "method":"peak_channel", "completeness":peak})
            rows.append({"test":"point_spectral_line", "parameter":width,
                         "parameter_name":"line_fwhm_channels", "input_snr":snr,
                         "method":"matched_spectral", "completeness":matched})
    ax.set_title("Point spectral lines: solid=peak, dashed=matched")
    ax.set_xlabel("Integrated matched-filter input S/N"); ax.set_ylabel("Completeness")
    ax.legend(fontsize=8, ncol=2)

    # 4. Extended lines: deliberately use a point+single-channel detector and
    # compare it with an ideal 3-D matched filter.
    ax = axes[1, 1]
    cases = ((0., 1), (1., 3), (1., 5), (2., 5))
    for (size, width), color in zip(cases, ("black", "#0072B2", "#E69F00", "#CC79A7")):
        spatial_factor = 1.0/(1.0 + size**2)
        spectral_factor = float(line_profile(width).max())
        simple = recovery_curve(spatial_factor*spectral_factor); matched = recovery_curve(1.0)
        ax.plot(SNRS, simple, "o-", ms=3, color=color, label=f"{size:g} beam, {width} ch")
        ax.plot(SNRS, matched, color=color, ls="--", alpha=.65)
        for snr, value, ideal in zip(SNRS, simple, matched):
            rows.append({"test":"extended_spectral_line", "parameter":f"{size:g}/{width}",
                         "parameter_name":"source_fwhm_beams/line_fwhm_channels", "input_snr":snr,
                         "method":"point_peak", "completeness":value})
            rows.append({"test":"extended_spectral_line", "parameter":f"{size:g}/{width}",
                         "parameter_name":"source_fwhm_beams/line_fwhm_channels", "input_snr":snr,
                         "method":"matched_3d", "completeness":ideal})
    ax.set_title("Extended lines: solid=point/peak, dashed=3-D matched")
    ax.set_xlabel("Integrated matched-filter input S/N"); ax.set_ylabel("Completeness")
    ax.legend(fontsize=8)

    for ax in axes.ravel():
        ax.grid(alpha=.2)
        if "polarization" not in ax.get_title().lower(): ax.set_ylim(-.02, 1.02)
    fig.suptitle("DICE extension smoke tests (idealized image/cube domain)\n"
                 "4σ threshold; 20,000 trials per grid point")
    png = FIG/"extension_smoke_tests.png"
    fig.savefig(png, dpi=200); fig.savefig(png.with_suffix(".pdf")); plt.close(fig)
    write_rows(OUT/"extension_smoke_test_results.csv", rows)
    summary = {"status":"PASS", "n_trials_per_grid_point":N, "threshold_sigma":THRESHOLD,
        "tests":["extended_continuum", "linear_polarization", "point_spectral_line", "extended_spectral_line"],
        "scope":"idealized image/cube-domain implementation and bias smoke tests",
        "not_included":["visibility-domain extended models", "cross-hand calibration", "frequency-dependent uv coverage"]}
    (OUT/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(png)


if __name__ == "__main__": main()
