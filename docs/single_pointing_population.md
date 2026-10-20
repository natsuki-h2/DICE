# Single-pointing source-population test

`examples/06_casa_single_pointing_population.py` is an end-to-end CASA test
for measuring source recovery across a single primary beam efficiently.

## Workflow

1. Place equal-intrinsic-flux point sources on a uniform Cartesian grid within
   one pointing.
2. Create populated and blank one-channel continuum sky models.
3. Run `simobserve` with an explicit single-pointing file and identical noise
   settings for both models.
4. Form one uncorrected dirty image for each measurement set.
5. Extract all positive islands from the populated image in one pass.
6. Extract positive and negative islands from the blank image with the same
   thresholds.
7. Match the populated catalogue to the known input positions and report
   completeness in primary-beam bins.

The populated sky model is passed through `simobserve`, whose voltage-pattern
model applies the primary-beam attenuation. A plain `ft` injection is not used
for this test because the goal is specifically to sample position-dependent
primary-beam response.

## Run

From the repository root using CASA 6:

```bash
casa --nogui --nologger --log2term \
  -c examples/06_casa_single_pointing_population.py
```

Large CASA products are written under
`data/casa_single_pointing_population/` and ignored by Git. The local recovery
catalogue is `source_recovery.csv`.

## Initial validation result

The first CASA 6.7.5.18 run used 21 sources with 1 mJy intrinsic flux, a 5-RMS
seed, a 3-RMS flood threshold, and a PB cutoff of 0.2. Six sources were
recovered, all six detections matched an input source, and no positive or
negative islands were found in the blank field.

Completeness increased with PB response:

| PB range | Recovered | Completeness |
|---|---:|---:|
| 0.2-0.4 | 0/12 | 0.00 |
| 0.4-0.6 | 1/4 | 0.25 |
| 0.6-0.8 | 4/4 | 1.00 |
| 0.8-1.0 | 1/1 | 1.00 |

This is a smoke test, not a final completeness measurement. The next version
should repeat the experiment over multiple noise realizations and intrinsic
fluxes, and should report the PB-corrected recovered-to-injected flux ratio in addition to detection
completeness.

## Ten-realization Monte Carlo smoke test

The first Monte Carlo runner is:

```bash
casa --nogui --nologger --log2term \
  -c examples/07_casa_single_pointing_monte_carlo.py
```

It runs the validated single-pointing experiment ten times with different
thermal-noise seeds. Each realization retains its source-recovery and
detection catalogues. It also creates three aggregate tables:

- `realization_summary.csv`;
- `source_recovery_all.csv`;
- `detections_all.csv`.

The detection table distinguishes populated, blank-positive, and
blank-negative detections and records PB and peak S/N. It therefore contains
the inputs needed for both blank-derived and injection-matched contamination
estimates.

The smoke test completed all ten realizations in CASA 6.7.5.18. Across 210
injected sources, 68 were recovered. The aggregate PB-dependent completeness
was:

| PB range | Recovered | Completeness |
|---|---:|---:|
| 0.2-0.4 | 1/120 | 0.008 |
| 0.4-0.6 | 19/40 | 0.475 |
| 0.6-0.8 | 38/40 | 0.950 |
| 0.8-1.0 | 10/10 | 1.000 |

There were 70 detections in the populated images: 68 matched and two
unmatched. No positive or negative blank-field detections occurred at the
current 5-RMS seed threshold. Thus the first smoke test already gives a
non-zero injection-matched contamination fraction while the blank-derived estimate
is zero. More realizations and a flux/S/N grid are required to determine
whether this difference is statistical or systematic.

## S/N-dependent diagnostic figure

`examples/08_plot_single_pointing_snr.py` reads only the aggregate CSV files;
it does not rerun CASA or create new measurement sets. It produces a four-panel
diagnostic containing differential and cumulative completeness and
contamination. PB bins are identified by colour, individual realizations are
shown faintly, the mean is shown by a heavier line, and the 16th-84th
percentile interval is shaded. Solid and dashed contamination curves represent
injection-matched and blank-derived estimates, respectively.

The heavy central curve is a pooled count ratio, not an unweighted mean of
per-realization ratios. For example, cumulative completeness is calculated as
`sum(N_recovered) / sum(N_injected)` above each threshold across all
realizations. The faint curves retain the per-realization ratios, and their
16th-84th percentiles describe realization-to-realization scatter.

The initial figure is intentionally a smoke-test diagnostic. Because all
injected sources currently have the same intrinsic flux, expected S/N and PB
are strongly correlated: each PB colour occupies a different S/N range. The
next flux-grid experiment must make the PB-coloured curves overlap in S/N so
that residual PB-dependent behaviour can be tested at fixed expected S/N.

The plotting script writes PNG, PDF, and metric CSV outputs below the ignored
local `data/` directory. It can be run with any current Python environment that
provides NumPy and Matplotlib.

## Flux-grid diagnostic

`examples/09_single_pointing_flux_grid.py` creates an inexpensive flux-grid
smoke test from the ten existing image realizations. It uses dirty-image
linearity to separate the simulated source response from the blank image and
rescales that response over 0.25-5 mJy. This makes the PB-coloured curves
overlap in expected S/N and is useful for diagnosing the plotting and
extraction logic without creating additional CASA products.

This is not the final visibility-domain validation. A small number of flux
points must be checked with direct injection into copies of the same blank MS,
and the final ALMA2040 analysis should use the GHOST blank MS plus
`componentlist`, `ft`, and `uvsub(reverse=True)`.

The primary figures should show differential completeness and differential
source statistics. To retain the Hayatsu et al. (2025) definition,
contamination Q should additionally be reported as the cumulative ratio
`N_mock(>=S/N) / N_real(>=S/N)`. The current injection-matched false fraction is a
separate validation quantity and should be labelled accordingly.
