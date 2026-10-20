# GHOST configuration/time matrix analysis design

## Fixed GHOST inputs

All matrix jobs use the same declination, frequency setup, PWV, blank FITS
model, single pointing, natural weighting, and dirty-image/MS products.  The
only planned changes are array configuration and integration time.

## Native-resolution field-of-view constraint

GHOST selects a cell size that samples the synthesized beam and caps the
image at 2048 pixels.  This covers the full primary-beam analysis field for
F1--F4, but not for F5 and F6:

| Configuration | Cell (arcsec) | Image size | Image FOV (arcsec) | Requested PB FOV (arcsec) | Covered fraction |
|---|---:|---:|---:|---:|---:|
| F1 | 0.124541 | 360 | 44.83 | 44.83 | 1.00 |
| F2 | 0.093125 | 480 | 44.70 | 44.70 | 1.00 |
| F3 | 0.034151 | 1250 | 42.69 | 42.69 | 1.00 |
| F4 | 0.020579 | 2048 | 42.15 | 42.15 | 1.00 |
| F5 | 0.007234 | 2048 | 14.82 | 42.44 | 0.35 |
| F6 | 0.003668 | 2048 | 7.51 | 42.44 | 0.18 |

Consequently, a native-resolution F1-style comparison across PB bins
0.2--1.0 is not defined for F5/F6.  Treating missing low-PB bins as zero
completeness would be incorrect: they are outside the image, not failed
recoveries.

## Adopted two-track comparison

1. **Native-resolution track.** Preserve each configuration's natural dirty
   beam.  Report beam, RMS, and completeness only over PB bins actually
   covered.  F5/F6 results are explicitly central-PB results.
2. **Common-resolution track.** Apply a documented uv taper and common image
   grid to F1--F6 so the full PB field can be compared at a matched angular
   resolution.  This track measures how configuration/uv coverage affects
   detectability after resolution matching; it is not a native-resolution
   comparison.

The two tracks must remain separate in figures and tables.  The common-grid
settings will be validated on one realization per configuration before the
ten-realization production run.

## Noise, weights, and flags

The audited GHOST MS has constant `WEIGHT=1`, constant `SIGMA=1`, and no
flags.  The downloaded matrix therefore defines an idealized thermal-noise
and uv-sampling baseline.  Structured flag and heteroscedastic-weight tests
are separate diagnostic comparisons and must not be described as empirical ALMA
calibration residuals.

## Completed common-resolution production run

The validated common-resolution workflow uses natural weighting, a
`0.7arcsec` uv taper, the full F1 primary-beam grid, and exact smoothing to a
`1.0arcsec` circular restoring beam.  Each scenario has ten independent
visibility-noise realizations and 44 injected point sources per realization.
Uncertainty bands are obtained by resampling whole realizations, not
individual sources.

| Configuration | Time | Realizations | Injected | Recovered | Mean RMS (microJy/beam) |
|---|---:|---:|---:|---:|---:|
| F1 | 5 min | 10 | 440 | 244 | 7.615 |
| F2 | 5 min | 10 | 440 | 256 | 9.273 |
| F3 | 5 min | 10 | 440 | 258 | 12.227 |
| F4 | 5 min | 10 | 440 | 255 | 15.971 |
| F5 | 5 min | 10 | 440 | 258 | 22.128 |
| F6 | 5 min | 10 | 440 | 256 | 82.725 |
| F1 | 1 h | 10 | 440 | 264 | 2.202 |
| F1 | 5 h | 10 | 440 | 255 | 0.977 |
| F6 | 1 h | 10 | 440 | 272 | 24.332 |
| F6 | 5 h | 10 | 440 | 282 | 10.528 |

This is 100 independent noise realizations and 4,400 injected sources in
total.  For both F1 and F6, the measured RMS time dependence is close to the
thermal-noise expectation proportional to the inverse square root of
integration time.  Small non-monotonic differences in differential
completeness are compatible with the finite ten-realization samples and
threshold selection; they should not be interpreted as a time-dependent
loss of sensitivity.

Primary outputs are:

- `data/ghost_matrix/common_resolution_runs/figures/common_resolution_matrix_metrics.png`
- `data/ghost_matrix/common_resolution_runs/figures/common_resolution_time_metrics.png`
- matching PDF and CSV files in the same directory

The native-resolution products remain a separate diagnostic track because
the native F5/F6 images do not cover the full primary beam.
