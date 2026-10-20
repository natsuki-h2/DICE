# Hayatsu et al. (2025) detectability reproduction

## Scope

This test asks whether the public DICE image/cube injection workflow recovers
the line-completeness behavior reported in Hayatsu et al. (2025,
arXiv:2510.02303). It is separate from the completed catalogue-level
reproduction of Hayatsu et al. (2017).

The published experiment used:

- Gaussian line FWHM values of 2, 4, ..., 22 channels;
- input peak S/N values of 2.0, 2.5, ..., 7.0;
- matched spectral boxcar smoothing;
- ten positions over PB = 0.81--0.99 and multiple channels in each field;
- 100 different fields and 1,016,400 injections in total;
- 100 blank mock mosaics for cumulative contamination.

The local archive does not contain all 100 injected fields or blank mock
mosaics. It does contain the observed, Gaussian-replaced, and simulated
Field10 measurement sets and dirty cubes. The reproduction therefore has two
explicitly separated levels:

1. independently recalculate completeness using the Gaussian-replaced
   Field10 dirty cube as the noise template;
2. reconstruct the published cumulative-contamination fits from Table 2,
   without claiming to regenerate the unavailable 100-mosaic clump counts.

## Reproduction command

First export the locally archived CASA image to FITS if needed:

```bash
casa --nogui --nologger --log2term -c \
  "from casatasks import exportfits; exportfits(\
imagename='data/hayatsu_reproduction/adf22_mockdata_comparison/\
Field10_Spw0_rep.ms.dirtycube.image', \
fitsimage='data/hayatsu_reproduction/adf22_mockdata_comparison/\
Field10_Spw0_rep.ms.dirtycube.image.fits', overwrite=True, dropdeg=True)"
```

Then run:

```bash
PYTHONPATH=src python \
  examples/64_reproduce_hayatsu2025_detectability.py
```

The script:

- normalizes every valid spectral plane by a robust RMS;
- draws 10,000 position/channel samples from the PB-limited Field10 cube;
- injects Gaussian profiles with `inject_source_image`;
- applies a boxcar equal to the input Gaussian FWHM;
- uses a 4.0 smoothed-S/N recovery threshold;
- repeats the published 11 by 11 width--input-S/N grid;
- groups samples into 100 subsets for a descriptive 16--84% range;
- fits an increasing error function to each DICE completeness curve;
- reconstructs the published completeness and contamination functions with
  signs consistent with the plotted and tabulated curve directions.

The 4.0 threshold best reproduces the historical recovery transition. A
4.5 threshold can be selected explicitly to test the catalogue threshold used
in Hayatsu et al. (2017):

```bash
PYTHONPATH=src python \
  examples/64_reproduce_hayatsu2025_detectability.py --threshold 4.5
```

## Result

Across the 103 Table-3 cells with published completeness values:

- median `DICE - published` completeness: 0.000;
- mean absolute difference: 0.031;
- 90th percentile absolute difference: 0.085.

The fitted 50%-completeness input-S/N values agree closely over most widths:

| FWHM (channels) | Published midpoint | DICE midpoint | DICE - published |
|---:|---:|---:|---:|
| 2 | 4.71 | 4.54 | -0.17 |
| 4 | 3.61 | 3.55 | -0.06 |
| 6 | 3.03 | 2.97 | -0.06 |
| 8 | 2.96 | 2.61 | -0.35 |
| 10 | 2.40 | 2.37 | -0.03 |
| 12 | 2.15 | 2.19 | +0.04 |
| 14 | 2.16 | 2.06 | -0.10 |
| 16 | 1.98 | 1.98 | -0.00 |
| 18 | 1.00 | 1.92 | +0.92 |
| 20 | 2.35 | 1.91 | -0.45 |
| 22 | 1.71 | 1.91 | +0.20 |

The 18-channel published midpoint is poorly constrained by Table 3 because
all reported completeness values are already 0.97--1.00 over the sampled
input-S/N range; its fitted transition lies outside the measured grid.

The largest raw differences occur for 16--22-channel lines at input S/N = 2.
Published completeness is 0.83--0.97, whereas the Field10 DICE reproduction
gives approximately 0.52--0.58. Plausible contributors are:

- 100 independent fields in the publication versus one reused Field10 noise
  template here;
- visibility-domain injection in the publication versus explicit
  image/cube-plane injection here;
- the treatment of TDM/Hanning spectral correlation in the RMS after
  smoothing;
- historical `clumpfind` behavior versus the explicit DICE recovery rule;
- even-width boxcar centering and edge conventions;
- field-to-field PB, flagging, and noise differences.

These effects have not yet been uniquely separated. The discrepancy is
confined mainly to broad, lowest-input-S/N cells; the scientifically important
overall dependence of completeness on input S/N and line width is reproduced.

## Independent blank-realization test

The archived Gaussian-replaced Field10 measurement set was used to generate
100 new visibility-noise realizations with
`examples/65_generate_hayatsu2025_blank_realizations.py`. Each realization
was imaged with natural weighting, `niter=0`, 0.1-arcsec pixels, and the
archived 15.625-MHz channel width. Working measurement sets and CASA images
were deleted after recording the clump counts. Uncertainties use one complete
visibility-noise realization as the sampling unit.

All IDs 0--99 completed, giving 24,200 unique combinations of realization,
line width, S/N threshold, and polarity. Positive and negative blank counts
are closely symmetric. Representative mean positive counts are:

| Width | 4.0 sigma | 4.5 sigma | 5.0 sigma |
|---:|---:|---:|---:|
| 2 channels | 60.30 | 8.52 | 0.84 |
| 6 channels | 39.01 | 5.74 | 0.73 |
| 14 channels | 23.63 | 3.23 | 0.37 |
| 22 channels | 16.64 | 2.57 | 0.29 |

The probability of at least one positive Field10 blank clump at 5 sigma is
0.57, 0.52, 0.32, and 0.25 for these widths, respectively. It is
approximately zero by 6 sigma.

## Contamination limitation

The paper defines cumulative contamination as

```text
Q(x) = N_mock(>x) / N_real(>x)
```

for matched peak-S/N and FWHM selections. The printed signs of Eqs. 4 and 5
are reversed relative to both the plotted curves and Table 3. The
reconstruction uses a decreasing error function for Q and an increasing one
for completeness.

The 100 new Field10 realizations robustly determine a representative-pointing
blank numerator. They do not supply the full-mosaic real-detection
denominator. At high thresholds the observed Field10 counts are zero or one,
so its ratio becomes undefined, zero, or one. It is therefore not presented
as the final mosaic-scale contamination curve.

Exact full-mosaic contamination reproduction still requires recovery of the
original catalogues or regeneration of a full-mosaic blank ensemble with
search volume and PB coverage matched to the real ADF22 search. The local
four-SPW continuum-subtracted image products occupy about 67 GB, while the
internal disk currently lacks enough free space for a duplicated full set.

## Products

- `examples/64_reproduce_hayatsu2025_detectability.py`
- `examples/65_generate_hayatsu2025_blank_realizations.py`
- `examples/66_analyze_hayatsu2025_blank_realizations.py`
- `data/hayatsu_reproduction/hayatsu2025_reproduction/`
  - `hayatsu2025_dice_completeness.csv`
  - `hayatsu2025_fit_comparison.csv`
  - `hayatsu2025_field10_contamination_diagnostic.csv`
  - `hayatsu2025_reproduction_metadata.json`
  - `hayatsu2025_detectability_reproduction.png`
  - `hayatsu2025_detectability_reproduction.pdf`
- `data/hayatsu_reproduction/hayatsu2025_blank_realizations/`
  - `blank_clump_summary.csv`
  - `hayatsu2025_blank_field10_contamination.png`
  - `hayatsu2025_blank_field10_contamination.pdf`
