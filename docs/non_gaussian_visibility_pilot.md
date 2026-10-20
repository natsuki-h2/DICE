# Controlled non-Gaussian visibility-noise pilot

## Purpose

This pilot tests whether replacing GHOST's independent Gaussian visibility
noise can change dirty-image detection statistics.  It is a controlled stress
test, not an empirical model of calibrated ALMA residuals.

## Setup

- F1 and F6 five-minute visibility sampling;
- native synthesized beam, natural weighting, no taper, `niter=0`;
- 1024 x 1024 pixels and four pixels per minor-axis beam FWHM;
- ten independent realizations per configuration and noise scenario;
- visibility marginal variance held fixed before imaging;
- image RMS, pixel histogram, blank-image local peaks, and isolated point
  source recovery measured separately;
- 200 isolated-source trials per input S/N and realization;
- four-sigma detection threshold.

The tested scenarios are:

1. independent complex Gaussian noise (baseline);
2. independent Student-t noise with three degrees of freedom, rescaled to
   unit variance;
3. independent Gaussian noise with 0.1% of samples multiplied by ten, then
   rescaled to unit variance;
4. a synthetic case in which 10% of the variance is a common Gaussian mode
   shared by baselines at the same timestamp.

The fourth scenario is labelled *time-coherent*. It is deliberately simple
and should not be called an ALMA calibration-error model. Although its random
variables are Gaussian, their coherence makes the image noise spatially
non-stationary; pooling pixels then produces a strongly non-Gaussian image
histogram and persistent high-S/N structures.

## Initial result

Independent non-Gaussian visibility samples are strongly Gaussianized by the
large number of visibilities contributing to each dirty-image pixel. The
Student-t and sparse-outlier scenarios have nearly the same image RMS,
pixel kurtosis, false-positive tail, and completeness as the independent
Gaussian baseline in this test.

The coherent scenario behaves differently. Relative to the Gaussian median,
the image RMS increases by approximately 21% for F1 and 3% for F6. It also
produces localized, persistent high-S/N image structures and a large pooled
pixel excess kurtosis. At fixed *input S/N defined using each image's measured
RMS*, the median completeness curve changes only modestly because the
coherent structures occupy a limited fraction of the image. The blank-image
peak tail, however, extends to much higher measured S/N.

The main conclusion is therefore not simply that a heavy-tailed visibility
PDF matters. Correlation and coherence across time, baselines, antennas, or
frequency are more capable of surviving Fourier averaging and changing
source-catalogue reliability. A realistic next stage should obtain residual
statistics from calibrated ALMA measurement sets and preserve those axes of
correlation when constructing the simulation template.

Primary outputs are under
`data/ghost_matrix/non_gaussian_visibility_pilot/figures/`.
