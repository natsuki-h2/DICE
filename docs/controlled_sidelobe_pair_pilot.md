# Controlled bright-neighbour pilot

This pilot measures how a neighbouring bright point source changes recovery
of a target point source in a native-resolution dirty image.  It is a
controlled PSF-response experiment, not yet a realistic sky-population
simulation.

## Fixed setup

- configurations: F1, F5, F6
- five-minute GHOST visibility templates
- natural weighting, no uv taper, `niter=0`
- native synthesized beam sampled at four pixels per minor-axis FWHM
- target S/N: 3, 4, 5
- neighbour S/N: 10, 30, 100
- separation: 2, 5, 10, 20 synthesized beams
- eight position angles
- ten paired independent visibility-noise realizations
- 2,880 paired trials per configuration; 8,640 total

The target-only and target-plus-neighbour cases use the same noise image.
Because dirty imaging is linear, the parameter grid is evaluated by shifting
the measured dirty PSF.  A direct visibility-domain injection of a 4-sigma
target and 100-sigma neighbour at ten beams validates this construction: the
image-domain versus visibility-domain difference RMS is 0.0018 sigma for F5
and 0.0013 sigma for F6.

## Detection accounting using the injected-source catalogue

Because this is a mock observation, the injected source catalogue is known.
The primary contamination statistic is therefore measured directly rather
than inferred from an inverted image:

`injection-matched contamination = N(unmatched detections) / N(all detections)`.

A detection is matched to an injected source within the adopted matching
radius; every other detection is counted as unmatched.  The unmatched count
is reported alongside the fraction because a two-source cutout is much
sparser than a realistic sky, making the fraction strongly dependent on the
assumed source surface density.

At the target position the paired target-only, bright-only, and combined
images additionally separate the outcomes into:

- robust target recovery: recovered with and without the bright neighbour;
- sidelobe-boosted recovery: missed without the neighbour but recovered in
  the combined image, without a bright-only false association;
- sidelobe-suppressed recovery: recovered without the neighbour but missed
  in the combined image;
- sidelobe-only false association: the bright source alone creates a
  threshold-crossing peak at the nominal target position;
- missed: none of the above and the target remains undetected.

This decomposition distinguishes a genuine target whose measured peak is
shifted by a sidelobe from a pure sidelobe peak incorrectly associated with
the target.  Inverted-image counts remain useful when applying DICE to real
data without a complete injected-source catalogue, but they are not required for this mock test.

## Initial result at target S/N = 4

F1 is weakly affected: neighbour S/N up to 30 changes completeness by zero in
this finite pilot, while a 100-sigma neighbour changes completeness by only a
few percentage points depending on separation and position angle.

F5 and F6 are strongly affected.  Examples of the median target peak change
and completeness change are:

| Config. | Neighbour S/N | Separation | Median peak change | Completeness change |
|---|---:|---:|---:|---:|
| F5 | 10 | 2 beams | +1.84 sigma | +0.40 |
| F5 | 100 | 10 beams | +0.46 sigma | +0.21 |
| F5 | 100 | 20 beams | +0.80 sigma | +0.24 |
| F6 | 10 | 2 beams | +1.32 sigma | +0.31 |
| F6 | 100 | 10 beams | +3.70 sigma | +0.39 |
| F6 | 100 | 20 beams | +0.32 sigma | +0.09 |

The distance dependence is non-monotonic because it samples the dirty PSF,
not a smooth Gaussian beam.  Positive completeness changes do not represent
improved intrinsic sensitivity: a neighbour's dirty-beam response boosts the
measured target peak and can move it above the four-sigma threshold.

## Limits

- target and neighbour are point sources near the field centre;
- only eight position angles and ten noise realizations are sampled;
- only one neighbour is present;
- PB dependence, extended structure, calibration residuals, source counts,
  and clustering are not included;
- the pilot quantifies controlled dirty-PSF coupling and should not yet be
  presented as a population-averaged contamination correction.

Primary figures and tables are under
`data/ghost_matrix/controlled_sidelobe_pair_test/figures/`.
