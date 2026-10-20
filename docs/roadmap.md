# Roadmap

## Planned time-domain and transient extension

DICE's injection--recovery design can be extended from a single image or cube
to a sequence of epochs or time bins. A validated implementation should:

- inject persistent, variable, appearing, and fading point or resolved sources
  with explicit light curves;
- support visibility-domain injection when time-resolved measurement sets are
  available, and image/difference-plane injection otherwise;
- run an unchanged transient finder on both positive and negative difference
  images;
- report completeness and contamination versus peak or fluence, duration,
  cadence, primary-beam response, and sky position;
- preserve whole-epoch or whole-execution-block resampling so calibration
  residuals and correlated artifacts are not treated as independent pixels;
- distinguish transient recovery from variability of persistent sources and
  artifacts caused by unequal uv coverage, beams, flags, or calibration
  between epochs.

This is a future extension, not a currently validated DICE capability.
The existing 3-D connected-component extractor is nevertheless axis-agnostic,
so an advanced user can experimentally replace the spectral axis with time.
This should be documented as user-validated exploratory use, not as evidence
that DICE currently provides a calibrated transient-search pipeline.

## 2026-07-22 image injection and absorption support

- Added CASA-independent 2-D image and 3-D cube injection, including additive
  emission and subtractive absorption models with edge-safe placement.
- Kept PSF convolution, primary-beam attenuation, and units explicit rather
  than silently modifying a supplied image-plane model.
- Added absorption-cube conversion for continuum residual depth, fractional
  depth, and optical depth.
- Added an absorption extraction wrapper using the existing 3-D connected-
  component finder in positive absorption-depth space.
- Visibility injection now imports CASA lazily, so image/cube-only workflows
  can use DICE in a standard Python environment.

## 2026-07-07

Restarted the old `GenerateMockData` repository as **BliSS-ALMA**.

Decided:
- project name: BliSS-ALMA
- repository name: bliss-alma
- default branch: main
- scope: visibility-, image-, and cube-domain injection--recovery measurements that connect source extraction to completeness, contamination, flux recovery, astrometric accuracy, and corrected population constraints; current applications focus on ALMA data and include untargeted line searches and mock-data generation
- target paper: Astronomy and Computing software/method paper

On 2026-07-20 the project was renamed **DICE: Discovery and Image/Cube
Evaluation** to avoid name collisions and to describe the broader image/cube
scope. The historical `bliss_alma` Python namespace is temporarily retained
for compatibility. The intended GitHub repository name is `DICE`.

Next:
- extract mock-data generation routines from `old/GenerateMockdata.py`
- extract the visibility-domain source-injection core into `src/bliss_alma/injection.py`
- keep the old scripts for historical/reference purposes

## 2026-07-12

Added an initial ALMA2040/GHOST test plan:

- `docs/alma2040_test_plan.md`

This document was created in response to Tom's suggestion that the first useful quantities would be:

- continuum point-source extraction
- flux recovery and astrometric accuracy
- completeness

The plan defines the first minimal test as:

GHOST / ALMA2040 output  
→ artificial continuum point-source injection  
→ point-source extraction  
→ completeness, flux-recovery, and positional-offset estimates

This is intended as an internal working plan for BliSS-ALMA and the ALMA2040 collaboration, rather than a final analysis document.

## 2026-07-13

Updated the repository documentation to make the ALMA2040/GHOST test plan easier to find from the main README.

Current near-term priority:

1. keep the source-injection module minimal and usable;
2. prepare helper functions for matching, completeness, flux recovery, and positional offsets;
3. test the workflow with GHOST outputs once the network/setup allows.
## 2026-07-16

Completed the first CASA 6.7 noisy-continuum end-to-end test.

Workflow:

1. create a one-channel zero-valued sky model;
2. generate a thermal-noise reference MS with `simobserve`;
3. copy the reference MS and inject a 10 mJy continuum point source using
   `ft` and `uvsub(reverse=True)`;
4. image the reference and injected datasets with `tclean` in MFS mode;
5. measure the reference RMS, recovered peak, recovered-to-injected flux ratio, and S/N.

First test result:

- reference RMS: approximately 0.112 mJy/beam;
- recovered peak: approximately 9.940 mJy/beam;
- recovered-to-injected flux ratio: approximately 0.994;
- recovered S/N: approximately 88.4;
- recovered position: the phase-centre pixel.

Next:

- replace the simulated reference MS with a GHOST-generated ALMA2040 MS;
- test a grid of injected continuum flux densities;
- implement automated source extraction and completeness calculation.

## 2026-07-17

Added the first PyCupid/ClumpFind source-extraction interface.

Implemented:

- a validated adapter for the PyCupid `clumpfind()` API;
- conversion of integer clump labels into a source catalogue;
- unit tests using a controlled backend;
- export of CASA reference/injected images and RMS as NumPy files;
- an external PyCupid continuum extraction example.

The updated CASA 6.7.5.18 end-to-end test passed again and produced a
contiguous 128 x 128 injected image with the source at pixel `(64, 64)`.

Compatibility finding:

- official PyCupid wheels support only Python 2.7 and 3.4-3.6;
- CASA 6.7.5.18 uses Python 3.12 on Apple Silicon;
- the published source bundles x86_64 Starlink libraries;
- therefore the PyCupid backend needs a separate compatible runtime or a
  maintained rebuild before it can be executed on this computer.

Next:

- run the external example in a compatible PyCupid environment;
- validate the ClumpFind catalogue against the known injected source;
- apply the same extraction to blank positive and negative images for the
  first contamination test.

### Source-finder decision

PyCupid was removed as a runtime dependency because its official binaries are
limited to old Python/x86_64 environments. BliSS-ALMA now uses a configurable
SciPy seed-and-flood connected-component finder instead. RMS method, masks,
relative or absolute thresholds, minimum island size, connectivity, and
positive/negative polarity are all selectable.

First SciPy extraction result using a 5-RMS seed, 3-RMS flood threshold,
five-pixel minimum, and border rejection:

- the injected phase-centre source was recovered at pixel `(64, 64)`;
- four positive islands were found in the injected dirty image;
- three of those additional islands are associated with dirty-beam sidelobe
  structure and demonstrate why PSF-aware validation is needed;
- zero positive and zero negative islands were found in the blank reference
  image at this threshold.

### Single-pointing primary-beam population test

Added a CASA 6.7.5.18 end-to-end test that places 21 equal-intrinsic-flux
point sources on a uniform grid across one pointing. `simobserve` applies the
ALMA voltage-pattern model, and the full dirty image is extracted in one pass.

Initial 1 mJy result with a 5-RMS seed and 3-RMS flood threshold:

- 21 input sources;
- 6 extracted islands, all matched to input sources;
- zero unmatched positive islands;
- zero positive and zero negative islands in the matched blank simulation;
- completeness 0/12 for `0.2 <= PB < 0.4`;
- completeness 1/4 for `0.4 <= PB < 0.6`;
- completeness 4/4 for `0.6 <= PB < 0.8`;
- completeness 1/1 for `0.8 <= PB <= 1.0`.

This establishes the efficient population workflow needed for later mapping
tests. Next, repeat over a flux grid and calculate the PB-corrected recovered-to-injected flux ratio.

### Ten-realization smoke test

Generalized the single-pointing population experiment to ten independent
thermal-noise seeds. All ten CASA simulations completed and produced
aggregate source-recovery and detection catalogues.

- 210 sources injected;
- 68 sources recovered (`C = 0.324` overall);
- 70 populated-image detections;
- two unmatched populated-image detections (`Q_truth = 0.029` overall);
- zero blank-positive and zero blank-negative detections at the current
  threshold (`Q_blank = 0` in this limited sample);
- mean blank RMS `0.1122 mJy/beam`.

The mismatch between the two preliminary contamination estimates confirms
that both must be retained. Next: implement differential and cumulative
metrics, then add a flux grid and expand to 100 realizations.

Added the first S/N-dependent diagnostic plot with PB-coloured differential
and cumulative panels. This exposes an important limitation of the one-flux
smoke test: PB and expected S/N are currently correlated. The flux-grid test
will deliberately overlap the PB bins in S/N before physical differences
between PB locations are interpreted.

## 2026-07-20 GHOST F1 end-to-end validation

Completed the first ALMA2040/GHOST F1 test using job
`ALMA-GHOST-CCK6FJ4E` with 255 12-m antennas, declination
`-27d48m00.00`, 5 minutes on source, PWV 1.0 mm, 230 GHz centre frequency,
7.5 GHz bandwidth, natural weighting, and a zero-valued 2D FITS sky model.

Downloaded MS inspection:

- 971,550 rows, one 7.5-GHz channel, and two correlations;
- 10-s integrations with zero flagged data;
- baseline range 15.0--388.0 m (11.50--297.68 klambda);
- GHOST phase centre `12h00m00s -27d48m00s`;
- dirty-image RMS `6.7173 microJy/beam`;
- beam `0.7361 x 0.7166 arcsec`, PA `-89.50 deg`.

Re-imaging the downloaded MS with the logged GHOST parameters reproduced the
distributed dirty FITS exactly: RMS ratio 1.0, identical beam, and zero pixel
difference over all 72,701 common finite pixels.

The visibility-injection smoke test copied the blank MS and injected five
100-microJy intrinsic point sources with `componentlist`, `ft`, and
`uvsub(reverse=True)`. Because plain `ft` does not apply the primary beam,
each component was explicitly scaled to its apparent flux using the sampled
PB value before injection. Results:

- 5/5 injected sources detected and matched at zero-pixel offset;
- zero unmatched populated-image detections;
- zero blank-positive and zero blank-negative detections;
- recovered S/N 12.3--14.5;
- noise-cancelled recovered-to-injected response ratio 0.988--1.005.

F1 therefore passes the complete GHOST-MS-to-BliSS workflow. F2--F6 may now
be generated with the identical setup, changing only the configuration file.
The final large Monte Carlo experiment remains deferred.

### F1 blank-noise threshold and S/N--PB grid

The F1 dirty-image S/N distribution is close to standard normal over the
72,701 finite pixels inside `PB >= 0.2` (mean `3.4e-5`, standard deviation
`1.000`, skewness `0.012`, excess kurtosis `-0.006`). It contains no positive
pixel above 4 sigma and two negative pixels below -4 sigma. With the production
island settings there are zero positive and one negative blank island at a
4-sigma seed / 3-sigma flood threshold. This threshold is adopted for the F1
continuum test.

The first controlled grid used 10 deterministic batches, each containing 44
simultaneous sources: 11 apparent S/N values from 2.0 to 7.0 in steps of 0.5,
crossed with four PB bins. Thus every S/N--PB cell contains exactly 10 sources
(440 total). Apparent injected flux was fixed by `S/N * blank RMS`; intrinsic
flux was recorded as `apparent / PB` because CASA `ft` bypasses primary-beam
attenuation. Sources were separated by at least 20 pixels.

- 273/440 sources were recovered;
- all 273 populated-image detections matched an injected source;
- there were no unmatched positive detections;
- completeness reaches 0.9--1.0 in every PB bin by input S/N 5.0--5.5;
- both differential and cumulative contamination are zero in this limited F1
  test;
- the blank catalogue repeats the same GHOST noise realization in all ten
  batches, so this is a pipeline/grid validation, not the final uncertainty
  estimate.

Do not interpret small PB-to-PB fluctuations from only 10 objects per cell as
physical. Independent noise realizations remain necessary for the final
Monte Carlo, which is still deferred.

For presentation, differential Q is clipped at 1.0 when the finite-sample
negative/positive clump ratio exceeds unity. The stored counts remain the
canonical placement record. Future injection batches use a deterministic outward spiral,
rotated between batches, rather than random-looking positions; every PB bin
still contains 11 sources and the verified minimum separation remains at least
20 pixels.

The contamination diagnostic was subsequently corrected to use the inverted
blank-image (negative-clump) catalogue as `N_mock`; blank-positive clumps are
not the numerator. A separate catalogue floor of 2 sigma / 1.5 sigma flood was
used only to trace the Q curve below the adopted detection threshold. It does
not replace the production 4-sigma seed / 3-sigma flood criterion. The
cumulative negative-clump fraction falls rapidly toward 4 sigma, is non-zero
in the PB bin containing the one >4-sigma negative island, and is zero above
4.5 sigma in this F1 realization.

The diagnostic figure also records input apparent S/N against recovered peak
S/N. Sources detected near the threshold are preferentially those on positive
noise excursions, so recovered S/N is visibly biased high at low input S/N.
This selection effect should be quantified with independent noise realizations
in the final analysis, even though it is secondary to completeness and Q.

### F1 independent-noise production test

The earlier ten injection batches reused one GHOST noise realization and must
be treated only as a pipeline smoke test. The production F1 test instead
regenerated the MS `DATA` and `CORRECTED_DATA` columns for every realization.
The original GHOST visibility data are zero-mean complex Gaussian; `SIGMA` and
`WEIGHT` are constant unity, so the real and imaginary noise standard
deviations were measured directly from the full original DATA column and
preserved per polarization. Seeds `20400000`--`20400009` are recorded for exact
reproduction.

Ten independent visibility-noise realizations completed successfully:

- 44 spiral-position sources per realization and 440 injections total;
- 257 sources recovered with the adopted 4-sigma seed / 3-sigma flood finder;
- blank RMS range 6.488--6.759 microJy/beam, mean 6.620 microJy/beam;
- negative diagnostic clumps per realization: 155--182 at the 2-sigma
  catalogue floor;
- recovered sources per realization: 24--27;
- differential Q is now based on independent blank-negative and populated
  positive catalogues in every realization, split by PB and clipped to one for
  presentation only.

Only CSV/JSON tables, figures, and the seeds are retained. The regenerated MS
and CASA image products were removed after the run.
