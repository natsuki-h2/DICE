# SciPy source extraction

DICE uses a configurable seed-and-flood connected-component finder
implemented with `scipy.ndimage`.

## Axis-agnostic 3-D use

The connected-component finder operates on array topology and does not inspect
the physical meaning of the third axis. A user may therefore supply a
`(time, y, x)` cube instead of a `(frequency, y, x)` cube for exploratory
transient or variability searches. Positive and negative extraction, masks,
seed/flood thresholds, and connected-island measurements will run unchanged.

This is currently an unsupported, unvalidated use rather than a dedicated
time-domain API. DICE does not yet account for cadence, irregular time
sampling, exposure-dependent noise, time-axis matched filtering, or changes
in uv coverage, PSF, flags, and calibration between epochs. Connectivity and
`min_pixels` are also measured in array pixels, not physical time. Users
adopting this approach must construct a statistically comparable time cube,
select appropriate scaling and connectivity, and validate completeness and
false-positive rates for their own data.

## Detection procedure

1. Estimate the RMS, or accept an externally measured RMS.
2. Identify seed pixels above the seed threshold.
3. Identify connected flood islands above the lower flood threshold.
4. Retain islands containing a seed and at least the requested number of
   pixels.
5. Measure the peak, peak pixel, intensity-weighted centroid, pixel count, and
   summed pixel value for each island.

The same procedure can be applied to `-image` using `polarity="negative"`.
This gives symmetric positive and negative catalogues for contamination tests.

## Configurable parameters

- `rms`: externally supplied RMS;
- `rms_method`: `mad`, `std`, or `rms` when estimating internally;
- `rms_mask`: pixels excluded from the RMS estimate;
- `seed_sigma` and `flood_sigma`: thresholds in units of RMS;
- `seed_threshold` and `flood_threshold`: absolute threshold overrides;
- `min_pixels`: smallest retained island;
- `connectivity`: 2D/3D neighbour connectivity;
- `polarity`: `positive` or `negative`.
- `detection_mask`: pixels excluded from source detection;
- `exclude_border`: reject islands touching an array edge.

The initial continuum demonstration uses a 5-RMS seed, a 3-RMS flood
threshold, and a five-pixel minimum. These values are provisional and should
be varied during validation.

## Running

First generate the noisy continuum products with CASA:

```bash
casa --nogui --nologger --log2term \
  -c examples/04_casa_continuum_end_to_end.py
```

Then run extraction with the same CASA Python or another current Python with
NumPy and SciPy:

```bash
PYTHONPATH=src python3 examples/05_scipy_continuum_extraction.py
```

The example extracts:

- positive sources from the injected image;
- positive sources from the blank reference image;
- negative sources from the blank reference image.

This is the starting point for the blank-field contamination workflow.

For an end-to-end population test that extracts many injected sources in one
pass, see:

```bash
casa --nogui --nologger --log2term \
  -c examples/06_casa_single_pointing_population.py
```

This example samples a range of primary-beam responses in one pointing and
uses the same extraction settings for the populated, blank-positive, and
blank-negative images.

The first dirty-image test also detects several sidelobe islands around the
bright injected source. This is expected and useful: quantifying how such PSF
structure changes the extracted catalogue is a central ALMA2040 test rather
than something to hide inside the finder.
