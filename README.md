# DICE

[![Tests](https://github.com/natsuki-h2/DICE/actions/workflows/tests.yml/badge.svg)](https://github.com/natsuki-h2/DICE/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23487915.svg)](https://doi.org/10.5281/zenodo.23487915)

**DICE** stands for **Discovery and Image/Cube Evaluation**.

DICE is an injection--recovery framework for quantifying how a specified data
set and analysis procedure recover and contaminate an underlying source
population. Rather than centring on the optimization of a single source
finder, it combines visibility-, image-, or cube-domain source injection with
source extraction and matching to measure completeness, contamination, flux
recovery, and astrometric accuracy. These measurements can then be propagated
into corrected source counts and luminosity-function constraints. Current
applications focus on ALMA, while the image- and cube-domain workflows are not
limited to interferometric data.

## Start here: evaluate your own ALMA data

If this is your first DICE test, begin with the student quick start in your
preferred language:

- **[English: Quick Start for Student Testers](docs/student_quick_start_en.md)**
- **[日本語：学生テスター向けクイックスタート](docs/student_quick_start_ja.md)**

It begins with Git installation, first-time cloning, and the safe `git pull`
procedure for later sessions. The CASA-independent smoke test then verifies
installation, writes a small JSON result, and records provenance before any
measurement set is modified.

For the practical measurement-set-to-catalogue workflow, start with:

- **[ALMA line-search and detectability guide](docs/alma_line_search_detectability_guide.md)**

The guide proceeds from continuum imaging and pixel-distribution quality
checks through masking decisions, continuum subtraction, positive and negative
line searches, candidate validation, source injection, completeness, and
contamination measurements. It also covers emission and absorption searches
and lists the settings and intermediate products needed for reproducibility.

Historical workflow validation is documented separately:

- **[Hayatsu et al. (2017) source-finding reproduction](docs/hayatsu2017_source_finding_reproduction.md)**
- **[Hayatsu et al. (2025) detectability reproduction](docs/hayatsu2025_reproduction.md)**

Documentation, figure labels, and new API names should follow the
**[DICE terminology guide](docs/terminology.md)**.

## Questions and bug reports

For usage questions, reproducible errors, and documentation problems, open a
[GitHub issue](https://github.com/natsuki-h2/DICE/issues). Include the command,
complete error message, OS/Python/CASA versions, configuration JSON, and the
generated provenance JSON when available. Do not upload measurement sets,
unpublished coordinates, candidate catalogues, credentials, or other private
research data. Contact the person who provided access through your approved
institutional channel when private data must be discussed.

The project grew out of an earlier repository, **GenerateMockData**, which was originally developed in 2018 for ALMA mock-data generation and artificial-source injection tests.

The original scripts are kept in the `old/` directory only for historical
reference. They contain legacy, destructive cleanup commands and must not be
run on current data; see [`old/README.md`](old/README.md). The existing Python
import namespace remains `bliss_alma` during the DICE rename to avoid breaking
current analysis scripts; a compatibility-preserving module migration will
follow separately.

## Main features

- Completeness and contamination measurements
- Flux-recovery and astrometric-accuracy measurements
- Visibility-domain artificial source injection
- CASA-independent image/cube-plane source injection
- Forward modelling of source populations
- Mock ALMA data generation
- Untargeted line searches in ALMA data cubes
- Continuum-normalized absorption-line searches
- Positive and negative source detection
- Source-recovery tests for simulated ALMA configurations

## Current status

DICE 0.1.0 is the first public alpha release. It provides visibility-domain
point-source injection, image- and cube-domain injection, configurable source
extraction, and completeness/contamination measurements. The package remains
under active development; workflows and configuration fields may change before
the first stable release.

For reproducible citation of this exact release, use
[DICE v0.1.0 on Zenodo](https://doi.org/10.5281/zenodo.23487916). The DOI badge
at the top of this page resolves to the archive containing all released
versions.

The first implemented module is `src/bliss_alma/injection.py`, which provides a minimal CASA-based visibility-domain point-source injection workflow.

## Planned extensions

The current release supports point-like emission-line injection and
completeness measurements, as well as positive and negative residual
extraction for emission- and absorption-line candidate searches. Planned
post-release extensions include resolved continuum and spectral-line sources,
full polarization, and time-domain/transient searches. The transient workflow
is expected to support multi-epoch or
time-binned visibility/image data, injected light curves, difference-image
source finding, and completeness/contamination as functions of flux, duration,
cadence, and position. These capabilities are roadmap items and are not yet
part of the validated public API. Advanced users may already apply the
axis-agnostic 3-D connected-component finder to a `(time, y, x)` cube, but
this exploratory use does not supply time-domain noise, cadence, PSF, or
calibration handling; see
[`docs/source_extraction.md`](docs/source_extraction.md#axis-agnostic-3-d-use).

## Implemented modules

### `src/bliss_alma/injection.py`

This module provides CASA-based visibility-domain point-source injection and
CASA-independent 2-D image/3-D cube injection. Image-plane models are supplied
in explicit compatible units and are not silently PSF-convolved or PB-scaled.

The first implemented function is designed for injecting a continuum point source into an ALMA measurement set using CASA tasks such as:

- `cl.addcomponent`
- `ft`
- `uvsub(reverse=True)`

The public tutorial entry point is:

- [`scripts/tutorials/inject_point_source.py`](scripts/tutorials/inject_point_source.py)

### `src/bliss_alma/recovery.py`

This module provides simple source-recovery metrics, including:

- flux recovery ratio, defined as recovered flux / injected flux;
- positional offset in arcsec;
- a simple recovered/not-recovered matching criterion;
- completeness, defined as recovered sources / injected sources.

Use `flux_recovery_ratio` in new analyses. The older function name
`flux_fidelity` remains available only so that existing scripts continue to
run; catalogue fidelity/reliability is a different quantity.

The public tutorial entry point is:

- [`scripts/tutorials/basic_metrics.py`](scripts/tutorials/basic_metrics.py)

### `src/bliss_alma/extraction.py`

This module provides configurable SciPy-based source extraction for 2D images
and 3D cubes. It uses separate seed and flood thresholds and supports
identical positive and negative extraction for contamination tests.

- [`scripts/tutorials/extract_continuum_sources.py`](scripts/tutorials/extract_continuum_sources.py)
- [`docs/source_extraction.md`](docs/source_extraction.md)

For a measurement-set-to-catalogue workflow, including visibility-domain
injection, completeness, contamination, absorption searches, validation, and
the required archived outputs, see:

- [`docs/alma_line_search_detectability_guide.md`](docs/alma_line_search_detectability_guide.md)

### Single-pointing population test

The CASA end-to-end population example places equal-intrinsic-flux point
sources uniformly across one primary beam, simulates one noisy pointing,
extracts the full field in one pass, and measures completeness as a function
of primary-beam response.

- [`scripts/tutorials/simulate_single_pointing_population.py`](scripts/tutorials/simulate_single_pointing_population.py)
- [`scripts/tutorials/run_single_pointing_monte_carlo.py`](scripts/tutorials/run_single_pointing_monte_carlo.py)
- [`scripts/tutorials/plot_single_pointing_metrics.py`](scripts/tutorials/plot_single_pointing_metrics.py)
- [`docs/single_pointing_population.md`](docs/single_pointing_population.md)

## Related work

The original implementation and its application are described in:

Hayatsu et al. (2025), *Publications of the Astronomical Society of Japan*,
77, 1163--1180, [doi:10.1093/pasj/psae094](https://doi.org/10.1093/pasj/psae094).

## Citation and license

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Cite the
versioned Zenodo record for the software release used in an analysis; the DICE
paper citation will be added after publication.

DICE is distributed under the [BSD 3-Clause License](LICENSE). Contributions
and bug reports are welcome; see [`CONTRIBUTING.md`](CONTRIBUTING.md) and
[`SECURITY.md`](SECURITY.md).

## Notes

Local data products such as GHOST outputs, FITS files, and CASA measurement sets are not stored in this repository; see [`docs/local_data_policy.md`](docs/local_data_policy.md).
