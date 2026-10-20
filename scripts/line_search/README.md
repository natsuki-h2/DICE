# Generic line-search workflows

This directory is target-neutral. It provides reusable cube QA, signed
positive/negative searches, width-dependent selection metrics, and candidate
plots. Target-specific coordinates, candidate identifiers, spectra,
frequencies, and catalogues must remain outside the public repository.

## First public entry point

Copy `configs/line_search/example.json`, set the paths to a continuum-
subtracted FITS cube and matching primary-beam FITS image, then run:

```bash
python scripts/line_search/search_fits_cube.py my_config.json
```

The command searches the declared channel widths with identical positive and
negative thresholds, merges overlapping detections, writes a signed candidate
catalogue, and records provenance. Inspect continuum quality and define masks
before treating the catalogue scientifically; see
`docs/alma_line_search_detectability_guide.md`.

To create a multi-page moment-0 and spectrum review PDF after the search, pass
the same configuration to:

```bash
python scripts/line_search/plot_candidates.py my_config.json
```

This writes a candidate-review PDF, a sky-distribution PDF/PNG, a three-panel
pixel/S/N/channel diagnostic PDF/PNG (including per-channel RMS), and their provenance JSON in the
configured output directory. These are visual checks, not reliability
estimates.

Use `--radius-pixels N` to change the moment-0 cutout radius and
`--max-candidates N` to plot only the N highest-|S/N| catalogue entries.

A successful search ends with a message such as
`Wrote .../signed_line_candidates.csv (3 positive, 2 negative)` and writes both
the catalogue and its provenance JSON. Zero positive and zero negative
candidates can still be a successful program run, but the cube and PB inputs
must then be checked. A successful plotting run ends with
`Wrote .../candidate_review.pdf (5 candidates)`.

`pb_fits` must be the dimensionless CASA `*.pb` response exported to FITS,
normally near one at the field centre and decreasing outwards. It must not be
an `impbcor`/`*.pbcor` science image. The width bank is applied internally as
multi-channel boxcar sums, so a separate pre-binning step is not normally
required. A header-only catalogue does not by itself prove that the field has
no lines; validate the PB mask and continuum-subtracted cube first.
