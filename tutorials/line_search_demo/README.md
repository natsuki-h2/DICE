# Line-search demo

This small deterministic example contains Gaussian noise, three injected
positive lines at different positions, frequencies, and strengths, and no
injected negative source. Noise can therefore produce both positive and
negative marginal candidates. Compare the untargeted-search catalogue with
`injected_truth.csv`; a candidate is not automatically a real source.

From the DICE repository root, run:

```bash
python scripts/line_search/search_fits_cube.py tutorials/line_search_demo/demo_config.json
python scripts/line_search/plot_candidates.py tutorials/line_search_demo/demo_config.json
```

A successful run writes the catalogue, provenance JSON files, and three plot
products under `data/tutorials/line_search_demo/`:

- `candidate_review.pdf`: moment-0 and spectrum page for every candidate
- `candidate_sky_distribution.pdf`: all candidate positions within the PB
- `line_search_diagnostics.pdf`: pixel S/N, candidate S/N, and channel-S/N QA

PNG versions of the last two summaries are also written. The repository
includes one reference run under `expected_output/` so that first-time users
can compare file names and plots. The deterministic reference result contains
6 positive and 3 negative candidates. Compare these with
`expected_candidate_interpretation.csv`: three positive candidates match the
injected lines, while three positive and all three negative candidates are
noise peaks.

![Expected spatial candidate distribution](expected_output/candidate_sky_distribution.png)

Circles show positive candidates and downward triangles show negative
candidates. The PB boundary makes the searched region explicit.

![Expected line-search diagnostics](expected_output/line_search_diagnostics.png)

These panels show the pixel S/N distribution, positive and negative candidate
counts, and candidate S/N versus channel with the channel RMS overlaid.

The input FITS files can be regenerated exactly with:

```bash
python scripts/tutorials/create_line_search_demo.py
```
