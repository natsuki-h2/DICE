# Tutorials

Start with the CASA-independent smoke test:

```bash
python scripts/tutorials/basic_metrics.py
```

Before changing a measurement set, inspect the requested injection:

```bash
python scripts/tutorials/inject_point_source.py --dry-run
```

The Gaussian-line point-source tutorial requires a working-copy MS selected to
one science field and one spectral window:

```bash
python scripts/tutorials/inject_line_source.py --dry-run
```

The complete introductory sequence has now been migrated from
`examples/01`--`09`:

1. `basic_metrics.py` -- CASA-independent installation check;
2. `inject_point_source.py` -- visibility injection with a safe dry run;
3. `inject_line_source.py` -- Gaussian-line visibility injection with a safe dry run;
4. `continuum_end_to_end.py` -- minimal noisy CASA continuum simulation;
5. `extract_continuum_sources.py` -- signed SciPy extraction;
6. `simulate_single_pointing_population.py` -- PB-dependent population test;
7. `run_single_pointing_monte_carlo.py` -- ten whole-noise realizations;
8. `plot_single_pointing_metrics.py` -- completeness/contamination figures;
9. `run_single_pointing_flux_grid.py` -- linear dirty-response flux grid.
10. `plot_line_injection_helix.py` -- visualize and export the standard
    fixed-FWHM, five-source placement used for line-completeness experiments.

The reference placement figure, CSV, and provenance JSON are kept
under `tutorials/line_injection_demo/`. The line-search demo and its expected
diagnostic figures are under `tutorials/line_search_demo/`.

The numerical arrays and shared catalogue fields from the CASA workflows and
all aggregate Monte Carlo, flux-grid, and metric CSVs have been checked
against the legacy products. Set `DICE_TUTORIAL_DATA_ROOT` to run an isolated
copy without overwriting existing products.
