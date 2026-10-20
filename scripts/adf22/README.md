# ADF22 / Hayatsu reproduction

This directory contains the public archive-to-result workflows for the
Hayatsu et al. 2017 source search and Hayatsu et al. 2025 detectability
analysis.  Measurement sets are not distributed; archive identifiers,
preparation commands, configuration, and expected numerical checks are.

## Recommended sequence

1. `inspect_measurement_sets.py`: inspect archive measurement-set structure.
2. `run_uvcontsub.py`: perform visibility-domain continuum subtraction.
3. `image_contsub_cubes.py`: make native dirty spectral cubes.
4. `search_contsub_cubes.py SPW`: search one spectral window (SPW 0--3).
5. `plot_exact_recovery.py` and `summarize_all_spw.py`: reproduce and summarize
   the Hayatsu et al. (2017) search.
6. `reproduce_hayatsu2025.py`: reproduce the image-domain injection branch.
7. `generate_hayatsu2025_blank_realizations.py` and
   `analyze_hayatsu2025_blank_realizations.py`: validate blank-field false
   positives using whole noise realizations as the uncertainty unit.

Measurement sets are intentionally not distributed. See
`configs/adf22/workflow_manifest.json` for inputs, stages, and expected checks.
Final figures from the migrated plotting scripts receive a `.provenance.json`
sidecar containing the command, inputs, outputs, parameters, and Git commit.

The numbered `examples/40*`--`51*` and `examples/64`--`66` are retained as
legacy reference implementations. The migrated numerical products have been
cross-checked against them.
