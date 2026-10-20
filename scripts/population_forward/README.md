# Forward modelling of source populations for facility comparison

This project covers the P0 pipeline identity test, P1 array-geometry
isolation, and P2 continuum/[CII] bandwidth comparisons.  Synthetic
catalogue definitions and facility assumptions are stored under
`configs/population_forward/`; final figures receive provenance sidecars.

The public entry points are:

- `generate_catalogues.py` and `validate_p0_identity.py`;
- `build_wideband_ms.py` and `image_wideband_ms.py`;
- `validate_p0_visibility.py` and `run_p1_geometry.py`;
- `plot_p1_geometry.py`, `plot_p1_lf.py`, `plot_p2_continuum.py`, and
  `plot_p2_cii.py`;
- `validate_p2_cii_cube.py`.

The CASA-independent products from the migrated catalogue, P0 identity,
P1 plotting/LF, and P2 continuum/[CII] scripts are byte-identical to the
legacy CSV and JSON references. The CASA visibility entry points remain to be
rerun under the migrated paths before the numbered `examples/67`--`77` are
retired.
