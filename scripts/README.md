# DICE analysis scripts

The scripts in this directory are the public, project-oriented entry points
for running the analyses described in the DICE documentation. Reusable algorithms belong in
`src/bliss_alma`; these scripts should mainly parse a configuration, call the
library, and record provenance.

Projects are organized as follows:

- `tutorials/`: small data-independent examples;
- `ghost/`: GHOST and ALMA2040 continuum experiments;
- `adf22/`: Hayatsu et al. 2017 and 2025 reproduction workflows;
- `population_forward/`: forward modelling of source populations for ALMA-generation comparisons;
- `line_search/`: generic line-search tools without unpublished target data.

Student testers should begin with `docs/student_quick_start_ja.md`. The first
generic line-search command is `scripts/line_search/search_fits_cube.py`.

The numbered files in `examples/` remain available during migration.  A file
is removed from that legacy location only after its numerical products have
been checked against the corresponding project script.

Data are not distributed with the repository.  Each project script reads a
configuration under `configs/` and writes products below its project
directory in `data/`.  Publishable figures are stored in `figures/` together
with a `.provenance.json` sidecar.
