# Extension smoke tests

Four fast, idealized image/cube-domain tests were completed before building
new visibility templates:

- Gaussian extended continuum sources;
- linear polarization, polarized-intensity bias, and fractional leakage;
- point-like spectral lines;
- spatially extended spectral lines.

All curves use a four-sigma threshold and 20,000 Monte Carlo trials per grid
point. The continuum calculation includes Gaussian beam dilution. The
polarization calculation forms `P = sqrt(Q^2 + U^2)` and applies the simple
`sqrt(max(P^2 - sigma^2, 0))` debiasing estimator. The line tests compare a
single-channel/point-source detector with correctly matched spectral or 3-D
filters.

These are implementation and selection-bias smoke tests, not ALMA2040
forecasts. The audited F1 GHOST MS contains only XX and YY parallel hands and
one spectral channel. It therefore cannot provide a visibility-domain Q/U
leakage test or a spectral cube. Those stages require a full-polarization,
multi-channel MS. The audit is saved as
`data/ghost_matrix/extension_smoke_tests/ms_axis_audit.json`.

Primary outputs are under
`data/ghost_matrix/extension_smoke_tests/`.
