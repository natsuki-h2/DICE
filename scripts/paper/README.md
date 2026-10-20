# Paper-specific post-processing

These scripts reproduce post-processing used in the DICE paper. They are
published for transparency, but they are **not part of the stable DICE API**.
Their command-line options and output schemas may change between paper
revisions.

The core DICE package performs injection, extraction, and source matching.
The scripts here add three paper-level operations:

1. fitting increasing completeness and decreasing contamination curves;
2. reconstructing the Hayatsu et al. Category I--III interpretation from the
   fitted transition locations;
3. calculating 16--84 per cent ranges after grouping measurements by complete
   noise realization.

## Fit and classify metric curves

The input is a tidy CSV containing one row per sampled `x` value and group.
It must contain `x`, `completeness`, and `contamination` columns. Additional
columns supplied with `--group-columns` define independent fits.

```bash
python scripts/paper/fit_and_classify_metrics.py \
  metrics.csv fit_summary.csv classification.csv \
  --group-columns linewidth
```

The classification reproduces the historical transition-based interpretation:
Category I lies above both fitted midpoints, Category II lies between the two
midpoints, and Category III lies below both. This is a paper-level catalogue
summary, not a required DICE detection step.

## Uncertainty from complete realizations

For an input catalogue with `realization`, grouping columns, and either a
precomputed `value` or binary/count `numerator` and `denominator` columns:

```bash
python scripts/paper/summarize_realization_uncertainty.py \
  measurements.csv uncertainty.csv \
  --group-columns configuration,input_snr \
  --numerator recovered --denominator injected
```

The script first calculates one statistic per complete realization and then
reports the median and 16th/84th percentiles across realizations. It never
treats sources from the same image as independent noise realizations.

## Hayatsu et al. (2025)

The version-locked ADF22 reproduction remains in
`scripts/adf22/reproduce_hayatsu2025.py`. It contains the published tabulated
values and the exact figure-specific reconstruction. The generic scripts in
this directory make the fit, classification, and uncertainty conventions
available without presenting them as core package functionality.
