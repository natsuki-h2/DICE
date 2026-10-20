# ALMA2040 / GHOST Test Plan

## Background

This document summarizes the first ALMA2040 test plan for DICE.

The near-term goal is to connect DICE with ALMA2040/GHOST outputs and test simple continuum point-source recovery.

Tom suggested that the first useful quantities would be:

- continuum point-source extraction
- flux recovery and astrometric accuracy
- completeness

GHOST can provide blank-field or noisy FITS images, default detected objects, and outputs that can be downloaded directly from the web interface.

## First-step goal

The first test should be intentionally simple.

The goal is to start with a blank or noisy ALMA2040/GHOST output, inject continuum point sources with different flux densities, and check whether they can be recovered.

## Planned workflow

1. Select one representative ALMA2040 configuration in GHOST.
2. Generate or download a blank/noisy field.
3. Inject artificial continuum point sources.
4. Image or analyze the resulting data.
5. Extract point sources.
6. Compare input and recovered source properties.
7. Estimate completeness, flux recovery, and astrometric accuracy.

## Quantities to measure

### Completeness

Completeness will be estimated as the fraction of injected sources recovered above a chosen detection threshold.

Initially, this can be measured as a function of input flux density.

### Flux recovery and astrometric accuracy

For the first test, report the recovered-to-injected flux ratio and positional
offset explicitly:

- recovered flux / input flux
- positional offset between injected and recovered sources
- number of spurious detections in blank or negative maps

The precise definition can be refined later after the first GHOST test.

## Initial source model

The first source model will be a continuum point source.

More realistic models can be added later, such as:

- Gaussian extended sources
- multiple point sources
- clumpy high-redshift galaxy models
- line emitters

## Relation to DICE

This test connects directly to the first implemented DICE module:

- `src/bliss_alma/injection.py`

The current source-injection function provides a minimal CASA-based visibility-domain point-source injection workflow.

## Notes

The first priority is not to build the full line-search pipeline, but to create a minimal end-to-end test connecting:

GHOST / ALMA2040 output  
→ artificial source injection  
→ point-source extraction  
→ completeness, flux-recovery, and positional-offset estimates

Once the first setup is working, the next steps can be divided among collaborators if useful.

## 2026-07-14: GHOST access status

I attempted to access GHOST for the first ALMA2040/DICE test, but the server was temporarily unavailable with a "Service Unavailable" message.

Tom noted that the issue is likely related to the recent warm weather in Manchester and that the service should be available again later this week.

While waiting for the service to become available again, the DICE side will be prepared for the first continuum point-source recovery test.

The first test will remain intentionally simple:

1. select one representative ALMA2040/GHOST setup;
2. download a blank or noisy field;
3. inject continuum point sources;
4. measure completeness and the recovered-to-injected flux ratio.

Local data products should be stored outside Git tracking, for example under:

```text
data/ghost_first_test/
```
