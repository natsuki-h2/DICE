# Contributing to DICE

Bug reports, documentation corrections, and focused code contributions are
welcome in English or Japanese.

Before opening an issue, check the README, the student quick start, and existing
issues. Never upload proprietary measurement sets, unpublished source
coordinates, private candidate catalogues, credentials, or personal data.

For a code contribution:

1. Create a branch from the current `main` branch.
2. Keep the change focused and document any new configuration fields.
3. Add or update tests when behaviour changes.
4. Run `PYTHONPATH=src python -m pytest -q`.
5. Open a pull request describing the scientific purpose, implementation, and
   validation performed.

CASA-dependent changes should state the CASA version used. Generated data,
measurement sets, FITS products, and local paths must not be committed.
