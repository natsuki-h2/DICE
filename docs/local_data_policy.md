# Local data policy

DICE should not store large ALMA, GHOST, CASA, FITS, or measurement-set products directly in the GitHub repository.

Local test data should be stored under a local directory such as:

```text
data/
```

This directory is ignored by Git.

## Examples of local data products

- GHOST blank-field outputs
- noisy FITS images
- CASA measurement sets
- CASA image products
- PSF / PB / residual images
- source extraction catalogs generated during tests

## Reproducibility

The repository should instead store:

- setup descriptions
- configuration names
- command examples
- scripts used to generate or analyze the data
- small tables or summary results when appropriate
