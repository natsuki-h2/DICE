# Hayatsu et al. (2017) ADF22 source-finding reproduction

## Scope

This is an image-domain reproduction using the archived four-SPW dirty FITS
cubes. The search is performed over the `PB >= 0.8` region with spectral
boxcars of 6, 9, 12, 15, 18, and 21 channels. Candidate extraction is blind:
the published coordinates are used only after extraction for catalogue
matching.

The archived FITS cubes retain continuum emission, while the published search
used `uvcontsub` with `fitorder=0`. A constant continuum is therefore estimated
and subtracted independently at every image pixel before spectral filtering.
This is a close image-domain approximation, but not identical to rerunning
`uvcontsub` on the full measurement set.

## Recovery

All four published positive candidates are recovered within 0.05--0.15 arcsec:

| Candidate | Published S/N | Recovered S/N |
|---|---:|---:|
| ADF22-LineA | 6.5 | 6.58 |
| ADF22-LineB | 6.2 | 5.68 |
| ADF22-LineC | 7.8 | 6.97 |
| ADF22-LineD | 6.5 | 6.35 |

The strongest recovered inverted clump in SPW 3 is -6.22 sigma, consistent
with the published -6.3 sigma clump. Another -6.06 sigma clump appears in SPW
2 in this image-domain implementation; its status should be checked using the
original visibility-domain continuum subtraction before claiming exact
catalogue-level reproduction.

## Visibility-domain reproduction

The original 58.1 GiB, 103-pointing `Core_region.ms` was also recovered and
processed with the historical CASA continuum-subtraction semantics:
`uvcontsub_old`, `fitorder=0`, `solint='int'`, and the target SPW excluded
from a fit combining the other SPWs. Native natural-weighted mosaic dirty
cubes were then made on the original 1600 by 2400 pixel, 0.1 arcsec grid.
The search used the same PB >= 0.8 region and 6, 9, 12, 15, 18, and 21-channel
boxcars. Published positions and frequencies were used only after blind
candidate extraction.

All four candidates are recovered:

| Candidate | Published S/N | Visibility-domain reproduction |
|---|---:|---:|
| ADF22-LineA | 6.5 | 6.37 |
| ADF22-LineB | 6.2 | 5.99 |
| ADF22-LineC | 7.8 | 8.13 |
| ADF22-LineD | 6.5 | 6.67 |

The close agreement confirms catalogue-level reproduction of the four
reported positive candidates. The earlier image-plane continuum approximation
remains useful as a small image-domain test, but the visibility-domain values are
the reference results.

The identical blind positive/negative search was subsequently run over all
four SPWs. The published SPW 3 inverted clump is recovered at -6.41 sigma,
consistent with the reported approximately -6.3 sigma value. A second
-6.06-sigma clump occurs in SPW 2.

One additional +6.22-sigma peak occurs in SPW 3 at RA=334.385028 deg,
Dec=0.295500 deg. Positive broad-width peaks recur at the same sky position
in every SPW (5.50--5.62 sigma outside SPW 3). It is therefore classified as
a multi-SPW positive residual candidate rather than a new isolated line.
All other positive candidates not matched to the four published sources are
below 6 sigma; the strongest is 5.86 sigma.

## Products

- `examples/41_adf22_blind_line_search.py`
- `examples/42_plot_adf22_source_recovery.py`
- `data/hayatsu_reproduction/results/adf22_blind_candidates.csv`
- `data/hayatsu_reproduction/results/adf22_published_candidate_recovery.csv`
- `data/hayatsu_reproduction/results/adf22_four_candidate_recovery.png`
- `examples/44_uvcontsub_adf22_spw.py`
- `examples/45_image_adf22_contsub_spw.py`
- `examples/46_search_adf22_exact_contsub_spw.py`
- `examples/50_plot_adf22_exact_recovery.py`
- `data/hayatsu_reproduction/results/adf22_exact_published_candidate_recovery.csv`
- `data/hayatsu_reproduction/results/adf22_exact_four_candidate_recovery.png`
- `examples/51_summarize_adf22_exact_all_spw.py`
- `data/hayatsu_reproduction/results/adf22_exact_all_spw_candidates.csv`
- `data/hayatsu_reproduction/results/adf22_exact_all_spw_candidate_summary.png`
