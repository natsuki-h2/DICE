# DICE Quick Start for Student Testers

This guide is the shortest path for students and other first-time DICE users who already know the basics of ALMA data analysis but may be unfamiliar with Git, GitHub, or DICE.

DICE uses artificial-source injection and recovery to measure how a specified
data set and analysis procedure affect an underlying source population. Source
extraction is an essential step, but the broader purpose is to propagate the
measured completeness, contamination, flux recovery, and astrometric accuracy
into corrected source counts and population constraints.

If anything is unclear, first ask the person who introduced you to DICE, such as your instructor. You can also consult the [FAQ for student testers](https://github.com/natsuki-h2/DICE/issues/1) and previous [DICE Issues](https://github.com/natsuki-h2/DICE/issues). Instructions for asking a question and the information to include are given at the end of this page.

Optional Tutorial C provides a complete example using simulated data. Its candidate tables and diagnostic plots are written to `data/tutorials/line_search_demo/`, where you can compare them with the supplied reference output.

## GitHub-related terminology

- **repository**: the GitHub location containing the DICE code and documentation
- **clone**: make the first local copy of the repository on your computer
- **pull**: update an existing local copy from GitHub
- **root**: the top-level DICE directory containing `README.md`
- **modified (`M`)**: a Git-tracked file has local changes
- **untracked (`??`)**: a new local file or directory is not tracked by Git

## Step 0. Get DICE and update it safely

### 0-1. Open a terminal and check Git

```bash
git --version
```

If this prints `git version ...`, Git is ready. On macOS, follow the displayed instructions to install Command Line Tools if Git is missing. On Windows, install Git for Windows from <https://git-scm.com/download/win>.

### 0-2. First time only: clone the repository

Run this only if you do not already have a DICE directory. The example places it under `Documents`. `git clone` creates the directory, so do not run `mkdir DICE` first.

```bash
cd ~/Documents
git clone https://github.com/natsuki-h2/DICE.git
cd DICE
```

The HTTPS command above does not require a GitHub account because DICE is a
public repository. If you already use a GitHub SSH key, you may instead run:

```bash
git clone git@github.com:natsuki-h2/DICE.git
cd DICE
```

If Git reports `fatal: destination path 'DICE' already exists`, do not clone again. Continue with the next subsection.

If the HTTPS command reports `Repository not found`, check the spelling of the
repository URL and your network connection. Send the complete error to the
person supporting you if the problem remains. Do not send the original
measurement set or unpublished data.

### 0-3. Later sessions: pull before starting work

At the beginning of every session, first run:

```bash
cd ~/Documents/DICE
git status --short
```

No output means that Git has found no local changes requiring attention. If file names appear, inspect them before pulling. If you cannot resolve the result using the example below, send the complete output to the person supporting you. Avoid running `git reset --hard` or deleting files when you do not understand why they are listed, because either action may discard your work. To test whether a file can be removed, copy or move it outside the DICE directory first so that it can be restored.

For example, this does not indicate a broken DICE installation:

```text
M configs/line_search/example.json
?? configs/line_search/my_line_search.json
```

`M` means that the distributed `example.json` was changed locally. `??` means that your own JSON is untracked. If all of your analysis settings are saved in `my_line_search.json`, restore the distributed template with:

```bash
git restore configs/line_search/example.json
```

Run `git status --short` again. The following output is normal; keep this JSON and use it for your line search.

```text
?? configs/line_search/my_line_search.json
```

Continue if the output is empty or contains only your own analysis JSON marked with `??`. Ask for help if other entries remain.

```bash
git switch main
git pull --ff-only
```

`Already up to date.` means that you have the current version. Otherwise Git lists the updated files.

### 0-4. Confirm that you are at the repository root

```bash
pwd
ls
```

You are in the correct place if `pwd` ends in `DICE` and `README.md`, `src`, `scripts`, and `configs` are visible. If you cloned elsewhere, replace `~/Documents/DICE` with the actual path. Run all remaining commands from this DICE root.

Finally check Python:

```bash
python --version
```

Python 3.10 or later is required. If `python` is not found, try `python3 --version` and replace `python` with `python3` in later commands. If neither works, do not install several Python distributions experimentally; send your operating system and the full error to the person supporting you.

The four-line routine for later sessions is:

```bash
cd ~/Documents/DICE
git status --short
git switch main
git pull --ff-only
```

## Step 1. Install DICE and run a basic check

From the repository root, create a virtual environment and install DICE in editable mode:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[all]'
python -m pytest -q
python scripts/tutorials/basic_metrics.py
```

The installation check succeeds if the last command prints a JSON result and writes a result and provenance record under `data/tutorials/products/`.

If the computer cannot access the network but the required libraries and `setuptools` already exist in its Python environment, you can reuse the system packages and disable downloads. The online procedure above is normally preferred.

```bash
deactivate 2>/dev/null || true
python -m venv --system-site-packages .venv
source .venv/bin/activate
python -m pip install -e . --no-deps --no-build-isolation
```

If this still reports `Cannot import 'setuptools.build_meta'`, skip installation and run only the basic tests with:

```bash
PYTHONPATH=src python -m pytest -q
```

You can create `.venv` from within a Conda environment, but the double prompt `(.venv) (another-environment)` is easy to misread. In general, run `conda deactivate` before creating `.venv`, or use only one of the two environments.

Once these commands succeed, the DICE installation check is complete. The optional tutorials demonstrate separate features; none is required preprocessing for a line search on your own cube.

## Optional Tutorial A. Inspect a JSON configuration

Open `configs/tutorials/basic_examples.json` in a text editor. JSON is a text record of analysis settings, not a program. Values and paths can be edited, but take care not to remove required braces, brackets, commas, or quotation marks.

This file belongs to the basic tutorial. For a real line search, directly use a separate JSON copied from `configs/line_search/example.json`, as described in Step 2.

## Optional Tutorial B. Inject an artificial point source into visibility data

The DICE functions used here handle one point source at a time. Continuum and Gaussian-line sources are described separately below. Visibility injection of extended sources or polarization is outside this tutorial.

### B-1. Continuum point source

DICE registers the artificial source in a CASA component list, calculates its model visibility with `ft(..., usescratch=True)`, and adds it to the MS with `uvsub(reverse=True)`.
`uvsub(reverse=True)` changes `CORRECTED_DATA`; it leaves `DATA` unchanged. Image the injected MS from `CORRECTED_DATA`. If the uninjected reference MS has no `CORRECTED_DATA` column, compare its `DATA` with the injected MS's `CORRECTED_DATA`.

The script that adds the artificial source to the visibility data modifies the specified MS in place. **Never specify the original MS.** First check that the destination has enough free space and make a complete working copy. Replace the example path below with the path to your MS.

```bash
mkdir -p data/tutorials/working
cp -R /path/to/original.ms data/tutorials/working/my_injected.ms
```

If the MS is too large for a complete local copy, stop and agree on a working location with the person supporting you. Do not overwrite the original MS to save space.

Copy the example JSON included with DICE to a new name. Do not edit `basic_examples.json` itself.

```bash
cp configs/tutorials/basic_examples.json \
  configs/tutorials/my_visibility_injection.json
```

Edit the `visibility_injection` section in the copied JSON:

```json
"visibility_injection": {
  "vis": "data/tutorials/working/my_injected.ms",
  "direction": "J2000 03h32m39.53769 -27d46m43.5673",
  "flux_jy": 0.0001,
  "spw": "",
  "complist": "data/tutorials/working/my_source.cl"
}
```

- `vis`: working-copy MS that will receive the artificial source; never the original MS
- `direction`: sky position of the artificial source
- `flux_jy`: apparent flux density of the artificial source in Jy
- `spw`: SPW/channel selection passed to CASA `ft`; an empty string selects all SPWs
- `complist`: destination for the CASA component list created by DICE

The CASA `ft` operation used in this procedure does not automatically apply primary-beam attenuation. Therefore, `flux_jy` is the apparent flux added to the MS. To place an artificial source with intrinsic flux \(S_{\rm intrinsic}\) at a position with PB response \(p\), set `flux_jy = p * S_intrinsic` yourself and record the PB value used. DICE replaces an existing component list at `complist`, so use a new name inside the working directory.

First confirm that CASA 6 Python tasks are available in the current environment:

```bash
python -c "from casatasks import ft; print('CASA tasks available')"
```

If this fails, do not run the real injection from an ordinary Python environment. Ask how CASA 6 should be started on your system. In an environment with CASA tasks, run:

```bash
python scripts/tutorials/inject_point_source.py \
  --config configs/tutorials/my_visibility_injection.json --dry-run
```

The dry run does not modify the MS. Confirm that the displayed `vis`, `direction`, `flux_jy`, `spw`, and `complist` are correct and that `vis` is the working copy. Then remove `--dry-run` to perform the injection:

```bash
python scripts/tutorials/inject_point_source.py \
  --config configs/tutorials/my_visibility_injection.json
```

On success, the last terminal line gives `data/tutorials/manifests/visibility_injection.json`, where DICE records the inputs, settings, and execution provenance.

A provenance file does not demonstrate that the artificial source was recovered. Inspect the original and injected MSs with identical `plotms` or `listobs` settings, and make dirty images with identical `tclean` parameters. In the injected-minus-original image, confirm a point-source response at the requested sky position with the specified apparent flux. Keep all imaging and detection settings unchanged between the original and injected cases when measuring source finding, completeness, or contamination.

### B-2. Gaussian-line point source

Do not modify the original MS. Use a CASA task such as `mstransform` to make a reference MS containing only one science field and one SPW, and then copy that reference MS to make the injection MS. For safety, DICE currently stops line injection when the MS contains multiple fields or SPWs.

For example, run the following from CASA 6, replacing the field, SPW, and paths with values for your data:

```python
from casatasks import mstransform

mstransform(
    vis="/path/to/original.ms",
    outputvis="data/tutorials/working/line_reference.ms",
    field="target_name",
    spw="0",
    datacolumn="corrected",
    reindex=True,
)
```

Use `datacolumn="data"` only if the input MS has no `CORRECTED_DATA` column. Then copy the reference MS from the terminal:

```bash
cp -R data/tutorials/working/line_reference.ms \
  data/tutorials/working/one_field_one_spw.ms
```

Copy the example JSON and edit `visibility_line_injection`:

```bash
cp configs/tutorials/basic_examples.json \
  configs/tutorials/my_line_injection.json
```

```json
"visibility_line_injection": {
  "vis": "data/tutorials/working/one_field_one_spw.ms",
  "direction": "J2000 03h32m39.53769 -27d46m43.5673",
  "center_frequency_hz": 230000000000.0,
  "fwhm_kms": 300.0,
  "peak_snr": 5.0,
  "rms_jy_per_beam": 0.00002,
  "frequency_frame": "TOPO",
  "complist": "data/tutorials/working/line_model.cl"
}
```

- `direction`: sky position of the artificial line source
- `center_frequency_hz`: observed line-centre frequency in Hz
- `fwhm_kms`: Gaussian-line FWHM in km s\(^{-1}\)
- `peak_snr`: requested peak S/N in the injected dirty cube
- `rms_jy_per_beam`: cube RMS that defines the peak S/N, in Jy beam\(^{-1}\)
- `frequency_frame`: frequency frame matching the MS channel frequencies; `TOPO` is common for ALMA MSs, but check the MS rather than assuming it

For a Gaussian line, peak flux density \(S_{\rm peak}\), FWHM
\(\Delta v_{\rm FWHM}\), and integrated line flux \(I_{\rm line}\) obey

\[
I_{\rm line}
= \int S(v)\,dv
= S_{\rm peak}\,\Delta v_{\rm FWHM}
  \sqrt{\frac{\pi}{4\ln 2}}
\simeq 1.0645\,S_{\rm peak}\,\Delta v_{\rm FWHM}.
\]

DICE therefore uses

\[
S_{\rm peak}
= \frac{I_{\rm line}}
       {1.0645\,\Delta v_{\rm FWHM}}
\]

and samples the resulting Gaussian at the channel centres. For example,
\(I_{\rm line}=1\) Jy km s\(^{-1}\) and
\(\Delta v_{\rm FWHM}=100\) km s\(^{-1}\) give
\(S_{\rm peak}\simeq9.39\) mJy. If the FWHM is comparable to or smaller than
the channel width, the profile is undersampled and the integral measured from
the cube can differ from the continuous value. Record the channel width.

In the standard completeness experiment, **keep the FWHM fixed within one
cube** and vary only peak S/N, sky position, and centre frequency. Use a
separate cube for each FWHM. This gives one source-finding smoothing width per
cube and avoids interactions between different line widths. The tutorial uses
five artificial sources per cube.

`plan_line_injection_positions` can separate their centre channels by four
FWHM. The multi-helix, mixed-width placement remains available as an advanced
option when a large population must be packed into fewer cubes.

![Multi-helix placement for artificial line sources](../tutorials/line_injection_demo/line_injection_helix_example.png)

**Figure: standard artificial-line placement.** The left panel shows the sky
positions over the PB, coloured by input peak S/N. The right panel shows five
Gaussian lines with a fixed FWHM of 10 channels and protected intervals of
four FWHM. A separate cube is used for each other FWHM.

For a point source, DICE converts `peak_snr * rms_jy_per_beam` to the apparent
peak flux and then uses the Gaussian relation above to derive integrated flux.
To inject a physical integrated flux from a luminosity-function model, specify
`integrated_flux_jy_kms` instead of `peak_snr` and `rms_jy_per_beam`. Do not
specify both modes.

Before batching sources, `line_injection_crosstalk` can sample the isolated
dirty response of every source at all target positions. Split the batch if
the summed absolute contribution from the other artificial sources exceeds
the default 0.2σ. Bright-source dirty-beam sidelobes can otherwise change the
recovery of sources close to a 4σ threshold.

First inspect the request with a dry run:

```bash
python scripts/tutorials/inject_line_source.py \
  --config configs/tutorials/my_line_injection.json --dry-run
```

After checking the MS, artificial-source sky position, centre frequency, FWHM, and integrated flux, remove `--dry-run` and run in an environment providing the CASA 6 Python tasks:

```bash
python scripts/tutorials/inject_line_source.py \
  --config configs/tutorials/my_line_injection.json
```

On success, DICE writes `data/tutorials/manifests/visibility_line_injection.json`. As for the continuum case, dirty-image the uninjected reference MS and injected MS with identical settings apart from the corresponding data-column choice. In their difference cube, confirm a Gaussian line at the requested position, frequency, and width.

These tutorials demonstrate visibility injection. They are not required line-search preprocessing.

## Optional Tutorial C. Inspect a simulated cube and run a line search

A useful line-search review includes more than a candidate table: inspect the cube noise and PB as well as the candidates' spatial and frequency distributions. Before using your own data, run the complete search and diagnostic workflow on the small simulated dirty cube and PB:

```bash
python scripts/line_search/search_fits_cube.py tutorials/line_search_demo/demo_config.json
python scripts/line_search/plot_candidates.py tutorials/line_search_demo/demo_config.json
```

A successful run writes a candidate CSV, provenance JSON files, and the following plots under `data/tutorials/line_search_demo/`:

- `candidate_review.pdf`: candidate moment-0 maps with positive/negative contours and spectra
- `candidate_sky_distribution.pdf`: positions of positive and negative candidates within the PB
- `line_search_diagnostics.pdf`: pixel S/N distribution, candidate S/N distribution, and candidate S/N versus channel with the channel noise level

Use the pixel S/N distribution to assess whether the noise is approximately Gaussian. Use the candidate S/N distribution to compare how positive and negative counts change with threshold. Use the channel-S/N panel to check whether candidates cluster in particular channels or coincide with high channel RMS. PNG versions of the final two summaries are also written.

Reference products are included under `tutorials/line_search_demo/expected_output/`. The demo injects only three positive lines, but noise produces additional positive and negative candidates. Compare the catalogue with `tutorials/line_search_demo/expected_candidate_interpretation.csv` to see why a detected candidate is not automatically a real source.

![Spatial distribution of line candidates](../tutorials/line_search_demo/expected_output/candidate_sky_distribution.png)

**Figure: spatial candidate distribution.** Circles mark positive candidates
and downward triangles mark negative candidates. Inspect the usable PB area
and check for clustering toward an image edge or a problematic region.

![Basic line-search diagnostics](../tutorials/line_search_demo/expected_output/line_search_diagnostics.png)

**Figure: basic line-search diagnostics.** The panels summarize the pixel S/N
distribution, candidate counts versus S/N, and candidate S/N versus channel.
The last panel overlays the channel RMS, making it easier to recognize
candidates confined to unusually noisy channels.

## Step 2. Search your own line cube

### 2-1. Check the input data

The current public command does not search an MS directly. It reads two FITS files:

1. **`image_fits`**: a continuum-subtracted, native-resolution dirty spectral cube
2. **`pb_fits`**: the corresponding primary-beam response on the same sky grid and pixel grid

Do not specify a PB-corrected cube (`*.pbcor`) as `image_fits`.

For `pb_fits`, export the CASA `*.pb` image with `exportfits`. Its values are dimensionless, normally close to 1 near the field centre and approaching 0 toward the edge.

```python
from casatasks import exportfits

exportfits(
    imagename="target.contsub.dirty.pb",
    fitsimage="target.contsub.dirty.pb.fits",
    overwrite=True,
)
```

The output of `impbcor` and `*.image.pbcor` are PB-corrected science images, not PB responses. Supplying either as `pb_fits` creates an incorrect search mask and can produce a false zero-candidate result. A FITS `BUNIT` of `Jy/beam` normally indicates that the file is not a PB response.

Also confirm that `image_fits` is continuum-subtracted. The command can run on a cube that retains continuum, but continuum sources may be detected repeatedly at several search widths and frequency positions, producing a cluttered catalogue. Local baselines and line S/N also become harder to interpret, reducing scientific reliability for a normal emission-line search.

When an MS is available, normally subtract the continuum in the visibility domain and make the dirty cube from the continuum-subtracted MS. DICE also supports continuum-normalized absorption-line searches, but the continuum model, normalization, and candidate interpretation differ and are outside this quick start. See `docs/alma_line_search_detectability_guide.md`.

### 2-2. Create the JSON that you will actually use

`configs/line_search/example.json` is a distributed template. Do not edit it directly. Copy it to a new name, edit only the copy, and pass that copied JSON directly to the line-search command. This prevents conflicts when DICE is updated.

```bash
cp configs/line_search/example.json configs/line_search/my_line_search.json
```

Open `configs/line_search/my_line_search.json` and update at least:

- `image_fits`: continuum-subtracted dirty-cube FITS
- `pb_fits`: FITS exported from the CASA `*.pb` image
- `output_directory`: destination for results
- `spw_label`: a label such as the SPW number to retain in the catalogue
- `pb_min`: minimum PB response included in the search, for example 0.5
- `noise_mode`: S/N normalization; normally start with `global`, and consider `local` only after a global/local comparison for a mosaic or clearly spatially varying noise
- `threshold_snr`: candidate-detection S/N threshold
- `width_channels`: searched line widths in channels
- `merge_radius_pixels`: spatial radius used to merge duplicate detections
- `candidate_prefix`: prefix used for candidate IDs

`spw_label` is recorded in the output catalogue but does not change the detection. For a cube from SPW 25, `25` is a useful label.

For `noise_mode: "local"`, also specify a local-noise window relative to the synthesized-beam minor-axis FWHM, for example `"noise_box_beams": 15`. DICE recalculates the local RMS map for each boxcar width. Local mode is not always better: for a uniform-noise single pointing it can add estimation noise. Therefore use `global` for the first search unless the data show a clear reason to use local normalization. Set the searched PB region independently with `pb_min`.

### 2-3. Run the line search

```bash
python scripts/line_search/search_fits_cube.py configs/line_search/my_line_search.json
```

A successful run ends with a line similar to:

```text
Wrote .../signed_line_candidates.csv (3 positive, 2 negative)
```

The command succeeded if this line appears and `output_directory` contains:

- `signed_line_candidates.csv`
- `signed_line_candidates.provenance.json`

`(0 positive, 0 negative)` is still a successful program run, but it does not by itself prove that the cube contains no candidates. Check the input and search region as described in Step 2-4. A Python `Traceback` or `Error` means that the run did not succeed.

When `width_channels` contains several values, DICE performs the boxcar sum and noise re-estimation for each width internally. Do not normally bin the same cube manually beforehand. If the cube is already spectrally smoothed or binned, account for channel correlation and the velocity width per channel when choosing the width list.

Positive emission candidates receive IDs such as `*-E##`; negative absorption candidates receive IDs such as `*-A##` and are searched at the same absolute threshold. This creates a candidate list but does not certify its reliability. Follow `docs/alma_line_search_detectability_guide.md` for continuum and pixel-distribution checks, masks, candidate validation, and injection/recovery.

Astropy may print a `FITSFixedWarning` about time or observatory-coordinate corrections. If processing reaches `Wrote ... signed_line_candidates.csv`, this warning alone normally does not mean failure. Retain the full warning when reporting a problem.

### 2-4. Check the output

At minimum, the output directory should contain:

- `signed_line_candidates.csv`: catalogue of positive and negative candidates
- `signed_line_candidates.provenance.json`: record of inputs, settings, Git commit, and execution environment

A header-only CSV with zero candidates still means that the program completed. Before drawing a scientific conclusion, verify that the data, settings, and searched region match your intention:

1. a sufficient region satisfies `pb >= pb_min`
2. `image_fits` is the intended continuum-subtracted dirty spectral cube
3. the image and PB have corresponding sky grids, shapes, and channel axes
4. the width list covers the expected line velocity widths

Only after checking the inputs should you compare `threshold_snr` values or add widths. Do not choose a favourable threshold after inspecting an untargeted search; record thresholds to be compared in advance.

### 2-5. Make candidate-review plots

Pass the same JSON used for the search to:

```bash
python scripts/line_search/plot_candidates.py configs/line_search/my_line_search.json
```

This writes:

- `candidate_review.pdf`: a moment-0 cutout and spectrum for each candidate
- `candidate_sky_distribution.pdf`: all candidate positions within the PB
- `line_search_diagnostics.pdf`: pixel S/N, candidate S/N, and candidate S/N versus channel with the channel noise level

PNG versions of the final two summaries are also written. These plots support visual review; they do not by themselves establish candidate reliability. No candidate plots are made when the catalogue is empty.

You can change the moment-0 cutout radius and maximum number of candidates. For a 12-pixel radius and the ten candidates with the largest absolute S/N:

```bash
python scripts/line_search/plot_candidates.py configs/line_search/my_line_search.json \
  --radius-pixels 12 --max-candidates 10
```

A successful run ends with a line similar to:

```text
Wrote .../candidate_review.pdf (5 candidates)
```

Plotting succeeded if the output directory contains the figures and `candidate_plots.provenance.json`, and the PDFs open correctly.

## What to send when asking for help

For usage questions, reproducible errors, or unclear documentation, select `New issue` on [DICE Issues](https://github.com/natsuki-h2/DICE/issues). Questions may be written in English or Japanese. Log in to GitHub and use a short title that identifies the problem. While the repository is private, only users with DICE access can use its Issues page.

Include:

- the exact command you ran
- the complete error message
- operating system, Python version, and CASA version
- the configuration JSON
- the generated `.provenance.json`

Do not attach an original MS, unpublished source coordinates, a private candidate catalogue, credentials, or other confidential research data to a GitHub issue. Discuss those through the email or chat channel approved by your laboratory, with your instructor or the person who introduced you to DICE. If unsure, ask first without attaching data.
