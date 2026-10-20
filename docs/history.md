# Project history

DICE grew out of an earlier repository, **GenerateMockData**, originally developed in 2018 for ALMA mock-data generation and artificial-source injection tests.

The original purpose of `GenerateMockData` was to compare different ways of injecting artificial sources into real or simulated ALMA image cubes and visibility data.

The code was originally run inside CASA using:

```python
execfile('GenerateMockdata.py')

GenerateMockdata/
├── README
├── GenerateMockdata.py
└── Template/
    ├── configurationfiles/
    ├── ptgfiles/
    └── MeasurementSet/

The configuration files were created using aU.buildConfigurationFile, and the pointing data were based on QA2 reports downloaded from the ALMA archive.

The code was shared with project members on 22 March 2018 by Natsuki Hayatsu.

In 2026, this old repository was first reorganized as BliSS-ALMA and then
renamed DICE: Discovery and Image/Cube Evaluation, with the aim of developing
a public injection--recovery framework for measuring completeness,
contamination, flux recovery, and astrometric accuracy in astronomical images
and cubes. Source extraction is a component of that framework; the broader aim
is to propagate the measured response into corrected source counts and
population constraints.
