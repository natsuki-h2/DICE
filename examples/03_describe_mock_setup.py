"""
Example: describe a first ALMA2040/GHOST mock setup.
"""

from bliss_alma.mockdata import describe_mock_setup


setup = describe_mock_setup(
    configuration="F4",
    output_type="FITS image or MS",
    field_type="blank/noisy",
    notes="First ALMA2040/GHOST setup for continuum point-source recovery tests.",
)

print(setup)
