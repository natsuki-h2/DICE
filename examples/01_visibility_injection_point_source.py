"""
Example: visibility-domain point-source injection with BliSS-ALMA.

This example is intended to be run inside CASA after preparing a measurement set.
"""

from bliss_alma.injection import inject_point_source_visibility


vis = "blank_field.ms"

direction = "J2000 03h32m39.53769 -27d46m43.5673"
flux_jy = 1.0e-4  # 0.1 mJy

inject_point_source_visibility(
    vis=vis,
    direction=direction,
    flux_jy=flux_jy,
    spw="",  # continuum injection
    complist="model.cl",
)