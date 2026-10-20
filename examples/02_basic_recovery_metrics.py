"""
Example: basic recovery metrics with BliSS-ALMA.

This example demonstrates flux-recovery and completeness calculations.
"""

from bliss_alma.recovery import (
    flux_recovery_ratio,
    positional_offset_arcsec,
    is_recovered,
    completeness,
)


input_flux = 1.0e-4      # Jy
recovered_flux = 0.9e-4  # Jy

input_position = (0.0, 0.0)       # arcsec
recovered_position = (0.05, 0.08) # arcsec

max_offset_arcsec = 0.2

print("Flux recovery ratio:", flux_recovery_ratio(input_flux, recovered_flux))
print("Offset [arcsec]:", positional_offset_arcsec(input_position, recovered_position))
print("Recovered:", is_recovered(input_position, recovered_position, max_offset_arcsec))
print("Completeness:", completeness(n_injected=100, n_recovered=82))
