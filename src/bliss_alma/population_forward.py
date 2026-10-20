"""Source-population forward-modelling utilities for facility comparisons.

The functions in this module separate an intrinsic population model from an
instrumental observing operator. Luminosity functions are
expressed per dex and catalogues are drawn as Poisson point processes in
``(redshift, log10 luminosity)``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import trapezoid


FIDUCIAL_COSMOLOGY = FlatLambdaCDM(H0=70.0, Om0=0.3, Tcmb0=2.725)


@dataclass(frozen=True)
class SaundersLF:
    """Modified-Schechter (Saunders) luminosity function parameters."""

    alpha: float
    log10_lstar_lsun: float
    log10_phistar_mpc3_dex: float
    sigma: float = 0.5


@dataclass(frozen=True)
class DoublePowerLawLF:
    """Double-power-law luminosity function parameters."""

    alpha_bright: float
    beta_faint: float
    log10_lstar_lsun: float
    log10_phistar_mpc3_dex: float


ALPINE_IR_Z45_60 = SaundersLF(
    alpha=1.28,
    log10_lstar_lsun=12.16,
    log10_phistar_mpc3_dex=-3.73,
    sigma=0.5,
)

# Hayatsu et al. (2025) adopts the Hemmati et al. local DPL shape and leaves
# the high-z amplitude to be constrained by the untargeted search. This object is
# therefore a benchmark normalization, not a claimed z~6 measurement.
HAYATSU_CII_BENCHMARK = DoublePowerLawLF(
    alpha_bright=2.36,
    beta_faint=0.42,
    log10_lstar_lsun=np.log10(2.173e8),
    log10_phistar_mpc3_dex=np.log10(0.003),
)


def saunders_phi(log10_l_lsun, parameters=ALPINE_IR_Z45_60):
    """Return ``dN / (dV dlog10L)`` in Mpc^-3 dex^-1."""

    ratio = 10.0 ** (np.asarray(log10_l_lsun, dtype=float) - parameters.log10_lstar_lsun)
    phistar = 10.0**parameters.log10_phistar_mpc3_dex
    return phistar * ratio ** (1.0 - parameters.alpha) * np.exp(
        -(np.log10(1.0 + ratio) ** 2) / (2.0 * parameters.sigma**2)
    )


def double_power_law_phi(log10_l_lsun, parameters=HAYATSU_CII_BENCHMARK):
    """Return ``dN / (dV dlog10L)`` in Mpc^-3 dex^-1."""

    ratio = 10.0 ** (np.asarray(log10_l_lsun, dtype=float) - parameters.log10_lstar_lsun)
    phistar = 10.0**parameters.log10_phistar_mpc3_dex
    return phistar / (
        ratio**parameters.alpha_bright + ratio**parameters.beta_faint
    )


def survey_volume_mpc3(area_arcmin2, z_min, z_max, cosmology=FIDUCIAL_COSMOLOGY):
    """Comoving volume in a pencil beam."""

    sky_fraction = (area_arcmin2 * u.arcmin**2).to_value(u.sr) / (4.0 * np.pi)
    shell = cosmology.comoving_volume(z_max) - cosmology.comoving_volume(z_min)
    return float(shell.to_value(u.Mpc**3) * sky_fraction)


def expected_count(
    phi_function,
    log10_l_min,
    log10_l_max,
    area_arcmin2,
    z_min,
    z_max,
    *,
    n_grid=4097,
    cosmology=FIDUCIAL_COSMOLOGY,
):
    """Integrate a redshift-independent per-dex LF over a survey volume."""

    grid = np.linspace(log10_l_min, log10_l_max, n_grid)
    density = trapezoid(phi_function(grid), grid)
    return float(density * survey_volume_mpc3(area_arcmin2, z_min, z_max, cosmology))


def _inverse_cdf_samples(grid, weights, size, rng):
    cdf = np.cumsum(np.asarray(weights, dtype=float))
    if cdf[-1] <= 0:
        raise ValueError("Sampling weights must have a positive integral")
    cdf /= cdf[-1]
    return np.interp(rng.random(size), cdf, grid)


def draw_poisson_catalogue(
    phi_function,
    log10_l_min,
    log10_l_max,
    area_arcmin2,
    z_min,
    z_max,
    *,
    seed,
    field_count=10,
    n_grid=4097,
    cosmology=FIDUCIAL_COSMOLOGY,
):
    """Draw one physical catalogue with uniform sky positions.

    Redshifts are sampled from the comoving-volume element and luminosities
    from the supplied, redshift-independent LF.
    """

    rng = np.random.default_rng(seed)
    mean_count = expected_count(
        phi_function,
        log10_l_min,
        log10_l_max,
        area_arcmin2,
        z_min,
        z_max,
        n_grid=n_grid,
        cosmology=cosmology,
    )
    count = int(rng.poisson(mean_count))

    log_grid = np.linspace(log10_l_min, log10_l_max, n_grid)
    log_l = _inverse_cdf_samples(log_grid, phi_function(log_grid), count, rng)
    z_grid = np.linspace(z_min, z_max, n_grid)
    dvolume = cosmology.differential_comoving_volume(z_grid).to_value(
        u.Mpc**3 / u.sr
    )
    redshift = _inverse_cdf_samples(z_grid, dvolume, count, rng)

    field_area = area_arcmin2 / field_count
    side_arcsec = np.sqrt(field_area) * 60.0
    field_id = rng.integers(0, field_count, size=count)
    offset_ra = rng.uniform(-side_arcsec / 2.0, side_arcsec / 2.0, size=count)
    offset_dec = rng.uniform(-side_arcsec / 2.0, side_arcsec / 2.0, size=count)

    return {
        "expected_count": mean_count,
        "field_id": field_id,
        "offset_ra_arcsec": offset_ra,
        "offset_dec_arcsec": offset_dec,
        "redshift": redshift,
        "log10_l_lsun": log_l,
    }


def line_luminosity_to_flux_jykms(
    luminosity_lsun,
    redshift,
    rest_frequency_ghz=1900.5369,
    cosmology=FIDUCIAL_COSMOLOGY,
):
    """Convert line luminosity to velocity-integrated flux."""

    distance_mpc = cosmology.luminosity_distance(redshift).to_value(u.Mpc)
    return (
        np.asarray(luminosity_lsun)
        * (1.0 + np.asarray(redshift))
        / (1.04e-3 * rest_frequency_ghz * distance_mpc**2)
    )


def gaussian_line_peak_jy(integrated_flux_jykms, fwhm_kms):
    """Peak flux density of a Gaussian line."""

    gaussian_area_factor = np.sqrt(np.pi / (4.0 * np.log(2.0)))
    return np.asarray(integrated_flux_jykms) / (
        gaussian_area_factor * np.asarray(fwhm_kms)
    )


def modified_blackbody_flux_jy(
    lir_lsun,
    redshift,
    observed_frequency_ghz,
    *,
    dust_temperature_k=35.0,
    beta=1.8,
    include_cmb=True,
    cosmology=FIDUCIAL_COSMOLOGY,
):
    """Convert total IR luminosity to flux using an optically thin MBB.

    The SED is normalized over rest-frame 8--1000 micron.  When requested,
    CMB heating and background contrast follow the common da Cunha-style
    modified-blackbody prescription.
    """

    redshift = np.asarray(redshift, dtype=float)
    lir_lsun = np.asarray(lir_lsun, dtype=float)
    power = 4.0 + beta
    tcmb0 = cosmology.Tcmb0.to_value(u.K)
    if include_cmb:
        dust_t = (
            dust_temperature_k**power
            + tcmb0**power * ((1.0 + redshift) ** power - 1.0)
        ) ** (1.0 / power)
        background_t = tcmb0 * (1.0 + redshift)
    else:
        dust_t = np.full_like(redshift, dust_temperature_k)
        background_t = np.zeros_like(redshift)

    nu_grid = np.geomspace(
        (const.c / (1000.0 * u.micron)).to_value(u.Hz),
        (const.c / (8.0 * u.micron)).to_value(u.Hz),
        8192,
    )

    def planck(nu_hz, temperature):
        x = const.h.value * nu_hz / (const.k_B.value * temperature)
        return 2.0 * const.h.value * nu_hz**3 / const.c.value**2 / np.expm1(x)

    output = np.empty(np.broadcast(lir_lsun, redshift).shape, dtype=float)
    luminosities, redshifts = np.broadcast_arrays(lir_lsun, redshift)
    dust_temperatures = np.broadcast_to(dust_t, output.shape)
    background_temperatures = np.broadcast_to(background_t, output.shape)
    for index in np.ndindex(output.shape):
        td = dust_temperatures[index]
        bg = background_temperatures[index]
        shape_grid = nu_grid**beta * planck(nu_grid, td)
        if include_cmb:
            shape_grid -= nu_grid**beta * planck(nu_grid, bg)
        normalization = trapezoid(shape_grid, nu_grid)
        nu_rest = observed_frequency_ghz * 1e9 * (1.0 + redshifts[index])
        shape_at_nu = nu_rest**beta * planck(nu_rest, td)
        if include_cmb:
            shape_at_nu -= nu_rest**beta * planck(nu_rest, bg)
        lnu_w_hz = luminosities[index] * const.L_sun.value * shape_at_nu / normalization
        distance_m = cosmology.luminosity_distance(redshifts[index]).to_value(u.m)
        flux_w_m2_hz = (1.0 + redshifts[index]) * lnu_w_hz / (
            4.0 * np.pi * distance_m**2
        )
        output[index] = flux_w_m2_hz / 1e-26
    return output
