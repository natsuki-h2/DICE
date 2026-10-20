"""Visibility- and image-domain source injection utilities for DICE."""

import os
import shutil
from collections import namedtuple

import numpy as np


SPEED_OF_LIGHT_KMS = 299792.458

LineInjectionPlacement = namedtuple(
    "LineInjectionPlacement",
    ["source_index", "channel", "y", "x", "fwhm_channels", "helix_index"],
)


def plan_line_injection_positions(
    fwhm_channels,
    cube_shape,
    *,
    spatial_mask=None,
    min_spatial_separation_pixels=5.0,
    spectral_extent_fwhm=3.0,
    spatial_edge_margin_pixels=0,
    seed=0,
    angle_step_degrees=137.507764,
    max_attempts_per_source=10000,
):
    """Place artificial lines on deterministic, interleaved 3-D helices.

    Along each helix, sky positions follow an outward spiral while channel
    centres increase monotonically. Broad lines are packed first and reserve
    a total spectral interval of
    ``spectral_extent_fwhm * fwhm_channels``. When the next line does not fit,
    a new helix starts at the low-channel end with a rotated sky phase. Two
    lines are rejected only when their protected spectral intervals overlap
    and their sky positions are closer than the requested separation. The
    returned order matches the input ``fwhm_channels`` order.
    """
    shape = tuple(int(value) for value in cube_shape)
    if len(shape) != 3 or any(value < 1 for value in shape):
        raise ValueError("cube_shape must be positive (n_channel, n_y, n_x).")
    widths = np.asarray(fwhm_channels, dtype=float)
    if widths.ndim != 1 or widths.size == 0:
        raise ValueError("fwhm_channels must be a non-empty one-dimensional sequence.")
    if not np.all(np.isfinite(widths)) or np.any(widths <= 0):
        raise ValueError("fwhm_channels must contain finite positive values.")
    if min_spatial_separation_pixels < 0:
        raise ValueError("min_spatial_separation_pixels must be non-negative.")
    if spectral_extent_fwhm <= 0:
        raise ValueError("spectral_extent_fwhm must be positive.")
    if spatial_edge_margin_pixels < 0:
        raise ValueError("spatial_edge_margin_pixels must be non-negative.")
    if not np.isfinite(angle_step_degrees) or angle_step_degrees <= 0:
        raise ValueError("angle_step_degrees must be finite and positive.")
    if max_attempts_per_source < 1:
        raise ValueError("max_attempts_per_source must be at least one.")

    n_channel, n_y, n_x = shape
    valid = np.ones((n_y, n_x), dtype=bool)
    if spatial_mask is not None:
        valid = np.asarray(spatial_mask, dtype=bool).copy()
        if valid.shape != (n_y, n_x):
            raise ValueError("spatial_mask must have shape (n_y, n_x).")
    margin = int(spatial_edge_margin_pixels)
    if margin:
        valid[:margin, :] = False
        valid[-margin:, :] = False
        valid[:, :margin] = False
        valid[:, -margin:] = False
    sky_positions = np.argwhere(valid)
    if not sky_positions.size:
        raise ValueError("No valid sky positions remain after masking and margins.")

    # First-fit decreasing packs broad protected intervals into monotonic
    # channel tracks. Each track becomes one helix.
    tracks = []
    for source_index in np.argsort(-widths):
        width = float(widths[source_index])
        half_extent = 0.5 * spectral_extent_fwhm * width
        channel_max = float((n_channel - 1) - half_extent)
        if half_extent > channel_max:
            raise ValueError(
                f"Line {source_index} is too broad for the available channel range."
            )
        assigned = False
        for track in tracks:
            channel = int(np.ceil(track["end"] + half_extent))
            if channel <= channel_max:
                track["sources"].append((int(source_index), width, channel))
                track["end"] = channel + half_extent
                assigned = True
                break
        if not assigned:
            channel = int(np.ceil(half_extent))
            tracks.append({
                "sources": [(int(source_index), width, channel)],
                "end": channel + half_extent,
            })

    centre_y = 0.5 * (n_y - 1)
    centre_x = 0.5 * (n_x - 1)
    radii = np.hypot(sky_positions[:, 0] - centre_y, sky_positions[:, 1] - centre_x)
    radial_order = np.argsort(radii)
    sky_positions = sky_positions[radial_order]
    radii = radii[radial_order]
    base_phase = np.random.default_rng(seed).uniform(0.0, 2.0 * np.pi)
    angle_step = np.deg2rad(angle_step_degrees)
    placed = []
    n_helix = len(tracks)
    for helix_index, track in enumerate(tracks):
        source_count = len(track["sources"])
        phase = base_phase + 2.0 * np.pi * helix_index / max(1, n_helix)
        for step_index, (source_index, width, channel) in enumerate(track["sources"]):
            # Equal-area radial sampling keeps the spiral from crowding its
            # centre. Multiple helices use the same radial run with rotated
            # phase, giving the requested overlapping/interleaved pattern.
            radial_fraction = np.sqrt((step_index + 0.5) / source_count)
            radius_target = float(np.quantile(radii, radial_fraction))
            angle = phase + step_index * angle_step
            target_y = centre_y + radius_target * np.sin(angle)
            target_x = centre_x + radius_target * np.cos(angle)
            candidate_order = np.argsort(
                (sky_positions[:, 0] - target_y) ** 2
                + (sky_positions[:, 1] - target_x) ** 2
            )
            half_extent = 0.5 * spectral_extent_fwhm * width
            for candidate_index in candidate_order[:max_attempts_per_source]:
                y, x = sky_positions[candidate_index]
                conflict = False
                for old in placed:
                    sky_distance = np.hypot(float(x - old.x), float(y - old.y))
                    old_half_extent = 0.5 * spectral_extent_fwhm * old.fwhm_channels
                    spectral_overlap = abs(channel - old.channel) < (
                        half_extent + old_half_extent
                    )
                    if (
                        sky_distance < min_spatial_separation_pixels
                        and spectral_overlap
                    ):
                        conflict = True
                        break
                if not conflict:
                    placed.append(
                        LineInjectionPlacement(
                            source_index=source_index, channel=channel,
                            y=int(y), x=int(x), fwhm_channels=width,
                            helix_index=helix_index,
                        )
                    )
                    break
            else:
                raise RuntimeError(
                    f"Could not place line {source_index} on helix {helix_index} "
                    "without overlap. Increase the cube area, reduce the number "
                    "of sources, or relax the separation settings."
                )
    return sorted(placed, key=lambda row: row.source_index)


def gaussian_line_flux_density(
    frequencies_hz,
    *,
    center_frequency_hz,
    fwhm_kms,
    integrated_flux_jy_kms,
):
    """Return a Gaussian line profile sampled at channel frequencies.

    The velocity convention is the first-order radio definition around the
    requested line centre.  The returned values are apparent flux densities
    in Jy; their continuous Gaussian integral is ``integrated_flux_jy_kms``.
    """
    frequencies = np.asarray(frequencies_hz, dtype=float)
    if frequencies.ndim != 1 or frequencies.size < 3:
        raise ValueError("frequencies_hz must contain at least three channels.")
    if not np.all(np.isfinite(frequencies)) or np.any(frequencies <= 0):
        raise ValueError("frequencies_hz must contain finite positive values.")
    if not np.isfinite(center_frequency_hz) or center_frequency_hz <= 0:
        raise ValueError("center_frequency_hz must be finite and positive.")
    if not np.isfinite(fwhm_kms) or fwhm_kms <= 0:
        raise ValueError("fwhm_kms must be finite and positive.")
    if not np.isfinite(integrated_flux_jy_kms) or integrated_flux_jy_kms < 0:
        raise ValueError("integrated_flux_jy_kms must be finite and non-negative.")
    if not frequencies.min() <= center_frequency_hz <= frequencies.max():
        raise ValueError("center_frequency_hz must lie inside the MS channel range.")

    gaussian_area = fwhm_kms * np.sqrt(np.pi / (4.0 * np.log(2.0)))
    peak_flux_jy = integrated_flux_jy_kms / gaussian_area
    velocity_kms = SPEED_OF_LIGHT_KMS * (
        frequencies - center_frequency_hz
    ) / center_frequency_hz
    profile = peak_flux_jy * np.exp(
        -4.0 * np.log(2.0) * (velocity_kms / fwhm_kms) ** 2
    )
    return profile


def gaussian_line_integrated_flux(
    *,
    peak_snr,
    rms_jy_per_beam,
    fwhm_kms,
):
    """Convert a requested Gaussian-line peak S/N to integrated flux.

    This is the high-level point-source convention used by the tutorial.
    ``rms_jy_per_beam`` must describe the image/cube in which the peak S/N is
    defined. The returned apparent integrated flux is in Jy km s-1.
    """
    values = {
        "peak_snr": peak_snr,
        "rms_jy_per_beam": rms_jy_per_beam,
        "fwhm_kms": fwhm_kms,
    }
    for name, value in values.items():
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive.")
    peak_flux_jy = float(peak_snr) * float(rms_jy_per_beam)
    gaussian_area = float(fwhm_kms) * np.sqrt(np.pi / (4.0 * np.log(2.0)))
    return peak_flux_jy * gaussian_area


def inject_source_image(data, model, position=None, operation="add"):
    """Inject an array model into a 2-D image or 3-D cube.

    Parameters
    ----------
    data : array-like
        Reference image or cube. The input is not modified.
    model : array-like
        Source model in the same units as ``data``. It may have the full data
        shape, or a smaller shape to be placed at ``position``.
    position : sequence of int, optional
        Pixel coordinate of the model centre in data-axis order. Required for
        a smaller model and ignored when both shapes are equal.
    operation : {"add", "subtract"}, optional
        ``subtract`` is convenient for controlled absorption injections.

    Notes
    -----
    DICE does not silently convolve the model or apply a primary beam. Callers
    must supply a PSF-convolved apparent-sky model when injecting into an
    observed image. This keeps units and instrumental assumptions explicit.
    """
    array = np.asarray(data)
    source = np.asarray(model)
    if array.ndim not in (2, 3):
        raise ValueError("data must be a 2-D image or 3-D cube.")
    if source.ndim != array.ndim:
        raise ValueError("model and data must have the same number of dimensions.")
    if operation not in ("add", "subtract"):
        raise ValueError("operation must be 'add' or 'subtract'.")
    if any(model_size > data_size for model_size, data_size in zip(source.shape, array.shape)):
        raise ValueError("model cannot be larger than data on any axis.")

    output = np.array(array, dtype=np.result_type(array, source, np.float64), copy=True)
    sign = 1.0 if operation == "add" else -1.0
    if source.shape == output.shape:
        output += sign * source
        return output
    if position is None or len(position) != output.ndim:
        raise ValueError("position must give one centre pixel per data axis.")

    data_slices = []
    model_slices = []
    for axis, (centre, model_size, data_size) in enumerate(zip(position, source.shape, output.shape)):
        if int(centre) != centre:
            raise ValueError(f"position[{axis}] must be an integer pixel coordinate.")
        centre = int(centre)
        start = centre - model_size // 2
        stop = start + model_size
        data_start = max(start, 0); data_stop = min(stop, data_size)
        if data_start >= data_stop:
            raise ValueError("model does not overlap data at the requested position.")
        data_slices.append(slice(data_start, data_stop))
        model_slices.append(slice(data_start-start, data_stop-start))
    output[tuple(data_slices)] += sign * source[tuple(model_slices)]
    return output

def inject_point_source_visibility(
    vis,
    direction,
    flux_jy,
    spw,
    complist="model.cl",
):
    """
    Inject a point source into ALMA visibility data using CASA tasks.

    Parameters
    ----------
    vis : str
        Input measurement set.
    direction : str
        Source direction, e.g. 'J2000 03h32m39.53769 -27d46m43.5673'.
    flux_jy : float
        Flux density of the injected point source in Jy.
    spw : str
        Spectral window / channel selection passed to CASA `ft`.
        Example: '0:0~0'.
    complist : str
        Name of the temporary CASA component list.
    """

    # Keep CASA optional for image/cube-only DICE installations.
    from casatasks import ft, uvsub
    from casatools import componentlist

    if not os.path.isdir(vis):
        raise FileNotFoundError(f"Measurement set not found: {vis}")

    if flux_jy < 0:
        raise ValueError("flux_jy must be non-negative.")

    if os.path.exists(complist):
        shutil.rmtree(complist)

    cl = componentlist()
    try:
        cl.addcomponent(
            dir=direction,
            flux=flux_jy,
            fluxunit="Jy",
            shape="point",
        )
        cl.rename(complist)
    finally:
        cl.done()

    ft(
        vis=vis,
        spw=spw,
        complist=complist,
        usescratch=True,
    )

    uvsub(vis=vis, reverse=True)

    return vis


def inject_gaussian_line_point_source_visibility(
    vis,
    direction,
    center_frequency_hz,
    fwhm_kms,
    integrated_flux_jy_kms=None,
    peak_snr=None,
    rms_jy_per_beam=None,
    complist="line_model.cl",
    frequency_frame="TOPO",
):
    """Inject one Gaussian-line point source into a one-field, one-SPW MS.

    Specify either the apparent, PB-attenuated ``integrated_flux_jy_kms`` or
    both ``peak_snr`` and ``rms_jy_per_beam``. In the latter high-level mode,
    DICE converts the requested image-plane peak S/N to the integrated flux of
    a point-source Gaussian line. The input MS is modified in place and must
    therefore be a working copy. Restricting this helper to one field and one
    spectral window avoids applying uninitialised model visibilities elsewhere
    when ``uvsub`` operates on the complete MS.
    """
    from casatasks import ft, uvsub
    from casatools import componentlist, table

    if not os.path.isdir(vis):
        raise FileNotFoundError(f"Measurement set not found: {vis}")
    if not frequency_frame or not isinstance(frequency_frame, str):
        raise ValueError("frequency_frame must be a non-empty CASA frame name.")
    uses_integrated_flux = integrated_flux_jy_kms is not None
    uses_peak_snr = peak_snr is not None or rms_jy_per_beam is not None
    if uses_integrated_flux == uses_peak_snr:
        raise ValueError(
            "Specify exactly one flux mode: integrated_flux_jy_kms, or both "
            "peak_snr and rms_jy_per_beam."
        )
    if uses_peak_snr:
        if peak_snr is None or rms_jy_per_beam is None:
            raise ValueError("peak_snr and rms_jy_per_beam must be given together.")
        integrated_flux_jy_kms = gaussian_line_integrated_flux(
            peak_snr=peak_snr,
            rms_jy_per_beam=rms_jy_per_beam,
            fwhm_kms=fwhm_kms,
        )

    tb = table()
    try:
        tb.open(os.path.join(vis, "FIELD"))
        n_fields = tb.nrows()
    finally:
        tb.done()
    if n_fields != 1:
        raise ValueError(
            "Line visibility injection currently requires a one-field MS; "
            "make a working copy selected to one science field first."
        )

    tb = table()
    try:
        tb.open(os.path.join(vis, "SPECTRAL_WINDOW"))
        n_spws = tb.nrows()
        if n_spws == 1:
            channel_frequencies = np.asarray(tb.getcell("CHAN_FREQ", 0), dtype=float)
    finally:
        tb.done()
    if n_spws != 1:
        raise ValueError(
            "Line visibility injection currently requires a one-SPW MS; "
            "make a working copy selected to one spectral window first."
        )

    order = np.argsort(channel_frequencies)
    tabular_frequencies = channel_frequencies[order]
    if np.unique(tabular_frequencies).size != tabular_frequencies.size:
        raise ValueError("MS channel frequencies must be unique.")
    tabular_flux = gaussian_line_flux_density(
        tabular_frequencies,
        center_frequency_hz=center_frequency_hz,
        fwhm_kms=fwhm_kms,
        integrated_flux_jy_kms=integrated_flux_jy_kms,
    )

    if os.path.exists(complist):
        shutil.rmtree(complist)

    cl = componentlist()
    try:
        cl.addcomponent(
            dir=direction,
            flux=float(tabular_flux.max()),
            fluxunit="Jy",
            shape="point",
            freq=f"{frequency_frame} {center_frequency_hz}Hz",
        )
        cl.setspectrum(
            which=0,
            type="tabular",
            tabularfreqs=tabular_frequencies.tolist(),
            tabularflux=tabular_flux.tolist(),
            tabularframe=frequency_frame,
        )
        cl.rename(complist)
    finally:
        cl.done()

    ft(
        vis=vis,
        spw="",
        complist=complist,
        incremental=False,
        usescratch=True,
    )
    uvsub(vis=vis, reverse=True)
    return vis
