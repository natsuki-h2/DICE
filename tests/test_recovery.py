"""Tests for source-recovery measurements and compatibility names."""

import pytest

from bliss_alma.recovery import flux_fidelity, flux_recovery_ratio


def test_flux_recovery_ratio():
    assert flux_recovery_ratio(2.0, 1.5) == pytest.approx(0.75)


def test_legacy_flux_fidelity_name_matches_flux_recovery_ratio():
    assert flux_fidelity(2.0, 1.5) == flux_recovery_ratio(2.0, 1.5)


def test_flux_recovery_ratio_rejects_zero_injected_flux():
    with pytest.raises(ValueError, match="non-zero"):
        flux_recovery_ratio(0.0, 1.0)
