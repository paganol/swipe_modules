"""Tests for SwipeSpinScanningStrategy."""

import numpy as np
import pytest
from astropy.time import Time
from litebird_sim import RotQuaternion

from swipe_modules import SwipeSpinScanningStrategy

START_TIME = Time("2023-01-01", scale="tdb")


def _balloon_trajectory():
    time = Time(
        [
            "2023-01-01T00:00",
            "2023-01-01T04:00",
            "2023-01-01T08:00",
            "2023-01-01T12:00",
            "2023-01-01T16:00",
        ]
    )
    latitude_deg = np.array([78.0, 78.5, 79.0, 79.5, 80.0])
    longitude_deg = np.array([10.0, 11.0, 12.0, 13.0, 14.0])
    return time, latitude_deg, longitude_deg


def test_site_mode_generates_unit_quaternions():
    strategy = SwipeSpinScanningStrategy()
    result = strategy.generate_spin2ecl_quaternions(START_TIME, time_span_s=3600, delta_time_s=60)

    assert isinstance(result, RotQuaternion)
    assert result.quats.shape == (61, 4)
    assert result.sampling_rate_hz == pytest.approx(1.0 / 60)
    np.testing.assert_allclose(np.linalg.norm(result.quats, axis=1), 1.0, atol=1e-10)


def test_generate_requires_astropy_time():
    strategy = SwipeSpinScanningStrategy()
    with pytest.raises(AssertionError):
        strategy.generate_spin2ecl_quaternions(0.0, time_span_s=3600, delta_time_s=60)


@pytest.mark.parametrize("balloon_time_type", [Time, list])
def test_balloon_mode_accepts_time_array_and_plain_list(balloon_time_type):
    """Regression test: balloon_time used to break as a plain list of Time objects."""
    time, latitude_deg, longitude_deg = _balloon_trajectory()
    balloon_time = time if balloon_time_type is Time else list(time)

    strategy = SwipeSpinScanningStrategy(
        balloon_latitude_deg=latitude_deg,
        balloon_longitude_deg=longitude_deg,
        balloon_time=balloon_time,
    )
    result = strategy.generate_spin2ecl_quaternions(
        Time("2023-01-01T02:00", scale="tdb"), time_span_s=3600, delta_time_s=60
    )

    assert result.quats.shape == (61, 4)
    np.testing.assert_allclose(np.linalg.norm(result.quats, axis=1), 1.0, atol=1e-10)
