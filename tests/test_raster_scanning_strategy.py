"""Tests for SwipeRasterScanningStrategy and its sawtooth scan pattern."""

import numpy as np
import pytest
from astropy.time import Time
from litebird_sim import RotQuaternion

from swipe_modules import SwipeRasterScanningStrategy
from swipe_modules.raster_scanning_strategy import _sawtooth

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
    strategy = SwipeRasterScanningStrategy()
    result = strategy.generate_spin2ecl_quaternions(START_TIME, time_span_s=3600, delta_time_s=60)

    assert isinstance(result, RotQuaternion)
    assert result.quats.shape == (61, 4)
    assert result.sampling_rate_hz == pytest.approx(1.0 / 60)
    np.testing.assert_allclose(np.linalg.norm(result.quats, axis=1), 1.0, atol=1e-10)


def test_generate_requires_astropy_time():
    strategy = SwipeRasterScanningStrategy()
    with pytest.raises(AssertionError):
        strategy.generate_spin2ecl_quaternions(0.0, time_span_s=3600, delta_time_s=60)


@pytest.mark.parametrize("balloon_time_type", [Time, list])
def test_balloon_mode_accepts_time_array_and_plain_list(balloon_time_type):
    """Regression test: balloon_time used to break as a vectorized Time array."""
    time, latitude_deg, longitude_deg = _balloon_trajectory()
    balloon_time = time if balloon_time_type is Time else list(time)

    strategy = SwipeRasterScanningStrategy(
        balloon_latitude_deg=latitude_deg,
        balloon_longitude_deg=longitude_deg,
        balloon_time=balloon_time,
    )
    result = strategy.generate_spin2ecl_quaternions(
        Time("2023-01-01T02:00", scale="tdb"), time_span_s=3600, delta_time_s=60
    )

    assert result.quats.shape == (61, 4)
    np.testing.assert_allclose(np.linalg.norm(result.quats, axis=1), 1.0, atol=1e-10)


def test_sawtooth_zero_amplitude_returns_offset():
    assert _sawtooth(123.0, offset=5.0, amplitude=0.0, speed=0.1) == 5.0


def test_sawtooth_stays_within_offset_and_amplitude_range():
    offset, amplitude, speed = 10.0, 80.0, 0.1
    period = amplitude / speed
    times = np.linspace(0, 3 * period, 500)
    values = np.array([_sawtooth(t, offset, amplitude, speed) for t in times])

    assert values.min() >= offset - 1e-9
    assert values.max() <= offset + amplitude + 1e-9


def test_sawtooth_is_periodic():
    offset, amplitude, speed = 10.0, 80.0, 0.1
    period = amplitude / speed
    t = 37.0
    assert _sawtooth(t, offset, amplitude, speed) == pytest.approx(
        _sawtooth(t + period, offset, amplitude, speed)
    )
