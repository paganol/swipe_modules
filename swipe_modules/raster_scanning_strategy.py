"""Raster scanning strategy for SWIPE."""

from typing import Final
from uuid import UUID

import astropy.time
import astropy.units
import numpy as np
from litebird_sim import RotQuaternion, ScanningStrategy
from litebird_sim.imo import Imo
from litebird_sim.quaternions import (
    quat_left_multiply,
    quat_rotation_x,
    quat_rotation_y,
    quat_rotation_z,
)
from numba import njit
from scipy import interpolate

from .common import (
    _ct_jd_to_lst_rad,
    _equator_ecliptic_angle_rad,
)

__all__ = ["SwipeRasterScanningStrategy"]

# Scan pattern constants
_DEFAULT_SITE_LATITUDE_DEG: Final[float] = 78.2232
_DEFAULT_SITE_LONGITUDE_DEG: Final[float] = 15.6267
_DEFAULT_AZIMUTH_START_DEG: Final[float] = 220.0
_DEFAULT_AZIMUTH_AMPLITUDE_DEG: Final[float] = 80.0
_DEFAULT_AZIMUTH_SCAN_SPEED_DEG_PER_S: Final[float] = 0.1
_EQUATORIAL_COLATITUDE_DEG: Final[float] = 90.0


@njit
def _sawtooth(
    time: float,
    offset: float,
    amplitude: float,
    speed: float,
) -> float:
    """Generate a sawtooth waveform.

    Args:
        time: Time value.
        offset: Offset of the sawtooth.
        amplitude: Amplitude of the sawtooth.
        speed: Speed parameter.

    Returns:
        Sawtooth waveform value.
    """

    if amplitude == 0:
        return offset
    else:
        r = time / amplitude * speed
        x = 4 * np.abs(r - np.floor(r + 0.5))  # - 1
        return offset + amplitude * x / 2


@njit
def _SWIPEraster_spin_to_ecliptic(
    result: np.ndarray,
    colatitude_rad: float,
    longitude_rad: float,
    azimuth_start_rad: float,
    azimuth_amplitude_rad: float,
    azimuth_scan_speed_rad_per_s: float,
    time_s: float,
    time_jd: float,
) -> None:
    """Compute spin-to-ecliptic quaternion for raster scanning.

    Args:
        result: Output quaternion array (4-element).
        colatitude_rad: Latitude in radians.
        longitude_rad: Longitude in radians.
        azimuth_start_rad: Starting azimuth in radians.
        azimuth_amplitude_rad: Azimuth amplitude in radians.
        azimuth_scan_speed_rad_per_s: Scanning speed in rad/s.
        time_s: Time in seconds.
        time_jd: Time in Julian Days.
    """

    result[:] = quat_rotation_z(
        np.pi
        - _sawtooth(
            time_s,
            azimuth_start_rad,
            azimuth_amplitude_rad,
            azimuth_scan_speed_rad_per_s,
        )
    )

    quat_left_multiply(result, *quat_rotation_y(colatitude_rad))

    lst = _ct_jd_to_lst_rad(time_jd, longitude_rad)
    quat_left_multiply(result, *quat_rotation_z(lst))

    obl = -_equator_ecliptic_angle_rad(time_jd)
    quat_left_multiply(result, *quat_rotation_x(obl))


@njit
def _SWIPEraster_all_spin_to_ecliptic(
    result_matrix: np.ndarray,
    colatitude_rad: np.ndarray,
    longitude_rad: np.ndarray,
    azimuth_start_rad: float,
    azimuth_amplitude_rad: float,
    azimuth_scan_speed_rad_per_s: float,
    time_vector_s: np.ndarray,
    time_vector_jd: np.ndarray,
) -> None:
    """Compute spin-to-ecliptic quaternions for all times in raster scanning.

    Args:
        result_matrix: Output quaternion matrix (N x 4).
        colatitude_rad: Colatitude values in radians.
        longitude_rad: Longitude values in radians.
        azimuth_start_rad: Starting azimuth in radians.
        azimuth_amplitude_rad: Azimuth amplitude in radians.
        azimuth_scan_speed_rad_per_s: Scanning speed in rad/s.
        time_vector_s: Time vector in seconds.
        time_vector_jd: Time vector in Julian Days.
    """

    for row in range(result_matrix.shape[0]):
        _SWIPEraster_spin_to_ecliptic(
            result=result_matrix[row, :],
            colatitude_rad=colatitude_rad[row],
            longitude_rad=longitude_rad[row],
            azimuth_start_rad=azimuth_start_rad,
            azimuth_amplitude_rad=azimuth_amplitude_rad,
            azimuth_scan_speed_rad_per_s=azimuth_scan_speed_rad_per_s,
            time_s=time_vector_s[row],
            time_jd=time_vector_jd[row],
        )


class SwipeRasterScanningStrategy(ScanningStrategy):
    """Sky scanning strategy for SWIPE with raster scanning pattern.

    This class defines the scanning parameters for the SWIPE balloon-borne
    instrument, supporting both fixed sites and balloon trajectories.

    Attributes:
        site_colatitude_rad (float): Colatitude of the site in radians.
        site_longitude_rad (float): Longitude of the site in radians.
        longitude_speed_rad_per_sec (float): Longitude drift speed in rad/s.
        azimuth_start_rad (float): Starting azimuth in radians.
        azimuth_amplitude_rad (float): Azimuth scan amplitude in radians.
        azimuth_scan_speed_rad_per_s (float): Azimuth scan speed in rad/s.
        start_time (astropy.time.Time | None): Start time of observation.
        balloon_colatitude_rad (ndarray | None): Balloon colatitude trajectory.
        balloon_longitude_rad (ndarray | None): Balloon longitude trajectory.
        balloon_time (list[astropy.time.Time] | None): Balloon time points.

    Example:
        Use :meth:`.from_imo` to create an instance from IMO parameters:

        >>> imo = Imo()
        >>> sstr = SwipeRasterScanningStrategy.from_imo(
        ...     imo=imo,
        ...     url="/releases/v0.0/balloon/scanning_parameters/",
        ... )
    """

    def __init__(
        self,
        site_latitude_deg: float = _DEFAULT_SITE_LATITUDE_DEG,
        site_longitude_deg: float = _DEFAULT_SITE_LONGITUDE_DEG,
        longitude_speed_deg_per_sec: float = 0.0,
        azimuth_start_deg: float = _DEFAULT_AZIMUTH_START_DEG,
        azimuth_amplitude_deg: float = _DEFAULT_AZIMUTH_AMPLITUDE_DEG,
        azimuth_scan_speed_deg_per_s: float = _DEFAULT_AZIMUTH_SCAN_SPEED_DEG_PER_S,
        start_time: astropy.time.Time | None = None,
        balloon_latitude_deg: np.ndarray | None = None,
        balloon_longitude_deg: np.ndarray | None = None,
        balloon_time: list[astropy.time.Time] | None = None,
    ) -> None:

        self.site_colatitude_rad = np.deg2rad(_EQUATORIAL_COLATITUDE_DEG - site_latitude_deg)
        self.site_longitude_rad = np.deg2rad(site_longitude_deg)
        self.longitude_speed_rad_per_sec = np.deg2rad(longitude_speed_deg_per_sec)

        self.azimuth_start_rad = np.deg2rad(azimuth_start_deg)
        self.azimuth_amplitude_rad = np.deg2rad(azimuth_amplitude_deg)
        self.azimuth_scan_speed_rad_per_s = np.deg2rad(azimuth_scan_speed_deg_per_s)

        self.start_time = start_time

        if balloon_latitude_deg is None:
            self.balloon_colatitude_rad = None
        else:
            self.balloon_colatitude_rad = np.deg2rad(
                _EQUATORIAL_COLATITUDE_DEG - balloon_latitude_deg
            )

        if balloon_longitude_deg is None:
            self.balloon_longitude_rad = None
        else:
            self.balloon_longitude_rad = np.deg2rad(balloon_longitude_deg)

        if (balloon_latitude_deg is None) and (balloon_longitude_deg is None):
            print(
                "site_latitude_deg, site_longitude_deg and longitude_speed_deg_per_sec used"
            )
        else:
            print(
                "site_latitude_deg, site_longitude_deg and longitude_speed_deg_per_sec ignored"
            )
            print("a tabulated trajectory will be used")

        self.balloon_time = balloon_time

    def __repr__(self) -> str:
        return (
            (
                f"SwipeRasterScanningStrategy(site_colatitude_rad={self.site_colatitude_rad}, "
                f"site_longitude_rad={self.site_longitude_rad},"
                f"longitude_speed_rad_per_sec={self.longitude_speed_rad_per_sec}, "
                f"azimuth_start_rad={self.azimuth_start_rad}, "
                f"azimuth_amplitude_rad={self.azimuth_amplitude_rad}, "
                f"azimuth_scan_speed_rad_per_s={self.azimuth_scan_speed_rad_per_s}, "
                f"start_time={self.start_time})"
            )
            if (
                (self.balloon_colatitude_rad is None)
                and (self.balloon_longitude_rad is None)
            )
            else (
                f"SwipeRasterScanningStrategy(colatitude_range_rad=[{self.balloon_colatitude_rad.min()},{self.balloon_colatitude_rad.max()}],"
                f"azimuth_start_rad={self.azimuth_start_rad}, "
                f"azimuth_amplitude_rad={self.azimuth_amplitude_rad}, "
                f"azimuth_scan_speed_rad_per_s={self.azimuth_scan_speed_rad_per_s}, "
                f"start_time={self.start_time})"
            )
        )

    def all_spin_to_ecliptic(
        self,
        result_matrix: np.ndarray,
        colatitude_rad: np.ndarray,
        longitude_rad: np.ndarray,
        azimuth_start_rad: float,
        azimuth_amplitude_rad: float,
        azimuth_scan_speed_rad_per_s: float,
        time_vector_s: np.ndarray,
        time_vector_jd: np.ndarray,
    ) -> None:
        """Compute spin-to-ecliptic quaternions for array of times.

        Args:
            result_matrix: Output quaternion matrix (N x 4).
            colatitude_rad: Colatitude values in radians.
            longitude_rad: Longitude values in radians.
            azimuth_start_rad: Starting azimuth in radians.
            azimuth_amplitude_rad: Azimuth amplitude in radians.
            azimuth_scan_speed_rad_per_s: Scanning speed in rad/s.
            time_vector_s: Time vector in seconds.
            time_vector_jd: Time vector in Julian Days.
        """
        assert result_matrix.shape == (len(time_vector_s), 4), "Result matrix size mismatch"
        assert len(time_vector_jd) == len(time_vector_s), "Time vector size mismatch"

        _SWIPEraster_all_spin_to_ecliptic(
            result_matrix=result_matrix,
            colatitude_rad=colatitude_rad,
            longitude_rad=longitude_rad,
            azimuth_start_rad=azimuth_start_rad,
            azimuth_amplitude_rad=azimuth_amplitude_rad,
            azimuth_scan_speed_rad_per_s=azimuth_scan_speed_rad_per_s,
            time_vector_s=time_vector_s,
            time_vector_jd=time_vector_jd,
        )

    @staticmethod
    def from_imo(imo: Imo, url: str | UUID) -> "SwipeRasterScanningStrategy":
        """Read scanning strategy parameters from the IMO database.

        Args:
            imo: An IMO database instance.
            url: IMO reference path or UUID for the scanning parameters.
                Example: "/releases/v0.0/balloon/scanning_parameters/"

        Returns:
            SwipeRasterScanningStrategy instance with IMO parameters.

        Example:
            >>> imo = Imo()
            >>> sstr = SwipeRasterScanningStrategy.from_imo(
            ...     imo=imo,
            ...     url="/releases/v0.0/balloon/scanning_parameters/",
            ... )
            >>> print(sstr)
        """
        obj = imo.query(url)
        return SwipeRasterScanningStrategy(
            site_latitude_deg=obj.metadata["site_latitude_deg"],
            site_longitude_deg=obj.metadata["site_longitude_deg"],
            longitude_speed_deg_per_sec=obj.metadata["longitude_speed_deg_per_sec"],
            azimuth_start_deg=obj.metadata["azimuth_start_deg"],
            azimuth_amplitude_deg=obj.metadata["azimuth_amplitude_deg"],
            azimuth_scan_speed_deg_per_s=obj.metadata["azimuth_scan_speed_deg_per_s"],
        )

    def generate_spin2ecl_quaternions(
        self,
        start_time: astropy.time.Time,
        time_span_s: float,
        delta_time_s: float,
    ) -> RotQuaternion:
        """Generate spin-to-ecliptic quaternions for a time span.

        Args:
            start_time: Start time of quaternion generation.
            time_span_s: Duration to cover in seconds.
            delta_time_s: Sampling interval in seconds.

        Returns:
            RotQuaternion: Quaternion time series.
        """

        assert isinstance(start_time, astropy.time.Time), "start_time must be astropy.time.Time"

        pointing_freq_hz = 1.0 / delta_time_s

        num_of_quaternions = ScanningStrategy.optimal_num_of_quaternions(
            time_span_s=time_span_s, delta_time_s=delta_time_s
        )

        spin2ecliptic_quats = np.empty((num_of_quaternions, 4))

        time, time_s = ScanningStrategy.get_times(
            start_time=start_time,
            delta_time_s=delta_time_s,
            num_of_quaternions=num_of_quaternions,
        )

        if self.balloon_time is None:

            colatitude_rad = np.repeat(self.site_colatitude_rad, num_of_quaternions)
            longitude_rad = np.mod(
                self.longitude_speed_rad_per_sec * time_s + self.site_longitude_rad,
                2 * np.pi,
            )
            time_jd = time.jd

        else:
            assert (
                len(self.balloon_colatitude_rad)
                == len(self.balloon_longitude_rad)
                == len(self.balloon_time)
            )

            assert self.balloon_time[0] <= start_time

            end_time = start_time + time_span_s * astropy.units.second

            assert self.balloon_time[-1] >= end_time

            time_jd = time.jd
            balloon_time_jd = astropy.time.Time(self.balloon_time).jd

            # interpolate
            fcolat = interpolate.interp1d(
                balloon_time_jd, self.balloon_colatitude_rad, kind="cubic"
            )
            colatitude_rad = fcolat(time_jd)

            flon = interpolate.interp1d(
                balloon_time_jd, np.unwrap(self.balloon_longitude_rad), kind="cubic"
            )
            longitude_rad = np.mod(flon(time_jd), 2 * np.pi)

        self.all_spin_to_ecliptic(
            result_matrix=spin2ecliptic_quats,
            colatitude_rad=colatitude_rad,
            longitude_rad=longitude_rad,
            azimuth_start_rad=self.azimuth_start_rad,
            azimuth_amplitude_rad=self.azimuth_amplitude_rad,
            azimuth_scan_speed_rad_per_s=self.azimuth_scan_speed_rad_per_s,
            time_vector_s=time_s,
            time_vector_jd=time_jd,
        )

        return RotQuaternion(
            start_time=start_time,
            sampling_rate_hz=pointing_freq_hz,
            quats=spin2ecliptic_quats,
        )
