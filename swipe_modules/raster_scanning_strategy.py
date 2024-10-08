# -*- encoding: utf-8 -*-

import numpy as np
import astropy
from scipy import interpolate
from typing import Union, List
from numba import njit
from litebird_sim import (
    ScanningStrategy,
    RotQuaternion,
)

from litebird_sim.quaternions import (
    quat_rotation_x,
    quat_rotation_y,
    quat_rotation_z,
    quat_left_multiply,
    rotate_x_vector,
    rotate_z_vector,
)
from litebird_sim.imo import Imo
from uuid import UUID

from .common import (
    _ct_jd_to_lst_rad,
    _equinox_precession_rad,
    _equator_ecliptic_angle_rad,
)


@njit
def _sawtooth(
    time,
    offset,
    amplitude,
    speed,
):

    if amplitude == 0:
        return offset
    else:
        r = time / amplitude * speed
        x = 4 * np.abs((r - np.floor(r + 0.5)))  # - 1
        return offset + amplitude * x / 2


@njit
def _SWIPEraster_spin_to_ecliptic(
    result,
    colatitude_rad,
    longitude_rad,
    azimuth_start_rad,
    azimuth_amplitude_rad,
    azimuth_scan_speed_rad_per_s,
    time_s,
    time_jd,
):

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
    eqx = _equinox_precession_rad(time_jd)
    quat_left_multiply(result, *quat_rotation_z(lst + eqx))

    obl = -_equator_ecliptic_angle_rad(time_jd)
    quat_left_multiply(result, *quat_rotation_x(obl))


@njit
def _SWIPEraster_all_spin_to_ecliptic(
    result_matrix,
    colatitude_rad,
    longitude_rad,
    azimuth_start_rad,
    azimuth_amplitude_rad,
    azimuth_scan_speed_rad_per_s,
    time_vector_s,
    time_vector_jd,
):

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
    """A class containing the parameters of the sky scanning strategy
    for SWIPE

    The constructor accepts the following parameters:

    - `site_latitude_deg`: latitude of the launching site in deg

    - `site_longitude_deg`: longitude of the launching site in deg

    - `longitude_speed_deg_per_sec`: longitude speed in deg/sec

    - `azimuth_start_deg`: start azimuth for constant elevation scan

    - `azimuth_amplitude_deg`: amplitude in azimuth for constant elevation scan

    - `azimuth_scan_speed_deg_per_s`: scan speed in azimuth

    - `start_time`: an ``astropy.time.Time`` object representing the
      start of the observation. It's currently unused, but it is meant
      to represent the time when the scan starts

    - `balloon_latitude_deg`: latitude of a tabulated trajectory

    - `balloon_longitude_deg`: longitude of a tabulated trajectory

    - `balloon_time`: list of astropy.time.Time

    These fields are available once the object has been initialized.

    You can create an instance of this class using the class method
    :meth:`.from_imo`, which reads the
    parameters from the IMO.

    """

    def __init__(
        self,
        site_latitude_deg=78.2232,
        site_longitude_deg=15.6267,
        longitude_speed_deg_per_sec=0,
        azimuth_start_deg=220.0,
        azimuth_amplitude_deg=80.0,
        azimuth_scan_speed_deg_per_s=0.1,
        start_time: Union[astropy.time.Time, None] = None,
        balloon_latitude_deg: Union[np.ndarray, None] = None,
        balloon_longitude_deg: Union[np.ndarray, None] = None,
        balloon_time: Union[List[astropy.time.Time], None] = None,
    ):

        self.site_colatitude_rad = np.deg2rad(90.0 - site_latitude_deg)
        self.site_longitude_rad = np.deg2rad(site_longitude_deg)
        self.longitude_speed_rad_per_sec = np.deg2rad(longitude_speed_deg_per_sec)

        self.azimuth_start_rad = np.deg2rad(azimuth_start_deg)
        self.azimuth_amplitude_rad = np.deg2rad(azimuth_amplitude_deg)
        self.azimuth_scan_speed_rad_per_s = np.deg2rad(azimuth_scan_speed_deg_per_s)

        self.start_time = start_time

        if balloon_latitude_deg is None:
            self.balloon_colatitude_rad = None
        else:
            self.balloon_colatitude_rad = np.deg2rad(90.0 - balloon_latitude_deg)

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

    def __repr__(self):
        return (
            (
                "SwipeRasterScanningStrategy(site_colatitude_rad={site_colatitude_rad}, "
                "site_longitude_rad={site_longitude_rad},"
                "longitude_speed_rad_per_sec={longitude_speed_rad_per_sec}, "
                "azimuth_start_rad={azimuth_start_rad}, "
                "azimuth_amplitude_rad={azimuth_amplitude_rad}, "
                "azimuth_scan_speed_rad_per_s={azimuth_scan_speed_rad_per_s}, "
                "start_time={start_time})".format(
                    site_colatitude_rad=self.site_colatitude_rad,
                    site_longitude_rad=self.site_longitude_rad,
                    longitude_speed_rad_per_sec=self.longitude_speed_rad_per_sec,
                    azimuth_start_rad=self.azimuth_start_rad,
                    azimuth_amplitude_rad=self.azimuth_amplitude_rad,
                    azimuth_scan_speed_rad_per_s=self.azimuth_scan_speed_rad_per_s,
                    start_time=self.start_time,
                )
            )
            if (
                (self.balloon_colatitude_rad is None)
                and (self.balloon_longitude_rad is None)
            )
            else (
                "SwipeRasterScanningStrategy(colatitude_range_rad=[{min_colatitude_rad},{max_colatitude_rad}],"
                "azimuth_start_rad={azimuth_start_rad}, "
                "azimuth_amplitude_rad={azimuth_amplitude_rad}, "
                "azimuth_scan_speed_rad_per_s={azimuth_scan_speed_rad_per_s}, "
                "start_time={start_time})".format(
                    min_colatitude_rad=self.balloon_colatitude_rad.min(),
                    max_colatitude_rad=self.balloon_colatitude_rad.max(),
                    azimuth_start_rad=self.azimuth_start_rad,
                    azimuth_amplitude_rad=self.azimuth_amplitude_rad,
                    azimuth_scan_speed_rad_per_s=self.azimuth_scan_speed_rad_per_s,
                    start_time=self.start_time,
                )
            )
        )

    def all_spin_to_ecliptic(
        self,
        result_matrix,
        colatitude_rad,
        longitude_rad,
        azimuth_start_rad,
        azimuth_amplitude_rad,
        azimuth_scan_speed_rad_per_s,
        time_vector_s,
        time_vector_jd,
    ):
        assert result_matrix.shape == (len(time_vector_s), 4)
        assert len(time_vector_jd) == len(time_vector_s)

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
    def from_imo(imo: Imo, url: Union[str, UUID]):
        """Read the definition of the scanning strategy from the IMO

        This function returns a :class:`.SwipeRasterScanningStrategy`
        object containing the set of parameters that define the
        scanning strategy of the balloon.

        Args:

            imo (:class:`.Imo`): an instance of the :class:`.Imo` class

            url (str or ``UUID``): a reference to the data file
                containing the definition of the scanning strategy. It can
                be either a string like
                ``/releases/v0.0/balloon/scanning_parameters/`` or a
                UUID.

        Example::

            imo = Imo()
            sstr = SwipeRasterScanningStrategy.from_imo(
                imo=imo,
                url="/releases/v0.0/balloon/scanning_parameters/",
            )
            print(sstr)

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

        assert type(start_time) == astropy.time.Time

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

            time_jd = [d.jd for d in time]
            balloon_time_jd = [d.jd for d in self.balloon_time]

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
