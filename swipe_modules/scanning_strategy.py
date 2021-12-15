# -*- encoding: utf-8 -*-

import numpy as np
import astropy
from typing import Union, List
from numba import njit
from litebird_sim import (
    ScanningStrategy,
    calculate_sun_earth_angles_rad,
    Spin2EclipticQuaternions,
)

from litebird_sim.quaternions import (
    quat_rotation_x,
    quat_rotation_y,
    quat_rotation_z,
    quat_left_multiply,
    rotate_x_vector,
    rotate_z_vector,
)


EQUATOR_ECLIPTIC_ANGLE_RAD = 0.408407045  # 23.4 deg in radians


@njit
def SWIPE_spin_to_ecliptic(
    result,
    sun_earth_angle_rad,
    site_colatitude_rad,
    site_longitude_rad,
    longitude_speed_rad_per_sec,
    spin_rate_hz,
    time_s,
):

    result[:] = quat_rotation_z(2 * np.pi * spin_rate_hz * time_s)
    quat_left_multiply(result, *quat_rotation_y(site_colatitude_rad))
    quat_left_multiply(
        result,
        *quat_rotation_z((longitude_speed_rad_per_sec + 2 * np.pi / 24 / 3600) * time_s + site_longitude_rad)
    )
    quat_left_multiply(result, *quat_rotation_x(-EQUATOR_ECLIPTIC_ANGLE_RAD))
    quat_left_multiply(result, *quat_rotation_z(sun_earth_angle_rad))


@njit
def SWIPE_all_spin_to_ecliptic(
    result_matrix,
    sun_earth_angles_rad,
    site_colatitude_rad,
    site_longitude_rad,
    longitude_speed_rad_per_sec,
    spin_rate_hz,
    time_vector_s,
):

    for row in range(result_matrix.shape[0]):
        SWIPE_spin_to_ecliptic(
            result=result_matrix[row, :],
            sun_earth_angle_rad=sun_earth_angles_rad[row],
            site_colatitude_rad=site_colatitude_rad,
            site_longitude_rad=site_longitude_rad,
            longitude_speed_rad_per_sec=longitude_speed_rad_per_sec,
            spin_rate_hz=spin_rate_hz,
            time_s=time_vector_s[row],
        )


class SwipeScanningStrategy(ScanningStrategy):
    """A class containing the parameters of the sky scanning strategy 
    for SWIPE

    """

    def __init__(
        self,
        site_latitude_deg=78.2232,
        site_longitude_deg=15.6267,
        longitude_speed_deg_per_sec=0,
        spin_rate_rmp=1.0,
        start_time=astropy.time.Time("2023-01-01", scale="tdb"),
    ):

        self.site_colatitude_rad = np.deg2rad(90.0 - site_latitude_deg)
        self.site_longitude_rad = np.deg2rad(site_longitude_deg)
        self.longitude_speed_rad_per_sec = np.deg2rad(longitude_speed_deg_per_sec)
        self.spin_rate_hz = spin_rate_rmp / 60.0
        self.start_time = start_time

    def __repr__(self):
        return (
            "SwipeScanningStrategy(site_colatitude_rad={site_colatitude_rad}, "
            "site_longitude_rad={site_longitude_rad},"
            "longitude_speed_rad_per_sec={longitude_speed_rad_per_sec}, "
            "spin_rate_hz={spin_rate_hz}, "
            "start_time={start_time})".format(
                site_colatitude_rad=self.site_colatitude_rad,
                site_longitude_rad=self.site_longitude_rad,
                longitude_speed_rad_per_sec=self.longitude_speed_rad_per_sec,
                spin_rate_hz=self.spin_rate_hz,
                start_time=self.start_time,
            )
        )

    def all_spin_to_ecliptic(
        self,
        result_matrix,
        sun_earth_angles_rad,
        site_colatitude_rad,
        site_longitude_rad,
        longitude_speed_rad_per_sec,
        spin_rate_hz,
        time_vector_s,
    ):
        assert result_matrix.shape == (len(time_vector_s), 4)
        assert len(sun_earth_angles_rad) == len(time_vector_s)

        SWIPE_all_spin_to_ecliptic(
            result_matrix=result_matrix,
            sun_earth_angles_rad=sun_earth_angles_rad,
            site_colatitude_rad=self.site_colatitude_rad,
            site_longitude_rad=self.site_longitude_rad,
            longitude_speed_rad_per_sec=self.longitude_speed_rad_per_sec,
            spin_rate_hz=self.spin_rate_hz,
            time_vector_s=time_vector_s,
        )

    def generate_spin2ecl_quaternions(
        self,
        start_time: Union[float, astropy.time.Time],
        time_span_s: float,
        delta_time_s: float,
    ) -> Spin2EclipticQuaternions:

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

        sun_earth_angles_rad = calculate_sun_earth_angles_rad(time)

        self.all_spin_to_ecliptic(
            result_matrix=spin2ecliptic_quats,
            sun_earth_angles_rad=sun_earth_angles_rad,
            site_colatitude_rad=self.site_colatitude_rad,
            site_longitude_rad=self.site_longitude_rad,
            longitude_speed_rad_per_sec=self.longitude_speed_rad_per_sec,
            spin_rate_hz=self.spin_rate_hz,            
            time_vector_s=time_s,
        )

        return Spin2EclipticQuaternions(
            start_time=start_time,
            pointing_freq_hz=pointing_freq_hz,
            quats=spin2ecliptic_quats,
        )
