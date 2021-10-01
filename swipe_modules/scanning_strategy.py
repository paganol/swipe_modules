from litebird_sim import ScanningStrategy, calculate_sun_earth_angles_rad
from numpy import np

EQUATOR_ECLIPTIC_ANGLE_RAD = 0.408407045 #23.4 deg in radians

@njit
def spin_to_ecliptic(
    result,
    site_colatitude_rad,
    site_longitude_rad,
    longitude_speed_rad_per_sec,    
    spin_rate_hz,
    time_s,
):

    result[:] = quat_rotation_z(2 * np.pi * spin_rate_hz * time_s)
    quat_left_multiply(result, *quat_rotation_y(site_colatitude_rad))
    quat_left_multiply(
        result, *quat_rotation_x(longitude_speed_rad_per_sec * time_s + site_longitude_rad)
    )
    quat_left_multiply(result, *quat_rotation_z(EQUATOR_ECLIPTIC_ANGLE_RAD))

@njit
def all_spin_to_ecliptic(
    result_matrix,
    site_colatitude_rad,
    site_longitude_rad,
    longitude_speed_rad_per_sec,    
    spin_rate_hz,
    time_vector_s,
):

    for row in range(result_matrix.shape[0]):
        spin_to_ecliptic(
            result=result_matrix[row, :],
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
        spin_sun_angle_rad,
        spin_rate_hz,
        start_time=astropy.time.Time("2023-01-01", scale="tdb"),
    ):

        self.spin_sun_angle_rad = spin_sun_angle_rad
        self.precession_rate_hz = precession_rate_hz
        self.spin_rate_hz = spin_rate_hz
        self.start_time = start_time

    def __repr__(self):
        return (
            "SwipeScanningStrategy(spin_sun_angle_rad={spin_sun_angle_rad}, "
            "precession_rate_hz={precession_rate_hz}, "
            "spin_rate_hz={spin_rate_hz}, "
            "start_time={start_time})".format(
                spin_sun_angle_rad=self.spin_sun_angle_rad,
                precession_rate_hz=self.precession_rate_hz,
                spin_rate_hz=self.spin_rate_hz,
                start_time=self.start_time,
            )
        )

    def __str__(self):
        return """SWIPE scanning strategy:
    angle between the Sun and the spin axis:       {spin_sun_angle_deg:.1f}°
    rotations around the precession angle:         {precession_rate_hr} rot/hr
    rotations around the spinning axis:            {spin_rate_hr} rot/hr
    start time of the simulation:                  {start_time}""".format(
            spin_sun_angle_deg=np.rad2deg(self.spin_sun_angle_rad),
            precession_rate_hr=3600.0 * self.precession_rate_hz,
            spin_rate_hr=3600.0 * self.spin_rate_hz,
            start_time=self.start_time,
        )

    def all_spin_to_ecliptic(self, result_matrix, sun_earth_angles_rad, time_vector_s):
        assert result_matrix.shape == (len(time_vector_s), 4)
        assert len(sun_earth_angles_rad) == len(time_vector_s)

        SWIPE_all_spin_to_ecliptic(
            result_matrix=result_matrix,
            sun_earth_angles_rad=sun_earth_angles_rad,
            spin_sun_angle_rad=self.spin_sun_angle_rad,
            precession_rate_hz=self.precession_rate_hz,
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
            time_vector_s=time_s,
        )

        return Spin2EclipticQuaternions(
            start_time=start_time,
            pointing_freq_hz=pointing_freq_hz,
            quats=spin2ecliptic_quats,
        )




