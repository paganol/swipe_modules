import pkg_resources
import numpy as np
from numba import njit

data_directory = pkg_resources.resource_filename("swipe_modules", "data")

EQUATOR_ECLIPTIC_ANGLE_RAD = -0.4090449142  # -23.4365472133 deg in radians

@njit
def _ct_jd_to_lst_rad(
    time_jd,
    longitude_rad,
):
    t0 = time_jd - 2451545.0
    t = t0 / 36525.0

    theta = 280.46061837 + 360.98564736629 * t0 + t**2 * (0.000387933 - t / 38710000)

    lst = np.deg2rad(theta) + longitude_rad
    if lst < 0:
        lst = 2 * np.pi + np.mod(lst, 2 * np.pi)
    else:
        lst = np.mod(lst, 2 * np.pi)

    return lst

