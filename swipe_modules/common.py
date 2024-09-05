import pkg_resources
import numpy as np
from numba import njit

data_directory = pkg_resources.resource_filename("swipe_modules", "data")


@njit
def _ct_jd_to_lst_rad(
    time_jd,
    longitude_rad,
):
    """Computes the local sidereal time function of the Julian Day
    Adapted from the astrolib routing ct2lst.pro
    """
    t0 = time_jd - 2451545.0
    t = t0 / 36525.0

    theta = 280.46061837 + 360.98564736629 * t0 + t**2 * (0.000387933 - t / 38710000)

    lst = np.deg2rad(theta) + longitude_rad
    if lst < 0:
        lst = 2 * np.pi + np.mod(lst, 2 * np.pi)
    else:
        lst = np.mod(lst, 2 * np.pi)

    return lst


@njit
def _equinox_precession_rad(
    time_jd,
):
    """Computes the precession of the equinox wrt J2000
    very rough estimation but sufficient for the SWIPE beam, error ~ 2 arcmin
    From https://en.wikipedia.org/wiki/Axial_precession#Values
    https://syrte.obspm.fr/iau2006/aa03_412_P03.pdf
    """
    t = (time_jd - 2451545.0) / 36525.0

    return -(0.02438029195 * t + 5.3592991461e-6 * t * t + 5.5608129223e-9 * t * t * t)


@njit
def _equator_ecliptic_angle_rad(
    time_jd,
):
    """
    Computes the obliquity of the ecliptic wrt J2000
    https://en.wikipedia.org/wiki/Ecliptic#Obliquity_of_the_ecliptic
    up to third order
    See also https://syrte.obspm.fr/iau2006/aa03_412_P03.pdf

    """

    t = (time_jd - 2451545.0) / 36525.0

    return (
        0.4090925998
        - 0.0002270710639 * t
        - 8.8769385011e-10 * t * t
        + 9.7127572873e-9 * t * t * t
    )
