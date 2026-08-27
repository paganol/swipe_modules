"""Common utility functions for SWIPE scanning strategies."""

from importlib.resources import files
from typing import Final

import numpy as np
from numba import njit

__all__ = [
    "data_directory",
]

data_directory: Final[str] = str(files("swipe_modules") / "data")

# Constants for astronomical calculations
_J2000_JD: Final[float] = 2451545.0
_SECONDS_PER_CENTURY: Final[float] = 36525.0


@njit
def _ct_jd_to_lst_rad(
    time_jd: float,
    longitude_rad: float,
) -> float:
    """Compute the local sidereal time as a function of Julian Day.

    Adapted from the astrolib routine ct2lst.pro.

    Args:
        time_jd: Time in Julian Days.
        longitude_rad: Longitude in radians.

    Returns:
        Local sidereal time in radians.
    """
    t0 = time_jd - _J2000_JD
    t = t0 / _SECONDS_PER_CENTURY

    theta = 280.46061837 + 360.98564736629 * t0 + t**2 * (0.000387933 - t / 38710000)

    lst = np.deg2rad(theta) + longitude_rad
    if lst < 0:
        lst = 2 * np.pi + np.mod(lst, 2 * np.pi)
    else:
        lst = np.mod(lst, 2 * np.pi)

    return lst


@njit
def _equinox_precession_rad(time_jd: float) -> float:
    """Compute the precession of the equinox with respect to J2000.

    Very rough estimation but sufficient for the SWIPE beam (~2 arcmin error).

    References:
        - https://en.wikipedia.org/wiki/Axial_precession#Values
        - https://syrte.obspm.fr/iau2006/aa03_412_P03.pdf

    Args:
        time_jd: Time in Julian Days.

    Returns:
        Precession angle in radians.
    """
    t = (time_jd - _J2000_JD) / _SECONDS_PER_CENTURY

    return -(0.02438029195 * t + 5.3592991461e-6 * t * t + 5.5608129223e-9 * t * t * t)


@njit
def _equator_ecliptic_angle_rad(time_jd: float) -> float:
    """Compute the obliquity of the ecliptic with respect to J2000.

    Third-order polynomial approximation.

    References:
        - https://en.wikipedia.org/wiki/Ecliptic#Obliquity_of_the_ecliptic
        - https://syrte.obspm.fr/iau2006/aa03_412_P03.pdf

    Args:
        time_jd: Time in Julian Days.

    Returns:
        Obliquity angle in radians.
    """

    t = (time_jd - _J2000_JD) / _SECONDS_PER_CENTURY

    return (
        0.4090925998
        - 0.0002270710639 * t
        - 8.8769385011e-10 * t * t
        + 9.7127572873e-9 * t * t * t
    )
