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
    t = (time_jd - 2451545.0) / 36525.0

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

    t = (time_jd - 2451545.0) / 36525.0

    return (
        0.4090925998
        - 0.0002270710639 * t
        - 8.8769385011e-10 * t * t
        + 9.7127572873e-9 * t * t * t
    )


@njit
def _chiA_rad(time_jd: float) -> float:
    """Compute chi_A precession component in radians.

    Args:
        time_jd: Time in Julian Days.

    Returns:
        chi_A angle in radians.
    """
    t = (time_jd - _J2000_JD) / _SECONDS_PER_CENTURY

    return (
        5.1160448512764897e-05 * t
        + 1.1541668417966057e-05 * t * t
        + 5.454153912482279e-09 * t * t * t
    )


@njit
def _zitaA_rad(time_jd: float) -> float:
    """Compute zeta_A precession component in radians.

    Args:
        time_jd: Time in Julian Days.

    Returns:
        zeta_A angle in radians.
    """
    t = (time_jd - _J2000_JD) / _SECONDS_PER_CENTURY

    return (
        1.2593605507709181e-05
        + 0.01118019594596964 * t
        + 1.4636573514065e-06 * t * t
        + 8.710308038918257e-08 * t * t * t
        + 1.5853407372281827e-10 * t * t * t * t
        + 9.696273622190719e-13 * t * t * t * t * t
    )

@njit
def _zetaA_rad(time_jd: float) -> float:
    """Compute zeta_A precession component in radians.

    Args:
        time_jd: Time in Julian Days.

    Returns:
        zeta_A angle in radians.
    """
    t = (time_jd - _J2000_JD) / _SECONDS_PER_CENTURY

    return (
        1.2593605507709181e-05
        + 0.011180192901339722 * t
        + 5.307638369914167e-06 * t * t
        + 8.836844409687845e-08 * t * t * t
        + 2.278624301214819e-10 * t * t * t * t
        + 1.4544410433286078e-12 * t * t * t * t * t
    )

@njit
def _thetaA_rad(time_jd: float) -> float:
    """Compute theta_A precession component in radians.

    Args:
        time_jd: Time in Julian Days.

    Returns:
        theta_A angle in radians.
    """
    t = (time_jd - _J2000_JD) / _SECONDS_PER_CENTURY

    return (
        + 0.0097165957880331 * t
        + 2.0698407438860407e-06 * t * t
        + 2.027738069377445e-07 * t * t * t
        + 2.913730223468311e-10 * t * t * t * t
        + 4.848136811095359e-13 * t * t * t * t * t
    )
