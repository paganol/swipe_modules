"""Cross-check the hand-rolled astronomy formulas in common.py against erfa.

erfa (bundled with astropy as a dependency) is the reference IAU SOFA
implementation. These tests pin the internal approximations to it so a
regression like the precession double-counting bug (lst + eqx, off by
19-33 arcmin for 2023-2040) cannot silently reappear.
"""

import erfa
import numpy as np
import pytest
from astropy.time import Time

from swipe_modules.common import _ct_jd_to_lst_rad, _equator_ecliptic_angle_rad

ARCSEC = 206264.80624709636  # radians -> arcsec
SWIPE_BEAM_TOLERANCE_ARCSEC = 120.0  # ~2 arcmin, per common.py's own docstrings

DATES = ["2000-01-01", "2020-01-01", "2026-08-27", "2040-01-01", "2050-01-01"]


@pytest.mark.parametrize("date_str", DATES)
def test_obliquity_matches_erfa(date_str):
    t = Time(date_str, scale="tdb")
    ours = _equator_ecliptic_angle_rad(t.jd)
    reference = erfa.obl06(t.jd1, t.jd2)
    assert abs(ours - reference) * ARCSEC < 1.0


@pytest.mark.parametrize("date_str", DATES)
def test_lst_matches_gmst06(date_str):
    """lst alone (no extra precession term) should track erfa's GMST."""
    t = Time(date_str, scale="utc")
    ours = _ct_jd_to_lst_rad(t.jd, 0.0)
    reference = erfa.gmst06(t.ut1.jd1, t.ut1.jd2, t.tt.jd1, t.tt.jd2)
    diff_arcsec = abs((ours - reference + np.pi) % (2 * np.pi) - np.pi) * ARCSEC
    assert diff_arcsec < SWIPE_BEAM_TOLERANCE_ARCSEC / 10


if __name__ == "__main__":
    for date_str in DATES:
        t = Time(date_str, scale="utc")
        ours = _ct_jd_to_lst_rad(t.jd, 0.0)
        reference = erfa.gmst06(t.ut1.jd1, t.ut1.jd2, t.tt.jd1, t.tt.jd2)
        diff_arcsec = abs((ours - reference + np.pi) % (2 * np.pi) - np.pi) * ARCSEC
        print(f"{date_str}: internal lst vs erfa gmst06 -> {diff_arcsec:.3f} arcsec")
