"""
Tests for SatelliteEphemerisToECEF — broadcast ephemeris to ECEF conversion.
"""
import sys
import os
import pytest
import numpy as np

project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(project_path, 'src'))

from gnssmultipath import SatelliteEphemerisToECEF
from gnssmultipath.Geodetic_functions import date2gpstime_vectorized

RINEX_NAV = os.path.join(project_path, "TestData", "NavigationFiles", "v3",
                         "BRDC00IGS_R_20220010000_01D_MN.rnx")

X_REC, Y_REC, Z_REC = 3149785.9652, 598260.8822, 5495348.4927

GREGORIAN_EPOCH = np.array([[2022, 1, 1, 2, 40, 0.0]])

# Approximate orbit radius per constellation, used as a sanity band [m]
ORBIT_RADIUS_RANGE = {
    'G': (25.0e6, 28.0e6),
    'R': (24.0e6, 27.0e6),
    'E': (28.0e6, 31.0e6),
    'C': (25.0e6, 44.0e6),   # BeiDou has MEO, IGSO and GEO satellites
}


@pytest.fixture(scope="module")
def tow_epoch():
    _, tow = date2gpstime_vectorized(GREGORIAN_EPOCH)
    return np.asarray(tow, dtype=float).ravel()


def _converter(systems):
    return SatelliteEphemerisToECEF(RINEX_NAV, X_REC, Y_REC, Z_REC, systems, data_rate=60)


def _first_position(result, system):
    return np.ravel(list(result[system]['position'].values())[0])[:3]


class TestGregorianTimeFormat:
    """`time_fmt='GREGORIAN'` must be equivalent to passing time-of-week.

    Regression test: the Gregorian branch used to build an ``[N, 2]``
    ``[week, tow]`` array, which crashed the GLONASS Runge-Kutta integrator
    with a broadcasting error and silently produced wrong Keplerian orbits
    because the week number was differenced against ``toe``.
    """

    @pytest.mark.parametrize("prn", ["G20", "E11", "C11", "R10"])
    def test_matches_tow_input(self, prn, tow_epoch):
        by_tow = _converter(prn[0]).get_sat_ecef_coordinates(tow_epoch, PRN=prn)
        by_greg = _converter(prn[0]).get_sat_ecef_coordinates(
            GREGORIAN_EPOCH, time_fmt="GREGORIAN", PRN=prn)

        pos_tow = np.array([np.ravel(by_tow[i])[0] for i in range(3)])
        pos_greg = np.array([np.ravel(by_greg[i])[0] for i in range(3)])
        np.testing.assert_allclose(pos_greg, pos_tow, rtol=0, atol=1e-9)

    def test_glonass_does_not_raise(self):
        """This used to raise ValueError from the Runge-Kutta integrator."""
        result = _converter(["G", "R", "E", "C"]).get_sat_ecef_coordinates(
            GREGORIAN_EPOCH, time_fmt="GREGORIAN")
        assert set(result) == {'G', 'R', 'E', 'C'}

    @pytest.mark.parametrize("system", ["G", "R", "E", "C"])
    def test_orbit_radius_is_physical(self, system):
        result = _converter(["G", "R", "E", "C"]).get_sat_ecef_coordinates(
            GREGORIAN_EPOCH, time_fmt="GREGORIAN")
        radius = np.linalg.norm(_first_position(result, system))
        low, high = ORBIT_RADIUS_RANGE[system]
        assert low < radius < high

    def test_accepts_a_flat_gregorian_row(self, tow_epoch):
        by_greg = _converter("G").get_sat_ecef_coordinates(
            GREGORIAN_EPOCH.ravel(), time_fmt="GREGORIAN", PRN="G20")
        by_tow = _converter("G").get_sat_ecef_coordinates(tow_epoch, PRN="G20")
        np.testing.assert_allclose(np.ravel(by_greg[0]), np.ravel(by_tow[0]), atol=1e-9)

    def test_multiple_epochs(self):
        epochs = np.array([[2022, 1, 1, 2, 40, 0.0],
                           [2022, 1, 1, 2, 50, 0.0]])
        xs, _, _, _ = _converter("G").get_sat_ecef_coordinates(
            epochs, time_fmt="GREGORIAN", PRN="G20")
        assert np.ravel(xs).size == 2


class TestPerSystemEarthConstants:
    """The propagator must use each constellation's own ICD constants."""

    def test_beidou_moves_when_the_rotation_rate_changes(self, tow_epoch):
        from gnssmultipath import constants

        correct = _converter("C").get_sat_ecef_coordinates(tow_epoch, PRN="C11")
        pos_correct = np.array([np.ravel(correct[i])[0] for i in range(3)])

        original = constants.EARTH_ROTATION_RATE['C']
        constants.EARTH_ROTATION_RATE['C'] = constants.EARTH_ROTATION_RATE['G']
        try:
            wrong = _converter("C").get_sat_ecef_coordinates(tow_epoch, PRN="C11")
        finally:
            constants.EARTH_ROTATION_RATE['C'] = original
        pos_wrong = np.array([np.ravel(wrong[i])[0] for i in range(3)])

        # Using the GPS rotation rate shifts BeiDou along-track by metres
        assert np.linalg.norm(pos_correct - pos_wrong) > 1.0

    def test_gps_is_unaffected_by_the_beidou_constants(self, tow_epoch):
        from gnssmultipath import constants

        first = _converter("G").get_sat_ecef_coordinates(tow_epoch, PRN="G20")
        pos_first = np.array([np.ravel(first[i])[0] for i in range(3)])

        original = constants.EARTH_ROTATION_RATE['C']
        constants.EARTH_ROTATION_RATE['C'] = 1.0
        try:
            second = _converter("G").get_sat_ecef_coordinates(tow_epoch, PRN="G20")
        finally:
            constants.EARTH_ROTATION_RATE['C'] = original
        pos_second = np.array([np.ravel(second[i])[0] for i in range(3)])

        np.testing.assert_allclose(pos_second, pos_first, atol=1e-9)


if __name__ == '__main__':
    pytest.main()
