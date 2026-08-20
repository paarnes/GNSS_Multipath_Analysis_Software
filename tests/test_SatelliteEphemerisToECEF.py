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


class TestDatetimeInput:
    """Datetime-like epochs must be accepted directly, e.g. ``RinexObsData.datetimes``."""

    DATETIMES = np.array(['2022-01-01T02:40:00', '2022-01-01T02:40:30'], dtype='datetime64[ns]')

    def test_matches_tow_input(self):
        _, tow = date2gpstime_vectorized(np.array([[2022, 1, 1, 2, 40, 0.0],
                                                   [2022, 1, 1, 2, 40, 30.0]]))
        by_tow = _converter(["G"]).get_sat_ecef_coordinates(np.asarray(tow, dtype=float))
        by_datetime = _converter(["G"]).get_sat_ecef_coordinates(self.DATETIMES)

        np.testing.assert_allclose(by_datetime['G']['position']['20'],
                                   by_tow['G']['position']['20'], atol=1e-9)

    def test_timestamps_are_kept(self):
        converter = _converter(["G"])
        converter.get_sat_ecef_coordinates(self.DATETIMES)
        np.testing.assert_array_equal(converter.epoch_timestamps, self.DATETIMES)

    def test_time_of_week_input_still_gets_timestamps(self, tow_epoch):
        converter = _converter(["G"])
        converter.get_sat_ecef_coordinates(tow_epoch)
        # The week is not in the input, so it is taken from the ephemerides.
        assert converter.epoch_timestamps[0] == np.datetime64('2022-01-01T02:40:00', 'ns')


class TestSystemTimeScales:
    """BeiDou ephemerides are referenced to BDT, which trails GPS time by 14 s."""

    @staticmethod
    def _ephemeris_module():
        # 'gnssmultipath.SatelliteEphemerisToECEF' resolves to the class on the package,
        # so the module itself has to be picked up from sys.modules.
        return sys.modules['gnssmultipath.SatelliteEphemerisToECEF']

    def test_beidou_uses_bdt(self, tow_epoch, monkeypatch):
        """Skipping the BDT offset moves a BeiDou MEO ~40 km along-track."""
        correct = _converter("C").get_sat_ecef_coordinates(tow_epoch, PRN="C11")
        pos_correct = np.array([np.ravel(correct[i])[0] for i in range(3)])

        monkeypatch.setattr(self._ephemeris_module(), 'BDT_GPST_OFFSET', 0.0)
        uncorrected = _converter("C").get_sat_ecef_coordinates(tow_epoch, PRN="C11")
        pos_uncorrected = np.array([np.ravel(uncorrected[i])[0] for i in range(3)])

        assert np.linalg.norm(pos_correct - pos_uncorrected) > 30_000.0

    def test_gps_is_not_shifted(self, tow_epoch, monkeypatch):
        """Only BeiDou gets the offset, so GPS must be unaffected by it."""
        first = _converter("G").get_sat_ecef_coordinates(tow_epoch, PRN="G20")
        pos_first = np.array([np.ravel(first[i])[0] for i in range(3)])

        monkeypatch.setattr(self._ephemeris_module(), 'BDT_GPST_OFFSET', 0.0)
        second = _converter("G").get_sat_ecef_coordinates(tow_epoch, PRN="G20")
        pos_second = np.array([np.ravel(second[i])[0] for i in range(3)])

        np.testing.assert_allclose(pos_second, pos_first, atol=1e-9)


class TestBeiDouGeoOrbits:
    """C01-C05 use the GEO branch of the BeiDou ICD, not the MEO/IGSO equations."""
    FULL_DAY = np.arange(0, 86400.0, 1800.0) + 172800.0  # tow for 2022-01-01, every 30 min

    @pytest.fixture(scope="class")
    @classmethod
    def beidou(cls):
        return _converter("C").get_sat_ecef_coordinates(cls.FULL_DAY)

    @pytest.mark.parametrize("prn", ["1", "2", "3", "4", "5"])
    def test_geo_orbit_radius(self, beidou, prn):
        radius = np.linalg.norm(beidou['C']['position'][prn], axis=1)
        # Geostationary radius is 42 164 km
        np.testing.assert_allclose(radius.mean(), 42.164e6, rtol=1e-3)

    @pytest.mark.parametrize("prn", ["1", "2", "3", "4", "5"])
    def test_geo_is_almost_stationary_in_ecef(self, beidou, prn):
        """A geostationary satellite barely moves in ECEF; the MEO equations would sweep
        it around the full orbit instead."""
        position = beidou['C']['position'][prn]
        span = np.linalg.norm(position.max(axis=0) - position.min(axis=0))
        assert span < 5.0e6

    def test_igso_and_meo_still_move(self, beidou):
        for prn in ('6', '11'):
            position = beidou['C']['position'][prn]
            span = np.linalg.norm(position.max(axis=0) - position.min(axis=0))
            assert span > 50.0e6


class TestGlonassEarthRotationCorrection:
    """GLONASS must end up in the same reference frame as the Keplerian systems."""

    def test_rotation_is_applied_when_the_receiver_is_known(self, tow_epoch):
        from gnssmultipath.SatelliteEphemerisToECEF import GLOStateVec2ECEF

        converter = _converter("R")
        ephemerides = converter.get_closest_ephemerides_for_PRN_at_time("R10", tow_epoch)

        rotated, _, _, _ = GLOStateVec2ECEF(X_REC, Y_REC, Z_REC).interpolate_glonass_coord_runge_kutta(
            ephemerides.copy(), tow_epoch)
        unrotated, _, _, _ = GLOStateVec2ECEF().interpolate_glonass_coord_runge_kutta(
            ephemerides.copy(), tow_epoch)

        offset = np.linalg.norm(rotated - unrotated)
        # omega_e * travel_time * orbit_radius is roughly 160 m for a GLONASS orbit
        assert 50.0 < offset < 400.0

    def test_rotation_preserves_the_orbit_radius(self, tow_epoch):
        from gnssmultipath.SatelliteEphemerisToECEF import GLOStateVec2ECEF

        converter = _converter("R")
        ephemerides = converter.get_closest_ephemerides_for_PRN_at_time("R10", tow_epoch)

        rotated, _, _, _ = GLOStateVec2ECEF(X_REC, Y_REC, Z_REC).interpolate_glonass_coord_runge_kutta(
            ephemerides.copy(), tow_epoch)
        unrotated, _, _, _ = GLOStateVec2ECEF().interpolate_glonass_coord_runge_kutta(
            ephemerides.copy(), tow_epoch)

        # A rotation about the Z-axis changes neither the radius nor the Z component
        np.testing.assert_allclose(np.linalg.norm(rotated, axis=1),
                                   np.linalg.norm(unrotated, axis=1), rtol=1e-12)
        np.testing.assert_allclose(rotated[:, 2], unrotated[:, 2], rtol=1e-12)


class TestDataFrameOutput:
    """``output_format='pd.DataFrame'`` gives a (timestamp, system, SV) indexed table."""

    @pytest.fixture(scope="class")
    @classmethod
    def frame(cls):
        converter = _converter(["G", "E"])
        return converter.get_sat_ecef_coordinates(
            np.array(['2022-01-01T02:40:00', '2022-01-01T02:40:30'], dtype='datetime64[ns]'),
            output_format="pd.DataFrame")

    def test_index_and_columns(self, frame):
        assert list(frame.index.names) == ['timestamp', 'system', 'SV']
        assert list(frame.columns) == ['X', 'Y', 'Z']

    def test_satellite_ids_are_zero_padded(self, frame):
        svs = frame.index.get_level_values('SV').unique()
        assert all(len(sv) == 3 and sv[0] in ('G', 'E') for sv in svs)
        assert 'G20' in svs

    def test_single_satellite_selection(self, frame):
        single = frame.xs(('G', 'G20'), level=('system', 'SV'))
        assert len(single) == 2
        assert list(single.columns) == ['X', 'Y', 'Z']

    def test_values_match_the_dict_output(self, frame):
        converter = _converter(["G", "E"])
        as_dict = converter.get_sat_ecef_coordinates(
            np.array(['2022-01-01T02:40:00', '2022-01-01T02:40:30'], dtype='datetime64[ns]'))
        np.testing.assert_allclose(frame.xs(('G', 'G20'), level=('system', 'SV')).to_numpy(),
                                   as_dict['G']['position']['20'], atol=1e-9)

    def test_to_dataframe_requires_computed_coordinates(self):
        with pytest.raises(ValueError):
            _converter(["G"]).to_dataframe()


if __name__ == '__main__':
    pytest.main()
