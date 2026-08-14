"""
Tests for gnssmultipath.constants — carrier frequencies and wavelengths.
"""
import sys
import os
import pytest
import numpy as np

project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(project_path, 'src'))

from gnssmultipath.constants import (
    CARRIER_FREQUENCIES,
    EARTH_ROTATION_RATE,
    GLONASS_FDMA,
    GM,
    SPEED_OF_LIGHT,
    SYSTEM_BANDS,
    band_label,
    build_frequency_overview,
    carrier_frequency,
    earth_gravitational_constant,
    earth_rotation_rate,
    wavelength,
)
from gnssmultipath.readers.GNSSObservationData import SystemObservations


# Legacy table that used to live in GNSS_MultipathAnalysis.py; the new
# constants module must reproduce it exactly (row index = band - 1).
_LEGACY_FREQ_OVERVIEW = {
    'G': np.array([[1.57542e+09], [1.22760e+09], [np.nan], [np.nan], [1.17645e+09],
                   [np.nan], [np.nan], [np.nan], [np.nan]]),
    'R': np.array([
        [1.602000e+09, 5.625000e+05], [1.246000e+09, 4.375000e+05],
        [1.202025e+09, 0.0], [1.600995e+09, 0.0], [np.nan, 0.0],
        [1.248060e+09, 0.0], [np.nan, 0.0], [np.nan, 0.0], [np.nan, 0.0],
    ]),
    'E': np.array([[1.575420e+09], [np.nan], [np.nan], [np.nan], [1.176450e+09],
                   [1.278750e+09], [1.207140e+09], [1.191795e+09], [np.nan]]),
    'C': np.array([[1.575420e+09], [1.561098e+09], [np.nan], [np.nan], [1.176450e+09],
                   [1.268520e+09], [1.207140e+09], [1.191795e+09], [np.nan]]),
}


def _legacy_glonass_expansion(slot2channel, max_glo_id=36):
    """Reproduce the original (buggy) GLONASS expansion from GNSS_MultipathAnalysis.

    The loop ran ``for j in range(max_GLO_ID)``, i.e. slots 0..35, so the
    highest GLONASS slot (36) was never assigned a frequency.
    """
    base = _LEGACY_FREQ_OVERVIEW['R']
    table = np.full((9, max_glo_id + 1), np.nan)
    for band_row in range(9):
        for slot in range(max_glo_id):          # off-by-one: excludes slot 36
            if slot in slot2channel:
                table[band_row, slot] = (base[band_row, 0]
                                         + slot2channel[slot] * base[band_row, 1])
    return table


class TestSpeedOfLight:

    def test_value(self):
        assert SPEED_OF_LIGHT == 299792458.0


class TestCarrierFrequency:

    @pytest.mark.parametrize("system, band, expected", [
        ('G', 1, 1575.42e6),
        ('G', 2, 1227.60e6),
        ('G', 5, 1176.45e6),
        ('E', 1, 1575.42e6),
        ('E', 5, 1176.45e6),
        ('E', 6, 1278.75e6),
        ('E', 7, 1207.14e6),
        ('E', 8, 1191.795e6),
        ('C', 1, 1575.42e6),
        ('C', 2, 1561.098e6),
        ('C', 5, 1176.45e6),
        ('C', 6, 1268.52e6),
        ('C', 7, 1207.14e6),
        ('C', 8, 1191.795e6),
        ('R', 3, 1202.025e6),
        ('R', 4, 1600.995e6),
        ('R', 6, 1248.06e6),
    ])
    def test_known_frequencies(self, system, band, expected):
        assert carrier_frequency(system, band) == pytest.approx(expected)

    def test_accepts_full_system_name(self):
        assert carrier_frequency('GPS', 1) == carrier_frequency('G', 1)

    def test_accepts_band_as_string(self):
        assert carrier_frequency('G', '1') == carrier_frequency('G', 1)

    @pytest.mark.parametrize("channel, expected", [
        (0, 1602.0e6),
        (1, 1602.0e6 + 562500.0),
        (-4, 1602.0e6 - 4 * 562500.0),
        (6, 1602.0e6 + 6 * 562500.0),
    ])
    def test_glonass_g1_fdma(self, channel, expected):
        assert carrier_frequency('R', 1, channel) == pytest.approx(expected)

    @pytest.mark.parametrize("channel, expected", [
        (0, 1246.0e6),
        (-7, 1246.0e6 - 7 * 437500.0),
    ])
    def test_glonass_g2_fdma(self, channel, expected):
        assert carrier_frequency('R', 2, channel) == pytest.approx(expected)

    def test_glonass_fdma_without_channel_raises(self):
        with pytest.raises(ValueError, match="FDMA"):
            carrier_frequency('R', 1)

    def test_glonass_cdma_ignores_channel(self):
        assert carrier_frequency('R', 3, 5) == carrier_frequency('R', 3)

    def test_unknown_system_raises(self):
        with pytest.raises(KeyError, match="Unknown GNSS system"):
            carrier_frequency('J', 1)

    def test_unsupported_band_raises(self):
        with pytest.raises(KeyError, match="not defined for GPS"):
            carrier_frequency('G', 7)


class TestWavelength:

    def test_gps_l1(self):
        assert wavelength('G', 1) == pytest.approx(0.190293672798, abs=1e-11)

    def test_matches_c_over_f(self):
        assert wavelength('E', 5) == pytest.approx(SPEED_OF_LIGHT / 1176.45e6)

    def test_glonass_channel_dependent(self):
        assert wavelength('R', 1, 0) != wavelength('R', 1, 1)


class TestBandLabel:

    def test_gps_l1(self):
        assert band_label('G', 1) == 'L1 (1575.42 MHz)'

    def test_galileo_e5b(self):
        assert band_label('E', 7) == 'E5b (1207.14 MHz)'

    def test_beidou_b1i_keeps_precision(self):
        assert band_label('C', 2) == 'B1I (1561.098 MHz)'

    def test_glonass_is_channel_aware(self):
        assert band_label('R', 1).startswith('G1 (1602 + k x')


class TestSystemBands:

    def test_all_declared_bands_have_a_frequency(self):
        for system, bands in SYSTEM_BANDS.items():
            for band in bands:
                assert band in CARRIER_FREQUENCIES[system]


class TestBuildFrequencyOverview:

    def test_non_glonass_shape_and_indexing(self):
        overview = build_frequency_overview({1: 'G', 2: 'E'})
        assert overview[1].shape == (9, 1)
        assert overview[1][0, 0] == pytest.approx(1575.42e6)
        assert overview[1][4, 0] == pytest.approx(1176.45e6)
        assert np.isnan(overview[1][6, 0])

    def test_glonass_expanded_per_slot(self):
        slot2channel = {1: 1, 2: -4, 36: 5}
        overview = build_frequency_overview({1: 'R'}, slot2channel)
        table = overview[1]
        assert table.shape == (9, 37)
        assert table[0, 1] == pytest.approx(1602.0e6 + 562500.0)
        assert table[0, 2] == pytest.approx(1602.0e6 - 4 * 562500.0)
        assert table[1, 2] == pytest.approx(1246.0e6 - 4 * 437500.0)
        # CDMA bands do not depend on the channel number
        assert table[2, 1] == pytest.approx(1202.025e6)
        # Slots without observations stay NaN
        assert np.isnan(table[0, 5])

    def test_glonass_without_channels_raises(self):
        with pytest.raises(ValueError, match="GLONASS k-numbers"):
            build_frequency_overview({1: 'R'})

    @pytest.mark.parametrize("system", ['G', 'E', 'C'])
    def test_parity_with_legacy_table(self, system):
        overview = build_frequency_overview({1: system})
        np.testing.assert_allclose(
            overview[1], _LEGACY_FREQ_OVERVIEW[system], equal_nan=True)

    def test_parity_with_legacy_glonass_expansion(self):
        """Slots 0-35 must match the previous implementation exactly."""
        slot2channel = {1: 1, 2: -4, 3: 5, 7: 5, 24: 2}
        legacy = _legacy_glonass_expansion(slot2channel)
        overview = build_frequency_overview({1: 'R'}, slot2channel)
        np.testing.assert_allclose(overview[1], legacy, equal_nan=True)


class TestGlonassFdmaTable:

    def test_step_sizes(self):
        assert GLONASS_FDMA[1] == (1602.0e6, 562500.0)
        assert GLONASS_FDMA[2] == (1246.0e6, 437500.0)


class TestEarthModelParameters:
    """Each constellation's ICD defines its own GM and Earth rotation rate.

    Reference values: Table 3.4, "Physical parameters of GNSS almanac and
    ephemeris models".
    """

    @pytest.mark.parametrize("system, expected", [
        ('C', 398600.4418e9),
        ('E', 398600.4418e9),
        ('R', 398600.4418e9),
        ('G', 398600.5e9),
        ('J', 398600.5e9),
    ])
    def test_gm(self, system, expected):
        assert GM[system] == expected
        assert earth_gravitational_constant(system) == expected

    @pytest.mark.parametrize("system, expected", [
        ('C', 7.2921150e-5),
        ('E', 7.2921151467e-5),
        ('R', 7.292115e-5),
        ('G', 7.2921151467e-5),
        ('J', 7.2921151467e-5),
    ])
    def test_earth_rotation_rate(self, system, expected):
        assert EARTH_ROTATION_RATE[system] == expected
        assert earth_rotation_rate(system) == expected

    def test_gps_and_galileo_share_the_rotation_rate(self):
        assert earth_rotation_rate('G') == earth_rotation_rate('E')

    def test_beidou_rotation_rate_differs_from_gps(self):
        assert earth_rotation_rate('C') != earth_rotation_rate('G')

    def test_gps_gm_differs_from_the_others(self):
        assert earth_gravitational_constant('G') != earth_gravitational_constant('E')

    def test_accepts_full_system_name(self):
        assert earth_gravitational_constant('BeiDou') == GM['C']
        assert earth_rotation_rate('GLONASS') == EARTH_ROTATION_RATE['R']

    def test_unknown_system_falls_back_to_gps(self):
        assert earth_gravitational_constant('X') == GM['G']
        assert earth_rotation_rate('X') == EARTH_ROTATION_RATE['G']

    def test_beidou_rotation_rate_error_is_metre_level(self):
        """Using the GPS rate for BeiDou shifts the orbit along-track.

        The -omega*toe term makes the error grow across the BDT week, so it
        reaches tens of metres by the end of the week for a MEO satellite.
        """
        d_omega = earth_rotation_rate('G') - earth_rotation_rate('C')
        r_meo = 27.9064e6
        end_of_week = 604800.0
        assert r_meo * d_omega * end_of_week == pytest.approx(24.8, abs=0.5)

    def test_gm_error_is_decimetre_level_over_an_hour(self):
        """Using the GPS GM for Galileo shifts the orbit along-track."""
        a = 29.5998e6
        n0 = np.sqrt(earth_gravitational_constant('E') / a**3)
        d_n0 = 0.5 * (GM['G'] - GM['E']) / GM['G'] * n0
        assert a * d_n0 * 3600 == pytest.approx(0.96, abs=0.05)


class TestEarthModelParametersInUse:
    """The orbit propagator must pick the constants up from the system code."""

    def test_kepler2ecef_uses_the_system_specific_constants(self):
        import inspect
        from gnssmultipath.SatelliteEphemerisToECEF import Kepler2ECEF

        source = inspect.getsource(Kepler2ECEF.kepler2ecef)
        assert 'earth_gravitational_constant(gnss_sys)' in source
        assert 'earth_rotation_rate(gnss_sys)' in source
        assert '3.986005e14' not in source
        assert '7.2921151467e-5' not in source

    def test_glonass_equations_use_the_glonass_constants(self):
        import inspect
        from gnssmultipath.SatelliteEphemerisToECEF import GLOStateVec2ECEF

        source = inspect.getsource(GLOStateVec2ECEF.glonass_diff_eq)
        assert "earth_gravitational_constant('R')" in source
        assert "earth_rotation_rate('R')" in source



class TestGlonassSlot36Regression:
    """GLONASS slot 36 must receive a carrier frequency.

    The previous inline implementation looped ``range(max_GLO_ID)`` and
    therefore left the whole column for slot 36 as NaN, which silently
    excluded R36 from any dual-frequency processing.
    """

    SLOT2CHANNEL = {1: 1, 2: -4, 8: 6, 24: 2, 36: 5}

    def test_legacy_implementation_drops_slot_36(self):
        """Guard the premise of this regression test."""
        legacy = _legacy_glonass_expansion(self.SLOT2CHANNEL)
        assert np.all(np.isnan(legacy[:, 36]))

    def test_slot_36_g1_frequency(self):
        table = build_frequency_overview({1: 'R'}, self.SLOT2CHANNEL)[1]
        assert table[0, 36] == pytest.approx(1602.0e6 + 5 * 562500.0)

    def test_slot_36_g2_frequency(self):
        table = build_frequency_overview({1: 'R'}, self.SLOT2CHANNEL)[1]
        assert table[1, 36] == pytest.approx(1246.0e6 + 5 * 437500.0)

    def test_slot_36_cdma_bands_are_channel_independent(self):
        table = build_frequency_overview({1: 'R'}, self.SLOT2CHANNEL)[1]
        assert table[2, 36] == pytest.approx(1202.025e6)   # G3
        assert table[3, 36] == pytest.approx(1600.995e6)   # G1a
        assert table[5, 36] == pytest.approx(1248.060e6)   # G2a

    def test_only_slot_36_differs_from_legacy(self):
        """Every other slot must be bit-identical to the legacy behaviour."""
        legacy = _legacy_glonass_expansion(self.SLOT2CHANNEL)
        table = build_frequency_overview({1: 'R'}, self.SLOT2CHANNEL)[1]
        np.testing.assert_allclose(table[:, :36], legacy[:, :36], equal_nan=True)
        assert np.any(np.isfinite(table[:, 36]))

    def test_slot_36_reaches_the_system_accessor(self):
        """The observation accessor must expose a frequency for R36 too."""
        sys_obs = SystemObservations(
            {1: np.zeros((37, 1))}, ['C1C'], system_code='R',
            glo_slot2channel=self.SLOT2CHANNEL,
        )
        assert sys_obs.frequency('C1C', prn=36) == pytest.approx(
            1602.0e6 + 5 * 562500.0)
        assert np.isfinite(sys_obs.frequency('C1C')[36])

    def test_slots_above_the_limit_are_ignored(self):
        table = build_frequency_overview({1: 'R'}, {36: 5, 40: 1})[1]
        assert table.shape == (9, 37)
        assert np.isfinite(table[0, 36])

    def test_custom_max_slot_is_respected(self):
        table = build_frequency_overview({1: 'R'}, {24: 2}, max_glo_id=24)[1]
        assert table.shape == (9, 25)
        assert table[0, 24] == pytest.approx(1602.0e6 + 2 * 562500.0)

