"""
Tests for ObsCode — the parsed RINEX observation code (signal) metadata.
"""
import sys
import os
import pytest

project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(os.path.join(project_path, 'src'))

from gnssmultipath.readers.ObsCode import ObsCode
from gnssmultipath.constants import SPEED_OF_LIGHT


class TestParsingRinex3:

    def test_fields(self):
        sig = ObsCode.parse('L5Q', 'E')
        assert sig.code == 'L5Q'
        assert sig.system == 'E'
        assert sig.obs_type == 'L'
        assert sig.band == 5
        assert sig.attribute == 'Q'

    def test_pseudorange(self):
        sig = ObsCode.parse('C1C', 'G')
        assert sig.is_pseudorange
        assert not sig.is_phase
        assert sig.type_name == 'pseudorange'

    def test_phase(self):
        sig = ObsCode.parse('L2W', 'G')
        assert sig.is_phase
        assert sig.type_name == 'phase'

    def test_doppler(self):
        assert ObsCode.parse('D1C', 'G').is_doppler

    def test_snr(self):
        sig = ObsCode.parse('S1C', 'G')
        assert sig.is_snr
        assert sig.type_name == 'snr'

    def test_system_name(self):
        assert ObsCode.parse('C2I', 'C').system_name == 'BeiDou'

    def test_band_description_is_system_specific(self):
        assert ObsCode.parse('L1C', 'G').band_description == 'L1 (1575.42 MHz)'
        assert ObsCode.parse('L1X', 'E').band_description == 'E1 (1575.42 MHz)'

    def test_strips_whitespace(self):
        assert ObsCode.parse(' C1C ', 'G').code == 'C1C'


class TestParsingRinex2:

    def test_two_character_code(self):
        sig = ObsCode.parse('C1', 'G')
        assert sig.obs_type == 'C'
        assert sig.band == 1
        assert sig.attribute == ''

    def test_p_code_is_pseudorange(self):
        sig = ObsCode.parse('P2', 'G')
        assert sig.is_pseudorange
        assert sig.type_name == 'pseudorange'
        assert sig.band == 2

    def test_phase(self):
        sig = ObsCode.parse('L1', 'R')
        assert sig.is_phase
        assert sig.attribute == ''


class TestFrequency:

    def test_gps_l1(self):
        assert ObsCode.parse('L1C', 'G').frequency() == pytest.approx(1575.42e6)

    def test_gps_l2(self):
        assert ObsCode.parse('C2W', 'G').frequency() == pytest.approx(1227.60e6)

    def test_galileo_e5b(self):
        assert ObsCode.parse('L7X', 'E').frequency() == pytest.approx(1207.14e6)

    def test_beidou_b1i(self):
        assert ObsCode.parse('C2I', 'C').frequency() == pytest.approx(1561.098e6)

    def test_wavelength(self):
        sig = ObsCode.parse('L1C', 'G')
        assert sig.wavelength() == pytest.approx(SPEED_OF_LIGHT / 1575.42e6)

    def test_glonass_needs_channel(self):
        with pytest.raises(ValueError, match="FDMA"):
            ObsCode.parse('C1C', 'R').frequency()

    def test_glonass_with_channel(self):
        sig = ObsCode.parse('C1C', 'R')
        assert sig.frequency(glonass_channel=-4) == pytest.approx(1602.0e6 - 4 * 562500.0)

    def test_rinex2_p_code_frequency(self):
        assert ObsCode.parse('P2', 'G').frequency() == pytest.approx(1227.60e6)


class TestErrors:

    @pytest.mark.parametrize("bad", ['', 'C', 'C1CC', 'CCCC'])
    def test_invalid_length_raises(self, bad):
        with pytest.raises(ValueError, match="not a valid RINEX observation code"):
            ObsCode.parse(bad, 'G')

    def test_unknown_type_raises(self):
        with pytest.raises(ValueError, match="unknown observation type"):
            ObsCode.parse('Z1C', 'G')

    def test_non_numeric_band_raises(self):
        with pytest.raises(ValueError, match="non-numeric frequency band"):
            ObsCode.parse('CXC', 'G')

    def test_non_string_raises(self):
        with pytest.raises(ValueError, match="must be a string"):
            ObsCode.parse(123, 'G')


class TestStringInterop:

    def test_str_returns_code(self):
        assert str(ObsCode.parse('C1C', 'G')) == 'C1C'

    def test_equals_string(self):
        assert ObsCode.parse('C1C', 'G') == 'C1C'
        assert ObsCode.parse('C1C', 'G') != 'L1C'

    def test_equals_other_obscode(self):
        assert ObsCode.parse('C1C', 'G') == ObsCode.parse('C1C', 'G')
        assert ObsCode.parse('C1C', 'G') != ObsCode.parse('C1C', 'E')

    def test_hashable_and_deduplicates(self):
        codes = {ObsCode.parse('C1C', 'G'), ObsCode.parse('C1C', 'G')}
        assert len(codes) == 1

    def test_is_immutable(self):
        sig = ObsCode.parse('C1C', 'G')
        with pytest.raises(Exception):
            sig.band = 2

    def test_repr(self):
        text = repr(ObsCode.parse('L5Q', 'E'))
        assert 'L5Q' in text and 'phase' in text and 'band=5' in text
