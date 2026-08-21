"""
Tests for reading RINEX navigation files.
Verifies that headers and ephemeris blocks are read correctly for
RINEX v2 (GPS), v3 (multi-constellation), and v4 formats.
"""
import sys
import os
import pytest
import numpy as np

project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(os.path.join(project_path, 'src'))

from gnssmultipath.readers.RinexNav import (
    RinexNav, RinexNavData,
    _read_header, _extract_v4_blocks,
    _ORBIT_LINES, _RINEX4_BODY_LEN, _SUPPORTED_V4_MESSAGES,
)


class TestRinexNavHeaderDetection:
    """Test that the header reader detects RINEX version correctly."""

    def test_v3_gps_header(self, nav_gps_file):
        version, header = _read_header(nav_gps_file)
        assert version == 3
        assert len(header) > 0

    def test_v3_galileo_header(self, nav_galileo_file):
        version, header = _read_header(nav_galileo_file)
        assert version == 3
        assert len(header) > 0

    def test_v3_glonass_header(self, nav_glonass_file):
        version, header = _read_header(nav_glonass_file)
        assert version == 3
        assert len(header) > 0

    def test_v3_beidou_header(self, nav_beidou_file):
        version, header = _read_header(nav_beidou_file)
        assert version == 3
        assert len(header) > 0

    def test_v3_mixed_header(self, nav_mixed_file):
        version, header = _read_header(nav_mixed_file)
        assert version == 3
        assert len(header) > 0

    def test_v4_mixed_header(self, nav_v4_mixed_file):
        version, header = _read_header(nav_v4_mixed_file)
        assert version == 4
        assert len(header) > 0

    def test_v2_gps_header(self, nav_v2_gps_file):
        version, header = _read_header(nav_v2_gps_file)
        assert version == 2
        assert len(header) > 0

    def test_v2_glonass_header(self, nav_v2_glonass_file):
        version, header = _read_header(nav_v2_glonass_file)
        assert version == 2
        header_text = "\n".join(header)
        assert "GLONASS NAV DATA" in header_text

    def test_header_contains_end_of_header(self, nav_gps_file):
        _, header = _read_header(nav_gps_file)
        header_text = "\n".join(header)
        assert "END OF HEADER" in header_text



class TestRinexV3GPSNav:
    """Test reading RINEX v3 GPS navigation file."""

    @pytest.fixture(autouse=True, scope="class")
    @classmethod
    def _read_nav(cls, nav_gps_file):
        cls.result = RinexNav.read_nav(nav_gps_file, desired_GNSS=["G"])

    def test_result_has_ephemerides(self):
        assert hasattr(self.result, "ephemerides")

    def test_result_has_header(self):
        assert hasattr(self.result, "header")

    def test_result_has_nepochs(self):
        assert hasattr(self.result, "nepochs")
        assert self.result.nepochs > 0

    def test_ephemerides_not_empty(self):
        eph = self.result.ephemerides
        assert len(eph) > 0

    def test_ephemerides_shape(self):
        eph = self.result.ephemerides
        # Should have 36 columns (PRN + 35 parameters)
        assert eph.shape[1] == 36

    def test_all_satellites_are_gps(self):
        eph = self.result.ephemerides
        prns = eph[:, 0]
        for prn in prns:
            assert str(prn).startswith("G"), f"Expected GPS PRN, got {prn}"

    def test_first_gps_satellite(self):
        eph = self.result.ephemerides
        first_prn = str(eph[0, 0])
        assert first_prn.startswith("G")

    def test_header_has_leap_seconds_info(self):
        header_text = "\n".join(self.result.header)
        assert "LEAP SECONDS" in header_text


class TestRinexV3GalileoNav:
    """Test reading RINEX v3 Galileo navigation file."""

    @pytest.fixture(autouse=True, scope="class")
    @classmethod
    def _read_nav(cls, nav_galileo_file):
        cls.result = RinexNav.read_nav(nav_galileo_file, desired_GNSS=["E"])

    def test_ephemerides_not_empty(self):
        assert len(self.result.ephemerides) > 0

    def test_all_satellites_are_galileo(self):
        eph = self.result.ephemerides
        prns = eph[:, 0]
        for prn in prns:
            assert str(prn).startswith("E"), f"Expected Galileo PRN, got {prn}"

    def test_ephemerides_shape(self):
        assert self.result.ephemerides.shape[1] == 36

    def test_first_satellite(self):
        first_prn = str(self.result.ephemerides[0, 0])
        assert "E" in first_prn



class TestRinexV3GLONASSNav:
    """Test reading RINEX v3 GLONASS navigation file."""

    @pytest.fixture(autouse=True, scope="class")
    @classmethod
    def _read_nav(cls, nav_glonass_file):
        cls.result = RinexNav.read_nav(nav_glonass_file, desired_GNSS=["R"])

    def test_ephemerides_not_empty(self):
        assert len(self.result.ephemerides) > 0

    def test_all_satellites_are_glonass(self):
        eph = self.result.ephemerides
        prns = eph[:, 0]
        for prn in prns:
            assert str(prn).startswith("R"), f"Expected GLONASS PRN, got {prn}"

    def test_ephemerides_shape(self):
        assert self.result.ephemerides.shape[1] == 36

    def test_glonass_fcn(self):
        assert hasattr(self.result, "glonass_fcn")
        if self.result.glonass_fcn is not None:
            assert len(self.result.glonass_fcn) > 0



class TestRinexV3BeiDouNav:
    """Test reading RINEX v3 BeiDou navigation file."""

    @pytest.fixture(autouse=True, scope="class")
    @classmethod
    def _read_nav(cls, nav_beidou_file):
        cls.result = RinexNav.read_nav(nav_beidou_file, desired_GNSS=["C"])

    def test_ephemerides_not_empty(self):
        assert len(self.result.ephemerides) > 0

    def test_all_satellites_are_beidou(self):
        eph = self.result.ephemerides
        prns = eph[:, 0]
        for prn in prns:
            assert str(prn).startswith("C"), f"Expected BeiDou PRN, got {prn}"

    def test_ephemerides_shape(self):
        assert self.result.ephemerides.shape[1] == 36



class TestRinexV3MixedNav:
    """Test reading RINEX v3 mixed navigation file with all systems."""

    @pytest.fixture(autouse=True, scope="class")
    @classmethod
    def _read_nav(cls, nav_mixed_file):
        cls.result = RinexNav.read_nav(
            nav_mixed_file,
            desired_GNSS=["G", "R", "E", "C"],
        )

    def test_ephemerides_not_empty(self):
        assert len(self.result.ephemerides) > 0

    def test_has_multiple_systems(self):
        eph = self.result.ephemerides
        prns = eph[:, 0]
        systems = set(str(p)[0] for p in prns)
        assert len(systems) >= 2, f"Expected multi-system, got {systems}"

    def test_has_gps(self):
        eph = self.result.ephemerides
        prns = [str(p) for p in eph[:, 0]]
        assert any(p.startswith("G") for p in prns)

    def test_has_galileo(self):
        eph = self.result.ephemerides
        prns = [str(p) for p in eph[:, 0]]
        assert any(p.startswith("E") for p in prns)



class TestRinexV3NavFiltering:
    """Test that filtering by system works correctly."""

    def test_filter_gps_only_from_mixed(self, nav_mixed_file):
        result = RinexNav.read_nav(nav_mixed_file, desired_GNSS=["G"])
        eph = result.ephemerides
        prns = [str(p) for p in eph[:, 0]]
        for prn in prns:
            assert prn.startswith("G"), f"Expected GPS only, got {prn}"

    def test_filter_galileo_only_from_mixed(self, nav_mixed_file):
        result = RinexNav.read_nav(nav_mixed_file, desired_GNSS=["E"])
        eph = result.ephemerides
        prns = [str(p) for p in eph[:, 0]]
        for prn in prns:
            assert prn.startswith("E"), f"Expected Galileo only, got {prn}"


class TestRinexV4NavHeaderDetection:
    """Test that the RINEX v4 header is correctly identified."""

    def test_v4_version_detected(self, nav_v4_mixed_file):
        version, header = _read_header(nav_v4_mixed_file)
        assert version == 4

    def test_v4_header_contains_end_of_header(self, nav_v4_mixed_file):
        _, header = _read_header(nav_v4_mixed_file)
        assert any("END OF HEADER" in line for line in header)

    def test_v4_header_contains_leap_seconds(self, nav_v4_mixed_file):
        _, header = _read_header(nav_v4_mixed_file)
        assert any("LEAP SECONDS" in line for line in header)


class TestRinexV4NavMessageFiltering:
    """Test the _extract_v4_blocks function directly.

    Uses BRD400DLR_S_20230710000_01D_MN_rin_v4.rnx which contains:
        EPH records: LNAV(428 GPS), CNAV(334 GPS), FDMA(1240 GLO),
                     INAV(2943 GAL), FNAV(2897 GAL), D1(893 BDS),
                     D2(168 BDS), CNV1(634 BDS), CNV2(1024 BDS),
                     SBAS(9072), QZSS LNAV, IRNSS
        Non-EPH records: STO(126), ION(130), EOP(32)
    Supported families: GPS LNAV, GLONASS FDMA, Galileo INAV/FNAV/IFNV, BeiDou D1/D2/D1D2.
    """

    @pytest.fixture(scope="class")
    @classmethod
    def v4_lines(cls, nav_v4_mixed_file):
        with open(nav_v4_mixed_file) as f:
            return f.readlines()

    def test_gps_lnav_block_count(self, v4_lines):
        blocks = _extract_v4_blocks(v4_lines, {"G"})
        assert len(blocks) == 428

    def test_glonass_fdma_block_count(self, v4_lines):
        blocks = _extract_v4_blocks(v4_lines, {"R"})
        assert len(blocks) == 1240

    def test_galileo_inav_fnav_block_count(self, v4_lines):
        blocks = _extract_v4_blocks(v4_lines, {"E"})
        assert len(blocks) == 2943 + 2897  # INAV + FNAV

    def test_beidou_d1_d2_block_count(self, v4_lines):
        blocks = _extract_v4_blocks(v4_lines, {"C"})
        assert len(blocks) == 893 + 168  # D1 + D2

    def test_all_systems_block_count(self, v4_lines):
        blocks = _extract_v4_blocks(v4_lines, {"G", "R", "E", "C"})
        expected = 428 + 1240 + 2943 + 2897 + 893 + 168
        assert len(blocks) == expected

    def test_gps_cnav_is_skipped(self, v4_lines):
        """GPS CNAV messages must not appear."""
        blocks = _extract_v4_blocks(v4_lines, {"G"})
        for block in blocks:
            assert block[0].strip().startswith("G"), "Non-GPS block leaked through"
        assert len(blocks) == 428

    def test_beidou_cnv1_cnv2_are_skipped(self, v4_lines):
        """BeiDou CNV1/CNV2 messages must not appear."""
        blocks = _extract_v4_blocks(v4_lines, {"C"})
        assert len(blocks) == 893 + 168

    def test_sbas_qzss_irnss_are_skipped(self, v4_lines):
        """SBAS, QZSS, and IRNSS are not in supported systems."""
        blocks = _extract_v4_blocks(v4_lines, {"G", "R", "E", "C"})
        for block in blocks:
            sys_char = block[0].strip()[0]
            assert sys_char in {"G", "R", "E", "C"}, f"Unexpected system: {sys_char}"

    def test_non_eph_records_are_skipped(self, v4_lines):
        """STO, ION, EOP records must not appear in output."""
        blocks = _extract_v4_blocks(v4_lines, {"G", "R", "E", "C"})
        for block in blocks:
            for line in block:
                assert not line.startswith("> STO")
                assert not line.startswith("> ION")
                assert not line.startswith("> EOP")

    def test_gps_block_body_length_is_8(self, v4_lines):
        blocks = _extract_v4_blocks(v4_lines, {"G"})
        for i, block in enumerate(blocks):
            assert len(block) == 8, f"GPS block {i} has {len(block)} lines, expected 8"

    def test_glonass_block_body_length_is_5(self, v4_lines):
        blocks = _extract_v4_blocks(v4_lines, {"R"})
        for i, block in enumerate(blocks):
            assert len(block) == 5, f"GLONASS block {i} has {len(block)} lines, expected 5"

    def test_galileo_block_body_length_is_8(self, v4_lines):
        blocks = _extract_v4_blocks(v4_lines, {"E"})
        for i, block in enumerate(blocks):
            assert len(block) == 8, f"Galileo block {i} has {len(block)} lines, expected 8"

    def test_beidou_block_body_length_is_8(self, v4_lines):
        blocks = _extract_v4_blocks(v4_lines, {"C"})
        for i, block in enumerate(blocks):
            assert len(block) == 8, f"BeiDou block {i} has {len(block)} lines, expected 8"

    def test_empty_for_unsupported_system_only(self, v4_lines):
        blocks = _extract_v4_blocks(v4_lines, {"S"})
        assert len(blocks) == 0


class TestRinexV4NavFullRead:
    """Test the full RinexNav.read_nav pipeline on the RINEX v4 file."""

    @pytest.fixture(autouse=True, scope="class")
    @classmethod
    def _read_all(cls, nav_v4_mixed_file):
        cls.result_all = RinexNav.read_nav(
            nav_v4_mixed_file, desired_GNSS=["G", "R", "E", "C"], data_rate=30
        )

    def test_result_keys(self):
        assert hasattr(self.result_all, "ephemerides")
        assert hasattr(self.result_all, "header")
        assert hasattr(self.result_all, "nepochs")
        assert hasattr(self.result_all, "glonass_fcn")

    def test_ephemerides_shape_36_columns(self):
        assert self.result_all.ephemerides.shape[1] == 36

    def test_all_four_systems_present(self):
        eph = self.result_all.ephemerides
        systems = {str(prn)[0] for prn in eph[:, 0]}
        assert systems == {"G", "R", "E", "C"}

    def test_total_ephemerides_count(self):
        """After data_rate=30 filtering, total count should be 3737."""
        assert len(self.result_all.ephemerides) == 3737

    def test_nepochs_matches_ephemerides_length(self):
        assert self.result_all.nepochs == len(self.result_all.ephemerides)


class TestRinexV4NavPerSystemRead:
    """Test per-system reads on a RINEX v4 file match expected counts and content."""

    def test_gps_only_count_and_purity(self, nav_v4_mixed_file):
        result = RinexNav.read_nav(nav_v4_mixed_file, desired_GNSS=["G"], data_rate=30)
        eph = result.ephemerides
        assert len(eph) == 393
        assert eph.shape[1] == 36
        assert all(str(p).startswith("G") for p in eph[:, 0])

    def test_gps_unique_prns(self, nav_v4_mixed_file):
        result = RinexNav.read_nav(nav_v4_mixed_file, desired_GNSS=["G"], data_rate=30)
        unique_prns = set(str(p) for p in result.ephemerides[:, 0])
        assert len(unique_prns) == 32

    def test_glonass_only_count_and_purity(self, nav_v4_mixed_file):
        result = RinexNav.read_nav(nav_v4_mixed_file, desired_GNSS=["R"], data_rate=30)
        eph = result.ephemerides
        assert len(eph) == 1240
        assert eph.shape[1] == 36
        assert all(str(p).startswith("R") for p in eph[:, 0])

    def test_glonass_fcn_populated(self, nav_v4_mixed_file):
        result = RinexNav.read_nav(nav_v4_mixed_file, desired_GNSS=["R"], data_rate=30)
        fcn = result.glonass_fcn
        assert fcn is not None
        assert len(fcn) == 26
        assert fcn[1] == 1    # R01 -> FCN +1
        assert fcn[2] == -4   # R02 -> FCN -4
        assert fcn[3] == 5    # R03 -> FCN +5

    def test_galileo_only_count_and_purity(self, nav_v4_mixed_file):
        result = RinexNav.read_nav(nav_v4_mixed_file, desired_GNSS=["E"], data_rate=30)
        eph = result.ephemerides
        assert len(eph) == 1049
        assert eph.shape[1] == 36
        assert all(str(p).startswith("E") for p in eph[:, 0])

    def test_galileo_unique_prns(self, nav_v4_mixed_file):
        result = RinexNav.read_nav(nav_v4_mixed_file, desired_GNSS=["E"], data_rate=30)
        unique_prns = set(str(p) for p in result.ephemerides[:, 0])
        assert len(unique_prns) == 26

    def test_beidou_only_count_and_purity(self, nav_v4_mixed_file):
        result = RinexNav.read_nav(nav_v4_mixed_file, desired_GNSS=["C"], data_rate=30)
        eph = result.ephemerides
        assert len(eph) == 1055
        assert eph.shape[1] == 36
        assert all(str(p).startswith("C") for p in eph[:, 0])

    def test_beidou_unique_prns(self, nav_v4_mixed_file):
        result = RinexNav.read_nav(nav_v4_mixed_file, desired_GNSS=["C"], data_rate=30)
        unique_prns = set(str(p) for p in result.ephemerides[:, 0])
        assert len(unique_prns) == 44


class TestRinexV4NavEphemerisValues:
    """Verify parsed ephemeris values against known values from the raw file."""

    @pytest.fixture(autouse=True, scope="class")
    @classmethod
    def _read_all(cls, nav_v4_mixed_file):
        cls.result = RinexNav.read_nav(
            nav_v4_mixed_file, desired_GNSS=["G", "R", "E", "C"], data_rate=30
        )

    def _first_entry_for(self, sys_char):
        eph = self.result.ephemerides
        for row in eph:
            if str(row[0]).startswith(sys_char):
                return row
        pytest.fail(f"No {sys_char} entry found")

    def test_gps_g01_first_epoch(self):
        """First G01 LNAV block: 2023 03 12 00 00 00, af0=2.037500962615e-04."""
        row = self._first_entry_for("G")
        assert str(row[0]) == "G01"
        assert float(row[1]) == 2023
        assert float(row[2]) == 3
        assert float(row[3]) == 12
        assert float(row[4]) == 0
        assert float(row[5]) == 0
        assert float(row[6]) == 0
        assert abs(float(row[7]) - 2.037500962615e-04) < 1e-16

    def test_glonass_r01_first_epoch(self):
        """First R01 FDMA block: 2023 03 12 00 15 00, af0=2.458319067955e-05."""
        row = self._first_entry_for("R")
        assert str(row[0]) == "R01"
        assert float(row[1]) == 2023
        assert float(row[2]) == 3
        assert float(row[3]) == 12
        assert abs(float(row[7]) - 2.458319067955e-05) < 1e-17

    def test_galileo_e01_first_epoch(self):
        """First E01 INAV/FNAV block: 2023 03 12 00 00 00, af0=-1.709302887321e-05."""
        row = self._first_entry_for("E")
        assert str(row[0]) == "E01"
        assert float(row[1]) == 2023
        assert abs(float(row[7]) - (-1.709302887321e-05)) < 1e-17

    def test_beidou_c01_first_epoch(self):
        """First C01 D1/D2 block: 2023 03 12 00 00 00, af0=9.050882654265e-04."""
        row = self._first_entry_for("C")
        assert str(row[0]) == "C01"
        assert float(row[1]) == 2023
        assert abs(float(row[7]) - 9.050882654265e-04) < 1e-16

    def test_gps_g01_sqrtA(self):
        """G01 LNAV sqrtA should be 5.153656053543e+03 (column index 17)."""
        row = self._first_entry_for("G")
        assert abs(float(row[17]) - 5.153656053543e+03) < 1e-6


class TestRinexV4NavDataRateFiltering:
    """Test that data_rate filtering works correctly on v4 parsed data."""

    def test_higher_data_rate_gives_fewer_ephemerides(self, nav_v4_mixed_file):
        result_30 = RinexNav.read_nav(nav_v4_mixed_file, desired_GNSS=["G"], data_rate=30)
        result_120 = RinexNav.read_nav(nav_v4_mixed_file, desired_GNSS=["G"], data_rate=120)
        assert len(result_120.ephemerides) < len(result_30.ephemerides)

    def test_data_rate_0_returns_all_blocks(self, nav_v4_mixed_file):
        result = RinexNav.read_nav(nav_v4_mixed_file, desired_GNSS=["G"], data_rate=0)
        assert len(result.ephemerides) >= 428


class TestRinexV2GPSNav:
    """Test reading RINEX v2 GPS navigation file."""

    @pytest.fixture(autouse=True, scope="class")
    @classmethod
    def _read_nav(cls, nav_v2_gps_file):
        cls.result = RinexNav.read_nav(nav_v2_gps_file)

    def test_result_has_ephemerides(self):
        assert hasattr(self.result, "ephemerides")

    def test_result_has_header(self):
        assert hasattr(self.result, "header")

    def test_ephemerides_not_empty(self):
        eph = self.result.ephemerides
        assert len(eph) > 0

    def test_ephemerides_shape(self):
        eph = self.result.ephemerides
        assert eph.shape[1] == 36

    def test_all_satellites_are_gps(self):
        eph = self.result.ephemerides
        prns = eph[:, 0]
        for prn in prns:
            prn_str = str(prn)
            assert prn_str.startswith("G") or prn_str.replace(".", "").isdigit(), \
                f"Unexpected PRN format: {prn}"

    def test_header_contains_version(self):
        header_text = "\n".join(self.result.header)
        assert "2.11" in header_text or "2" in header_text

    def test_nepochs(self):
        assert self.result.nepochs > 0

    def test_v2_dataframe_output(self, nav_v2_gps_file):
        result = RinexNav.read_nav(nav_v2_gps_file, dataframe=True)
        import pandas as pd
        assert isinstance(result.ephemerides, pd.DataFrame)
        assert len(result.ephemerides) > 0


class TestRinexV2GLONASSNav:
    """Test reading a RINEX v2.10 GLONASS navigation file."""

    @pytest.fixture(autouse=True, scope="class")
    @classmethod
    def _read_nav(cls, nav_v2_glonass_file):
        cls.result = RinexNav.read_nav(nav_v2_glonass_file)

    def test_result_has_ephemerides(self):
        assert hasattr(self.result, "ephemerides")

    def test_ephemerides_not_empty(self):
        assert len(self.result.ephemerides) > 0

    def test_ephemerides_shape(self):
        # GLONASS rows are zero-padded out to the unified 36-column layout.
        assert self.result.ephemerides.shape[1] == 36

    def test_all_satellites_are_glonass(self):
        prns = self.result.ephemerides[:, 0]
        for prn in prns:
            assert str(prn).startswith("R"), f"Expected GLONASS PRN, got {prn}"

    def test_prn_two_digit_format(self):
        # All PRNs should be R + 2-digit slot number (e.g. R01, R24).
        for prn in self.result.ephemerides[:, 0]:
            prn_str = str(prn)
            assert len(prn_str) == 3, f"Expected 3-char PRN, got {prn_str!r}"
            assert prn_str[1:].isdigit()

    def test_glonass_fcn_populated(self):
        assert hasattr(self.result, "glonass_fcn")
        assert self.result.glonass_fcn is not None
        assert len(self.result.glonass_fcn) > 0
        # Slot numbers must be ints; FCN values must be ints in [-7, 6].
        for slot, fcn in self.result.glonass_fcn.items():
            assert isinstance(slot, (int, np.integer))
            assert isinstance(fcn, (int, np.integer))
            assert -7 <= int(fcn) <= 6

    def test_desired_gnss_filter_includes_glonass(self, nav_v2_glonass_file):
        # ``desired_GNSS=['R']`` must keep all rows; ``['G']`` must drop them.
        keep = RinexNav.read_nav(nav_v2_glonass_file, desired_GNSS=["R"])
        drop = RinexNav.read_nav(nav_v2_glonass_file, desired_GNSS=["G"])
        assert len(keep.ephemerides) == len(self.result.ephemerides)
        assert len(drop.ephemerides) == 0

    def test_v2_glonass_dataframe_output(self, nav_v2_glonass_file):
        import pandas as pd
        result = RinexNav.read_nav(nav_v2_glonass_file, dataframe=True)
        assert isinstance(result.ephemerides, pd.DataFrame)
        assert len(result.ephemerides) > 0


class TestRinexV3NavDataFrame:
    """Test that RINEX v3 reader can return a DataFrame."""

    def test_dataframe_output(self, nav_gps_file):
        result = RinexNav.read_nav(nav_gps_file, desired_GNSS=["G"], dataframe=True)
        import pandas as pd
        assert isinstance(result.ephemerides, pd.DataFrame)
        assert len(result.ephemerides) > 0
