"""
Tests for reading RINEX observation files.
Verifies that headers and observation blocks are read correctly for
both RINEX 3.04 and RINEX 2.11 formats.
"""
import sys
import os
import pytest
import numpy as np
from numpy.testing import assert_almost_equal

project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(os.path.join(project_path, 'src'))

from gnssmultipath.readers.readRinexObs import (
    readRinexObs,
    readRinexObs304,
    rinexReadObsFileHeader304,
    RinexObsData,
)


class TestRinex304Header:
    """Test header parsing for RINEX 3.04 observation files."""

    @pytest.fixture(autouse=True, scope="class")
    def _read_header(self, rinex304_obs_file):
        """Read header once for the class."""
        cls = type(self)
        (cls.success, cls.rinexVersion, cls.gnssType, cls.markerName,
         cls.recType, cls.antDelta, cls.GNSSsystems, cls.numOfObsCodes,
         cls.obsCodes, cls.obsCodeIndex, cls.tFirstObs, cls.tLastObs,
         cls.tInterval, cls.timeSystem, cls.numHeaderLines,
         cls.clockOffsetsON, cls.rinexProgr, cls.rinexDate,
         cls.leapSec, cls.approxPosition, cls.GLO_Slot2ChannelMap,
         cls.eof, cls.fid
        ) = rinexReadObsFileHeader304(
            rinex304_obs_file,
            includeAllGNSSsystems=1,
            includeAllObsCodes=1,
            desiredGNSSsystems=None,
            desiredObsCodes=None,
            desiredObsBands=None,
        )

    def test_success(self):
        assert self.success == 1

    def test_rinex_version(self):
        assert "3.04" in str(self.rinexVersion).strip()

    def test_gnss_type_mixed(self):
        assert self.gnssType.strip() == "M"

    def test_approx_position(self):
        assert_almost_equal(self.approxPosition[0], 3149785.9652, decimal=2)
        assert_almost_equal(self.approxPosition[1], 598260.8822, decimal=2)
        assert_almost_equal(self.approxPosition[2], 5495348.4927, decimal=2)

    def test_interval(self):
        assert_almost_equal(self.tInterval, 30.0, decimal=1)

    def test_time_of_first_obs(self):
        assert self.tFirstObs[0] == 2022
        assert self.tFirstObs[1] == 1
        assert self.tFirstObs[2] == 1
        assert self.tFirstObs[3] == 0  # hour
        assert self.tFirstObs[4] == 0  # min
        assert_almost_equal(self.tFirstObs[5], 0.0, decimal=3)

    def test_time_of_last_obs(self):
        assert self.tLastObs[0] == 2022
        assert self.tLastObs[1] == 1
        assert self.tLastObs[2] == 1
        assert self.tLastObs[3] == 3   # hour
        assert self.tLastObs[4] == 39  # min
        assert_almost_equal(self.tLastObs[5], 30.0, decimal=3)

    def test_time_system(self):
        assert self.timeSystem.strip() == "GPS"

    def test_receiver_type(self):
        assert "TRIMBLE_NETR9" in self.recType

    def test_pgm_name(self):
        assert "CONVBIN" in self.rinexProgr

    def test_gnss_systems_present(self):
        systems = list(self.GNSSsystems.values())
        assert "G" in systems
        assert "R" in systems
        assert "E" in systems
        assert "C" in systems

    def test_num_gnss_systems(self):
        assert len(self.GNSSsystems) == 4

    def test_obs_codes_gps(self):
        # GPS should have 9 obs types: C1C L1C C1P C2W L2W C2X L2X C5X L5X
        gps_idx = [k for k, v in self.GNSSsystems.items() if v == "G"][0]
        gps_codes = self.obsCodes[gps_idx]["G"]
        assert len(gps_codes) == 9
        assert "C1C" in gps_codes
        assert "L1C" in gps_codes
        assert "L5X" in gps_codes

    def test_obs_codes_galileo(self):
        # Galileo should have 8 obs types: C1X L1X C7X L7X C5X L5X C8X L8X
        gal_idx = [k for k, v in self.GNSSsystems.items() if v == "E"][0]
        gal_codes = self.obsCodes[gal_idx]["E"]
        assert len(gal_codes) == 8
        assert "C1X" in gal_codes
        assert "L8X" in gal_codes

    def test_obs_codes_glonass(self):
        # GLONASS should have 8 obs types: C1C L1C C1P L1P C2P L2P C2C L2C
        glo_idx = [k for k, v in self.GNSSsystems.items() if v == "R"][0]
        glo_codes = self.obsCodes[glo_idx]["R"]
        assert len(glo_codes) == 8
        assert "C1C" in glo_codes
        assert "L2C" in glo_codes

    def test_obs_codes_beidou(self):
        # BeiDou should have 6 obs types: C2X L2X C7X L7X C6X L6X
        bds_idx = [k for k, v in self.GNSSsystems.items() if v == "C"][0]
        bds_codes = self.obsCodes[bds_idx]["C"]
        assert len(bds_codes) == 6
        assert "C2X" in bds_codes
        assert "L6X" in bds_codes

    def test_antenna_delta(self):
        # antDelta is a list of strings from header
        assert_almost_equal(float(self.antDelta[0]), 0.0, decimal=3)
        assert_almost_equal(float(self.antDelta[1]), 0.0, decimal=3)
        assert_almost_equal(float(self.antDelta[2]), 0.0, decimal=3)

    def test_leap_seconds(self):
        # Leap seconds not explicitly in header -> may be NaN
        # If present, should be numeric
        if not np.isnan(self.leapSec):
            assert self.leapSec >= 0

    def test_glonass_slot_channel_map(self):
        # Header has 22 GLONASS slots
        assert len(self.GLO_Slot2ChannelMap) > 0
        # R01 has slot 1
        assert self.GLO_Slot2ChannelMap[1] == 1
        # R02 has slot -4
        assert self.GLO_Slot2ChannelMap[2] == -4



class TestRinex304ObservationBlocks:
    """Test reading observation blocks for RINEX 3.04 file."""

    @pytest.fixture(autouse=True, scope="class")
    def _read_obs(self, rinex304_obs_file):
        """Read the cropped observation file once for the class."""
        cls = type(self)
        rinex_data = readRinexObs(rinex304_obs_file)
        cls.GNSS_obs = rinex_data.GNSS_obs
        cls.GNSS_LLI = rinex_data.GNSS_LLI
        cls.GNSS_SS = rinex_data.GNSS_SS
        cls.GNSS_SVs = rinex_data.GNSS_SVs
        cls.time_epochs = rinex_data.time_epochs
        cls.nepochs = rinex_data.nepochs
        cls.GNSSsystems = rinex_data.GNSSsystems
        cls.obsCodes = rinex_data.obsCodes
        cls.approxPosition = rinex_data.approxPosition
        cls.max_sat = rinex_data.max_sat
        cls.tInterval = rinex_data.tInterval
        cls.markerName = rinex_data.markerName
        cls.rinexVersion = rinex_data.rinexVersion
        cls.recType = rinex_data.recType
        cls.timeSystem = rinex_data.timeSystem
        cls.leapSec = rinex_data.leapSec
        cls.gnssType = rinex_data.gnssType
        cls.rinexProgr = rinex_data.rinexProgr
        cls.rinexDate = rinex_data.rinexDate
        cls.antDelta = rinex_data.antDelta
        cls.tFirstObs = rinex_data.tFirstObs
        cls.tLastObs = rinex_data.tLastObs
        cls.clockOffsetsON = rinex_data.clockOffsetsON
        cls.GLO_Slot2ChannelMap = rinex_data.GLO_Slot2ChannelMap
        cls.success = rinex_data.success

    def test_success(self):
        assert self.success == 1

    def test_nepochs_positive(self):
        assert self.nepochs > 0

    def test_nepochs_consistent_with_time(self):
        # 3h39m30s at 30s interval = 439 epochs
        expected = int((3 * 3600 + 39 * 60 + 30) / 30) + 1
        assert self.nepochs == expected

    def test_time_epochs_shape(self):
        assert self.time_epochs.shape[0] == self.nepochs
        assert self.time_epochs.shape[1] == 2  # GPS week, time of week

    def test_time_epochs_gps_week(self):
        # 2022-01-01 corresponds to GPS week 2190
        assert self.time_epochs[0, 0] == 2190

    def test_gnss_systems(self):
        systems = list(self.GNSSsystems.values())
        assert "G" in systems
        assert "R" in systems
        assert "E" in systems
        assert "C" in systems

    def test_gnss_obs_has_systems(self):
        for sys_code in ["G", "R", "E", "C"]:
            assert sys_code in self.GNSS_obs

    def test_gnss_obs_gps_epoch1_has_data(self):
        # Epochs are 1-indexed. First epoch should have GPS data
        gps_obs = self.GNSS_obs["G"]
        assert 1 in gps_obs  # epoch 1
        epoch1 = gps_obs[1]  # ndarray shape [max_PRN+1, num_obs_types]
        assert epoch1.shape[0] > 0

    def test_gnss_obs_gps_first_epoch_prns(self):
        # Epoch 1 data is an array indexed by PRN. Non-zero rows indicate observed PRNs.
        gps_obs = self.GNSS_obs["G"]
        epoch1 = gps_obs[1]
        expected_prns = [30, 15, 16, 18, 1, 8, 27, 14, 21, 10, 23]
        for prn in expected_prns:
            assert np.any(epoch1[prn] != 0), f"GPS PRN {prn} has no data in epoch 1"

    def test_gnss_obs_values_nonzero(self):
        # Check that specific observation values match expected data
        gps_obs = self.GNSS_obs["G"]
        epoch1 = gps_obs[1]
        # G30 C1C should be ~24850337.312
        # obsCodes from readRinexObs is list-based (not nested dict like header)
        gps_idx = [k for k, v in self.GNSSsystems.items() if v == "G"][0]
        obs_list = self.obsCodes[gps_idx]
        if isinstance(obs_list, dict):
            obs_list = list(obs_list.values())[0]
        c1c_idx = list(obs_list).index("C1C")
        assert epoch1[30, c1c_idx] == pytest.approx(24850337.312, abs=1.0)

    def test_gnss_obs_galileo_first_epoch(self):
        gal_obs = self.GNSS_obs["E"]
        epoch1 = gal_obs[1]
        expected_prns = [31, 3, 26, 1, 8, 14, 13, 33]
        for prn in expected_prns:
            assert np.any(epoch1[prn] != 0), f"Galileo PRN {prn} has no data in epoch 1"

    def test_gnss_obs_glonass_first_epoch(self):
        glo_obs = self.GNSS_obs["R"]
        epoch1 = glo_obs[1]
        expected_prns = [8, 15, 24, 7, 17, 1, 14, 23]
        for prn in expected_prns:
            assert np.any(epoch1[prn] != 0), f"GLONASS PRN {prn} has no data in epoch 1"

    def test_gnss_obs_beidou_first_epoch(self):
        bds_obs = self.GNSS_obs["C"]
        epoch1 = bds_obs[1]
        expected_prns = [26, 13, 6, 16, 27, 5, 29, 30, 9]
        for prn in expected_prns:
            assert np.any(epoch1[prn] != 0), f"BeiDou PRN {prn} has no data in epoch 1"

    def test_gnss_svs_epoch0_satellite_count(self):
        # First epoch should have 36 total satellites across all systems
        total_sats = 0
        for sys_code in ["G", "R", "E", "C"]:
            svs = self.GNSS_SVs[sys_code]
            total_sats += svs[0, 0]  # count at epoch 0
        assert total_sats == 36

    def test_tinterval(self):
        assert_almost_equal(self.tInterval, 30.0, decimal=1)

    def test_approx_position(self):
        assert_almost_equal(self.approxPosition[0], 3149785.9652, decimal=2)

    def test_max_sat_per_system(self):
        # max_sat should be an array with max PRN for each system
        assert len(self.max_sat) == len(self.GNSSsystems)
        for ms in self.max_sat:
            assert ms > 0

    def test_rinex_version_string(self):
        assert "3.04" in str(self.rinexVersion).strip()


class TestRinex304FilteredSystems:
    """Test that filtering GNSS systems works correctly."""

    @pytest.fixture(autouse=True, scope="class")
    def _read_obs_gps_only(self, rinex304_obs_file):
        cls = type(self)
        (cls.GNSS_obs, _, _, cls.GNSS_SVs, _, cls.nepochs,
         cls.GNSSsystems, cls.obsCodes, *_, cls.success
        ) = readRinexObs304(
            rinex304_obs_file,
            includeAllGNSSsystems=0,
            desiredGNSSsystems=["G"],
        )

    def test_success(self):
        assert self.success == 1

    def test_only_gps_system(self):
        systems = list(self.GNSSsystems.values())
        assert "G" in systems
        assert len(systems) == 1

    def test_gps_obs_has_data(self):
        assert "G" in self.GNSS_obs
        assert len(self.GNSS_obs["G"]) > 0



class TestRinex211ObservationFile:
    """Test reading a RINEX 2.11 observation file."""

    @pytest.fixture(autouse=True, scope="class")
    def _read_obs(self, rinex211_obs_file):
        cls = type(self)
        rinex_data = readRinexObs(rinex211_obs_file)
        cls.GNSS_obs = rinex_data.GNSS_obs
        cls.GNSS_LLI = rinex_data.GNSS_LLI
        cls.GNSS_SS = rinex_data.GNSS_SS
        cls.GNSS_SVs = rinex_data.GNSS_SVs
        cls.time_epochs = rinex_data.time_epochs
        cls.nepochs = rinex_data.nepochs
        cls.GNSSsystems = rinex_data.GNSSsystems
        cls.obsCodes = rinex_data.obsCodes
        cls.approxPosition = rinex_data.approxPosition
        cls.max_sat = rinex_data.max_sat
        cls.tInterval = rinex_data.tInterval
        cls.markerName = rinex_data.markerName
        cls.rinexVersion = rinex_data.rinexVersion
        cls.recType = rinex_data.recType
        cls.timeSystem = rinex_data.timeSystem
        cls.leapSec = rinex_data.leapSec
        cls.gnssType = rinex_data.gnssType
        cls.rinexProgr = rinex_data.rinexProgr
        cls.rinexDate = rinex_data.rinexDate
        cls.antDelta = rinex_data.antDelta
        cls.tFirstObs = rinex_data.tFirstObs
        cls.tLastObs = rinex_data.tLastObs
        cls.clockOffsetsON = rinex_data.clockOffsetsON
        cls.GLO_Slot2ChannelMap = rinex_data.GLO_Slot2ChannelMap
        cls.success = rinex_data.success

    def test_success(self):
        assert self.success == 1

    def test_rinex_version(self):
        assert "2" in str(self.rinexVersion).strip().split(".")[0]

    def test_interval(self):
        assert_almost_equal(self.tInterval, 1.0, decimal=1)

    def test_time_of_first_obs(self):
        assert self.tFirstObs[0] == 2020
        assert self.tFirstObs[1] == 11
        assert self.tFirstObs[2] == 5

    def test_time_system(self):
        assert "GPS" in str(self.timeSystem).strip()

    def test_approx_position(self):
        assert_almost_equal(self.approxPosition[0], 3172305.6262, decimal=2)
        assert_almost_equal(self.approxPosition[1], 603531.5081, decimal=2)
        assert_almost_equal(self.approxPosition[2], 5481983.8910, decimal=2)

    def test_marker_name(self):
        assert "GMGD" in str(self.markerName).strip()

    def test_nepochs_positive(self):
        assert self.nepochs > 0

    def test_obs_codes(self):
        # Should have 7 obs types: C1, L1, D1, P1, P2, L2, D2
        # For v2, obsCodes is {sys_idx: {system_code: [codes]}}
        for sys_idx in self.obsCodes:
            for sys_code in self.obsCodes[sys_idx]:
                codes = self.obsCodes[sys_idx][sys_code]
                assert len(codes) == 7
                assert "C1" in codes
                assert "L1" in codes
                assert "L2" in codes
                break  # just check first system
            break

    def test_gnss_obs_has_data(self):
        # Should have at least GPS data
        has_data = False
        for sys_key in self.GNSS_obs:
            if len(self.GNSS_obs[sys_key]) > 0:
                has_data = True
                break
        assert has_data

    def test_time_epochs_shape(self):
        assert self.time_epochs.shape[0] == self.nepochs


class TestRinex210EventFlagsObservationFile:
    """

    The file ``OPEC0010_truncated.22o`` is a truncated RINEX 2.10 mixed
    observation file (GPS + GLONASS) with ``nObsCodes == 10`` (an exact
    multiple of 5) and event-flag blocks (flag 2 with 6 records, flag 5
    with 1 record) embedded between epochs.

    Before the fix, ``rinexReadObsBlock211`` over-read one continuation
    line per satellite when ``nObsCodes`` was a multiple of 5, leaving
    leftover observation lines in the stream. The next call to
    ``rinexReadObsBlockHead211`` then misinterpreted a data line as an
    epoch header, raising
    ``ValueError: invalid literal for int() with base 10: '  '``.
    """

    def test_parses_without_exception(self, rinex210_obs_file_event_flags):
        data = readRinexObs(rinex210_obs_file_event_flags)
        assert data.success == 1

    def test_expected_systems_and_epochs(self, rinex210_obs_file_event_flags):
        data = readRinexObs(rinex210_obs_file_event_flags)
        # Truncated to the first 5 epochs.
        assert data.nepochs == 5
        assert set(data.GNSS_obs.keys()) == {"G", "R"}

    def test_obs_code_count_is_multiple_of_five(
        self, rinex210_obs_file_event_flags
    ):
        # Triggering condition: exactly 10 obs codes per system.
        data = readRinexObs(rinex210_obs_file_event_flags)
        for sys_idx in data.obsCodes:
            for sys_code in data.obsCodes[sys_idx]:
                assert len(data.obsCodes[sys_idx][sys_code]) == 10

