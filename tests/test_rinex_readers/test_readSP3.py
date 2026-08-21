"""
Tests for reading SP3 precise orbit files.
Verifies that both the modern SP3Reader (pandas-based) and the legacy
readSP3Nav (numpy-based) correctly parse headers and position data.
"""
import sys
import os
import pytest
import numpy as np
from numpy.testing import assert_almost_equal

project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(os.path.join(project_path, 'src'))

from gnssmultipath.readers.SP3Reader import SP3Reader
from gnssmultipath.readers.read_SP3Nav import readSP3Nav



class TestSP3Reader2022:
    """Test SP3Reader with 2022 multi-GNSS SP3 file."""

    @pytest.fixture(autouse=True, scope="class")
    @classmethod
    def _read_sp3(cls, sp3_file_2022):
        cls.reader = SP3Reader(filepaths=sp3_file_2022)
        cls.df = cls.reader.read()
        cls.metadata = cls.reader.get_metadata()

    def test_dataframe_not_empty(self):
        assert len(self.df) > 0

    def test_dataframe_columns(self):
        expected_cols = {"Epoch", "Satellite", "X", "Y", "Z", "Clock Bias"}
        assert expected_cols.issubset(set(self.df.columns))

    def test_epoch_interval(self):
        assert_almost_equal(self.metadata["epoch_interval_sec"], 300.0, decimal=1)

    def test_num_epochs(self):
        # n_epochs in metadata comes from SP3 header (GPS week number 2190)
        # Actual unique epochs in data is 289
        unique_epochs = len(self.df["Epoch"].unique())
        assert unique_epochs == 289

    def test_gnss_systems_detected(self):
        systems = self.metadata["gnss_systems_in_file"]
        assert "G" in systems
        assert "R" in systems
        assert "E" in systems
        assert "C" in systems

    def test_satellite_names_format(self):
        sats = self.df["Satellite"].unique()
        for sat in sats:
            assert len(sat) == 3, f"Satellite name should be 3 chars, got '{sat}'"
            assert sat[0] in "GREJC", f"Unexpected system char in '{sat}'"

    def test_gps_satellite_present(self):
        gps_sats = self.df[self.df["Satellite"].str.startswith("G")]
        assert len(gps_sats) > 0

    def test_galileo_satellite_present(self):
        gal_sats = self.df[self.df["Satellite"].str.startswith("E")]
        assert len(gal_sats) > 0

    def test_coordinates_in_meters(self):
        # SP3 stores in km, reader converts to meters
        # Typical GNSS orbit radius ~26000 km = 26e6 m
        g01 = self.df[self.df["Satellite"] == "G01"].iloc[0]
        distance = np.sqrt(g01["X"]**2 + g01["Y"]**2 + g01["Z"]**2)
        # Should be roughly 20,000-30,000 km = 2e7 to 3e7 meters
        assert 2e7 < distance < 3e7, f"Unexpected orbit distance: {distance}"

    def test_first_epoch_gps_g01_position(self):
        # From test data: PG01  13882.271938 -21710.006216   5357.125462
        g01_first = self.df[self.df["Satellite"] == "G01"].iloc[0]
        assert_almost_equal(g01_first["X"], 13882.271938 * 1000, decimal=0)
        assert_almost_equal(g01_first["Y"], -21710.006216 * 1000, decimal=0)
        assert_almost_equal(g01_first["Z"], 5357.125462 * 1000, decimal=0)

    def test_clock_bias_in_seconds(self):
        # Clock biases should be in seconds (converted from microseconds)
        g01_first = self.df[self.df["Satellite"] == "G01"].iloc[0]
        # From data: 469.127768 microseconds = 4.69127768e-4 seconds
        assert_almost_equal(g01_first["Clock Bias"], 469.127768e-6, decimal=9)

    def test_epochs_sorted(self):
        epochs = self.df["Epoch"].values
        assert all(epochs[i] <= epochs[i+1] for i in range(len(epochs)-1))


class TestSP3Reader2020:
    """Test SP3Reader with 2020 multi-GNSS SP3 file."""

    @pytest.fixture(autouse=True, scope="class")
    @classmethod
    def _read_sp3(cls, sp3_file_2020):
        cls.reader = SP3Reader(filepaths=sp3_file_2020)
        cls.df = cls.reader.read()
        cls.metadata = cls.reader.get_metadata()

    def test_dataframe_not_empty(self):
        assert len(self.df) > 0

    def test_epoch_interval(self):
        assert_almost_equal(self.metadata["epoch_interval_sec"], 300.0, decimal=1)

    def test_num_epochs(self):
        unique_epochs = len(self.df["Epoch"].unique())
        assert unique_epochs == 576

    def test_gnss_systems(self):
        systems = self.metadata["gnss_systems_in_file"]
        assert "G" in systems
        assert "R" in systems
        assert "E" in systems
        assert "C" in systems

    def test_first_epoch_g01(self):
        # PG01  22209.803938  14683.385019   2254.797886
        g01_first = self.df[self.df["Satellite"] == "G01"].iloc[0]
        assert_almost_equal(g01_first["X"], 22209.803938 * 1000, decimal=0)
        assert_almost_equal(g01_first["Y"], 14683.385019 * 1000, decimal=0)
        assert_almost_equal(g01_first["Z"], 2254.797886 * 1000, decimal=0)


class TestSP3ReaderFiltering:
    """Test that GNSS system filtering works in SP3Reader."""

    def test_filter_gps_only(self, sp3_file_2022):
        reader = SP3Reader(filepaths=sp3_file_2022, desiredGNSSsystems=["G"])
        df = reader.read()
        sats = df["Satellite"].unique()
        for sat in sats:
            assert sat.startswith("G"), f"Expected GPS only, got {sat}"

    def test_filter_galileo_only(self, sp3_file_2022):
        reader = SP3Reader(filepaths=sp3_file_2022, desiredGNSSsystems=["E"])
        df = reader.read()
        sats = df["Satellite"].unique()
        for sat in sats:
            assert sat.startswith("E"), f"Expected Galileo only, got {sat}"

    def test_filter_glonass_only(self, sp3_file_2022):
        reader = SP3Reader(filepaths=sp3_file_2022, desiredGNSSsystems=["R"])
        df = reader.read()
        sats = df["Satellite"].unique()
        for sat in sats:
            assert sat.startswith("R"), f"Expected GLONASS only, got {sat}"


class TestSP3ReaderMultipleFiles:
    """Test reading multiple SP3 files at once."""

    def test_read_two_files(self, sp3_file_2020, sp3_file_2022):
        reader = SP3Reader(filepaths=[sp3_file_2020, sp3_file_2022])
        df = reader.read()
        assert len(df) > 0
        # Should have data from both files
        epochs = df["Epoch"].unique()
        assert len(epochs) > 289  # more than single file


class TestSP3ReaderUnitConversion:
    """Test coordinate and clock bias unit conversions."""

    def test_coords_in_km_when_disabled(self, sp3_file_2022):
        reader = SP3Reader(filepaths=sp3_file_2022, coords_in_meter=False)
        df = reader.read()
        g01 = df[df["Satellite"] == "G01"].iloc[0]
        # Should be in km: ~13882 km
        assert_almost_equal(g01["X"], 13882.271938, decimal=2)

    def test_clock_in_microseconds_when_disabled(self, sp3_file_2022):
        reader = SP3Reader(filepaths=sp3_file_2022, clock_bias_in_sec=False)
        df = reader.read()
        g01 = df[df["Satellite"] == "G01"].iloc[0]
        # Should be in microseconds: ~469.127768
        assert_almost_equal(g01["Clock Bias"], 469.127768, decimal=2)



class TestReadSP3Nav2022:
    """Test legacy readSP3Nav function with 2022 file."""

    @pytest.fixture(autouse=True, scope="class")
    @classmethod
    def _read_sp3(cls, sp3_file_2022):
        (cls.sat_pos, cls.epoch_dates, cls.navGNSSsystems,
         cls.nEpochs, cls.epochInterval, cls.success
        ) = readSP3Nav(sp3_file_2022)

    def test_success(self):
        assert self.success == 1

    def test_nepochs(self):
        assert self.nEpochs > 0

    def test_epoch_interval(self):
        assert_almost_equal(self.epochInterval, 300.0, decimal=1)

    def test_gnss_systems(self):
        assert "G" in self.navGNSSsystems
        assert "R" in self.navGNSSsystems
        assert "E" in self.navGNSSsystems
        assert "C" in self.navGNSSsystems

    def test_sat_pos_has_gps(self):
        assert "G" in self.sat_pos
        assert len(self.sat_pos["G"]) > 0

    def test_sat_pos_gps_epoch0(self):
        # GPS PRN 1 at epoch 0 should have position
        gps = self.sat_pos["G"]
        assert 0 in gps
        assert 1 in gps[0]  # PRN 1
        pos = gps[0][1]  # shape (1, 3)
        distance = np.sqrt(pos[0, 0]**2 + pos[0, 1]**2 + pos[0, 2]**2)
        assert 2e7 < distance < 3e7

    def test_sat_pos_gps_g01_values(self):
        # PG01: 13882.271938 -21710.006216 5357.125462 (km) -> meters
        pos = self.sat_pos["G"][0][1]  # shape (1, 3)
        assert_almost_equal(pos[0, 0], 13882.271938 * 1000, decimal=0)
        assert_almost_equal(pos[0, 1], -21710.006216 * 1000, decimal=0)
        assert_almost_equal(pos[0, 2], 5357.125462 * 1000, decimal=0)

    def test_epoch_dates_shape(self):
        assert self.epoch_dates.shape[0] == self.nEpochs
        assert self.epoch_dates.shape[1] == 6  # YYYY, MM, DD, hh, mm, ss

    def test_first_epoch_date(self):
        # First epoch: 2022 1 1 0 0 0
        assert str(self.epoch_dates[0, 0]).strip() == "2022"
        assert str(self.epoch_dates[0, 1]).strip() == "1"
        assert str(self.epoch_dates[0, 2]).strip() == "1"

    def test_sat_pos_galileo(self):
        assert "E" in self.sat_pos
        assert len(self.sat_pos["E"]) > 0

    def test_sat_pos_glonass(self):
        assert "R" in self.sat_pos
        assert len(self.sat_pos["R"]) > 0

    def test_sat_pos_beidou(self):
        assert "C" in self.sat_pos
        assert len(self.sat_pos["C"]) > 0


class TestReadSP3Nav2020:
    """Test legacy readSP3Nav function with 2020 file."""

    @pytest.fixture(autouse=True, scope="class")
    @classmethod
    def _read_sp3(cls, sp3_file_2020):
        (cls.sat_pos, cls.epoch_dates, cls.navGNSSsystems,
         cls.nEpochs, cls.epochInterval, cls.success
        ) = readSP3Nav(sp3_file_2020)

    def test_success(self):
        assert self.success == 1

    def test_nepochs(self):
        assert self.nEpochs > 0

    def test_first_epoch_date(self):
        assert str(self.epoch_dates[0, 0]).strip() == "2020"
        assert str(self.epoch_dates[0, 1]).strip() == "11"
        assert str(self.epoch_dates[0, 2]).strip() == "5"


class TestReadSP3NavFiltering:
    """Test GNSS system filtering in legacy readSP3Nav."""

    def test_gps_only(self, sp3_file_2022):
        sat_pos, _, navGNSSsystems, _, _, success = readSP3Nav(
            sp3_file_2022, desiredGNSSsystems=["G"]
        )
        assert success == 1
        assert "G" in sat_pos
        assert len(sat_pos["G"]) > 0

    def test_galileo_only(self, sp3_file_2022):
        sat_pos, _, navGNSSsystems, _, _, success = readSP3Nav(
            sp3_file_2022, desiredGNSSsystems=["E"]
        )
        assert success == 1
        assert "E" in sat_pos
        assert len(sat_pos["E"]) > 0
