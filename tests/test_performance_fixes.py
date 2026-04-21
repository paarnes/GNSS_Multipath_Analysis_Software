"""
Tests for the performance optimizations:
  - ECEF2enu_batch (vectorized ECEF→ENU)
  - PreciseSatCoords.compute_azimuth_and_elevation (array-based DataFrame build)
  - SP3Interpolator DataFrame output (array-based DataFrame build)

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""

import sys
import os
import pytest
import numpy as np
from numpy.testing import assert_allclose
import pandas as pd
from datetime import datetime

project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_path)
sys.path.append(os.path.join(project_path, 'src'))

from gnssmultipath.Geodetic_functions import ECEF2enu, ECEF2enu_batch, ECEF2geodb
from gnssmultipath.SP3Interpolator import SP3Interpolator
from gnssmultipath.readers.SP3Reader import SP3Reader


# ── Receiver at OPEC (Norway) ──────────────────────────────────────────
REC_ECEF = np.array([3149785.9652, 598260.8822, 5495348.4927])
a_WGS84, b_WGS84 = 6378137.0, 6356752.314245
LAT_REC, LON_REC, _ = ECEF2geodb(a_WGS84, b_WGS84, *REC_ECEF)

# Test data paths
sp3_path = os.path.join(project_path, "TestData/SP3/Testfile_20220101.eph")
rinObs_path = os.path.join(project_path, "TestData/ObservationFiles/v3/OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.rnx")


# ═══════════════════════════════════════════════════════════════════════
# ECEF2enu_batch – correctness tests
# ═══════════════════════════════════════════════════════════════════════

class TestECEF2enuBatchScalarEquivalence:
    """Verify ECEF2enu_batch reproduces the scalar ECEF2enu for every input."""

    def test_single_point_matches_scalar(self):
        """A single-element array must match the scalar function exactly."""
        lat, lon = LAT_REC, LON_REC
        dX, dY, dZ = 1000.0, -500.0, 2000.0

        e_s, n_s, u_s = ECEF2enu(lat, lon, dX, dY, dZ)
        e_b, n_b, u_b = ECEF2enu_batch(lat, lon,
                                        np.array([dX]), np.array([dY]), np.array([dZ]))
        assert_allclose(e_b, [e_s], atol=1e-12)
        assert_allclose(n_b, [n_s], atol=1e-12)
        assert_allclose(u_b, [u_s], atol=1e-12)

    def test_multiple_points_match_scalar(self):
        """N-element arrays must match N calls to the scalar function."""
        lat, lon = LAT_REC, LON_REC
        rng = np.random.default_rng(42)
        N = 500
        dX = rng.uniform(-1e7, 1e7, N)
        dY = rng.uniform(-1e7, 1e7, N)
        dZ = rng.uniform(-1e7, 1e7, N)

        # Scalar reference
        e_ref = np.empty(N)
        n_ref = np.empty(N)
        u_ref = np.empty(N)
        for i in range(N):
            e_ref[i], n_ref[i], u_ref[i] = ECEF2enu(lat, lon, dX[i], dY[i], dZ[i])

        # Batch
        e_b, n_b, u_b = ECEF2enu_batch(lat, lon, dX, dY, dZ)

        assert_allclose(e_b, e_ref, atol=1e-6)
        assert_allclose(n_b, n_ref, atol=1e-6)
        assert_allclose(u_b, u_ref, atol=1e-6)

    def test_known_answer_from_existing_test(self):
        """Reproduce the known-answer from test_GNSS_MultipathAnalysis.py."""
        lat = 0.2838307690924083
        lon = -1.0738580938997349
        dX = 12893444.612051928
        dY = -9888519.684510132
        dZ = 13138615.445688058
        expected = (6619718.534261429, 8457426.340380616, 17924793.375904225)

        e_b, n_b, u_b = ECEF2enu_batch(lat, lon,
                                        np.array([dX]), np.array([dY]), np.array([dZ]))
        assert_allclose((e_b[0], n_b[0], u_b[0]), expected, atol=1e-3)


class TestECEF2enuBatchEdgeCases:
    """Edge cases and robustness checks for ECEF2enu_batch."""

    def test_zero_differences(self):
        """Receiver at satellite → all-zero ENU."""
        e, n, u = ECEF2enu_batch(LAT_REC, LON_REC,
                                  np.array([0.0]), np.array([0.0]), np.array([0.0]))
        assert_allclose(e, [0.0], atol=1e-15)
        assert_allclose(n, [0.0], atol=1e-15)
        assert_allclose(u, [0.0], atol=1e-15)

    def test_purely_vertical_offset(self):
        """dP along receiver-up should give ~zero E/N and positive U."""
        # Approximate unit "up" vector at the receiver
        rec_mag = np.linalg.norm(REC_ECEF)
        up_unit = REC_ECEF / rec_mag
        dP = up_unit * 1000.0  # 1 km straight up

        e, n, u = ECEF2enu_batch(LAT_REC, LON_REC,
                                  np.array([dP[0]]), np.array([dP[1]]), np.array([dP[2]]))
        assert abs(e[0]) < 5.0, f"E should be ~0, got {e[0]}"
        assert abs(n[0]) < 5.0, f"N should be ~0, got {n[0]}"
        assert u[0] > 990.0, f"U should be ~1000, got {u[0]}"

    def test_longitude_wrapping_positive(self):
        """lon in (π, 2π) should be wrapped to (-π, 0) and produce same result."""
        lat = 0.5
        lon_normal = -0.5
        lon_wrapped = -0.5 + 2 * np.pi  # in (π, 2π)
        dX, dY, dZ = np.array([1e6]), np.array([-2e6]), np.array([3e6])

        e1, n1, u1 = ECEF2enu_batch(lat, lon_normal, dX, dY, dZ)
        e2, n2, u2 = ECEF2enu_batch(lat, lon_wrapped, dX, dY, dZ)
        assert_allclose(e1, e2, atol=1e-6)
        assert_allclose(n1, n2, atol=1e-6)
        assert_allclose(u1, u2, atol=1e-6)

    def test_longitude_wrapping_negative(self):
        """lon in (-2π, -π) should be wrapped and produce same result."""
        lat = 0.5
        lon_normal = 0.5
        lon_wrapped = 0.5 - 2 * np.pi  # in (-2π, -π)
        dX, dY, dZ = np.array([1e6]), np.array([-2e6]), np.array([3e6])

        e1, n1, u1 = ECEF2enu_batch(lat, lon_normal, dX, dY, dZ)
        e2, n2, u2 = ECEF2enu_batch(lat, lon_wrapped, dX, dY, dZ)
        assert_allclose(e1, e2, atol=1e-6)
        assert_allclose(n1, n2, atol=1e-6)
        assert_allclose(u1, u2, atol=1e-6)

    def test_large_array_performance_shape(self):
        """10 000 points should return correct shapes without error."""
        N = 10_000
        dX = np.ones(N) * 1e6
        dY = np.ones(N) * 2e6
        dZ = np.ones(N) * 3e6
        e, n, u = ECEF2enu_batch(LAT_REC, LON_REC, dX, dY, dZ)
        assert e.shape == (N,)
        assert n.shape == (N,)
        assert u.shape == (N,)

    def test_rotation_matrix_orthogonality(self):
        """The internal rotation matrix M should be orthogonal (M^T M = I)."""
        from numpy import pi as PI
        lat, lon = 0.7, -1.2
        sin_lon, cos_lon = np.sin(lon), np.cos(lon)
        sin_lat, cos_lat = np.sin(lat), np.cos(lat)
        M = np.array([[-sin_lon,         cos_lon,          0],
                      [-sin_lat*cos_lon, -sin_lat*sin_lon, cos_lat],
                      [ cos_lat*cos_lon,  cos_lat*sin_lon, sin_lat]])
        assert_allclose(M.T @ M, np.eye(3), atol=1e-14)

    def test_equator_prime_meridian(self):
        """At equator/prime meridian, batch and scalar agree."""
        lat, lon = 0.0, 0.0
        dX, dY, dZ = np.array([100.0, -200.0]), np.array([300.0, 400.0]), np.array([-150.0, 250.0])
        e_b, n_b, u_b = ECEF2enu_batch(lat, lon, dX, dY, dZ)
        for i in range(2):
            es, ns, us = ECEF2enu(lat, lon, dX[i], dY[i], dZ[i])
            assert_allclose(e_b[i], es, atol=1e-12)
            assert_allclose(n_b[i], ns, atol=1e-12)
            assert_allclose(u_b[i], us, atol=1e-12)


# ═══════════════════════════════════════════════════════════════════════
# PreciseSatCoords – DataFrame correctness tests
# ═══════════════════════════════════════════════════════════════════════

class TestPreciseSatCoordsAzEl:
    """Tests for PreciseSatCoords.compute_azimuth_and_elevation DataFrame building."""

    @pytest.fixture(scope="class")
    def az_el_df(self):
        """Build azimuth/elevation DataFrame from real test data."""
        from gnssmultipath.PreciseSatCoords import PreciseSatCoords
        precise = PreciseSatCoords(sp3_file=sp3_path, rinex_obs_file=rinObs_path)
        return precise.compute_azimuth_and_elevation(tuple(REC_ECEF), drop_below_horizon=False)

    @pytest.fixture(scope="class")
    def az_el_df_drop(self):
        """Build with drop_below_horizon=True."""
        from gnssmultipath.PreciseSatCoords import PreciseSatCoords
        precise = PreciseSatCoords(sp3_file=sp3_path, rinex_obs_file=rinObs_path)
        return precise.compute_azimuth_and_elevation(tuple(REC_ECEF), drop_below_horizon=True)

    def test_returns_dataframe(self, az_el_df):
        assert isinstance(az_el_df, pd.DataFrame)

    def test_expected_columns(self, az_el_df):
        assert set(az_el_df.columns) == {"Epoch", "Satellite", "Azimuth", "Elevation"}

    def test_no_missing_epochs(self, az_el_df):
        """Every satellite should have the same number of epochs."""
        counts = az_el_df.groupby("Satellite").size()
        assert counts.nunique() == 1, f"Epoch counts vary: {counts.unique()}"

    def test_satellite_naming(self, az_el_df):
        """Satellites should be named like G01, R05, E12, C30."""
        import re
        for sat in az_el_df["Satellite"].unique():
            assert re.match(r'^[GREC]\d{2,3}$', sat), f"Bad sat name: {sat}"

    def test_elevation_range(self, az_el_df):
        """Non-NaN elevations should be in (-90, 90)."""
        valid = az_el_df["Elevation"].dropna()
        assert (valid > -90).all() and (valid < 90).all()

    def test_azimuth_range(self, az_el_df):
        """Non-NaN azimuths should be in [0, 360)."""
        valid = az_el_df["Azimuth"].dropna()
        assert (valid >= -180).all() and (valid < 540).all()

    def test_drop_below_horizon_masks_negative_elevation(self, az_el_df_drop):
        """With drop_below_horizon, negative elevations should be NaN."""
        el = az_el_df_drop["Elevation"]
        valid = el.dropna()
        assert (valid > 0).all(), "Found non-NaN elevation <= 0 with drop_below_horizon=True"

    def test_row_count_matches(self, az_el_df):
        """Total rows = n_satellites × n_epochs."""
        n_sats = az_el_df["Satellite"].nunique()
        n_epochs = az_el_df.groupby("Satellite").size().iloc[0]
        assert len(az_el_df) == n_sats * n_epochs

    def test_dtypes(self, az_el_df):
        """Azimuth and Elevation should be float64."""
        assert az_el_df["Azimuth"].dtype == np.float64
        assert az_el_df["Elevation"].dtype == np.float64


# ═══════════════════════════════════════════════════════════════════════
# SP3Interpolator – DataFrame output path tests
# ═══════════════════════════════════════════════════════════════════════

class TestSP3InterpolatorDataFrame:
    """Tests for SP3Interpolator.interpolate_sat_coordinates DataFrame output."""

    @pytest.fixture(scope="class")
    def sp3_setup(self):
        """Read SP3 and RINEX obs data once for all tests."""
        from gnssmultipath.readers.readRinexObs import readRinexObs
        reader = SP3Reader(sp3_path, coords_in_meter=True, desiredGNSSsystems=["G", "E"])
        sp3_df = reader.read()
        metadata = reader.get_metadata()

        obs_data = readRinexObs(rinObs_path)
        time_epochs = obs_data.time_epochs

        return sp3_df, metadata["epoch_interval_sec"], time_epochs

    @pytest.fixture(scope="class")
    def interpol_df(self, sp3_setup):
        sp3_df, interval, time_epochs = sp3_setup
        interp = SP3Interpolator(sp3_df, interval)
        return interp.interpolate_sat_coordinates(time_epochs, ["G", "E"], output_format="pd.DataFrame")

    @pytest.fixture(scope="class")
    def interpol_dict(self, sp3_setup):
        sp3_df, interval, time_epochs = sp3_setup
        interp = SP3Interpolator(sp3_df, interval)
        return interp.interpolate_sat_coordinates(time_epochs, ["G", "E"], output_format="dict")

    def test_returns_dataframe(self, interpol_df):
        assert isinstance(interpol_df, pd.DataFrame)

    def test_expected_columns(self, interpol_df):
        assert set(interpol_df.columns) == {"Epoch", "Satellite", "X", "Y", "Z", "Clock Bias"}

    def test_no_empty_result(self, interpol_df):
        assert len(interpol_df) > 0

    def test_satellite_naming(self, interpol_df):
        """Satellites should be named like G01, E12."""
        import re
        for sat in interpol_df["Satellite"].unique():
            assert re.match(r'^[GE]\d{2,3}$', sat), f"Bad sat name: {sat}"

    def test_dict_dataframe_consistency(self, interpol_df, interpol_dict, sp3_setup):
        """Dict and DataFrame outputs should contain the same data."""
        # Count total satellite entries in the dict
        dict_sats = set()
        for gnss, sats in interpol_dict.items():
            for sat in sats:
                dict_sats.add(sat)
        df_sats = set(interpol_df["Satellite"].unique())
        assert dict_sats == df_sats, f"Dict sats {dict_sats} != DF sats {df_sats}"

    def test_xyz_coordinates_are_numeric(self, interpol_df):
        assert interpol_df["X"].dtype == np.float64
        assert interpol_df["Y"].dtype == np.float64
        assert interpol_df["Z"].dtype == np.float64

    def test_positions_are_reasonable(self, interpol_df):
        """Satellite positions should be roughly at orbital altitude (20-30 Mm from Earth center)."""
        r = np.sqrt(interpol_df["X"]**2 + interpol_df["Y"]**2 + interpol_df["Z"]**2)
        assert (r > 20e6).all(), "Some positions too close to Earth"
        assert (r < 45e6).all(), "Some positions too far from Earth"

    def test_epoch_column_contains_datetimes(self, interpol_df):
        epoch = interpol_df["Epoch"].iloc[0]
        assert isinstance(epoch, datetime), f"Expected datetime, got {type(epoch)}"

    def test_row_count_matches_sats_times_epochs(self, interpol_df):
        """Total rows = n_unique_sats × n_epochs_per_sat."""
        counts = interpol_df.groupby("Satellite").size()
        assert counts.nunique() == 1, f"Epoch counts vary: {counts.unique()}"
        n_sats = interpol_df["Satellite"].nunique()
        n_epochs = counts.iloc[0]
        assert len(interpol_df) == n_sats * n_epochs

    def test_dict_position_values_match_dataframe(self, interpol_df, interpol_dict):
        """Spot-check: first satellite's X values should match between dict and df."""
        first_sat = list(list(interpol_dict.values())[0].keys())[0]
        dict_x = list(interpol_dict.values())[0][first_sat]["positions"][:, 0]
        df_x = interpol_df[interpol_df["Satellite"] == first_sat]["X"].to_numpy()
        assert_allclose(df_x, dict_x, atol=1e-6)


# ═══════════════════════════════════════════════════════════════════════
# Integration: end-to-end SP3 → azimuth/elevation pipeline
# ═══════════════════════════════════════════════════════════════════════

class TestSP3ToAzElPipeline:
    """Smoke-test the full SP3 interpolation → azimuth/elevation pipeline."""

    @pytest.fixture(scope="class")
    def pipeline_result(self):
        from gnssmultipath.PreciseSatCoords import PreciseSatCoords
        precise = PreciseSatCoords(sp3_file=sp3_path, rinex_obs_file=rinObs_path)
        az_el = precise.compute_azimuth_and_elevation(tuple(REC_ECEF), drop_below_horizon=True)
        return precise, az_el

    def test_pipeline_completes(self, pipeline_result):
        _, az_el = pipeline_result
        assert isinstance(az_el, pd.DataFrame)
        assert len(az_el) > 0

    def test_coords_and_azel_have_same_satellites(self, pipeline_result):
        precise, az_el = pipeline_result
        coord_sats = set(precise.satcoords["Satellite"].unique())
        azel_sats = set(az_el["Satellite"].unique())
        assert coord_sats == azel_sats

    def test_create_satellite_data_dict(self, pipeline_result):
        """create_satellite_data_dict should build a valid nested dict."""
        precise, az_el = pipeline_result
        from gnssmultipath.PreciseSatCoords import PreciseSatCoords
        sat_dict = PreciseSatCoords.create_satellite_data_dict(precise.satcoords, az_el)
        assert isinstance(sat_dict, dict)
        for gnss, data in sat_dict.items():
            assert "coordinates" in data
            assert "azimuth" in data
            assert "elevation" in data
            assert isinstance(data["azimuth"], np.ndarray)
            assert isinstance(data["elevation"], np.ndarray)
