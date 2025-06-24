"""
This module is using pytest to run tests on the software.
Each pull request to the master branch will trigger
a workflow that runs these tests.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com

"""

import sys
import os
import pytest
import numpy as np
from numpy.testing import assert_almost_equal
import warnings


# Define the project path and append it to the system path
project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_path)  # REMOVE AFTER FILE IS MOVED INTO GNSS MULTIPATH
sys.path.append(os.path.join(project_path, 'src'))

from gnssmultipath.SP3PositionEstimator import SP3PositionEstimator
from gnssmultipath.SatelliteEphemerisToECEF import SatelliteEphemerisToECEF
from gnssmultipath.readRinexObs import readRinexObs


# Test data paths
rinObs = os.path.join(project_path,"TestData/ObservationFiles/OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.rnx")
sp3 = os.path.join(project_path,"TestData/SP3/Testfile_20220101.eph")


# Initial coordinates for the receiver
initial_coordinates = np.array([3149785.9652, 598260.8822, 5495348.4927])  # Exact coordinates for receiver
x_rec, y_rec, z_rec = initial_coordinates.T
desired_time = np.array([2022, 1, 1, 1, 5, 30.0000000])


def test_with_initial_coordinates_GPS():
    desired_system = "G"  # Desired system for GNSS

    # Initialize the SP3PositionEstimator object
    GNSSPos = SP3PositionEstimator(sp3_data=sp3,
                                   rinex_obs_file=rinObs,
                                   desired_time=desired_time,
                                   desired_system=desired_system,
                                   x_rec_approx=x_rec,
                                   y_rec_approx=y_rec,
                                   z_rec_approx=z_rec,
                                   elevation_cut_off_angle=10
                                   )
    # Estimate position and extract the results
    estimated_position, stats = GNSSPos.estimate_position()
    X, Y, Z, dT = estimated_position.T
    computed_pos = np.array([X, Y, Z])
    computed_clock_error = np.array([dT])

    # Expected coordinates for the receiver
    expected_coords = np.array([3149786.9820, 598260.5953, 5495355.5783])
    expected_clock_error = np.array([2.33716370e-08])

    # Use assert_almost_equal to compare the computed and expected positions
    assert_almost_equal(computed_pos, expected_coords, decimal=3)
    assert_almost_equal(computed_clock_error, expected_clock_error, decimal=8)


def test_with_initial_coordinates_GLONASS():
    desired_system = "R"  # Desired system for GNSS

    # Initialize the SP3PositionEstimator object
    GNSSPos = SP3PositionEstimator(sp3_data=sp3,
                                   rinex_obs_file=rinObs,
                                   desired_time=desired_time,
                                   desired_system=desired_system,
                                   x_rec_approx=x_rec,
                                   y_rec_approx=y_rec,
                                   z_rec_approx=z_rec,
                                   elevation_cut_off_angle=10
                                   )
    # Estimate position and extract the results
    estimated_position, stats = GNSSPos.estimate_position()
    X, Y, Z, dT = estimated_position.T
    computed_pos = np.array([X, Y, Z])
    computed_clock_error = np.array([dT])

    # Expected coordinates for the receiver
    expected_coords = np.array([3149784.1907, 598261.5753, 5495359.2244])
    expected_clock_error = np.array([0.00000003])

    # Use assert_almost_equal to compare the computed and expected positions
    assert_almost_equal(computed_pos, expected_coords, decimal=3)
    assert_almost_equal(computed_clock_error, expected_clock_error, decimal=8)


def test_with_initial_coordinates_BeiDou():
    desired_system = "C"  # Desired system for GNSS

    # Initialize the SP3PositionEstimator object
    GNSSPos = SP3PositionEstimator(sp3_data=sp3,
                                   rinex_obs_file=rinObs,
                                   desired_time=desired_time,
                                   desired_system=desired_system,
                                   x_rec_approx=x_rec,
                                   y_rec_approx=y_rec,
                                   z_rec_approx=z_rec,
                                   elevation_cut_off_angle=10
                                   )
    # Estimate position and extract the results
    estimated_position, stats = GNSSPos.estimate_position()
    X, Y, Z, dT = estimated_position.T
    computed_pos = np.array([X, Y, Z])
    computed_clock_error = np.array([dT])

    # Expected coordinates for the receiver
    expected_coords = np.array([3149788.6315,  598263.6383, 5495336.5707 ])
    expected_clock_error = np.array([-0.00000003])

    # Use assert_almost_equal to compare the computed and expected positions
    assert_almost_equal(computed_pos, expected_coords, decimal=3)
    assert_almost_equal(computed_clock_error, expected_clock_error, decimal=8)


def test_with_initial_coordinates_Galileo():
    desired_system = "E"  # Desired system for GNSS

    # Initialize the SP3PositionEstimator object
    GNSSPos = SP3PositionEstimator(sp3_data=sp3,
                                   rinex_obs_file=rinObs,
                                   desired_time=desired_time,
                                   desired_system=desired_system,
                                   x_rec_approx=x_rec,
                                   y_rec_approx=y_rec,
                                   z_rec_approx=z_rec,
                                   elevation_cut_off_angle=10
                                   )
    # Estimate position and extract the results
    estimated_position, stats = GNSSPos.estimate_position()
    X, Y, Z, dT = estimated_position.T
    computed_pos = np.array([X, Y, Z])
    computed_clock_error = np.array([dT])

    # Expected coordinates for the receiver
    expected_coords = np.array([3149789.3839, 598261.6005, 5495362.9073])
    expected_clock_error = np.array([0.00000003])

    # Use assert_almost_equal to compare the computed and expected positions
    assert_almost_equal(computed_pos, expected_coords, decimal=3)
    assert_almost_equal(computed_clock_error, expected_clock_error, decimal=8)


def test_with_initial_coordinates_Galileo_with_20_el_cutoff_and_no_approx_pos():
    desired_system = "E"  # Desired system for GNSS

    # Initialize the SP3PositionEstimator object
    GNSSPos = SP3PositionEstimator(sp3_data=sp3,
                                   rinex_obs_file=rinObs,
                                   desired_time=desired_time,
                                   desired_system=desired_system,
                                   elevation_cut_off_angle=20,
                                   )
    # Estimate position and extract the results
    estimated_position, stats = GNSSPos.estimate_position()
    X, Y, Z, dT = estimated_position.T
    computed_pos = np.array([X, Y, Z])
    computed_clock_error = np.array([dT])

    # Expected coordinates for the receiver
    expected_coords = np.array([3149787.3134, 598262.1634, 5495359.4333])
    expected_clock_error = np.array([0.00000002])

    # Use assert_almost_equal to compare the computed and expected positions
    assert_almost_equal(computed_pos, expected_coords, decimal=3)
    assert_almost_equal(computed_clock_error, expected_clock_error, decimal=8)


if __name__ == '__main__':
    pytest.main()
