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

from gnssmultipath.BroadNavPositionEstimator import BroadNavPositionEstimator
from gnssmultipath.SatelliteEphemerisToECEF import SatelliteEphemerisToECEF
from gnssmultipath.readRinexObs import readRinexObs


# Test data paths
rinObs = os.path.join(project_path,"TestData/ObservationFiles/OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.rnx")
rinNav = os.path.join(project_path,"TestData/NavigationFiles/BRDC00IGS_R_20220010000_01D_MN.rnx")


# Initial coordinates for the receiver
initial_coordinates = np.array([3149785.9652, 598260.8822, 5495348.4927])  # Exact coordinates for receiver
x_rec, y_rec, z_rec = initial_coordinates.T
desired_time = np.array([2022, 1, 1, 0, 0, 30.0000000])  # Desired time for the estimate
desired_system = "G"  # Desired system for GNSS



def test_with_initial_coordinates():
    # Initialize the BroadNavPositionEstimator object
    GNSSPos = BroadNavPositionEstimator(
        rinObs, rinNav, desired_time, desired_system, x_rec_approx=x_rec, y_rec_approx=y_rec, z_rec_approx=z_rec, elevation_cut_off_angle=15
    )
    # Estimate position and extract the results
    estimated_position, stats = GNSSPos.estimate_position()
    X, Y, Z, dT = estimated_position.T
    computed_pos = np.array([X, Y, Z])
    computed_clock_error = np.array([dT])

    # Expected coordinates for the receiver
    expected_coords = np.array([3149788.2203, 598262.2791, 5495355.6211])
    expected_clock_error = np.array([2.7228175783920506e-08])

    # Use assert_almost_equal to compare the computed and expected positions
    assert_almost_equal(computed_pos, expected_coords, decimal=3)
    assert_almost_equal(computed_clock_error, expected_clock_error, decimal=8)


def test_without_initial_coordinates():
    # Initialize the BroadNavPositionEstimator object
    GNSSPos = BroadNavPositionEstimator(
        rinObs, rinNav, desired_time, desired_system, elevation_cut_off_angle=15)
    # Estimate position and extract the results
    estimated_position, stats = GNSSPos.estimate_position()
    X, Y, Z, dT = estimated_position.T
    computed_pos = np.array([X, Y, Z])
    computed_clock_error = np.array([dT])

    # Expected coordinates for the receiver
    expected_coords = np.array([3149788.2203, 598262.2791, 5495355.6211])
    expected_clock_error = np.array([2.7228175783920506e-08])

    # Use assert_almost_equal to compare the computed and expected positions
    assert_almost_equal(computed_pos, expected_coords, decimal=3)
    assert_almost_equal(computed_clock_error, expected_clock_error, decimal=8)


def test_GLONASS_without_initial_coordinates():
    # Initialize the BroadNavPositionEstimator object
    desired_system ="R"
    desired_time = np.array([2022, 1, 1, 0, 0, 30.0000000])
    GNSSPos = BroadNavPositionEstimator(
        rinObs, rinNav, desired_time, desired_system)
    # Estimate position and extract the results
    estimated_position, stats = GNSSPos.estimate_position()
    X, Y, Z, dT = estimated_position.T
    computed_pos = np.array([X, Y, Z])
    computed_clock_error = np.array([dT])

    # Expected coordinates for the receiver
    expected_coords = np.array([3149782.9046, 598279.0639, 5495355.8332])
    expected_clock_error = np.array([5.083982295608008e-08])

    # Use assert_almost_equal to compare the computed and expected positions
    assert_almost_equal(computed_pos, expected_coords, decimal=3)
    assert_almost_equal(computed_clock_error, expected_clock_error, decimal=8)


def test_GLONASS_with_initial_coordinates():
    # Initialize the BroadNavPositionEstimator object
    desired_system ="R"
    desired_time = np.array([2022, 1, 1, 0, 0, 30.0000000])
    GNSSPos = BroadNavPositionEstimator(
        rinObs, rinNav, desired_time, desired_system,
        x_rec_approx=x_rec, y_rec_approx=y_rec, z_rec_approx=z_rec)
    # Estimate position and extract the results
    estimated_position, stats = GNSSPos.estimate_position()
    X, Y, Z, dT = estimated_position.T
    computed_pos = np.array([X, Y, Z])
    computed_clock_error = np.array([dT])

    # Expected coordinates for the receiver
    expected_coords = np.array([3149782.9046, 598279.0639, 5495355.8332])
    expected_clock_error = np.array([5.083982295608008e-08])

    # Use assert_almost_equal to compare the computed and expected positions
    assert_almost_equal(computed_pos, expected_coords, decimal=3)
    assert_almost_equal(computed_clock_error, expected_clock_error, decimal=8)



def test_with_arrays_as_input_galileo():
    fasit = np.array([3149785.9652, 598260.8822, 5495348.4927])
    x_rec, y_rec, z_rec = fasit.T
    desired_sys = "E"
    desired_time = np.array([2022, 1, 1, 0, 0, 30.0000000])
    navObj = SatelliteEphemerisToECEF(rinNav,  x_rec, y_rec, z_rec, desired_sys)
    GNSS_obs, _, _, _, time_epochs, _, GNSSsystems, \
        obsCodes, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ = readRinexObs(rinObs)

    GNSSPos = BroadNavPositionEstimator(desired_time=desired_time,
                                    desired_system=desired_sys,
                                    navdata=navObj,
                                    GNSS_obs=GNSS_obs,
                                    time_epochs=time_epochs,
                                    GNSSsystems=GNSSsystems,
                                    obsCodes=obsCodes,
                                    elevation_cut_off_angle=15)
    estimated_position, stats = GNSSPos.estimate_position()
    X,Y,Z,dT = estimated_position.T

    computed_pos = np.array([X, Y, Z])
    computed_clock_error = np.array([dT])

    expected_coords = np.array([3149787.01967, 598262.0575, 5495355.1109 ])
    expected_clock_error = np.array([1.7165079447804373e-08])

    # Use assert_almost_equal to compare the computed and expected positions
    assert_almost_equal(computed_pos, expected_coords, decimal=3)
    assert_almost_equal(computed_clock_error, expected_clock_error, decimal=8)



if __name__ == '__main__':
    pytest.main()
