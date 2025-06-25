"""
This module is using pytest to run tests on the GNSSPositionEstimator Python class.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com

"""

import sys
import os
import numpy as np
from numpy.testing import assert_almost_equal

# Define the project path and append it to the system path
project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_path)
sys.path.append(os.path.join(project_path, 'src'))

from gnssmultipath.GNSSPositionEstimator import GNSSPositionEstimator



# Test data paths
rinObs = os.path.join(project_path,"TestData/ObservationFiles/OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.rnx")
rinNav = os.path.join(project_path,"TestData/NavigationFiles/BRDC00IGS_R_20220010000_01D_MN.rnx")


# Initial coordinates for the receiver
initial_coordinates = np.array([3149785.9652, 598260.8822, 5495348.4927])  # Exact coordinates for receiver
x_rec, y_rec, z_rec = initial_coordinates.T
desired_time = np.array([2022, 1, 1, 0, 0, 30.0000000])  # Desired time for the estimate
desired_system = "G"  # Desired system for GNSS



def test_with_initial_coordinates():
    # Initialize the GNSSPositionEstimator object
    GNSSPos = GNSSPositionEstimator(
        rinex_obs_file = rinObs,
        rinex_nav_file = rinNav,
        desired_time = desired_time,
        desired_system = desired_system,
        x_rec_approx=x_rec,
        y_rec_approx=y_rec,
        z_rec_approx=z_rec,
        elevation_cut_off_angle=15
    )
    # Estimate position and extract the results
    estimated_position, _ = GNSSPos.estimate_position()
    X, Y, Z, dT = estimated_position.T
    computed_pos = np.array([X, Y, Z])
    computed_clock_error = np.array([dT])

    # Expected coordinates for the receiver
    expected_coords = np.array([3149788.2203, 598262.2791, 5495355.6211])
    expected_clock_error = np.array([2.7228175783920506e-08])

    # Use assert_almost_equal to compare the computed and expected positions
    assert_almost_equal(computed_pos, expected_coords, decimal=3)
    assert_almost_equal(computed_clock_error, expected_clock_error, decimal=8)



def test_with_WGS84_UTM32_as_output():
    # Initialize the GNSSPositionEstimator object
    GNSSPos = GNSSPositionEstimator(
        rinex_obs_file = rinObs,
        rinex_nav_file = rinNav,
        desired_time = desired_time,
        desired_system = desired_system,
        elevation_cut_off_angle=15,
        crs="EPSG:32632"
    )
    # Estimate position and extract the results
    estimated_position, _ = GNSSPos.estimate_position()
    X, Y, Z, dT = estimated_position.T
    computed_pos = np.array([X, Y, Z])
    computed_clock_error = np.array([dT])

    # Expected coordinates for the receiver
    expected_coords = np.array([598128.600, 6642363.650, 71.237])
    expected_clock_error = np.array([2.7228175783920506e-08])

    # Use assert_almost_equal to compare the computed and expected positions
    assert_almost_equal(computed_pos, expected_coords, decimal=2)
    assert_almost_equal(computed_clock_error, expected_clock_error, decimal=8)



