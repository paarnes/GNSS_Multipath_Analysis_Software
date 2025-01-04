"""
Module for estimating the receiver position based on pseudoranges and intepolated satellites coordinates based
on data from SP3 files
"""

import warnings
import logging
import numpy as np
import pandas as pd
from datetime import timedelta, datetime
from typing import Literal, Tuple, Dict, Optional, Union
from gnssmultipath.Geodetic_functions import date2gpstime, date2gpstime_vectorized
from gnssmultipath.SP3Reader import SP3Reader
from gnssmultipath.SP3Interpolator import SP3Interpolator
from gnssmultipath.readRinexObs import readRinexObs
from gnssmultipath.StatisticalAnalysis import StatisticalAnalysis

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)

c = 299792458  # Speed of light [m/s]


class SP3PositionEstimator:
    """
    A class to estimate the position of a GNSS receiver using pseudoranges and SP3 satellite position data.

    Parameter:
    ----------
    - sp3_data: str. Path to the SP3 file containing precise satellite coordinates or a pandas dataframe with the satellite coordinates.
    - rinex_obs_file: str. Path to the RINEX observation file containing pseudorange measurements.
    - desired_time: np.ndarray. Desired time for position estimation as a 1D array in the format
                  [year, month, day, hour, minute, second].
    - desired_system: Literal["G", "E", "R"], optional. GNSS system to use for position estimation.
                    Default is "G" (GPS).
    - x_rec_approx: float, optional. Approximate receiver x-coordinate in ECEF (meters). Default is 0.0.
    - y_rec_approx: float, optional. Approximate receiver y-coordinate in ECEF (meters). Default is 0.0.
    - z_rec_approx: float, optional. Approximate receiver z-coordinate in ECEF (meters). Default is 0.0.
    - GNSS_obs: dict. Pseudorange observations for all GNSS systems.
    - time_epochs: np.ndarray. Array of observation epochs (GPS time).
    - GNSSsystems: dict. Mapping of GNSS systems to indices.
    - obsCodes: dict. Observation codes for each GNSS system.

    """

    def __init__(
        self,
        sp3_data: Union[str, pd.DataFrame],
        desired_time: np.ndarray,
        rinex_obs_file: Optional[str] = None,
        desired_system: Literal["G", "E", "R"] = "G",
        x_rec_approx: float = 0.0,
        y_rec_approx: float = 0.0,
        z_rec_approx: float = 0.0,
        GNSS_obs: Optional[dict] = None,
        time_epochs: Optional[np.ndarray] = None,
        GNSSsystems: Optional[dict] = None,
        obsCodes: Optional[dict] = None,
        sp3_metadata_dict: Optional[dict] = None
    ):
        # Desired time in GPS format
        self.desired_time = np.atleast_2d(date2gpstime_vectorized(desired_time)).T
        self.desired_sys = desired_system

        # Load SP3 data
        if isinstance(sp3_data, str):
            sp3_reader = SP3Reader(sp3_data, coords_in_meter=True, desiredGNSSsystems=[desired_system])
            self.sp3_df = sp3_reader.read()
            self.sp3_metadata_dict = sp3_reader.get_metadata()
            self.sp3_epoch_interval = self.sp3_metadata_dict["epoch_interval_sec"]
        else:
            self.sp3_df = sp3_data
            self.sp3_metadata_dict = sp3_metadata_dict
            self.sp3_epoch_interval = self.sp3_metadata_dict["epoch_interval_sec"]

        # Load RINEX observation data
        if rinex_obs_file:
            self.GNSS_obs, _, _, _, self.time_epochs, _, self.GNSSsystems, \
                self.obsCodes, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ = readRinexObs(rinex_obs_file)
        else:
            self.GNSS_obs = GNSS_obs
            self.time_epochs = time_epochs
            self.GNSSsystems = GNSSsystems
            self.obsCodes = obsCodes

        # Receiver approximate coordinates
        self.x_rec_approx = x_rec_approx
        self.y_rec_approx = y_rec_approx
        self.z_rec_approx = z_rec_approx

        # GNSS system and signal indexing
        self.__gnss_index_mapper_dict = self.__gnss_index_mapper()
        self.sys_idx = self.__gnss_index_mapper_dict[self.desired_sys]
        self.sys_obs_codes = self.obsCodes[self.sys_idx][self.desired_sys]
        self.signal_idx = next((i for i, code in enumerate(self.sys_obs_codes) if code.startswith("C")), -1)

    def __gnss_index_mapper(self) -> dict:
        """Map GNSS system identifiers to their corresponding indices."""
        return {self.GNSSsystems[k]: k for k in self.GNSSsystems}

    @staticmethod
    def __find_idx_of_epoch_closest_in_time(desired_time: np.ndarray, obs_time_epochs: np.ndarray) -> int:
        """Find the index of the epoch in `obs_time_epochs` closest to `desired_time`."""
        total_time_diff = np.abs(obs_time_epochs[:, 0] - desired_time[0, 0]) * 604800 + \
                          np.abs(obs_time_epochs[:, 1] - desired_time[1, 0])
        return np.argmin(total_time_diff) + 1  # Observations are 1-indexed

    def __extract_pseudoranges_and_satellites(self, idx_of_closest_eph) -> Tuple[np.array, np.array]:
        """
        Extract valid pseudoranges and corresponding satellite indices.

        Filters out satellites that are not present in the SP3 file.

        Parameter:
        ----------
        idx_of_closest_eph: int. Index of the closest epoch in the observation file.

        Returns:
        -------
        prn_lst: np.ndarray. Valid satellite PRNs.
        rji: np.ndarray. Corresponding pseudoranges for the valid satellites.
        """
        pseudoranges = self.GNSS_obs[self.desired_sys][idx_of_closest_eph]
        sat_nr = np.where(pseudoranges[:, self.signal_idx] != 0.0)[0]
        rji = pseudoranges[sat_nr, self.signal_idx]

        # Construct PRN list for satellites in RINEX
        prn_lst = [f"{self.desired_sys}{str(prn).zfill(2)}" for prn in sat_nr]

        # Filter out satellites not present in SP3 data
        available_prns = self.sp3_df['Satellite'].unique()
        valid_indices = [i for i, prn in enumerate(prn_lst) if prn in available_prns]

        # Update sat_nr and rji with valid indices only
        sat_nr = sat_nr[valid_indices]
        rji = rji[valid_indices]

        return sat_nr, rji

    def __compute_relativistic_correction(self, sp3_interpolator, prn_lst, signal_transmission_time):
        """
        Compute relativistic corrections for the satellite clock bias.

        NOT IN USE.

        Parameters:
        ----------
        sp3_interpolator: SP3Interpolator. Interpolator object for satellite positions.
        prn_lst: List[str]. List of satellite PRNs.
        signal_transmission_time: float. Signal transmission time in GPS seconds.

        Returns:
        -------
        dT_rel: np.ndarray. Relativistic corrections for each satellite.
        """
        delta_t = self.sp3_epoch_interval  # Time step for velocity approximation (seconds)

        t_delta_plus = np.array([self.desired_time[0], self.desired_time[1] + delta_t])
        t_delta_minus = np.array([self.desired_time[0], self.desired_time[1] - delta_t])

        # Interpolate satellite positions at t+delta_t, t, and t-delta_t
        pos_plus = sp3_interpolator.interpolate_sat_coordinates(
            t_delta_plus, self.desired_sys
        )
        pos_minus = sp3_interpolator.interpolate_sat_coordinates(
            t_delta_minus, self.desired_sys
        )
        pos_now = sp3_interpolator.interpolate_sat_coordinates(
            signal_transmission_time, self.desired_sys
        )

        # Extract satellite positions for PRNs in prn_lst
        pos_plus = sp3_interpolator.filter_by_prn(pos_plus, prn_lst)
        pos_minus = sp3_interpolator.filter_by_prn(pos_minus, prn_lst)
        pos_now = sp3_interpolator.filter_by_prn(pos_now, prn_lst)

        # Compute velocities
        velocities = (pos_plus[['X', 'Y', 'Z']].to_numpy() - pos_minus[['X', 'Y', 'Z']].to_numpy()) / (2 * delta_t)
        positions = pos_now[['X', 'Y', 'Z']].to_numpy()

        # Compute relativistic corrections
        dT_rel = -2 * np.sum(positions * velocities, axis=1) / c**2
        return dT_rel

    def estimate_position(self, max_iterations: int = 15, tol: float = 1e-8) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
        """
        Estimate the receiver's position using a vectorized least-squares algorithm.

        Parameters:
        ----------
        max_iterations: int, optional. Maximum number of iterations for convergence (default is 15).
        tol: float, optional. Convergence tolerance for position updates (default is 1e-8).

        Returns:
        -------
        position: numpy.ndarray. Estimated receiver position in ECEF coordinates (x, y, z) and clock bias (dt).
        """
        x, y, z = self.x_rec_approx, self.y_rec_approx, self.z_rec_approx
        idx_of_closest_eph = self.__find_idx_of_epoch_closest_in_time(self.desired_time, self.time_epochs)
        sat_nr, rji = self.__extract_pseudoranges_and_satellites(idx_of_closest_eph)
        prn_lst = [f"{self.desired_sys}{str(prn).zfill(2)}" for prn in sat_nr]

        num_sats = len(sat_nr)
        if num_sats < 4:
            raise ValueError("Insufficient satellites for position estimation (minimum 4 required).")

        dTi0 = 0.0  # Initial clock bias
        improvement = 10
        iteration = 0
        sp3_interpolator = SP3Interpolator(self.sp3_df, self.sp3_epoch_interval)

        while improvement > tol and iteration < max_iterations:
            # Interpolate satellite positions

            signal_transmission_time = np.array([self.desired_time[0], self.desired_time[1] - rji[0] / c])
            # dT_rel = self.__compute_relativistic_correction(sp3_interpolator, prn_lst, signal_transmission_time)
            sat_positions = sp3_interpolator.interpolate_sat_coordinates(signal_transmission_time, self.desired_sys)
            # sat_positions = sp3_interpolator.interpolate_sat_coordinates(self.desired_time, self.desired_sys) # gir mest nøyaktig hvis ikke mircosekundene tas med


            sat_positions_tmp = sp3_interpolator.filter_by_prn(sat_positions, prn_lst)
            ordered_positions = sat_positions_tmp.set_index('Satellite').loc[prn_lst, ['X', 'Y', 'Z', 'Clock Bias']].reset_index()
            sat_positions = np.array(ordered_positions.iloc[:, 1:4])
            clock_biases = np.array(ordered_positions['Clock Bias'])
            X, Y, Z = sat_positions[:, 0], sat_positions[:, 1], sat_positions[:, 2]

            # Correct pseudoranges using satellite clock bias
            rji_corrected = rji + c * clock_biases

            # Compute pseudorange corrections
            diff = np.column_stack([X - x, Y - y, Z - z])
            rho = np.linalg.norm(diff, axis=1)
            l = rji_corrected - (rho + c * dTi0)


            A = np.column_stack([-diff[:, 0] / rho, -diff[:, 1] / rho, -diff[:, 2] / rho, np.ones_like(rho)])

            # Solve for updates
            N = A.T @ A
            h = A.T @ l
            dx = np.linalg.solve(N, h)

            x += dx[0]
            y += dx[1]
            z += dx[2]
            dTi0 += dx[3] / c
            improvement = np.max(np.abs(dx))
            iteration += 1

        stats_report = StatisticalAnalysis(A, l, N, h).run_statistical_analysis()
        return np.array([x, y, z, dTi0]), stats_report



if __name__ == "__main__":
    # Turn off scientific notation
    np.set_printoptions(suppress=True)

    fasit = np.array([3149785.9652, 598260.8822, 5495348.4927])

    sp3_file = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\SP3\Testfile_20220101.eph"
    rinex_obs_file = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.rnx"
    desired_time = np.array([2022, 1, 1, 2, 0, 30.0000000])
    # desired_time = np.array([2022, 1, 1, 2, 3, 30.0000000])
    desired_time = np.array([2022, 1, 1, 0, 0, 30.0000000])
    desired_system = "G"  # GPS

    position_estimator = SP3PositionEstimator(sp3_file, rinex_obs_file, desired_time, desired_system)
    estimated_position, stats = position_estimator.estimate_position()

    print("Estimated Position:\n", estimated_position)
    print("Difference:\n", fasit- estimated_position[:3])
    # print("Statistics Report:", stats)

