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
from gnssmultipath.Geodetic_functions import date2gpstime, date2gpstime_vectorized, compute_satellite_azimut_and_elevation_angle, gpstime2date_arrays_with_microsec
from gnssmultipath.SP3Reader import SP3Reader
from gnssmultipath.SP3Interpolator import SP3Interpolator
from gnssmultipath.readRinexObs import readRinexObs
from gnssmultipath.StatisticalAnalysis import StatisticalAnalysis

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)

c = 299792458  # Speed of light [m/s]
OMEGA_EARTH = 7.2921159e-5 # Earth's rotational speed in rad/s



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
        sp3_metadata_dict: Optional[dict] = None,
        elevation_cut_off_angle: Optional[int] = 10
    ):
        # Desired time in GPS format
        self.desired_time = np.atleast_2d(date2gpstime_vectorized(desired_time)).T
        self.desired_sys = desired_system
        self.elevation_cut_off_angle = elevation_cut_off_angle

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

    def __estimate_without_low_satellites(
        self,
        X: np.ndarray,
        Y: np.ndarray,
        Z: np.ndarray,
        Rji: np.ndarray,
        sat_nr: list,
        x: float,
        y: float,
        z: float,
        t: float,
        dTi0: float
    )  -> Tuple[float, float, float, float, np.ndarray, np.ndarray, np.ndarray, list]:
        """
        Filter out satellites with an elevation angle below the specified threshold and recompute the position and clock bias.

        Parameters:
        ----------
        X, Y, Z: Arrays of satellite ECEF coordinates.
        dT_rel: Array of relativistic clock corrections.
        Rji: Array of corrected pseudoranges.
        sat_nr: List of satellite numbers.
        x, y, z: Current receiver position.
        t: Current time.
        dTi0: Current clock bias.

        Returns:
        -------
        Tuple containing the updated receiver position (x, y, z), clock bias (dTi0),
        the design matrix (A), the vector of observations (l), the normal matrix (N),
        the observation vector (h), and the list of active satellite numbers.
        """
        # Compute azimuth and elevation for all satellites
        azimuths, elevations = compute_satellite_azimut_and_elevation_angle(X, Y, Z, x, y, z)

        # Filter satellites with elevation < 15 degrees for example
        active_sat_indices = np.where(elevations >= self.elevation_cut_off_angle)[0]
        low_sat_indices = np.where(elevations < self.elevation_cut_off_angle)[0]
        for i in low_sat_indices:
            PRN = f"{self.desired_sys}{str(sat_nr[i]).zfill(2)}"
            logger.info(f"\nINFO(SP3PositionEstimator): {PRN} is exluded due to low elevation: {np.round(elevations[i],3)}°")


        # Update satellite-related data based on active satellites
        X = X[active_sat_indices]
        Y = Y[active_sat_indices]
        Z = Z[active_sat_indices]
        Rji = Rji[active_sat_indices]
        sat_nr = [sat_nr[i] for i in active_sat_indices]

        # Recompute design matrix and observation corrections for active satellites
        diff = np.column_stack([X - x, Y - y, Z - z])  # Shape (num_sats, 3)
        rho = np.linalg.norm(diff, axis=1)
        dax = -diff[:, 0] / rho
        day = -diff[:, 1] / rho
        daz = -diff[:, 2] / rho
        dadT = np.ones_like(rho)

        l = Rji - (rho + c * dTi0)

        A = np.column_stack([dax, day, daz, dadT])

        # Solve for corrections
        N = A.T @ A
        h = A.T @ l
        dx = np.linalg.solve(N, h)  # Solve for parameter updates

        # Update receiver position and clock bias
        x += dx[0]
        y += dx[1]
        z += dx[2]
        dTi0 += dx[3] / c

        # return X, Y, Z, sat_nr, x, y, z, dTi0, A,l,N,h
        return x, y, z, dTi0, A,l,N,h, sat_nr

    def compute_relativistic_correction(self, sp3_interpolator, prn_lst, signal_transmission_time):
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


    def update_receiver_coordinates(self, x, y, z):
        """Update receiver coordinates in the class attributes."""
        self.x_rec_approx = x
        self.y_rec_approx = y
        self.z_rec_approx = z


    def apply_sagnac_correction(self, X_sat: float, Y_sat: float, Z_sat: float,
                                X_rec: float, Y_rec: float, Z_rec: float) -> Tuple[float, float, float]:
        """
        Applies the Sagnac (Earth rotation) correction by rotating the satellite's
        ECEF position from transmit-time frame to the receiver's receive-time frame.

        Parameters
        ----------
        X_sat, Y_sat, Z_sat : float
            Satellite ECEF coordinates at signal-transmission time (meters).
        X_rec, Y_rec, Z_rec : float
            Receiver ECEF coordinates at signal-reception time (meters).

        Returns
        -------
        X_corr, Y_corr, Z_corr : float
            The satellite ECEF coordinates after applying the Earth-rotation correction
            to align them with Earth's orientation at receiver time.
        """
        # 1) Compute geometric range and travel time
        dx = X_sat - X_rec
        dy = Y_sat - Y_rec
        dz = Z_sat - Z_rec
        R = np.sqrt(dx*dx + dy*dy + dz*dz)  # 3D distance
        dt = R / c                          # signal travel time

        # 2) Rotation angle around the Z-axis
        # Typically, we rotate by +Omega*dt if we want to bring satellite coords
        # forward in time from transmit to receive epoch.
        # Adjust sign if you see an opposite sign offset in your final solutions.
        alpha = OMEGA_EARTH * dt

        cos_alpha = np.cos(alpha)
        sin_alpha = np.sin(alpha)

        # 3) Apply rotation around Earth's Z-axis
        # This rotates the point (X_sat, Y_sat) about Z by 'alpha'.
        # X' = cos(alpha)*X - sin(alpha)*Y
        # Y' = sin(alpha)*X + cos(alpha)*Y
        # Z' = Z
        # Here, we treat the satellite’s position as if Earth has rotated alpha from transmit to receive.
        X_rot =  cos_alpha * X_sat + sin_alpha * Y_sat
        Y_rot = -sin_alpha * X_sat + cos_alpha * Y_sat
        Z_rot =  Z_sat  # no change in the Z coordinate for rotation about Z-axis

        return X_rot, Y_rot, Z_rot


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

        # Preallocate arrays for satellite data
        t = self.desired_time[1]
        X, Y, Z, dT_rel = np.zeros(num_sats), np.zeros(num_sats), np.zeros(num_sats), np.zeros(num_sats)
        dTj, Tj_GPS, Rji = np.zeros(num_sats), np.full((num_sats, 1), t), np.zeros(num_sats)

        while improvement > tol and iteration < max_iterations:
            for i, PRN in enumerate(prn_lst):
                signal_transmission_time = np.array([self.desired_time[0], Tj_GPS[i]])
                rel_corr = sp3_interpolator.compute_relativistic_correction_single_sat(PRN, signal_transmission_time)
                dT_rel[i] = rel_corr.item() if hasattr(rel_corr, "item") else rel_corr
                interpolated_positions, interpolated_clock_bias = sp3_interpolator.interpolate_single_satellite(PRN, signal_transmission_time)
                # Ensure interpolated_positions is a 1D array before unpacking
                interpolated_positions = np.asarray(interpolated_positions).flatten()
                X[i], Y[i], Z[i] = interpolated_positions[0], interpolated_positions[1], interpolated_positions[2]
                dTj[i] = interpolated_clock_bias.item() if hasattr(interpolated_clock_bias, "item") else interpolated_clock_bias

                # Apply Sagnac correction
                X_sagnac, Y_sagnac, Z_sagnac = self.apply_sagnac_correction(X[i], Y[i], Z[i], x, y, z)
                X[i], Y[i], Z[i] = X_sagnac, Y_sagnac, Z_sagnac

                # Correct observed pseudorange
                Rji[i] = rji[i] + c * (dTj[i] + dT_rel[i])

                # Compute time when signal left the satellite
                Tj_GPS[i] = t - Rji[i] / c + dTj[i]

            # Compute the vector differences between receiver and satellite positions
            diff = np.column_stack([X - x, Y - y, Z - z])  # Shape (num_sats, 3)

            # Compute distances and partial derivatives
            rho = np.linalg.norm(diff, axis=1)
            dax = -diff[:, 0] / rho
            day = -diff[:, 1] / rho
            daz = -diff[:, 2] / rho
            dadT = np.ones_like(rho)

            # Observation corrections
            # l = Rji + c * dT_rel - (rho + c * dTi0)
            l = Rji - (rho + c * dTi0)
            A = np.column_stack([dax, day, daz, dadT])

            # Solve for corrections
            N = A.T @ A
            h = A.T @ l
            dx = np.linalg.solve(N, h)  # Solve for parameter updates

            # Update receiver position and clock bias
            x += dx[0]
            y += dx[1]
            z += dx[2]
            dTi0 += dx[3] / c

            self.update_receiver_coordinates(x, y, z)

            # Update iteration info
            iteration += 1
            logger.info(f"Iteration {iteration}: dx = {dx}")
            improvement = np.max(np.abs(dx))

        # Try to filter low satellites and update position
        try:
            x, y, z, dTi0, A, l, N, h, sat_nr = self.__estimate_without_low_satellites(X, Y, Z, Rji, sat_nr, x, y, z, signal_transmission_time, dTi0)
        except Exception as e:
            logger.warning(
                f"WARNING (GNSS_MultipathAnalysis): Failed to perform the final position estimation after removing low-elevation satellites. "
                f"Error: {str(e)}"
            )
            pass

        stats_report = StatisticalAnalysis(A, l, N, h).run_statistical_analysis()
        return np.array([x, y, z, dTi0]), stats_report






if __name__ == "__main__":
    pass

