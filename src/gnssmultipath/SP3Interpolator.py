import pandas as pd
import numpy as np
from datetime import datetime
from typing import Literal, Tuple
from gnssmultipath import PickleHandler
from gnssmultipath.readRinexObs import readRinexObs
from gnssmultipath.SP3Reader import SP3Reader
from gnssmultipath.SatelliteEphemerisToECEF import SatelliteEphemerisToECEF, Kepler2ECEF
from gnssmultipath.Geodetic_functions import date2gpstime, date2gpstime_vectorized, gpstime2date_arrays, gpstime2date_arrays_with_microsec
from tqdm import tqdm


c = 299792458  # Speed of light [m/s]

class SP3Interpolator:
    """
    SP3 files already provide satellite positions in the Earth-Centered Earth-Fixed (ECEF) frame.
    These positions are valid for the provided epoch timestamp and do not require
    an additional correction for Earth's rotation when interpolated directly.

    Example:
    -------
    .. code-block:: python
        interpolator = SP3Interpolator(sp3_df, epochInterval_sp3)
        interpolated_positions = interpolator.interpolate_sat_coordinates(time_epochs, gnss_systems)
    """

    def __init__(self, sp3_dataframe, epoch_interval, receiver_position: tuple = None):
        """
        Initializes the SP3 Interpolator with the provided SP3 DataFrame.

        Parameter:
        ----------
        - sp3_dataframe: Pandas DataFrame containing SP3 data (columns: ['Epoch', 'Satellite', 'X', 'Y', 'Z', 'Clock Bias']).
        - epoch_interval: Interval between each epoch in seconds.
        - receiver_position: Tuple of receiver ECEF coordinates (x, y, z) in meters. Defaults to None.
        """
        self.sp3_dataframe = sp3_dataframe
        self.epoch_interval = epoch_interval
        self.receiver_position = receiver_position  # Receiver ECEF coordinates

    @staticmethod
    def epoch_to_seconds(epoch):
        """
        Convert an epoch (datetime object) into seconds since 2000-01-01.

        Parameter:
        ----------
        - epoch: Datetime object representing the epoch.

        Return:
        ------
        Total seconds since the reference epoch (January 1st, 2000).
        """
        base_time = datetime(2000, 1, 1)
        delta = epoch - base_time
        return delta.total_seconds()

    @staticmethod
    def interppol(x, y, n):
        """
        Polynomial interpolation using Neville's algorithm.

        Parameter:
        ----------
        - x: Array of x values (time differences from the target epoch)
        - y: Array of y values (positions or velocities to interpolate)
        - n: Number of data points for interpolation

        Return:
        ------
        - Interpolated value
        """
        y_copy = y.copy()  # Avoid modifying the original array
        for j in range(1, n):
            for i in range(n - j):
                y_copy[i] = (x[i + j] * y_copy[i] - x[i] * y_copy[i + 1]) / (x[i + j] - x[i])
        return y_copy[0]

    def compute_relativistic_correction_single_sat(self, prn, time_epochs):
        """
        Compute the relativistic clock correction for a single satellite.

        Parameter:
        ----------
        - prn: PRN of the satellite (e.g., 'G12').
        - time_epochs: Array of observation times in GPS time format (week, TOW).

        Return:
        ------
        Relativistic clock correction values as a NumPy array.
        """
        delta_t = self.epoch_interval  # Time step for velocity approximation (seconds)

        # Times for velocity computation
        t_delta_plus = time_epochs + np.array([[0, delta_t]]).T
        t_delta_minus = time_epochs - np.array([[0, delta_t]]).T

        # Interpolate satellite positions at t+delta_t, t, and t-delta_t
        pos_plus, _ = self.interpolate_single_satellite(prn, t_delta_plus)
        pos_minus, _ = self.interpolate_single_satellite(prn, t_delta_minus)
        pos_now, _ = self.interpolate_single_satellite(prn, time_epochs)

        # Compute satellite velocities
        velocities = (pos_plus - pos_minus) / (2 * delta_t)

        # Compute relativistic corrections
        corrections = -2 * np.sum(pos_now * velocities, axis=1) / c**2

        return corrections

    def interpolate_single_satellite(self, prn, time_epochs, n_interpol_points=7):
        """
        Interpolates satellite positions and clock biases for a single satellite specified by its PRN.

        Parameter:
        ----------
        - prn: PRN of the satellite to interpolate (e.g., 'G12').
        - time_epochs: Array of observation times in GPS time format (week, TOW).
        - n_interpol_points: Number of nearest points for interpolation.

        Return:
        ------
        - Interpolated positions and clock biases as a dictionary.
        """
        # Convert GPS time to datetime objects
        if len(time_epochs) > 2:
            observation_times = gpstime2date_arrays_with_microsec(time_epochs[:, 0], time_epochs[:, 1])
        else:
            observation_times = gpstime2date_arrays_with_microsec(time_epochs[0], time_epochs[1])

        # Convert observation times to seconds since the reference epoch
        observation_seconds = np.array([self.epoch_to_seconds(datetime(*obs)) for obs in observation_times])

        # Filter the SP3 DataFrame to include only the specified satellite
        satellite_data = self.sp3_dataframe[self.sp3_dataframe['Satellite'] == prn].copy()

        if satellite_data.empty:
            raise ValueError(f"No data found for satellite {prn}.")

        # Convert the 'Epoch' column to seconds since the reference epoch
        satellite_data['Epoch_Seconds'] = satellite_data['Epoch'].apply(self.epoch_to_seconds)

        # Sort by epoch to ensure consistent ordering
        satellite_data = satellite_data.sort_values(by='Epoch_Seconds')

        # Extract satellite data for vectorized processing
        satellite_seconds = satellite_data['Epoch_Seconds'].to_numpy()
        satellite_positions = satellite_data[['X', 'Y', 'Z']].to_numpy()
        satellite_clock_bias = satellite_data['Clock Bias'].to_numpy()

        # Compute time differences for all observation times
        time_diffs = np.abs(satellite_seconds[:, None] - observation_seconds)

        # Find the indices of the nearest points for each observation time
        nearest_indices = np.argsort(time_diffs, axis=0)[:n_interpol_points, :]

        # Prepare arrays for interpolation
        interpolated_positions = np.zeros((len(observation_seconds), 3))  # (epochs, X/Y/Z)
        interpolated_clock_bias = np.zeros(len(observation_seconds))  # (epochs, Clock Bias)

        for obs_idx in range(len(observation_seconds)):
            # Select nearest points for the current observation time
            idx = nearest_indices[:, obs_idx]
            nearest_times = satellite_seconds[idx]
            nearest_positions = satellite_positions[idx]
            nearest_clock_biases = satellite_clock_bias[idx]

            # Interpolate positions using Neville's algorithm
            time_diff = nearest_times - observation_seconds[obs_idx]
            for i in range(3):  # X, Y, Z
                interpolated_positions[obs_idx, i] = self.interppol(
                    time_diff, nearest_positions[:, i], len(nearest_times)
                )

            # Interpolate clock bias
            interpolated_clock_bias[obs_idx] = self.interppol(
                time_diff, nearest_clock_biases, len(nearest_times)
            )

        return interpolated_positions, interpolated_clock_bias


    def interpolate_sat_coordinates(self, time_epochs, gnss_systems, n_interpol_points=7, output_format: Literal["pd.DataFrame", "dict"] = "pd.DataFrame"):
        """
        Interpolates satellite positions and clock biases for all systems and satellites for given time epochs.

        Parameter:
        ----------
        - time_epochs: Array of observation times in GPS time format (week, TOW).
        - gnss_systems: List of GNSS systems to include (e.g., ['G', 'R', 'E']).
        - n_interpol_points: Number of nearest points for interpolation.
        - output_format: Desired output format. Options are 'dict' or 'dataframe'.


        Return:
        -------
        - Interpolated positions and clock biases in the specified output format.
        """

        # Convert GPS time to datetime objects
        if len(time_epochs) > 2:
            # observation_times = gpstime2date_arrays(time_epochs[:, 0], time_epochs[:, 1])
            observation_times = gpstime2date_arrays_with_microsec(time_epochs[:, 0], time_epochs[:, 1])
        else:
            # observation_times = gpstime2date_arrays(time_epochs[0], time_epochs[1])
            observation_times = gpstime2date_arrays_with_microsec(time_epochs[0], time_epochs[1])

        # Convert observation times to seconds since the reference epoch
        observation_seconds = np.array([self.epoch_to_seconds(datetime(*obs)) for obs in observation_times])

        # Create a dictionary to store results for all GNSS systems
        interpolated_positions = {}

        # Progress bar setup
        total_satellites = sum(
            len(self.sp3_dataframe[self.sp3_dataframe['Satellite'].str[0] == gnss]['Satellite'].unique())
            for gnss in gnss_systems
        )
        bar_format = '{desc}: {percentage:3.0f}%|{bar}| ({n_fmt}/{total_fmt})'
        pbar = tqdm(total=total_satellites, desc="Interpolating satellite coordinates", position=0, leave=True, bar_format=bar_format)

        # Loop through each GNSS system
        for gnss in gnss_systems:
            # Filter the SP3 DataFrame to include only the current GNSS system
            gnss_data = self.sp3_dataframe[self.sp3_dataframe['Satellite'].str[0] == gnss]

            # Initialize dictionary for the current GNSS system
            interpolated_positions[gnss] = {}

            # Interpolate for each satellite in the current GNSS system
            for satellite in gnss_data['Satellite'].unique():
                # Filter data for the current satellite
                satellite_data = gnss_data[gnss_data['Satellite'] == satellite].copy()

                # Convert the 'Epoch' column to seconds since the reference epoch
                satellite_data['Epoch_Seconds'] = satellite_data['Epoch'].apply(self.epoch_to_seconds)

                # Sort by epoch to ensure consistent ordering
                satellite_data = satellite_data.sort_values(by='Epoch_Seconds')

                # Extract satellite data for vectorized processing
                satellite_seconds = satellite_data['Epoch_Seconds'].to_numpy()
                satellite_positions = satellite_data[['X', 'Y', 'Z']].to_numpy()
                satellite_clock_bias = satellite_data['Clock Bias'].to_numpy()

                # Compute time differences for all observation times
                time_diffs = np.abs(satellite_seconds[:, None] - observation_seconds)

                # Find the indices of the nearest points for each observation time
                nearest_indices = np.argsort(time_diffs, axis=0)[:n_interpol_points, :]

                # Prepare arrays for interpolation
                interpolated_positions_sat = np.zeros((len(observation_seconds), 3))  # (epochs, X/Y/Z)
                interpolated_clock_bias = np.zeros(len(observation_seconds))  # (epochs, Clock Bias)

                for obs_idx in range(len(observation_seconds)):
                    # Select nearest points for the current observation time
                    idx = nearest_indices[:, obs_idx]
                    nearest_times = satellite_seconds[idx]
                    nearest_positions = satellite_positions[idx]
                    nearest_clock_biases = satellite_clock_bias[idx]

                    # Interpolate positions using Neville's algorithm
                    time_diff = nearest_times - observation_seconds[obs_idx]
                    for i in range(3):  # X, Y, Z
                        interpolated_positions_sat[obs_idx, i] = self.interppol(
                            time_diff, nearest_positions[:, i], len(nearest_times)
                        )
                    # Interpolate clock bias
                    interpolated_clock_bias[obs_idx] = self.interppol(
                        time_diff, nearest_clock_biases, len(nearest_times)
                    )

                # Add to the results dictionary under the current GNSS system
                interpolated_positions[gnss][satellite] = {
                    "positions": interpolated_positions_sat,
                    "clock_bias": interpolated_clock_bias
                }
                pbar.update(1)

        pbar.close()

        # If output_format is 'dataframe', convert the dictionary to a DataFrame
        if output_format == "pd.DataFrame":
            rows = []
            for gnss, satellites in interpolated_positions.items():
                for satellite, data in satellites.items():
                    positions = data["positions"]
                    clock_biases = data["clock_bias"]
                    for idx, (time, position, clock_bias) in enumerate(zip(observation_times, positions, clock_biases)):
                        rows.append({
                            "Epoch": datetime(*time),  # Convert to datetime object
                            "Satellite": satellite,
                            "X": position[0],
                            "Y": position[1],
                            "Z": position[2],
                            "Clock Bias": clock_bias
                        })

            return pd.DataFrame(rows)

        return interpolated_positions


    @staticmethod
    def filter_by_prn(interpolated_data, prn_list:list):
        """
        Filters the interpolated data for specific PRN numbers.

        Parameter:
        ----------
        - interpolated_data: Interpolated data as a dictionary or DataFrame.
        - prn_list: List of PRN numbers to include (e.g., ['G01', 'E02', 'R03']).

        Return:
        ------
        - Filtered data as the same type as input (dict or DataFrame).
        """
        if isinstance(interpolated_data, pd.DataFrame):
            return interpolated_data[interpolated_data['Satellite'].str[0:].isin(prn_list)]
        elif isinstance(interpolated_data, dict):
            filtered_data = {}
            for gnss, satellites in interpolated_data.items():
                filtered_data[gnss] = {
                    sat: positions
                    for sat, positions in satellites.items()
                    if sat[0:] in prn_list
                }
            return filtered_data
        else:
            raise TypeError("Unsupported data type for filtering. Expected dict or pd.DataFrame.")

    @staticmethod
    def filter_by_system(interpolated_data, gnss_systems: list):
        """
        Filters the interpolated data for specific GNSS systems.

        Parameter:
        ----------
        - interpolated_data: Interpolated data as a dictionary or DataFrame.
        - gnss_systems: List of GNSS systems to include (e.g., ['G', 'R', 'E']).

        Return:
        ------
        - Filtered data as the same type as input (dict or DataFrame).
        """
        if isinstance(interpolated_data, pd.DataFrame):
            return interpolated_data[interpolated_data['Satellite'].str[0].isin(gnss_systems)]
        elif isinstance(interpolated_data, dict):
            return {gnss: satellites for gnss, satellites in interpolated_data.items() if gnss in gnss_systems}
        else:
            raise TypeError("Unsupported data type for filtering. Expected dict or pd.DataFrame.")





if __name__ == "__main__":
    gnss_systems = ["G","R","E","C"]
    rinObs = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.rnx"
    sp3 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\SP3\Testfile_20220101.eph"
    rinNav = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\BRDC00IGS_R_20220010000_01D_MN.rnx"
    results_rnav = PickleHandler.read_zstd_pickle(r"C:\Users\perhe\Desktop\TEST BROADCAST\analysisResults.pkl")

    x_rec_approx, y_rec_approx, z_rec_approx = [3149785.9652, 598260.8822, 5495348.4927]
    navdata = SatelliteEphemerisToECEF(rinNav, x_rec_approx, y_rec_approx, z_rec_approx, gnss_systems)
    desired_time = [2022, 1, 1, 0, 3, 0]
    des_time = date2gpstime_vectorized(np.array([desired_time]))[-1]
    xyz_nav = navdata.get_sat_ecef_coordinates(desired_time=des_time, PRN="G01")

    GNSS_obs,_, _, _, time_epochs, nepochs, GNSSsystems,\
    obsCodes, _, _, tInterval, _, _, _, _, _, _,\
    _, _, _, _, _, _, _, _ = readRinexObs(rinObs)

    observation_times = gpstime2date_arrays(time_epochs[:,0],time_epochs[:,1])

    # sat_coord_sp3, epoch_dates_sp3, navGNSSsystems_sp3, nEpochs_sp3, epochInterval_sp3, _= readSP3Nav(sp3)
    sp3_reader = SP3Reader(sp3, coords_in_meter=True, desiredGNSSsystems=gnss_systems)
    sp3_df = sp3_reader.read()

    # Print metadata
    metadata = sp3_reader.get_metadata()
    nEpochs_sp3 = metadata["n_epochs"]
    epochInterval_sp3 = metadata["epoch_interval_sec"]
    navGNSSsystems_sp3 = ["G"]

    ### Interpolate
    # interpolator = SP3Interpolator(sp3_df, epochInterval_sp3)
    interpolator = SP3Interpolator(sp3_df, epochInterval_sp3)




    ## Interpolate single sataellite
    # target_epoch = datetime(*desired_time) # Target epoch for interpolation
    # Interpolate the position of satellite 'G01'
    # interpolated_position = interpolator.interpolate('G01', target_epoch)
    # print(f"Interpolated Position: {interpolated_position}")
    # print(f"Difference:\n{np.array(xyz_nav)[0:3] - np.atleast_2d(interpolated_position).T}")


    # Interpolate all satellites for all systems
    interpolated_positions = interpolator.interpolate_sat_coordinates(time_epochs, gnss_systems)


    # interpolated_positions = interpolator.interpolate_sat_coordinates(time_epochs, gnss_systems, output_format="dict")
    # diff = interpolated_positions["G"]["G01"]["positions"] - results_rnav["Sat_position"]["G"]["position"]["1"][0:440]
    # diff_E = interpolated_positions["E"]["E01"]["positions"] - results_rnav["Sat_position"]["E"]["position"]["1"][0:440]
    # diff_R = interpolated_positions["R"]["R01"]["positions"] - results_rnav["Sat_position"]["R"]["position"]["1"][0:440]
    # print(diff_R)





