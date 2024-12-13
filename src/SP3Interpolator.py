import pandas as pd
import numpy as np
from datetime import datetime
from typing import Literal, Tuple
from gnssmultipath import PickleHandler
from gnssmultipath.readRinexObs import readRinexObs
from gnssmultipath.SP3Reader import SP3Reader
from gnssmultipath.SatelliteEphemerisToECEF import SatelliteEphemerisToECEF, Kepler2ECEF
from gnssmultipath.Geodetic_functions import date2gpstime, date2gpstime_vectorized, gpstime2date_arrays
from tqdm import tqdm


class SP3Interpolator:
    """
    SP3 files already provide satellite positions in the Earth-Centered Earth-Fixed (ECEF) frame.
    These positions are valid for the provided epoch timestamp and do not require
    an additional correction for Earth's rotation when interpolated directly.
    """
    
    def __init__(self, sp3_dataframe, epoch_interval, receiver_position: tuple = None):
        """
        Initializes the SP3 Interpolator with the provided SP3 DataFrame.

        :param sp3_dataframe: Pandas DataFrame containing SP3 data (columns: ['Epoch', 'Satellite', 'X', 'Y', 'Z', 'Clock Bias']).
        :param epoch_interval: Interval between each epoch in seconds.
        :param receiver_position: Tuple of receiver ECEF coordinates (x, y, z) in meters. Defaults to None.
        """
        self.sp3_dataframe = sp3_dataframe
        self.epoch_interval = epoch_interval
        self.receiver_position = receiver_position  # Receiver ECEF coordinates

    @staticmethod
    def epoch_to_seconds(epoch):
        """
        Convert an epoch (datetime object) into seconds since 2000-01-01.

        :param epoch: Datetime object representing the epoch.
        :return: Total seconds since the reference epoch (January 1st, 2000).
        """
        base_time = datetime(2000, 1, 1)
        delta = epoch - base_time
        return delta.total_seconds()

    @staticmethod
    def interppol(x, y, n):
        """
        Polynomial interpolation using Neville's algorithm.

        :param x: Array of x values (time differences from the target epoch)
        :param y: Array of y values (positions or velocities to interpolate)
        :param n: Number of data points for interpolation
        :return: Interpolated value
        """
        y_copy = y.copy()  # Avoid modifying the original array
        for j in range(1, n):
            for i in range(n - j):
                y_copy[i] = (x[i + j] * y_copy[i] - x[i] * y_copy[i + 1]) / (x[i + j] - x[i])
        return y_copy[0]
    
    

    def interpolate_sat_coordinates(self, time_epochs, gnss_systems, n_interpol_points=7, output_format:Literal["pd.DataFrame","dict"]="pd.DataFrame"):
        """
        Interpolates satellite positions for all systems and satellites for given time epochs.
    
        :param time_epochs: Array of observation times in GPS time format (week, TOW).
        :param gnss_systems: List of GNSS systems to include (e.g., ['G', 'R', 'E']).
        :param n_interpol_points: Number of nearest points for interpolation.
        :param output_format: Desired output format. Options are 'dict' (default) or 'dataframe'.
        :return: Interpolated positions in the specified output format.
        """
        from gnssmultipath.Geodetic_functions import gpstime2date_arrays
    
        # Convert GPS time to datetime objects
        observation_times = gpstime2date_arrays(time_epochs[:, 0], time_epochs[:, 1])
    
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
    
                # Compute time differences for all observation times
                time_diffs = np.abs(satellite_seconds[:, None] - observation_seconds)
    
                # Find the indices of the nearest points for each observation time
                nearest_indices = np.argsort(time_diffs, axis=0)[:n_interpol_points, :]
    
                # Prepare arrays for interpolation
                interpolated_positions_sat = np.zeros((len(observation_seconds), 3))  # (epochs, X/Y/Z)
    
                for obs_idx in range(len(observation_seconds)):
                    # Select nearest points for the current observation time
                    idx = nearest_indices[:, obs_idx]
                    nearest_times = satellite_seconds[idx]
                    nearest_positions = satellite_positions[idx]
    
                    # Interpolate positions using Neville's algorithm
                    time_diff = nearest_times - observation_seconds[obs_idx]
                    for i in range(3):  # X, Y, Z
                        interpolated_positions_sat[obs_idx, i] = self.interppol(
                            time_diff, nearest_positions[:, i], len(nearest_times)
                        )
    
                # Add to the results dictionary under the current GNSS system
                interpolated_positions[gnss][satellite] = interpolated_positions_sat
                pbar.update(1)
    
        pbar.close()
    
        # If output_format is 'dataframe', convert the dictionary to a DataFrame
        if output_format == "pd.DataFrame":
            rows = []
            for gnss, satellites in interpolated_positions.items():
                for satellite, positions in satellites.items():
                    for idx, (time, position) in enumerate(zip(observation_times, positions)):
                        rows.append({
                            "Time": datetime(*time),  # Convert to datetime object
                            "Satellite": satellite,
                            "X": position[0],
                            "Y": position[1],
                            "Z": position[2]
                        })
    
            return pd.DataFrame(rows)
    
        return interpolated_positions
    
    
    
    @staticmethod
    def filter_by_prn(interpolated_data, prn_list:list):
        """
        Filters the interpolated data for specific PRN numbers.

        :param interpolated_data: Interpolated data as a dictionary or DataFrame.
        :param prn_list: List of PRN numbers to include (e.g., ['G01', 'E02', 'R03']).
        :return: Filtered data as the same type as input (dict or DataFrame).
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
    
        :param interpolated_data: Interpolated data as a dictionary or DataFrame.
        :param gnss_systems: List of GNSS systems to include (e.g., ['G', 'R', 'E']).
        :return: Filtered data as the same type as input (dict or DataFrame).
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
    sp3_df = sp3_reader.read_file()
    
    # Print metadata
    metadata = sp3_reader.get_metadata()
    nEpochs_sp3 = metadata["Number of Epochs"]
    epochInterval_sp3 = metadata["Epoch Interval (s)"]
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

    
    
    # diff = interpolated_positions["G"]["G01"] - results_rnav["Sat_position"]["G"]["position"]["1"][0:440]
    # diff_E = interpolated_positions["E"]["E01"] - results_rnav["Sat_position"]["E"]["position"]["1"][0:440]
    # diff_R = interpolated_positions["R"]["R01"] - results_rnav["Sat_position"]["R"]["position"]["1"][0:440]
    # print(diff_R)
    

    
    
    
