import pandas as pd
import numpy as np
from datetime import datetime
from gnssmultipath import PickleHandler
from gnssmultipath.readRinexObs import readRinexObs
from gnssmultipath.read_SP3Nav_NEW import SP3Reader
from gnssmultipath.SatelliteEphemerisToECEF import SatelliteEphemerisToECEF, Kepler2ECEF
from gnssmultipath.Geodetic_functions import date2gpstime, date2gpstime_vectorized, gpstime2date_arrays
from tqdm import tqdm


class SP3Interpolator:
    """
    SP3 files already provide satellite positions in the Earth-Centered Earth-Fixed (ECEF) frame.
    These positions are valid for the provided epoch timestamp and do not require
    an additional correction for Earth's rotation when interpolated directly.
    """
    
    
    def __init__(self, sp3_dataframe, epoch_interval, receiver_position:tuple=None):
        """
        Initializes the SP3 Interpolator with the provided SP3 DataFrame.

        :param sp3_dataframe: Pandas DataFrame containing SP3 data (columns: ['Epoch', 'Satellite', 'X', 'Y', 'Z', 'Clock Bias']).
        :param epoch_interval: Interval between each epoch in seconds.
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

    def interpolate(self, satellite, target_epoch, n_interpol_points:int=7):
        """
        Interpolates the satellite position at the target_epoch.

        :param satellite: Satellite identifier (e.g., 'G01').
        :param target_epoch: Target epoch as a datetime object.
        :return: Interpolated satellite position (x, y, z).
        """
        # Convert target_epoch into seconds since the reference epoch
        target_epoch_seconds = self.epoch_to_seconds(target_epoch)

        # Filter the DataFrame for the specific satellite
        satellite_data = self.sp3_dataframe[self.sp3_dataframe['Satellite'] == satellite].copy()

        # Convert the 'Epoch' column to seconds since the reference epoch
        satellite_data['Epoch_Seconds'] = satellite_data['Epoch'].apply(self.epoch_to_seconds)

        # Sort by epoch to ensure consistent ordering
        satellite_data = satellite_data.sort_values(by='Epoch_Seconds')

        # Select the nearest epochs
        satellite_data['Time_Diff'] = np.abs(satellite_data['Epoch_Seconds'] - target_epoch_seconds)
        nearest_points = satellite_data.nsmallest(n_interpol_points, 'Time_Diff')  # Use n_interpol_points points for cubic interpolation

        # Ensure sufficient points are available
        if len(nearest_points) < 4:
            raise ValueError(f"Insufficient data points for interpolation of {satellite} at {target_epoch}.")

        # Prepare data for interpolation
        time_diffs = nearest_points['Epoch_Seconds'].to_numpy() - target_epoch_seconds
        positions = nearest_points[['X', 'Y', 'Z']].to_numpy()

        # Interpolate positions using Neville's algorithm
        interpolated_positions = np.zeros(3)
        for i in range(3):
            interpolated_positions[i] = self.interppol(time_diffs, positions[:, i], len(time_diffs))
            
    
        return interpolated_positions





def interpolate_all_satellites(sp3_dataframe, time_epochs, gnss_systems, n_interpol_points=7):
    """
    Interpolates satellite positions for all systems and satellites for given time epochs.

    :param sp3_dataframe: Pandas DataFrame containing SP3 data (columns: ['Epoch', 'Satellite', 'X', 'Y', 'Z', 'Clock Bias']).
    :param time_epochs: Array of observation times in GPS time format (week, TOW).
    :param gnss_systems: List of GNSS systems to include (e.g., ['G', 'R', 'E']).
    :param n_interpol_points: Number of nearest points for interpolation.
    :return: Dictionary containing interpolated positions for each GNSS system and satellite at each epoch.
    """

    # Convert GPS time to datetime objects
    observation_times = gpstime2date_arrays(time_epochs[:, 0], time_epochs[:, 1])

    # Create a dictionary to store results for all GNSS systems
    interpolated_positions = {}

    # Loop through each GNSS system
    for gnss in gnss_systems:
        # Filter the SP3 DataFrame to include only the current GNSS system
        gnss_data = sp3_dataframe[sp3_dataframe['Satellite'].str[0] == gnss]

        # Initialize dictionary for the current GNSS system
        interpolated_positions[gnss] = {}

        # Interpolate for each satellite in the current GNSS system
        for satellite in gnss_data['Satellite'].unique():
            # Filter data for the current satellite
            satellite_data = gnss_data[gnss_data['Satellite'] == satellite].copy()

            # Convert the 'Epoch' column to seconds since the reference epoch
            satellite_data['Epoch_Seconds'] = satellite_data['Epoch'].apply(SP3Interpolator.epoch_to_seconds)

            # Sort by epoch to ensure consistent ordering
            satellite_data = satellite_data.sort_values(by='Epoch_Seconds')

            # Initialize list for this satellite's interpolated positions
            satellite_positions = []

            # Process each observation time
            for obs_time in observation_times:
                obs_time = datetime(*obs_time)
                # Convert observation time to seconds since the reference epoch
                obs_time_seconds = SP3Interpolator.epoch_to_seconds(obs_time)

                # Calculate time differences and find nearest points
                satellite_data['Time_Diff'] = np.abs(satellite_data['Epoch_Seconds'] - obs_time_seconds)
                nearest_points = satellite_data.nsmallest(n_interpol_points, 'Time_Diff')

                # Ensure sufficient points for interpolation
                if len(nearest_points) < 4:
                    satellite_positions.append([np.nan, np.nan, np.nan])  # Insufficient data
                    continue

                # Prepare data for interpolation
                time_diffs = nearest_points['Epoch_Seconds'].to_numpy() - obs_time_seconds
                positions = nearest_points[['X', 'Y', 'Z']].to_numpy()

                # Interpolate positions using Neville's algorithm
                interpolated_position = np.zeros(3)
                for i in range(3):
                    interpolated_position[i] = SP3Interpolator.interppol(time_diffs, positions[:, i], len(time_diffs))

                satellite_positions.append(interpolated_position)

            # Add to the results dictionary under the current GNSS system
            interpolated_positions[gnss][satellite] = np.array(satellite_positions)

    return interpolated_positions




def interpolate_all_satellites_with_tqdm(sp3_dataframe, time_epochs, gnss_systems, n_interpol_points=7):
    """
    Interpolates satellite positions for all systems and satellites for given time epochs.

    :param sp3_dataframe: Pandas DataFrame containing SP3 data (columns: ['Epoch', 'Satellite', 'X', 'Y', 'Z', 'Clock Bias']).
    :param time_epochs: Array of observation times in GPS time format (week, TOW).
    :param gnss_systems: List of GNSS systems to include (e.g., ['G', 'R', 'E']).
    :param n_interpol_points: Number of nearest points for interpolation.
    :return: Dictionary containing interpolated positions for each GNSS system and satellite at each epoch.
    """

    # Convert GPS time to datetime objects
    observation_times = gpstime2date_arrays(time_epochs[:, 0], time_epochs[:, 1])

    # Create a dictionary to store results for all GNSS systems
    interpolated_positions = {}

    # Progress bar setup
    total_satellites = sum(
        len(sp3_dataframe[sp3_dataframe['Satellite'].str[0] == gnss]['Satellite'].unique())
        for gnss in gnss_systems
    )
    total_epochs = len(observation_times)
    total_steps = total_satellites * total_epochs
    bar_format = '{desc}: {percentage:3.0f}%|{bar}|'
    pbar = tqdm(total=total_steps, desc="Interpolating satellite coordinates", position=0, leave=True, bar_format=bar_format)

    # Loop through each GNSS system
    for gnss in gnss_systems:
        # Filter the SP3 DataFrame to include only the current GNSS system
        gnss_data = sp3_dataframe[sp3_dataframe['Satellite'].str[0] == gnss]

        # Initialize dictionary for the current GNSS system
        interpolated_positions[gnss] = {}

        # Interpolate for each satellite in the current GNSS system
        for satellite in gnss_data['Satellite'].unique():
            # Filter data for the current satellite
            satellite_data = gnss_data[gnss_data['Satellite'] == satellite].copy()

            # Convert the 'Epoch' column to seconds since the reference epoch
            satellite_data['Epoch_Seconds'] = satellite_data['Epoch'].apply(SP3Interpolator.epoch_to_seconds)

            # Sort by epoch to ensure consistent ordering
            satellite_data = satellite_data.sort_values(by='Epoch_Seconds')

            # Initialize list for this satellite's interpolated positions
            satellite_positions = []

            # Process each observation time
            for obs_time in observation_times:
                obs_time = datetime(*obs_time)
                # Convert observation time to seconds since the reference epoch
                obs_time_seconds = SP3Interpolator.epoch_to_seconds(obs_time)

                # Calculate time differences and find nearest points
                satellite_data['Time_Diff'] = np.abs(satellite_data['Epoch_Seconds'] - obs_time_seconds)
                nearest_points = satellite_data.nsmallest(n_interpol_points, 'Time_Diff')

                # Ensure sufficient points for interpolation
                if len(nearest_points) < 4:
                    satellite_positions.append([np.nan, np.nan, np.nan])  # Insufficient data
                    pbar.update(1)
                    continue

                # Prepare data for interpolation
                time_diffs = nearest_points['Epoch_Seconds'].to_numpy() - obs_time_seconds
                positions = nearest_points[['X', 'Y', 'Z']].to_numpy()

                # Interpolate positions using Neville's algorithm
                interpolated_position = np.zeros(3)
                for i in range(3):
                    interpolated_position[i] = SP3Interpolator.interppol(time_diffs, positions[:, i], len(time_diffs))

                satellite_positions.append(interpolated_position)
                pbar.update(1)

            # Add to the results dictionary under the current GNSS system
            interpolated_positions[gnss][satellite] = np.array(satellite_positions)

    pbar.close()
    return interpolated_positions





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
    
    # Target epoch for interpolation
    # target_epoch = datetime(2022, 1, 1, 0, 0, 0)
    target_epoch = datetime(*desired_time)
    
    # Interpolate the position of satellite 'G01'
    interpolated_position = interpolator.interpolate('G01', target_epoch)
    print(f"Interpolated Position: {interpolated_position}")
    
    
    print(f"Difference:\n{np.array(xyz_nav)[0:3] - np.atleast_2d(interpolated_position).T}")
    
    
    #%%    
    # Interpolate all satellites
    interpolated_positions = interpolate_all_satellites(sp3_df, time_epochs, gnss_systems)
    
    # Print results for satellite 'G01'
    satellite_id = 'G01'
    print(f"Interpolated positions for {satellite_id}:")
    print(interpolated_positions["G"][satellite_id])

    diff = interpolated_positions["G"]["G01"] - results_rnav["Sat_position"]["G"]["position"]["1"][0:440]
    print(diff)
    
    
    
