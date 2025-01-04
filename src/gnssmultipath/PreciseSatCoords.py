import numpy as np
import warnings
from tqdm import tqdm
import pandas as pd
from datetime import datetime
from typing import Optional, List, Tuple, Dict
from gnssmultipath.readRinexObs import readRinexObs
from gnssmultipath.SP3Reader import SP3Reader
from gnssmultipath.SP3Interpolator import SP3Interpolator
import gnssmultipath.Geodetic_functions as geodf
from gnssmultipath import PickleHandler

warnings.filterwarnings("ignore")





class PreciseSatCoords:
    """
    Class to interpolate precise satellite coordinates for each GNSS system.

    Parameters:
    ----------
    rinex_obs_file: str.Path to the RINEX observation file.

    sp3_file: str. Path to the SP3 file containing precise satellite coordinates.

    navGNSSsystems: List[str]. List of GNSS systems to process (default is ["G", "R", "E", "C"]).

    Example:
    --------
    .. code-block:: python
        precise_orbits = PreciseOrbits2ECEF(rinex_obs_file, sp3_file, ["G", "E"])
        positions = precise_orbits.interpolate_to_epoch(epoch_time=[2024, 1, 1, 12, 0, 0], system="G", PRN=5)
    """

    def __init__(self, sp3_file: str, rinex_obs_file: str=None, time_epochs: np.ndarray=None, GNSSsystems: Dict=None):

        if rinex_obs_file:
            _,_, _, _, self.time_epochs, _, self.GNSSsystems,\
            _, _, _, _, _, _, _, _, _, _,\
            _, _, _, _, _, _, _, _ = readRinexObs(rinex_obs_file)
        else:
            self.time_epochs = time_epochs
            self.GNSSsystems = GNSSsystems

        self.gnss_systems = list(self.GNSSsystems.values())

        # Read SP3
        sp3_reader = SP3Reader(sp3_file, coords_in_meter=True, desiredGNSSsystems=self.gnss_systems)
        self.sp3_df = sp3_reader.read()
        self.sp3_metadata_dict = sp3_reader.get_metadata()
        self.sp3_epoch_interval = self.sp3_metadata_dict["epoch_interval_sec"]
        self.satcoords = self.interpolate_sp3()

    def interpolate_sp3(self):
        """
        Interpolate the precise satellite coordinates to the
        GNSS observation time stamps.

        """
        sp3_interpol = SP3Interpolator(self.sp3_df, self.sp3_epoch_interval)
        sat_coords = sp3_interpol.interpolate_sat_coordinates(self.time_epochs, self.gnss_systems)
        return sat_coords


    def compute_azimuth_and_elevation(self, receiver_position: Tuple[float, float, float], drop_below_horizon: bool = False):
        """
        Computes the azimuth and elevation angles for all satellites and systems based on satellite and
        receiver ECEF coordinates. Utilizes vectorization for better performance.

        :param receiver_position: Tuple containing receiver ECEF coordinates (x_rec, y_rec, z_rec).
        :param drop_below_horizon: Boolean to drop satellites below the horizon.
        :return: DataFrame with azimuth and elevation angles for all satellites.
        """
        x_rec, y_rec, z_rec = receiver_position

        # WGS 84 ellipsoid parameters
        a = 6378137.0  # Semi-major axis
        b = 6356752.314245  # Semi-minor axis

        # Convert receiver position to geodetic coordinates
        lat_rec, lon_rec, _ = geodf.ECEF2geodb(a, b, x_rec, y_rec, z_rec)

        # Initialize results
        results = []


        # Progress bar setup
        bar_fmt = '{desc}: {percentage:3.0f}%|{bar}| ({n_fmt}/{total_fmt})'
        total_satellites = self.satcoords['Satellite'].nunique()
        with tqdm(total=total_satellites, desc="Computing azimuth and elevation", unit="satellite", bar_format=bar_fmt) as pbar:
            for satellite in self.satcoords['Satellite'].unique():
                sat_data = self.satcoords[self.satcoords['Satellite'] == satellite]

                # Satellite ECEF coordinates
                X = sat_data['X'].to_numpy()
                Y = sat_data['Y'].to_numpy()
                Z = sat_data['Z'].to_numpy()

                # Compute coordinate differences
                dX = X - x_rec
                dY = Y - y_rec
                dZ = Z - z_rec

                # Convert from ECEF to ENU (east, north, up)
                east, north, up = np.vectorize(geodf.ECEF2enu)(lat_rec, lon_rec, dX, dY, dZ)

                # Calculate azimuth angle and correct for quadrants
                azimuth = np.rad2deg(np.arctan(east/north))
                azimuth = np.where((east > 0) & (north < 0) | ((east < 0) & (north < 0)), azimuth + 180, azimuth)
                azimuth = np.where((east < 0) & (north > 0), azimuth + 360, azimuth)


                # Calculate elevation angle
                elevation = np.rad2deg(np.arctan(up / np.sqrt(east**2 + north**2)))

                # Optionally drop satellites below the horizon
                if drop_below_horizon:
                    mask = elevation > 0
                    azimuth = np.where(mask, azimuth, np.nan)
                    elevation = np.where(mask, elevation, np.nan)

                # Store results
                for epoch, az, el in zip(sat_data['Epoch'], azimuth, elevation):
                    results.append({
                        "Epoch": epoch,
                        "Satellite": satellite,
                        "Azimuth": az,
                        "Elevation": el
                    })

                pbar.update(1)

        return pd.DataFrame(results)


    @staticmethod
    def create_satellite_data_dict(interpol_coords_df: pd.DataFrame, azimuth_elevation_df: pd.DataFrame) -> Dict[str, Dict[str, np.ndarray]]:
        """
        Creates a dictionary with interpolated satellite coordinates, azimuth, and elevation data.

        Parameter:
        ---------
        - interpolated_coords_df: DataFrame containing satellite coordinates with columns ['Epoch', 'Satellite', 'X', 'Y', 'Z'].
        - az_el_df: DataFrame containing azimuth and elevation angles with columns ['Epoch', 'Satellite', 'Azimuth', 'Elevation'].

        Return:
        ------

        - Dictionary with structure:
                 {
                     'G': {
                         'coordinates': {prn: numpy_array},
                         'azimuth': numpy_array,
                         'elevation': numpy_array
                     },
                     'R': {...},
                     ...
                 }
        """
        # Create copies of the input DataFrames to prevent modifications to the original data
        interpol_coords_df = interpol_coords_df.copy()
        azimuth_elevation_df = azimuth_elevation_df.copy()

        # Ensure Satellite column is consistent
        interpol_coords_df['PRN'] = interpol_coords_df['Satellite'].str[1:].astype(int)
        interpol_coords_df['GNSS'] = interpol_coords_df['Satellite'].str[0]

        azimuth_elevation_df['PRN'] = azimuth_elevation_df['Satellite'].str[1:].astype(int)
        azimuth_elevation_df['GNSS'] = azimuth_elevation_df['Satellite'].str[0]

        # Normalize Epoch values to range [0, num_epochs - 1]
        unique_epochs = sorted(azimuth_elevation_df['Epoch'].unique())
        epoch_map = {epoch: idx for idx, epoch in enumerate(unique_epochs)}
        azimuth_elevation_df['Normalized_Epoch'] = azimuth_elevation_df['Epoch'].map(epoch_map)

        # Initialize dictionary structure
        satellite_data = {}

        # Get GNSS systems
        gnss_systems = interpol_coords_df['GNSS'].unique()

        # Initialize azimuth and elevation arrays for each GNSS
        num_epochs = len(unique_epochs)

        for gnss in gnss_systems:
            max_prn = 36 if gnss != "C" else 100
            satellite_data[gnss] = {
                'coordinates': {},
                'azimuth': np.full((num_epochs, max_prn + 1), np.nan),
                'elevation': np.full((num_epochs, max_prn + 1), np.nan)
            }

            # Filter GNSS-specific data
            gnss_coords = interpol_coords_df[interpol_coords_df['GNSS'] == gnss]
            gnss_az_el = azimuth_elevation_df[azimuth_elevation_df['GNSS'] == gnss]

            for prn in gnss_coords['PRN'].unique():
                # Get PRN-specific coordinates
                prn_coords = gnss_coords[gnss_coords['PRN'] == prn][['X', 'Y', 'Z']].to_numpy()
                satellite_data[gnss]['coordinates'][str(prn)] = prn_coords

                # Populate azimuth and elevation arrays
                prn_az_el = gnss_az_el[gnss_az_el['PRN'] == prn]
                normalized_epochs = prn_az_el['Normalized_Epoch'].to_numpy()
                satellite_data[gnss]['azimuth'][normalized_epochs, prn] = prn_az_el['Azimuth'].to_numpy()
                satellite_data[gnss]['elevation'][normalized_epochs, prn] = prn_az_el['Elevation'].to_numpy()

        return satellite_data







if __name__=="__main__":
    # RES = PickleHandler.read_zstd_pickle(r"C:\Users\perhe\Desktop\TEST\analysisResults.pkl")
    # rinObs = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx"
    rinObs = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.rnx"
    sp3 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\SP3\Testfile_20220101.eph"
    results_rnav_PRN1 = PickleHandler.read_zstd_pickle(r"C:\Users\perhe\Desktop\TEST BROADCAST\analysisResults.pkl")
    x_rec_approx, y_rec_approx, z_rec_approx = [3149785.9652, 598260.8822, 5495348.4927]

    sats_obj = PreciseSatCoords(sp3, rinObs)
    df_satcoords = sats_obj.satcoords
    df_az_el = sats_obj.compute_azimuth_and_elevation(receiver_position=(x_rec_approx, y_rec_approx, z_rec_approx))
    sat_dict = sats_obj.create_satellite_data_dict(df_satcoords, df_az_el)





