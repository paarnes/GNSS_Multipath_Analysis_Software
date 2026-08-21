import numpy as np
import warnings
from tqdm import tqdm
import pandas as pd
from datetime import datetime
from typing import Optional, List, Tuple, Dict, Union
from gnssmultipath.readers.readRinexObs import readRinexObs, RinexObsData
from gnssmultipath.readers.SP3Reader import SP3Reader
from gnssmultipath.SP3Interpolator import SP3Interpolator
import gnssmultipath.Geodetic_functions as geodf
from gnssmultipath.constants import WGS84_SEMI_MAJOR_AXIS, WGS84_SEMI_MINOR_AXIS
from gnssmultipath.utils.PickleHandler import PickleHandler

warnings.filterwarnings("ignore")





class PreciseSatCoords:
    """
    Class to interpolate precise satellite coordinates (SP3) for each GNSS system.

    Parameters:
    ----------
    sp3_file: str, Path or list of these. Path(s) to the SP3 file(s) with precise satellite coordinates.

    rinex_obs_file: str, Path or RinexObsData. Either a path to a RINEX observation file, or an
                    already parsed ``RinexObsData`` object (from ``readRinexObs``). Passing the parsed
                    object avoids reading the same observation file twice.

    time_epochs: np.ndarray. Observation epochs as [GPS week, time-of-week]. Used when no RINEX
                 observation file/object is given.

    GNSSsystems: List or dict. Systems to interpolate, e.g. ``["G", "E"]`` or ``{1: "G", 2: "E"}``.
                 Defaults to the systems in the RINEX file, or ``["G", "R", "E", "C"]``.

    The receiver position is not needed to interpolate the orbits, only to compute azimuth and
    elevation angles, and is therefore an argument to those methods.

    Example:
    --------
    .. code-block:: python

        rinex = readRinexObs(rinex_obs_file)
        precise = PreciseSatCoords(sp3_file, rinex_obs_file=rinex)
        df_coords = precise.satcoords                       # DataFrame with X, Y, Z and clock bias
        sat_data = precise.compute_satellite_azimut_and_elevation_angle(receiver_position)
    """

    def __init__(self, sp3_file, rinex_obs_file: Union[str, RinexObsData]=None, time_epochs: np.ndarray=None,
                 GNSSsystems=None):

        if rinex_obs_file is not None:
            obs_data = rinex_obs_file if isinstance(rinex_obs_file, RinexObsData) else readRinexObs(rinex_obs_file)
            self.time_epochs = obs_data.time_epochs
            self.GNSSsystems = GNSSsystems if GNSSsystems is not None else obs_data.GNSSsystems
        else:
            if time_epochs is None:
                raise ValueError("Either 'rinex_obs_file' or 'time_epochs' must be provided.")
            self.time_epochs = time_epochs
            self.GNSSsystems = GNSSsystems

        self.gnss_systems = self._normalize_gnss_systems(self.GNSSsystems)

        # Read SP3
        sp3_reader = SP3Reader(sp3_file, coords_in_meter=True, desiredGNSSsystems=self.gnss_systems)
        self.sp3_df = sp3_reader.read()
        self.sp3_metadata_dict = sp3_reader.get_metadata()
        self.sp3_epoch_interval = self.sp3_metadata_dict["epoch_interval_sec"]
        self.satcoords = self.interpolate_sp3()

    @staticmethod
    def _normalize_gnss_systems(gnss_systems) -> List[str]:
        """Accept a dict ({1: 'G'}), a list/tuple/set of system codes, a single code or None."""
        if gnss_systems is None:
            return ["G", "R", "E", "C"]
        if isinstance(gnss_systems, dict):
            gnss_systems = list(gnss_systems.values())
        elif isinstance(gnss_systems, str):
            gnss_systems = [gnss_systems]
        return [str(sys_code).strip().upper()[0] for sys_code in gnss_systems]

    def interpolate_sp3(self):
        """
        Interpolate the precise satellite coordinates to the
        GNSS observation time stamps.

        """
        sp3_interpol = SP3Interpolator(self.sp3_df, self.sp3_epoch_interval)
        sat_coords = sp3_interpol.interpolate_sat_coordinates(self.time_epochs, self.gnss_systems)
        return sat_coords


    def compute_satellite_azimut_and_elevation_angle(self,
                                                     receiver_position: Tuple[float, float, float],
                                                     drop_below_horizon: bool = False) -> Dict[str, Dict[str, np.ndarray]]:
        """
        Computes azimuth and elevation angles and returns them per GNSS system, using the same
        structure as ``SatelliteEphemerisToECEF.compute_satellite_azimut_and_elevation_angle``.
        Both the azimuth and elevation arrays have shape (n_epochs, max_PRN + 1) and are indexed
        by PRN number, which makes them directly usable in ``make_skyplot``.

        :param receiver_position: Receiver ECEF coordinates (x_rec, y_rec, z_rec).
        :param drop_below_horizon: Boolean to drop satellites below the horizon.
        :return: Dict on the form {'G': {'position': {PRN: array}, 'azimuth': array, 'elevation': array}, ...}
        """
        az_el_df = self.compute_azimuth_and_elevation(receiver_position, drop_below_horizon=drop_below_horizon)
        return self.create_satellite_data_dict(self.satcoords, az_el_df)


    def compute_azimuth_and_elevation(self, receiver_position: Tuple[float, float, float], drop_below_horizon: bool = False):
        """
        Computes the azimuth and elevation angles for all satellites and systems based on satellite and
        receiver ECEF coordinates. Utilizes vectorization for better performance.

        :param receiver_position: Tuple containing receiver ECEF coordinates (x_rec, y_rec, z_rec).
        :param drop_below_horizon: Boolean to drop satellites below the horizon.
        :return: DataFrame with azimuth and elevation angles for all satellites.
        """
        x_rec, y_rec, z_rec = receiver_position

        # Convert receiver position to geodetic coordinates
        lat_rec, lon_rec, _ = geodf.ECEF2geodb(WGS84_SEMI_MAJOR_AXIS, WGS84_SEMI_MINOR_AXIS, x_rec, y_rec, z_rec)

        # Initialize results as column arrays
        all_epochs = []
        all_sats = []
        all_az = []
        all_el = []


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
                east, north, up = geodf.ECEF2enu_batch(lat_rec, lon_rec, dX, dY, dZ)

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

                # Store results as arrays
                all_epochs.append(sat_data['Epoch'].to_numpy())
                all_sats.append(np.full(len(sat_data), satellite))
                all_az.append(azimuth)
                all_el.append(elevation)

                pbar.update(1)

        return pd.DataFrame({
            "Epoch": np.concatenate(all_epochs),
            "Satellite": np.concatenate(all_sats),
            "Azimuth": np.concatenate(all_az),
            "Elevation": np.concatenate(all_el),
        })


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
                         'position': {prn: numpy_array},
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
                'position': {},
                'azimuth': np.full((num_epochs, max_prn + 1), np.nan),
                'elevation': np.full((num_epochs, max_prn + 1), np.nan)
            }

            # Filter GNSS-specific data
            gnss_coords = interpol_coords_df[interpol_coords_df['GNSS'] == gnss]
            gnss_az_el = azimuth_elevation_df[azimuth_elevation_df['GNSS'] == gnss]

            for prn in gnss_coords['PRN'].unique():
                # Get PRN-specific coordinates
                prn_coords = gnss_coords[gnss_coords['PRN'] == prn][['X', 'Y', 'Z']].to_numpy()
                satellite_data[gnss]['position'][str(prn)] = prn_coords

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





