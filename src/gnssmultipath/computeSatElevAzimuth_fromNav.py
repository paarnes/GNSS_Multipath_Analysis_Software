"""
Module for satellite elevations and azimut angles based on broadcasted ephemerides.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""
from typing import Dict, Tuple, Union
from gnssmultipath.SatelliteEphemerisToECEF import SatelliteEphemerisToECEF
from gnssmultipath.BroadNavPositionEstimator import BroadNavPositionEstimator
from gnssmultipath.Geodetic_functions import gpstime2date
import numpy as np

def computeSatElevAzimuth_fromNav(rinex_nav_file: str, approxPosition: np.ndarray,
                                  GNSS_SVs: Dict, time_epochs:  np.ndarray,
                                  nav_data_rate: int, GNSS_obs: Dict,
                                  GNSSsystems: Dict, obsCodes: Dict) -> Tuple[Dict, Dict, Union[np.ndarray, None], Union[Dict, None]]:
    """
    A function for computing satellite elevations and azimut angles based on
    broadcasted ephemerides. Support all global navigation systems (GPS,GLONASS,Galileo & BeiDou).

    Parameter:
    -----
    - navigationFile: path to navigation file or list of paths to navigations files
    - approxPosition : Reciever position in ECEF from RINEX obseration file
    - GNSS_SVs : dictionary containing overview of observed satellites for each epoch and for each systems
                {'G': array([[11., 30., 15., ...,  0.,  0.,  0.],
                             [11., 30., 15., ...,  0.,  0.,  0.],
                             [ 9., 31., 19., ...,  0.,  0.,  0.]]),
                 'R': array([[ 8.,  8., 15., ...,  0.,  0.,  0.],
                             [ 8.,  8., 15., ...,  0.,  0.,  0.],
                             [ 8.,  8., 15., ...,  0.,  0.,  0.],
                             [ 8.,  9.,  2., ...,  0.,  0.,  0.]])}
    - time_epochs : numpy array with GPS time for GNSS observations (week, time-of-week)
                  np.array([[  2190.        , 518399.99999988],
                            [  2190.        , 518429.99999988],
                            [  2190.        , 518459.99999988]])
    - nav_data_rate : int. desired data rate for the broadcasted ephemerides (when reading the rinex nav file)
    - GNSS_obs: Dict with GNSS observation from RINEX observation file.
    - GNSSsystems: Dict Index/System code mapper (ex {1: 'G', 2: 'R', 3: 'E', 4: 'C'}).
    - obsCodes: Dict with all the observation codes for each system.

    Return:
    ------
    - sat_pos : Dictionary containing satellite posistion in ECEF, in addition to elavation and azimuth angels.
    - glo_fcn : Dictionary containing GLONASS frequency channel numbers extracted from RINEX navigation file
    - estimated_position: np.ndarray. Estimated approximate position of the receiver based on pseudoranges.
    - stats: Dictionary with statistical measures for the estimated approximate position of the receiver.
    """

    sat_pos = {}
    systems = GNSS_SVs.keys()
    estimated_position, stats = None, None

    x_rec, y_rec, z_rec = np.atleast_2d(approxPosition).flatten() # Extracting approx postion from RINEX obs-file
    navObj = SatelliteEphemerisToECEF(rinex_nav_file, x_rec, y_rec, z_rec, desired_systems=systems, data_rate=nav_data_rate)
    if all(coord == 0 for coord in [x_rec, y_rec, z_rec]):
        desired_time = np.array(gpstime2date(time_epochs[0,0], round(time_epochs[0,1],6)))
        GNSSPos = BroadNavPositionEstimator(desired_time=desired_time,
                                        navdata=navObj,
                                        GNSS_obs=GNSS_obs,
                                        time_epochs=time_epochs,
                                        GNSSsystems=GNSSsystems,
                                        obsCodes=obsCodes)

        estimated_position, stats = GNSSPos.estimate_position()

    # Initialize a converter object to convert the broadcast ephemerides to ECEF and interpolate to "time_epochs"
    sat_pos = navObj.get_sat_ecef_coordinates(time_epochs[:,1])
    sat_pos = navObj.compute_satellite_azimut_and_elevation_angle(drop_below_horizon = True)
    glo_fcn = navObj.glo_fcn

    return sat_pos, glo_fcn, estimated_position, stats
