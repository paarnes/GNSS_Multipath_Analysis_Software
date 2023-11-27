"""
Module for satellite elevations and azimut angles based on broadcasted ephemerides.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""

from gnssmultipath.SatelliteEphemerisToECEF import SatelliteEphemerisToECEF

def computeSatElevAzimuth_fromNav(rinex_nav_file, approxPosition, GNSS_SVs, time_epochs, nav_data_rate):
    """
    A function for computing satellite elevations and azimut angles based on
    broadcasted ephemerides. Support all global navigation systems (GPS,GLONASS,Galileo & BeiDou).

    Input:
    -----
    navigationFile: path to navigation file or list of paths to navigations files
    approxPosition : Reciever position in ECEF from RINEX obseration file
    GNSS_SVs : dictionary containing overview of observed satellites for each epoch and for each systems
                {'G': array([[11., 30., 15., ...,  0.,  0.,  0.],
                             [11., 30., 15., ...,  0.,  0.,  0.],
                             [ 9., 31., 19., ...,  0.,  0.,  0.]]),
                 'R': array([[ 8.,  8., 15., ...,  0.,  0.,  0.],
                             [ 8.,  8., 15., ...,  0.,  0.,  0.],
                             [ 8.,  8., 15., ...,  0.,  0.,  0.],
                             [ 8.,  9.,  2., ...,  0.,  0.,  0.]])}
    time_epochs : numpy array with GPS time for GNSS observations (week, time-of-week)  
                  np.array([[  2190.        , 518399.99999988],
                            [  2190.        , 518429.99999988],
                            [  2190.        , 518459.99999988]])      
    nav_data_rate : desired data rate for the broadcasted ephemerides (when reading the rinex nav file)

    Return:
    ------
    sat_pos : Dictionary containing satellite posistion in ECEF, in addition to elavation and azimuth angels.
    glo_fcn : Dictionary containing GLONASS frequency channel numbers extracted from RINEX navigation file
    """

    sat_pos = {}
    systems = GNSS_SVs.keys()
    x_rec, y_rec, z_rec = approxPosition.flatten() # Extracting approx postion from RINEX obs-file
    # Initialize a converter object to convert the broadcast ephemerides to ECEF and interpolate to "time_epochs"
    CONVERTER = SatelliteEphemerisToECEF(rinex_nav_file, x_rec, y_rec, z_rec, desired_systems=systems, data_rate=nav_data_rate)
    sat_pos = CONVERTER.get_sat_ecef_coordinates(time_epochs[:,1])
    sat_pos = CONVERTER.compute_satellite_azimut_and_elevation_angle(drop_below_horizon = True)
    glo_fcn = CONVERTER.glo_fcn

    return sat_pos, glo_fcn
