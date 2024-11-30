import warnings
import logging
import numpy as np
from typing import Literal
from gnssmultipath.Geodetic_functions import date2gpstime, date2gpstime_vectorized
from gnssmultipath.SatelliteEphemerisToECEF import SatelliteEphemerisToECEF, Kepler2ECEF
from gnssmultipath.readRinexObs import readRinexObs
import gnssmultipath.Geodetic_functions as geof
from pyproj import Transformer
from StatisticalAnalysis import StatisticalAnalysis
from typing import Tuple, Dict


warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)

c = 299792458  # Speed of light [m/s]





class GNSSPositionEstimator:
    """
    A class to estimate the position of a GNSS receiver using pseudoranges and ephemeris data.

    Parameters:
    ----------
    rinex_obs_file: str. Path to the RINEX observation file containing pseudorange measurements.
    rinex_nav_file: str. Path to the RINEX navigation file containing satellite ephemeris data.
    desired_time: np.ndarray. Desired time for position estimation as a 1D array in the format
                  [year, month, day, hour, minute, second].
    desired_system: Literal["G", "E", "C"], optional. GNSS system to use for position estimation.
                    Default is "G" (GPS).
                    - "G" = GPS
                    - "E" = Galileo
                    - "C" = BeiDou
    x_rec_approx: float, optional. Approximate receiver x-coordinate in ECEF (meters). Default is 0.0.
    y_rec_approx: float, optional. Approximate receiver y-coordinate in ECEF (meters). Default is 0.0.
    z_rec_approx: float, optional. Approximate receiver z-coordinate in ECEF (meters). Default is 0.0.

    Attributes:
    ----------
    desired_time: np.ndarray. Desired GPS time for position estimation, converted to GPS week and seconds.
    desired_sys: str. Selected GNSS system for position estimation.
    navdata: SatelliteEphemerisToECEF. Object containing satellite ephemeris data for the chosen system.
    GNSS_obs: dict. Pseudorange observations for all GNSS systems.
    time_epochs: np.ndarray. Array of observation epochs (GPS time).
    GNSSsystems: dict. Mapping of GNSS systems to indices.
    obsCodes: dict. Observation codes for each GNSS system.
    x_rec_approx: float. Approximate receiver x-coordinate in ECEF (meters).
    y_rec_approx: float. Approximate receiver y-coordinate in ECEF (meters).
    z_rec_approx: float. Approximate receiver z-coordinate in ECEF (meters).
    gnss_index_mapper: dict. Mapping of GNSS system identifiers to their corresponding indices.
    sys_idx: int. Index of the selected GNSS system in the observation data.
    sys_obs_codes: list. Observation codes for the selected GNSS system.
    signal_idx: int. Index of the first valid pseudorange signal in the observation codes.
    """

    def __init__(self, rinex_obs_file:str, rinex_nav_file:str, desired_time:np.array,
                 desired_system:Literal["G","E","C"]="G",
                 x_rec_approx:float=0.0, y_rec_approx:float=0.0, z_rec_approx:float=0.0):


        self.desired_time = np.atleast_2d(date2gpstime_vectorized(desired_time)).T
        self.desired_sys = desired_system
        self.navdata = SatelliteEphemerisToECEF(rinex_nav_file, x_rec_approx, y_rec_approx, z_rec_approx, self.desired_sys)
        self.GNSS_obs, _, _, _, self.time_epochs,_, self.GNSSsystems,\
              self.obsCodes, _, _, _, _, _, _, _, _, _,\
              _, _, _, _, _, _, _, _ = readRinexObs(rinex_obs_file)

        self.x_rec_approx = x_rec_approx
        self.y_rec_approx = y_rec_approx
        self.z_rec_approx = z_rec_approx

        self.__gnss_index_mapper_dict = self.__gnss_index_mapper()
        self.sys_idx = self.__gnss_index_mapper_dict[self.desired_sys]
        self.sys_obs_codes = self.obsCodes[self.sys_idx][self.desired_sys]
        self.signal_idx = next((i for i, code in enumerate(self.sys_obs_codes) if code.startswith("C")), -1)



    def __gnss_index_mapper(self) -> dict:
        """
        Map GNSS system identifiers to their corresponding indices.

        Returns:
        -------
        mapper: dict. Dictionary mapping GNSS system identifiers to indices.
        """
        mapper = {}
        for i in np.arange(0, len(self.GNSSsystems)):
            GNSSsystemIndex = list(self.GNSSsystems.keys())[i]
            mapper[self.GNSSsystems[GNSSsystemIndex]]= GNSSsystemIndex
        return mapper


    @staticmethod
    def __find_idx_of_epoch_closest_in_time(desired_time: np.ndarray, obs_time_epochs: np.ndarray) -> int:
        """
        Find the index of the epoch in `obs_time_epochs` where the difference in time is the smallest.

        Parameters:
        ----------
        desired_time: (numpy.ndarray). Desired GPS time as a 2D array [week, seconds].
        obs_time_epochs: (numpy.ndarray). Array of observed time epochs as a 2D array [week, seconds].

        Returns:
        -------
        closest_index: int: Index of the epoch with the smallest time difference.
        """
        # Compute the absolute difference between desired time and observation epochs
        week_diff = np.abs(obs_time_epochs[:, 0] - desired_time[0, 0])
        sec_diff = np.abs(obs_time_epochs[:, 1] - desired_time[1, 0])

        # Convert time differences to total seconds
        total_time_diff = week_diff * 604800 + sec_diff  # 604800 seconds in a GPS week

        # Find the index of the smallest difference
        closest_index = np.argmin(total_time_diff) + 1 # the key of the dict that contains pseudoranges starts on 1
        return closest_index


    @staticmethod
    def __compute_toc_for_single_ephemeris(ephemeris: np.ndarray):
        """
        Compute the time of clock (toc) for a single satellite.

        Parameters:
        -----------
        ephemeris : (numpy.ndarray). Ephemeris data for a single satellite.
                    The input must include the year, month, day, hour, minute
                    and second at positions [1, 2, 3, 4, 5, 6] respectively.

        Returns:
        --------
        toc: float. The computed `toc` value (time of clock correction).
        """

        year = int(ephemeris[1])  # Year from ephemeris
        month = int(ephemeris[2])   # Month
        day = int(ephemeris[3])       # Day
        hour = int(ephemeris[4])      # Hour
        minute = int(ephemeris[5])    # Minute
        second = float(ephemeris[6])  # Second

        # Convert date and time to GPS time
        _, toc = date2gpstime(year, month, day, hour, minute, second)

        return toc


    def __extract_pseudoranges_and_satellites(self, idx_of_closest_eph):
        """
        Extract valid pseudoranges and corresponding satellite indices.

        Parameters:
        ----------
        idx_of_closest_eph: int. Index of the closest ephemeris.

        Returns:
        -------
        sat_nr: numpy.ndarray. Array of satellite indices with valid pseudoranges.
        rji: numpy.ndarray. Array of valid pseudoranges for the satellites.
        """
        pseudoranges = self.GNSS_obs[self.desired_sys][idx_of_closest_eph]
        sat_nr = np.where(pseudoranges[:, self.signal_idx] != 0.0)[0]
        rji = pseudoranges[sat_nr, self.signal_idx]
        return sat_nr, rji



    def __estimate_without_low_satellites(
        self,
        X: np.ndarray,
        Y: np.ndarray,
        Z: np.ndarray,
        dT_rel: np.ndarray,
        Rji: np.ndarray,
        sat_nr: list,
        x: float,
        y: float,
        z: float,
        t: float,
        dTi0: float
    )  -> Tuple[float, float, float, float, np.ndarray, np.ndarray, np.ndarray, list]:
        """
        Filter out satellites with an elevation angle below 15 degrees and recompute the position and clock bias.
        
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
        azimuths, elevations = self.navdata.compute_azimuth_and_elevation_static(X, Y, Z, x, y, z)
        
        # Filter satellites with elevation < 15 degrees
        active_sat_indices = np.where(elevations >= 15)[0]
        low_sat_indices = np.where(elevations < 15)[0]
        for i in low_sat_indices:
            PRN = f"{self.desired_sys}{str(sat_nr[i]).zfill(2)}"
            logger.info(f"\nINFO(GNSSPositionEstimator): {PRN} is exluded due to low elevation: {np.round(elevations[i],3)}°")
            
    
        # Update satellite-related data based on active satellites
        X = X[active_sat_indices]
        Y = Y[active_sat_indices]
        Z = Z[active_sat_indices]
        dT_rel = dT_rel[active_sat_indices]
        Rji = Rji[active_sat_indices]
        sat_nr = [sat_nr[i] for i in active_sat_indices]
        
        # Recompute design matrix and observation corrections for active satellites
        diff = np.column_stack([X - x, Y - y, Z - z])  # Shape (num_sats, 3)
        rho = np.linalg.norm(diff, axis=1)
        dax = -diff[:, 0] / rho
        day = -diff[:, 1] / rho
        daz = -diff[:, 2] / rho
        dadT = np.ones_like(rho)
    
        l = Rji + c * dT_rel - (rho + c * dTi0)
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
    
        # return X, Y, Z, dT_rel, sat_nr, x, y, z, dTi0, A,l,N,h
        return x, y, z, dTi0, A,l,N,h, sat_nr


    def estimate_position(self, max_iterations:int = 15, tol:float = 1e-8) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
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

        # sat_nr = [32, 6, 25, 12, 14, 17, 19, 24]
        # rji = np.array([23033856.605, 24571877.822, 23668552.454, 20673456.591,
        #                   24594925.874, 22780048.894, 21428473.569, 20562322.798])
        
    
        num_sats = len(sat_nr)
        if num_sats < 4:
            raise ValueError("Insufficient satellites for position estimation (minimum 4 required).")

        # Initialize variables
        dTi0 = 0.0 # Initial clock bias
        improvement = 10
        iteration = 0
        t = self.desired_time[1][0]


        # Preallocate arrays for satellite data
        X, Y, Z, dT_rel = np.zeros(num_sats), np.zeros(num_sats), np.zeros(num_sats), np.zeros(num_sats)
        dTj, Tj_GPS, Rji = np.zeros(num_sats), np.zeros(num_sats), np.zeros(num_sats)


        while improvement > tol and iteration < max_iterations:
            for i, sat_num in enumerate(sat_nr):
                PRN = f"{self.desired_sys}{str(sat_num).zfill(2)}"
                efemeris = self.navdata.get_closest_ephemerides_for_PRN_at_time(PRN, self.desired_time[1])
                toc = self.__compute_toc_for_single_ephemeris(efemeris[0])
                # Extract satellite clock parameters
                a0 = float(efemeris[0][7]) # SV clock bias (seconds)
                a1 = float(efemeris[0][8]) # SV clock drift (sec/sec)
                a2 = float(efemeris[0][9]) # SV clock drift rate (sec/sec2)
                T_GD = float(efemeris[0][32]) # TGD (seconds)

                # Calculate satellite clock error
                dTj[i] = a0 + a1 * (t - toc) + a2 * (t - toc) ** 2 - T_GD

                # Correct observed pseudorange
                Rji[i] = rji[i] + dTj[i] * c

                # Compute time when signal left the satellite
                Tj_GPS[i] = t - Rji[i] / c

                # Compute satellite ECEF position and relativistic clock correction
                X_, Y_, Z_, dT_rel_ = Kepler2ECEF(x, y, z).kepler2ecef(efemeris, Tj_GPS[i])
                X[i], Y[i], Z[i], dT_rel[i]  = X_[0], Y_[0], Z_[0], dT_rel_[0]
                
                # GLONASS
                # glo_1 = self.navdata.get_sat_ecef_coordinates(desired_time=np.array([Tj_GPS[i]]), PRN=PRN)
                


            # Compute the vector differences between receiver and satellite positions
            diff = np.column_stack([X - x, Y - y, Z - z])  # Shape (num_sats, 3)

            # Compute distances and partial derivatives
            rho = np.linalg.norm(diff, axis=1)
            dax = -diff[:, 0] / rho
            day = -diff[:, 1] / rho
            daz = -diff[:, 2] / rho
            dadT = np.ones_like(rho)

            # Observation corrections
            l = Rji + c * dT_rel - (rho + c * dTi0)
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

            # Update iteration info
            iteration += 1
            logger.info(f"Iteration {iteration}: dx = {dx}")
            improvement = np.max(np.abs(dx))
            
            
        
        # Try to filter low satellites and update position
        try:
            x, y, z, dTi0, A, l, N, h, sat_nr = self.__estimate_without_low_satellites(X, Y, Z, dT_rel, Rji, sat_nr, x, y, z, t, dTi0)
        except Exception as e:
            logger.warning(
                f"WARNING (GNSS_MultipathAnalysis): Failed to perform the final position estimation after removing low-elevation satellites. "
                f"Error: {str(e)}"
            )
            pass
           
        stats_report = StatisticalAnalysis(A,l,N,h).run_statistical_analysis()

        return np.array([x, y, z, dTi0]), stats_report










if __name__=="__main__":

    # rinObs = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx"
    rinObs = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.rnx"
    rinObs = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped – Kopi.rnx"
    rinNav = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\BRDC00IGS_R_20220010000_01D_MN.rnx"


    fasit = np.array([3149785.9652, 598260.8822, 5495348.4927])
    x_rec, y_rec, z_rec = fasit.T
    # desired_time = np.array([2022, 1, 1, 2, 20, 0.0000000])
    desired_time = np.array([2022, 1, 1, 0, 0, 30.0000000])
    desired_system = "G"


    ### TEST
    # rinNav = r"C:\Users\perhe\OneDrive\Documents\Python_skript\Approx_rec_pos\BRDM00DLR_S_20191550000_01D_MN.rnx"
    # fasit = np.array([3173400.0, 605500.0, 5481000.0])
    # # fasit = np.array([3173289, 605742, 5481217])
    # x_rec, y_rec, z_rec = fasit.T
    # desired_time = np.array([2019, 6, 4, 7, 30, 0.0000000])

    GNSSPos = GNSSPositionEstimator(rinObs, rinNav, desired_time, desired_system, x_rec_approx=x_rec, y_rec_approx=y_rec, z_rec_approx=z_rec)
    # GNSSPos = GNSSPositionEstimator(rinObs, rinNav, desired_time, desired_system, x_rec_approx=0, y_rec_approx=0, z_rec_approx=0)
    estimated_position, stats = GNSSPos.estimate_position()
    X,Y,Z,dT = estimated_position.T
    # fasit = np.array([3173289, 605742, 5481217])


    # lon, lat, h = Transformer.from_crs("EPSG:4936", "EPSG:4258", always_xy=True).transform(X,Y,Z)
    # estimated__etrs89_nn2000 = np.array(Transformer.from_crs("EPSG:4936", "EPSG:5972").transform(X,Y,Z))
    # fasit_etrs89_nn2000 = Transformer.from_crs("EPSG:4936", "EPSG:5972").transform(*fasit.T)

    print(X,Y,Z,dT)
    print(f"\nDifference: {np.round(fasit- np.array([X,Y,Z]),3)}")
    # print(f"\nDifference in UTM Z32: {np.round(fasit_etrs89_nn2000 - estimated__etrs89_nn2000,3)}")
    # print(stats)


    # geof.ECEF2geodb(a=6378137, b=6356752.314245, X=X,Y=Y, Z=Z)
    # my_func = np.array(geof.ECEF2geodb(a=6378137, b=6356752.314245, X=X,Y=Y,Z=Z))
