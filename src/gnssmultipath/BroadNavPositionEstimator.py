import warnings
import logging
import numpy as np
from typing import Literal
from gnssmultipath.Geodetic_functions import date2gpstime, date2gpstime_vectorized
from gnssmultipath.SatelliteEphemerisToECEF import SatelliteEphemerisToECEF, Kepler2ECEF
from gnssmultipath.readRinexObs import readRinexObs
from gnssmultipath.StatisticalAnalysis import StatisticalAnalysis
from typing import Tuple, Dict, Optional, Union


warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)

c = 299792458  # Speed of light [m/s]



class BroadNavPositionEstimator:
    """
    A class to estimate the position of a GNSS receiver using pseudoranges and ephemeris data.

    Example:
    --------
    .. code-block:: python

        rinObs = "path/to/rinex_obs_file.rnx"
        rinNav = "path/to/rinex_nav_file.rnx"
        desired_time = np.array([2024, 1, 1, 12, 0, 0])  # Example date and time
        desired_system = "G"  # Use GPS
        GNSSPos = BroadNavPositionEstimator(rinObs, rinNav, desired_time, desired_system)
        estimated_position, stats = GNSSPos.estimate_position()



    Parameter:
    ----------
    - rinex_obs_file: str. Path to the RINEX observation file containing pseudorange measurements.
    - rinex_nav_file: str. Path to the RINEX navigation file containing satellite ephemeris data.
    - desired_time: np.ndarray. Desired time for position estimation as a 1D array in the format
                  [year, month, day, hour, minute, second].
    - desired_system: Literal["G", "E", "R"], optional. GNSS system to use for position estimation.
                    Default is "G" (GPS).
                    - "G" = GPS
                    - "E" = Galileo
                    - "R" = GLONASS
    - x_rec_approx: float, optional. Approximate receiver x-coordinate in ECEF (meters). Default is 0.0.
    - y_rec_approx: float, optional. Approximate receiver y-coordinate in ECEF (meters). Default is 0.0.
    - z_rec_approx: float, optional. Approximate receiver z-coordinate in ECEF (meters). Default is 0.0.

    Attribute:
    ----------
    - desired_time: np.ndarray. Desired GPS time for position estimation, converted to GPS week and seconds.
    - desired_sys: str. Selected GNSS system for position estimation.
    - navdata: SatelliteEphemerisToECEF. Object containing satellite ephemeris data for the chosen system.
    - GNSS_obs: dict. Pseudorange observations for all GNSS systems.
    - time_epochs: np.ndarray. Array of observation epochs (GPS time).
    - GNSSsystems: dict. Mapping of GNSS systems to indices.
    - obsCodes: dict. Observation codes for each GNSS system.
    - x_rec_approx: float. Approximate receiver x-coordinate in ECEF (meters).
    - y_rec_approx: float. Approximate receiver y-coordinate in ECEF (meters).
    - z_rec_approx: float. Approximate receiver z-coordinate in ECEF (meters).
    - gnss_index_mapper: dict. Mapping of GNSS system identifiers to their corresponding indices.
    - sys_idx: int. Index of the selected GNSS system in the observation data.
    - sys_obs_codes: list. Observation codes for the selected GNSS system.
    - signal_idx: int. Index of the first valid pseudorange signal in the observation codes.
    """

    def __init__(
        self,
        rinex_obs_file: Union[str, None] = None,
        rinex_nav_file: Union[str, None] = None,
        desired_time: np.ndarray = None,
        desired_system: Literal["G", "E", "R"] = "G",
        x_rec_approx: float = 0.0,
        y_rec_approx: float = 0.0,
        z_rec_approx: float = 0.0,
        navdata: Optional[SatelliteEphemerisToECEF] = None,
        GNSS_obs: Optional[dict] = None,
        time_epochs: Optional[np.ndarray] = None,
        GNSSsystems: Optional[dict] = None,
        obsCodes: Optional[dict] = None,
        elevation_cut_off_angle: Optional[int] = 10,
    ):
        if "C" in desired_system:
            raise ValueError(
                "Using BeiDou to compute approximate position is not supported yet. Please choose another system.")

        # Desired time in GPS format
        self.desired_time = np.atleast_2d(
            date2gpstime_vectorized(desired_time)).T
        self.desired_sys = desired_system
        self.elevation_cut_off_angle = elevation_cut_off_angle

        # Handle satellite ephemeris data
        if navdata:
            self.navdata = navdata
        elif rinex_nav_file:
            self.navdata = SatelliteEphemerisToECEF(
                rinex_nav_file, x_rec_approx, y_rec_approx, z_rec_approx, self.desired_sys)
        else:
            raise ValueError(
                "Either 'navdata' or 'rinex_nav_file' must be provided.")

        # Handle GNSS observation data
        if all(var is not None for var in [GNSS_obs, time_epochs, GNSSsystems, obsCodes]):
            self.GNSS_obs = GNSS_obs
            self.time_epochs = time_epochs
            self.GNSSsystems = GNSSsystems
            self.obsCodes = obsCodes
        elif rinex_obs_file:
            self.GNSS_obs, _, _, _, self.time_epochs, _, self.GNSSsystems, \
                self.obsCodes, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ = readRinexObs(
                    rinex_obs_file)
        else:
            raise ValueError(
                "Either GNSS observation arrays or 'rinex_obs_file' must be provided.")

        # Receiver approximate coordinates
        self.x_rec_approx = x_rec_approx
        self.y_rec_approx = y_rec_approx
        self.z_rec_approx = z_rec_approx

        # GNSS system and signal indexing
        self.__gnss_index_mapper_dict = self.__gnss_index_mapper()
        self.sys_idx = self.__gnss_index_mapper_dict[self.desired_sys]
        self.sys_obs_codes = self.obsCodes[self.sys_idx][self.desired_sys]
        self.signal_idx = next((i for i, code in enumerate(
            self.sys_obs_codes) if code.startswith("C")), -1)

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
            mapper[self.GNSSsystems[GNSSsystemIndex]] = GNSSsystemIndex
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
        # the key of the dict that contains pseudoranges starts on 1
        closest_index = np.argmin(total_time_diff) + 1
        return closest_index

    @staticmethod
    def __get_gpstime_for_single_ephemeris(ephemeris: np.ndarray):
        """
        Compute the gps time (gps week and time of week) from
        the ephemerid for a single satellite.

        Parameters:
        -----------
        ephemeris : (numpy.ndarray). Ephemeris data for a single satellite.
                    The input must include the year, month, day, hour, minute
                    and second at positions [1, 2, 3, 4, 5, 6] respectively.

        Returns:
        --------
        tow: float. The computed `week` value (GPS week).
        tow: float. The computed `tow` value (time of week).
        """

        year = int(ephemeris[1])  # Year from ephemeris
        month = int(ephemeris[2])   # Month
        day = int(ephemeris[3])       # Day
        hour = int(ephemeris[4])      # Hour
        minute = int(ephemeris[5])    # Minute
        second = float(ephemeris[6])  # Second

        # Convert date and time to GPS time
        week, tow = date2gpstime(year, month, day, hour, minute, second)

        return week, tow

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
    ) -> Tuple[float, float, float, float, np.ndarray, np.ndarray, np.ndarray, list]:
        """
        Filter out satellites with an elevation angle below the specified threshold (like 10 degrees for intance) and recompute the position and clock bias.

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
        azimuths, elevations = self.navdata.compute_azimuth_and_elevation_static(
            X, Y, Z, x, y, z)

        # Filter satellites with elevation < 15 degrees for example
        active_sat_indices = np.where(
            elevations >= self.elevation_cut_off_angle)[0]
        low_sat_indices = np.where(
            elevations < self.elevation_cut_off_angle)[0]
        for i in low_sat_indices:
            PRN = f"{self.desired_sys}{str(sat_nr[i]).zfill(2)}"
            logger.info(
                f"\nINFO(BroadNavPositionEstimator): {PRN} is exluded due to low elevation: {np.round(elevations[i], 3)}°")

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
        return x, y, z, dTi0, A, l, N, h, sat_nr

    def __compute_satellite_clock_error_polynomial(self, efemeris, t, toc) -> float:
        """
        Compute the satellite clock error for GPS, GALILEO, and BeiDou based on polynomial parameters.

        Parameters:
        -----------
        efemeris: Ephemeris data for the satellite.
        t: Current time.
        toc: Ephemeris time of clock.

        Returns:
        --------
        dTj: Satellite clock error (seconds).
        """
        a0 = float(efemeris[0][7])  # SV clock bias (seconds)
        a1 = float(efemeris[0][8])  # SV clock drift (sec/sec)
        a2 = float(efemeris[0][9])  # SV clock drift rate (sec/sec2)
        T_GD = float(efemeris[0][32])  # TGD (seconds)

        # Calculate satellite clock error
        dTj = a0 + a1 * (t - toc) + a2 * (t - toc) ** 2 - T_GD
        return dTj

    def update_receiver_coordinates(self, x, y, z):
        """Update receiver coordinates in the class attributes."""
        self.navdata.x_rec = x
        self.navdata.y_rec = y
        self.navdata.z_rec = z

    def estimate_position(self, max_iterations: int = 15, tol: float = 1e-8) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
        """
        Estimate the receiver's position using a vectorized least-squares algorithm.

        Parameter:
        ----------
        - max_iterations: int, optional. Maximum number of iterations for convergence (default is 15).
        - tol: float, optional. Convergence tolerance for position updates (default is 1e-8).

        Return:
        -------
        - position: numpy.ndarray. Estimated receiver position in ECEF coordinates (x, y, z) and clock bias (dt).
        """
        x, y, z = self.x_rec_approx, self.y_rec_approx, self.z_rec_approx
        idx_of_closest_eph = self.__find_idx_of_epoch_closest_in_time(
            self.desired_time, self.time_epochs)
        sat_nr, rji = self.__extract_pseudoranges_and_satellites(
            idx_of_closest_eph)

        num_sats = len(sat_nr)
        if num_sats < 4:
            raise ValueError(
                "Insufficient satellites for position estimation (minimum 4 required).")

        # Initialize variables
        dTi0 = 0.0  # Initial clock bias
        improvement = 10
        iteration = 0
        t = self.desired_time[1][0]

        # Preallocate arrays for satellite data
        X, Y, Z, dT_rel = np.zeros(num_sats), np.zeros(
            num_sats), np.zeros(num_sats), np.zeros(num_sats)
        dTj, Tj_GPS, Rji = np.zeros(num_sats), np.zeros(
            num_sats), np.zeros(num_sats)

        while improvement > tol and iteration < max_iterations:
            for i, sat_num in enumerate(sat_nr):
                PRN = f"{self.desired_sys}{str(sat_num).zfill(2)}"
                efemeris = self.navdata.get_closest_ephemerides_for_PRN_at_time(
                    PRN, self.desired_time[1])
                _, tow = self.__get_gpstime_for_single_ephemeris(efemeris[0])

                # Select the correct clock error computation method based on satellite system
                if self.desired_sys in ['G', 'E']:  # GPS, GALILEO, BeiDou
                    dTj[i] = self.__compute_satellite_clock_error_polynomial(
                        efemeris, t, tow)

                    # Correct observed pseudorange
                    Rji[i] = rji[i] + dTj[i] * c

                    # Compute time when signal left the satellite
                    Tj_GPS[i] = t - Rji[i] / c

                    # Compute satellite ECEF position and relativistic clock correction
                    # X[i], Y[i], Z[i], dT_rel[i] = self.navdata.get_sat_ecef_coordinates(
                    #     desired_time=np.array([Tj_GPS[i]]), PRN=PRN)
                    desired_time = np.atleast_1d(Tj_GPS[i]).astype(float)
                    coords = self.navdata.get_sat_ecef_coordinates(
                        desired_time=desired_time, PRN=PRN)
                    X[i] = np.asarray(coords[0]).item()
                    Y[i] = np.asarray(coords[1]).item()
                    Z[i] = np.asarray(coords[2]).item()
                    dT_rel[i] = np.asarray(coords[3]).item()
                else:
                    # Correct observed pseudorange
                    Rji[i] = rji[i] + dTj[i] * c

                    # Compute time when signal left the satellite
                    Tj_GPS[i] = t - Rji[i] / c

                    # Compute satellite ECEF position and relativistic clock correction (dT_rel in this case is not only relativitic correction (GLONASS efemerids))
                    desired_time = np.atleast_1d(Tj_GPS[i]).astype(float)
                    coords = self.navdata.get_sat_ecef_coordinates(
                        desired_time=desired_time, PRN=PRN)
                    X[i] = np.asarray(coords[0]).item()
                    Y[i] = np.asarray(coords[1]).item()
                    Z[i] = np.asarray(coords[2]).item()
                    dT_rel[i] = np.asarray(coords[3]).item()

            # Compute the vector differences between receiver and satellite positions
            # Shape (num_sats, 3)
            diff = np.column_stack([X - x, Y - y, Z - z])

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

            # Update reciever coordiantes (class attributes)
            self.update_receiver_coordinates(x, y, z)

            # Update iteration info
            iteration += 1
            logger.info(f"Iteration {iteration}: dx = {dx}")
            # print(f"Iteration {iteration}: dx = {dx}") ### DELETE PRINT
            improvement = np.max(np.abs(dx))

        # Try to filter low satellites and update position
        try:
            x, y, z, dTi0, A, l, N, h, sat_nr = self.__estimate_without_low_satellites(
                X, Y, Z, dT_rel, Rji, sat_nr, x, y, z, t, dTi0)
        except Exception as e:
            logger.warning(
                f"WARNING (GNSS_MultipathAnalysis): Failed to perform the final position estimation after removing low-elevation satellites. "
                f"Error: {str(e)}"
            )
            pass

        stats_report = StatisticalAnalysis(
            A, l, N, h).run_statistical_analysis()

        return np.array([x, y, z, dTi0]), stats_report


if __name__ == "__main__":

    # rinObs = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx"
    rinObs = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.rnx"
    rinObs = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped – Kopi.rnx"
    rinNav = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\BRDC00IGS_R_20220010000_01D_MN.rnx"

    fasit = np.array([3149785.9652, 598260.8822, 5495348.4927])
    x_rec, y_rec, z_rec = fasit.T
    # desired_time = np.array([2022, 1, 1, 2, 20, 0.0000000])
    desired_time = np.array([2022, 1, 1, 0, 0, 30.0000000])
    # desired_time = np.array([2022, 1, 1, 0, 5, 0.0000000])
    desired_system = "R"

    # Example of not using file as input
    fasit = np.array([3149785.9652, 598260.8822, 5495348.4927])
    x_rec, y_rec, z_rec = fasit.T
    desired_sys = "E"
    desired_time = np.array([2022, 1, 1, 0, 0, 30.0000000])
    navObj = SatelliteEphemerisToECEF(
        rinNav,  x_rec, y_rec, z_rec, desired_sys)
    GNSS_obs, _, _, _, time_epochs, _, GNSSsystems, \
        obsCodes, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ = readRinexObs(
            rinObs)

    GNSSPos = BroadNavPositionEstimator(desired_time=desired_time,
                                        desired_system=desired_sys,
                                        navdata=navObj,
                                        GNSS_obs=GNSS_obs,
                                        time_epochs=time_epochs,
                                        GNSSsystems=GNSSsystems,
                                        obsCodes=obsCodes)
    estimated_position, stats = GNSSPos.estimate_position()
    X, Y, Z, dT = estimated_position.T
    print(X, Y, Z, dT)
    print(f"\nDifference: {np.round(fasit - np.array([X, Y, Z]), 3)}")
    print(f"Stats:\n{stats}")
