import numpy as np
from typing import Literal
from gnssmultipath.Geodetic_functions import date2gpstime, date2gpstime_vectorized
from gnssmultipath.SatelliteEphemerisToECEF import SatelliteEphemerisToECEF, Kepler2ECEF
from gnssmultipath.readRinexObs import readRinexObs
import gnssmultipath.Geodetic_functions as geof
from pyproj import Transformer


c = 299792458  # Speed of light [m/s]





class GNSSPositionEstimator:
    """
    Parameters:

    pseudoranges (numpy array): Array of pseudorange observations (meters), shape (N,).
    timestamps (numpy array): Array of timestamps for each observation, shape (N,).
    satellite_positions (numpy array): Satellite coordinates in ECEF (meters), shape (N, 3).
    desired_time: Time to estimate position for ( Ex: np.array([2019, 6, 4, 7, 30, 0.0000000]))

    """

    def __init__(self, rinex_obs_file:str, rinex_nav_file:str, desired_time:np.array,
                 desired_system:Literal["G","E","C"]="G",
                 x_rec_approx:float=0.0, y_rec_approx:float=0.0, z_rec_approx:float=0.0):


        self.desired_time = np.atleast_2d(date2gpstime_vectorized(desired_time)).T
        self.desired_sys = desired_system
        self.navdata = SatelliteEphemerisToECEF(rinex_nav_file, x_rec_approx, y_rec_approx, z_rec_approx, self.desired_sys)
        self.GNSS_obs, _, _, self.GNSS_SVs, self.time_epochs,_, self.GNSSsystems,\
              self.obsCodes, approxPosition, _, _, _, _, _, _, _, _,\
              _, _, _, _, _, _, _, _ = readRinexObs(rinex_obs_file)
              
        self.x_rec_approx = x_rec_approx
        self.y_rec_approx = y_rec_approx
        self.z_rec_approx = z_rec_approx
        
        self.gnss_index_mapper = self.gnss_index_mapper()
        self.sys_idx = self.gnss_index_mapper[self.desired_sys]
        self.sys_obs_codes = self.obsCodes[self.sys_idx][self.desired_sys]
        self.signal_idx = next((i for i, code in enumerate(self.sys_obs_codes) if code.startswith("C")), -1)
        


    def gnss_index_mapper(self) -> dict:
        mapper = {}
        for i in np.arange(0, len(self.GNSSsystems)):
            GNSSsystemIndex = list(self.GNSSsystems.keys())[i]
            mapper[self.GNSSsystems[GNSSsystemIndex]]= GNSSsystemIndex 
        return mapper
        
    
    @staticmethod
    def find_idx_of_epoch_closest_in_time(desired_time: np.ndarray, obs_time_epochs: np.ndarray) -> int:
        """
        Find the index of the epoch in `obs_time_epochs` where the difference in time is the smallest.
    
        Parameters:
        ----------
        desired_time (numpy.ndarray): Desired GPS time as a 2D array [week, seconds].
        obs_time_epochs (numpy.ndarray): Array of observed time epochs as a 2D array [week, seconds].
    
        Returns:
        -------
        int: Index of the epoch with the smallest time difference.
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
    def compute_toc_for_single_ephemeris(ephemeris):
        """
        Compute the time of clock correction (toc) for a single satellite.
    
        Parameters:
        -----------
        ephemeris : list or numpy array
            Ephemeris data for a single satellite. The input must include the year, 
            month, day, hour, minute, and second at positions [1, 2, 3, 4, 5, 6] respectively.
    
        Returns:
        --------
        float
            The computed `toc` value (time of clock correction).
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
        
    
    
    def preprocess_satellites(self, sat_nr, rji, x, y, z, desired_tow):
        """
        Preprocess satellites by filtering based on elevation.
        
        Parameters:
        -----------
        sat_nr : list
            List of satellite PRNs.
        rji : numpy array
            Observed pseudoranges for each satellite.
        x, y, z : float
            Receiver's approximate ECEF coordinates.
        desired_tow : float
            Desired time of week (seconds).
        
        Returns:
        --------
        dict
            Preprocessed satellite data including positions, clock corrections, and valid PRNs.
        """
        valid_X, valid_Y, valid_Z = [], [], []
        valid_dT_rel, valid_rji = [], []
        valid_sat_nr, valid_ephemeris = [], []
        toc_lst = []

        for i, sat_num in enumerate(sat_nr):
            ephemeris = self.navdata.get_closest_ephemerides_for_PRN_at_time(f"{self.desired_sys}{str(sat_num).zfill(2)}", desired_tow)
            if not isinstance(ephemeris, np.ndarray):
                print(f"Satellite {str(sat_num).zfill(2)} skipped due to missing ephemeris.")
                continue

            a0 = float(ephemeris[0][7])
            a1 = float(ephemeris[0][8])
            a2 = float(ephemeris[0][9])
            T_GD = float(ephemeris[0][32])

            # Compute satellite clock correction
            toc = self.compute_toc_for_single_ephemeris(ephemeris[0])
            dTj = a0 + a1 * (desired_tow[0] - toc) + a2 * (desired_tow[0] - toc) ** 2 - T_GD
            Rji = rji[i] + dTj * c
            Tj_GPS = desired_tow[0] - Rji / c

            # Compute satellite position
            X_, Y_, Z_, dT_rel_ = Kepler2ECEF(x, y, z).kepler2ecef(ephemeris, desired_tow)
            az, el = geof.compute_azimut_elev(X_[0], Y_[0], Z_[0], x, y, z)

            # if el > 15:
            if True:
                valid_X.append(X_[0])
                valid_Y.append(Y_[0])
                valid_Z.append(Z_[0])
                valid_dT_rel.append(dT_rel_[0])
                toc_lst.append(toc)
                valid_rji.append(rji[i])
                valid_sat_nr.append(sat_num)
                valid_ephemeris.append(ephemeris)
            else:
                print(f"Satellite {sat_num} excluded due to low elevation: {np.round(el, 2)}°")

        return {
            "X": np.array(valid_X),
            "Y": np.array(valid_Y),
            "Z": np.array(valid_Z),
            "dT_rel": np.array(valid_dT_rel),
            "toc": np.array(toc_lst),
            "rji": np.array(valid_rji),
            "sat_nr": valid_sat_nr,
            "ephemeris": valid_ephemeris,
        }

    def estimate_position(self, max_iterations=10, tol=1e-6):
        """
        Estimate the receiver position using vectorized least-squares.
        """
        x, y, z = self.x_rec_approx, self.y_rec_approx, self.z_rec_approx
        dTi0 = 0.0  # Initial clock bias

        idx_of_closest_eph = self.find_idx_of_epoch_closest_in_time(self.desired_time, self.time_epochs)
        pseudoranges = self.GNSS_obs[self.desired_sys][idx_of_closest_eph]
        sat_nr = list(set(np.where(pseudoranges != 0.0)[0]))
        rji = pseudoranges[:, self.signal_idx][sat_nr]
        

        # Preprocess satellite data
        processed_data = self.preprocess_satellites(sat_nr, rji, x, y, z, self.desired_time[1])
        X, Y, Z = processed_data["X"], processed_data["Y"], processed_data["Z"]
        dT_rel, rji = processed_data["dT_rel"], processed_data["rji"]
        efemeris, sat_nr = processed_data["ephemeris"], processed_data["sat_nr"]
        toc = processed_data["toc"]
        t = self.desired_time[1][0]
        
        # Preallocate
        num_sats = len(sat_nr)
        dTj = np.zeros(num_sats)
        Tj_GPS = np.zeros(num_sats)
        Rji = np.zeros(num_sats)
        
        iteration = 0
        improvement = 10
        while improvement > tol and iteration < max_iterations:
        # while improvement > 1:
            for i, sat_num in enumerate(sat_nr):
                # Extract satellite clock parameters
                a0 = float(efemeris[i][0][7])
                a1 = float(efemeris[i][0][8])
                a2 = float(efemeris[i][0][9])
                T_GD = float(efemeris[i][0][32])
                
                # Calculate satellite clock error
                dTj[i] = a0 + a1 * (t - toc[i]) + a2 * (t - toc[i]) ** 2 - T_GD
        
                # Correct observed pseudorange
                Rji[i] = rji[i] + dTj[i] * c
        
                # Compute time when signal left the satellite
                Tj_GPS[i] = t - Rji[i] / c
        
                # Compute satellite ECEF position and relativistic clock correction
                X_, Y_, Z_, dT_rel_ = Kepler2ECEF(x, y, z).kepler2ecef(efemeris[i], Tj_GPS[i])
                X[i], Y[i], Z[i], dT_rel[i]  = X_[0], Y_[0], Z_[0], dT_rel_[0]
            
        
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
            print(f"Iteration {iteration}: dx = {dx}")
            improvement = np.max(np.abs(dx))
            
        return np.array([x, y, z, dTi0])
            
        
    
    
    
    




if __name__=="__main__":
    
    # rinObs = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx"
    rinObs = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.rnx"
    rinNav = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\BRDC00IGS_R_20220010000_01D_MN.rnx"

    
    fasit = np.array([3149785.9652, 598260.8822, 5495348.4927])
    x_rec, y_rec, z_rec = fasit.T
    desired_time = np.array([2022, 1, 1, 2, 20, 0.0000000])
    desired_time = np.array([2022, 1, 1, 0, 0, 30.0000000])
    desired_system = "G"
    
    
    ### TEST
    # rinNav = r"C:\Users\perhe\OneDrive\Documents\Python_skript\Approx_rec_pos\BRDM00DLR_S_20191550000_01D_MN.rnx"
    # fasit = np.array([3173400.0, 605500.0, 5481000.0])
    # desired_time = np.array([2019, 6, 4, 7, 30, 0.0000000])
    
    # GNSSPos = GNSSPositionEstimator(rinObs, rinNav, desired_time, desired_system, x_rec_approx=x_rec, y_rec_approx=y_rec, z_rec_approx=z_rec)
    GNSSPos = GNSSPositionEstimator(rinObs, rinNav, desired_time, desired_system, x_rec_approx=0, y_rec_approx=0, z_rec_approx=0)
    estimated_position = GNSSPos.estimate_position()
    X,Y,Z,dT = estimated_position.T
    
    lon, lat, h = Transformer.from_crs("EPSG:4936", "EPSG:4258").transform(X,Y,Z)
    estimated__etrs89_nn2000 = np.array(Transformer.from_crs("EPSG:4936", "EPSG:5972").transform(X,Y,Z))
    fasit_etrs89_nn2000 = Transformer.from_crs("EPSG:4936", "EPSG:5972").transform(*fasit.T)
    
    print(X,Y,Z,dT)
    print(f"\nDifference: {np.round(fasit- np.array([X,Y,Z]),3)}")
    print(f"\nDifference in UTM Z32: {np.round(fasit_etrs89_nn2000- estimated__etrs89_nn2000,3)}")
    
