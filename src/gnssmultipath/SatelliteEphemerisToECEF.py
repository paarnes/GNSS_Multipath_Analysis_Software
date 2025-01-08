"""
Module for converting from broadcasted ephemerides to ECEF and interpolated to the desired epoch.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""


from typing import Literal, Optional, List, Union, Tuple
import numpy as np
from numpy import ndarray
from tqdm import tqdm
from gnssmultipath.Geodetic_functions import date2gpstime_vectorized, get_leap_seconds, gpstime2date_arrays, ECEF2enu, ECEF2geodb
from gnssmultipath.RinexNav import Rinex_v3_Reader



class GLOStateVec2ECEF:

    """
    Class for interpolating a GLONASS state vector using 4th order Runge Kutta
    """
    def __init__(self):
        self.sys_code = "R" # glonass system code


    def interpolate_glonass_coord_runge_kutta(self, filtered_eph_data, time_epochs):
        """
        Function that use broadcast ephemerides (from rinex nav file) and interpolates to the
        desired time defined by "time_epochs" using 4th order Runge-kutta.

        Parameter:
        ------
        ephemerides: A array with ephemerides for current satellite
        time_epochs: An n_obs X 2 sized array with weeks and tow read from the rinex observation file.

        Return:
        -------
        pos            : The interpolated satellite coordiantes for the desired time defined by "time_epochs" [m]
        vel            : The satellite velocity for the desired time defined by "time_epochs" [m/s]
        clock_err      : The satellite absolute clock error for the desired time defined by "time_epochs" [seconds]
        clock_rate_err : The satellite clock rate error for the desired time defined by "time_epochs"

        """
        filtered_eph_data[:,0] = np.nan #reming char from string (ex G01 -> 9999)
        filtered_eph_data = filtered_eph_data.astype(float)

        ## Read in data:
        week,toc = date2gpstime_vectorized(filtered_eph_data[:, 1:7].astype(int))
        tauN = filtered_eph_data[:,7]   # SV clock bias (sec) (-TauN)
        gammaN = filtered_eph_data[:,8] # SV relative frequency bias (+GammaN)

        tow_rec = np.round(time_epochs, 6)  # extracting tow
        x_te = filtered_eph_data[:,10]  # X-coordinates at t_e in PZ-90 [km]
        y_te = filtered_eph_data[:,14]  # Y-coordinates at t_e in PZ-90 [km]
        z_te = filtered_eph_data[:,18]  # Z-coordinates at t_e in PZ-90 [km]

        vx_te = filtered_eph_data[:,11] # Velocity component at t_e in PZ-90 (v_x) [km/s]
        vy_te = filtered_eph_data[:,15] # Velocity component at t_e in PZ-90 (v_y) [km/s]
        vz_te = filtered_eph_data[:,19] # Velocity component at t_e in PZ-90 (v_z) [km/s]

        J_x = filtered_eph_data[:,12]   # Moon and sun acceleration at t_e [km/sec**2]
        J_y = filtered_eph_data[:,16]   # Moon and sun acceleration at t_e [km/sec**2]
        J_z = filtered_eph_data[:,20]   # Moon and sun acceleration at t_e [km/sec**2]

        ## -- Convert from UTC to GPST by adding leap seconds
        leap_sec = get_leap_seconds(week[0],toc[0]) # Get correct leap sec based on date (using first epoch)
        toc_gps_time = toc + leap_sec # convert from UTC to GPST

        ## -- Find time difference
        tdiff = tow_rec - toc_gps_time

        ## -- Clock correction (except for general relativity which is applied later)
        clock_err = tauN + tdiff * (gammaN)
        clock_rate_err = gammaN

        # Vectorized initialization
        init_state = np.array([x_te, y_te, z_te, vx_te, vy_te, vz_te]).T  # Transpose to make it a 2D array
        init_state = 1000*init_state        # converting from km to meters
        acc = 1000*np.array([J_x, J_y,J_z]) # converting from km to meters
        state = init_state
        tstep = 90
        tt = np.where(tdiff < 0, -tstep, tstep)
        while np.any(np.abs(tdiff) > 1e-9):
            tt = np.where(np.abs(tdiff) < tstep, tdiff, tt)
            k1 = self.glonass_diff_eq(state, acc)
            k2 = self.glonass_diff_eq(state + k1 * tt[:, None] / 2, -acc)
            k3 = self.glonass_diff_eq(state + k2 * tt[:, None] / 2, -acc)
            k4 = self.glonass_diff_eq(state + k3 * tt[:, None], -acc)
            state += (k1 + 2 * k2 + 2 * k3 + k4) * tt[:, None] / 6.0
            tdiff -= tt

        pos = state[:, :3]
        vel = state[:, 3:6]
        return pos, vel, clock_err, clock_rate_err


    def glonass_diff_eq(self, state:ndarray , acc:ndarray):
        """
        Calculate the derivatives of Glonass satellite coordinates using the Runge-Kutta interpolation method.

        INPUT:
        ------
        state : numpy.ndarray. State vector containing the coordinates and velocities of Glonass satellites.
                Shape: (n, 6), where n is the number of epochs.

        acc   : numpy.ndarray. Acceleration vector. Shape: (3, n), where n is the number of epochs.

        Returns:
        --------
        der_state : numpy.ndarray. Derivatives of the state vector representing the change in satellite coordinates and velocities.
                    Shape: (n, 6), where n is the number of epochs.

        Constants:
        -----------
        J2    : float. Second zonal coefficient of the spherical harmonic expression.
        mu    : float. Gravitational constant [m^3/s^2] (product of the mass of the Earth and the gravitational constant).
        omega : float. Earth rotation rate [rad/s].
        ae    : float. Semi-major axis of the PZ-90 ellipsoid [m].
        """

        J2 = 1.0826257e-3     # Second zonal coefficient of spherical harmonic expression.
        mu = 3.9860044e14     # Gravitational constant [m3/s2]   (product of the mass of the earth and and gravity constant)
        omega = 7.292115e-5   # Earth rotation rate    [rad/sek]
        ae = 6378136.0        # Semi-major axis PZ-90   [m]

        r = np.linalg.norm(state[:, :3], axis=1)  # Euclidean norm for the radius
        der_state = np.zeros((state.shape[0], 6))

        zero_indices = r**2 < 0
        der_state[zero_indices] = 0  # Set derivatives to zero for cases where r^2 < 0

        a = 1.5 * J2 * mu * (ae**2) / (r**5)
        b = 5 * (state[:, 2]**2) / (r**2)
        c = -mu / (r**3) - a * (1 - b)

        der_state[:, :3] = state[:, 3:]
        der_state[:, 3] = (c + omega**2) * state[:, 0] + 2 * omega * state[:, 4] + acc[0, :]
        der_state[:, 4] = (c + omega**2) * state[:, 1] - 2 * omega * state[:, 3] + acc[1, :]
        der_state[:, 5] = (c - 2 * a) * state[:, 2] + acc[2, :]

        return der_state


    def get_leap_seconds(self, week,tow):
        """
        Add leap seconds based on date. Input is week and tow for current obs.
        """
        year,month,day,_,_,_ = gpstime2date_arrays(week, tow) # convert to gregorian date
        ep_time = (year,month,day)
        if ep_time <= (2006, 1, 1):
            raise ValueError("Have no history on leap seconds before 2006")
        elif ep_time <= (2009, 1, 1):
            return 14
        elif ep_time <= (2012, 7, 1):
            return 15
        elif ep_time <= (2015, 7, 1):
            return 16
        elif ep_time <= (2017, 1, 1):
            return 17
        else:
            return 18










class Kepler2ECEF:
    """
    Class for converting from Kepler elements from a
    RINEX navigation file to Earth Centered Earth Fixed (ECEF).

    Initialize the class with a RINEX navigation file, the reciecer coordinates in ECEF
    (x_rec, y_rec, z_rec). data_rate is optional
    and set to 60 min as default. (60 min between each epoch)

    """

    def __init__(self, x_rec, y_rec, z_rec, data_rate=60):
        self.x_rec = x_rec
        self.y_rec = y_rec
        self.z_rec = z_rec
        self.data_rate = data_rate

    def kepler2ecef(self, filtered_eph_data, tow_rec):
        """
        The function calculates satellite coordinates and corrects for Earth's rotation. It transforms from Kepler elements to Earth-Centered, Earth-Fixed (ECEF) coordinates.
        Supports GPS,Galileo and BeiDou.

        Kepler-to-ECEF Earth Rotation Correction:
        The Earth rotation correction ensures the satellite's position aligns with the time of signal transmission.
        This is not necessary for SP3 files which already provide satellite positions in the Earth-Centered Earth-Fixed (ECEF) frame.
        """
        GM         = 3.986005e14      # Product of Earth's mass and the gravitational constant
        omega_e    = 7.2921151467e-5  # Earth's angular velocity [rad/second]
        c          = 299792458        # Speed of light [m/s]

        gnss_sys = filtered_eph_data[0,0][0]
        filtered_eph_data[:,0] = np.nan  # Remove cell containing a string representing PRN with system code (to be able to convert array to float)
        filtered_eph_data = filtered_eph_data.astype(float)

        # Extract data:
        M0         = filtered_eph_data[:,13] # Mean anomaly at reference epoch
        delta_n    = filtered_eph_data[:,12] # Mean motion difference from computed value [rad/second]
        e          = filtered_eph_data[:,15] # Eccentricity (e)
        A          = filtered_eph_data[:,17]**2  # Square root of semimajor axis in (sqrt(meter))
        OMEGA      = filtered_eph_data[:,20] # Longitude of ascending node (OMEGA0) in radians.
        i0         = filtered_eph_data[:,22] # Inclination angle at reference time (i0) [rad]
        omega      = filtered_eph_data[:,24] # Argument of perigee (omega) in [rad]
        OMEGA_dot  = filtered_eph_data[:,25] # Rate of right ascension (OMEGA_dot) [rad/second]
        i_dot      = filtered_eph_data[:,26] # Rate of inclination angle (i_dot) [rad/second]
        C_uc       = filtered_eph_data[:,14] # cosine harmonic correction coefficients (argument of perigee)
        C_us       = filtered_eph_data[:,16] # sine harminic Correction coefficients (argument of perigee)
        C_rc       = filtered_eph_data[:,23] # cosine harmonic correction term to the orbit radius (C_rc)
        C_rs       = filtered_eph_data[:,11] # sine harmonic correction term to the orbit radius (C_rs)
        C_ic       = filtered_eph_data[:,19] # cosine harmonic correction term to the angle of inclination (C_ic)
        C_is       = filtered_eph_data[:,21] # sine harmonic correction term to the angle of inclination (C_is)
        toe        = filtered_eph_data[:,18] # Reference "time of ephemeris" data (toe) (sec of GPST,BDT,GLONASST or GST) week)

        n0  = np.sqrt(GM/A**3) # Computed mean motion (n0) [rad/second]
        t_k = tow_rec - toe # Time difference between reception time and ephemeris reference time

        n_k = n0 + delta_n   # Corrected mean motion (n_k)
        M_k = M0 + n_k*t_k   # Mean anomaly (M_k) [rad]

        # Calculate eccentric anomaly
        E = M_k.copy()
        for _ in range(10):
            E = M_k + e * np.sin(E)
            dE = np.fmod(E - M_k, 2 * np.pi)
            mask = abs(dE) < 1.e-12
            if np.all(mask):
                break

        cosv = (np.cos(E) - e)/(1 - e*np.cos(E))
        sinv = (np.sqrt(1 - e**2)*np.sin(E))/(1-e*np.cos(E))
        tanv = sinv/cosv

        # Quadrant correction
        v = np.arctan(tanv)
        mask = (sinv > 0) & (cosv < 0) | (sinv < 0) & (cosv < 0)
        v[mask] += np.pi
        mask = sinv < 0
        v[mask] += 2 * np.pi

        theta   = v + omega
        # Orbital parameter corrections
        du_k    = C_uc*np.cos(2*theta) + C_us*np.sin(2*theta)
        dr_k    = C_rc*np.cos(2*theta) + C_rs*np.sin(2*theta)
        di_k    = C_ic*np.cos(2*theta) + C_is*np.sin(2*theta)

        # Corrected orbital parameters
        u_k     = theta + du_k
        r_k     = A*(1 - e*np.cos(E)) + dr_k
        i_k     = i0 + i_dot*t_k + di_k

        # Corrected longitude of ascending node
        OMEGA_k = OMEGA + (OMEGA_dot - omega_e)*t_k - omega_e*toe

        # Satellite position in the orbit
        x = r_k*np.cos(u_k)
        y = r_k*np.sin(u_k)

        # ECEF coordinates for the satellite
        X = x*np.cos(OMEGA_k) - y*np.sin(OMEGA_k)*np.cos(i_k)
        Y = x*np.sin(OMEGA_k) + y*np.cos(OMEGA_k)*np.cos(i_k)
        Z = y*np.sin(i_k)

        # Relativistic clock correction
        dT_rel = (-2/c**2)*np.sqrt(A*GM)*e*np.sin(E)

        if (abs(self.x_rec)) > 1.0 and (abs(self.y_rec)) > 1.0 and (abs(self.z_rec)) > 1.0:
            TRANS = 0
            TRANS0 = 0.075  # approximate signal travel time
            j = 0
            while np.any(abs(TRANS0 - TRANS) > 1e-10):
                j += 1
                if (j > 20):
                    print('Error: The travel time-rotation does not converge!')
                    break

                TRANS = TRANS0
                OMEGA_k = OMEGA + (OMEGA_dot - omega_e) * t_k - omega_e * (toe + TRANS)
                X = x * np.cos(OMEGA_k) - y * np.sin(OMEGA_k) * np.cos(i_k)
                Y = x * np.sin(OMEGA_k) + y * np.cos(OMEGA_k) * np.cos(i_k)
                Z = y * np.sin(i_k)
                dX = (X - self.x_rec)
                dY = (Y - self.y_rec)
                dZ = (Z - self.z_rec)
                DS = np.sqrt(dX**2 + dY**2 + dZ**2)
                TRANS0 = DS / c
        else:
            # Do not correct for the earth rotation
            OMEGA_k = OMEGA + (OMEGA_dot - omega_e)*t_k - omega_e*toe
            X = x*np.cos(OMEGA_k) - y*np.sin(OMEGA_k)*np.cos(i_k)
            Y = x*np.sin(OMEGA_k) + y*np.cos(OMEGA_k)*np.cos(i_k)
            Z = y*np.sin(i_k)

        return X, Y, Z, dT_rel







class SatelliteEphemerisToECEF:
    """
    Parent class for converting from broadcasted ephemerides to Earth Centered Earth Fixed (ECEF).
    Make use of the classes "GLOStateVec2ECEF" and "Kepler2ECEF".

    Input:
    -----

    - rinex_nav_file  : List or string. Takes in both a single RINEX navfile or a list of RINEX navigation file. If a list is provided, the data will
                      be merged in a single array (merged on to one file).
    - x_rec           : Receiver X-coordinates in ECEF
    - y_rec           : Receiver Y-coordinates in ECEF
    - z_rec           : Receiver Z-coordinates in ECEF
    - desired_systems : List of desired system codes ["G","R","E"]
    - data_rate       : Update rate of the ephemerides. Higher values correspond to less data and faster file reading (in minutes)


    Examples on how to use the class:

    - Compute satellite coordinates for all systems and satellites where the time is defined GPS time "Time-of-week" (seconds):
    .. code-block:: python
            sat_pos = SatelliteEphemerisToECEF(rin_nav_file, x_rec, y_rec, z_rec).get_sat_ecef_coordinates(tow_epochs) # where "tow_epochs" is an 1D array with time-of-week


    - Compute satellite coordinates for a specific satellite (in this case "R10") and where the time is defined in GPS time and "Time-of-week" (seconds):
    .. code-block:: python
            sat_pos = SatelliteEphemerisToECEF(rin_nav_file, x_rec, y_rec, z_rec).get_sat_ecef_coordinates(tow_epochs, PRN = 'R10')


    - Compute satellite coordinates for a specific satellite (in this case "G20") and where the time is defined in Gregorian time (year, month,day,hour,minute,seconds):
    .. code-block:: python
            time_epochs = np.array([[2020,   10,   30,   13,   22,   14],
                                    [2020,   10,   30,   13,   22,   15],
                                    [2020,   10,   30,   13,   22,   16]])

            sat_pos = SatelliteEphemerisToECEF(rin_nav_file, x_rec, y_rec, z_rec).get_sat_ecef_coordinates(time_epochs, time_fmt='GREGORIAN', PRN = 'G20')


    - Compute satellite azimuth and elevation angles
    .. code-block:: python
            CONVERTER = SatelliteEphemerisToECEF(rinnav, x_rec, y_rec, z_rec, data_rate=120) # create an object
            CONVERTER.get_sat_ecef_coordinates(time_epochs)                                  # compute satellite coordinates
            sat_pos = CONVERTER.compute_satellite_azimut_and_elevation_angle()               # compute satellite elevation and azimuth angles


    """

    def __init__(self, rinex_nav_file:Union[str, List[str]], x_rec, y_rec, z_rec, desired_systems: Optional[List[str]] = None, data_rate=60):

        if desired_systems is None:
            desired_systems = ["G", "R", "E", "C"]

        if isinstance(rinex_nav_file, list):
            self.ephemerides, self.glo_fcn = self.read_a_list_of_nav_files(rinex_nav_file, data_rate=data_rate)
        else:
            self.nav_data = Rinex_v3_Reader().read_rinex_nav(rinex_nav_file, desired_systems, data_rate=data_rate)
            self.ephemerides = self.nav_data['ephemerides']
            self.glo_fcn = self.nav_data['glonass_fcn']

        self.x_rec = x_rec
        self.y_rec = y_rec
        self.z_rec = z_rec
        self.max_sat_per_sys = {"G" : 36, "R" : 36, "E" : 36, "C" : 100}
        self.available_systems = self.find_available_systems_in_eph_data()
        self.available_systems = list(set(self.available_systems).intersection(desired_systems))
        self.sat_coord = {sys: {"position": {str(PRN): None for PRN in range(1, max_sat + 1)}} for sys, max_sat in self.max_sat_per_sys.items() if sys in self.available_systems}
        self.sat_coord_computed = False
        self.prn_overview = self.get_availible_satellites_for_a_system()
        self.system_code_mapper = {"G" : "GPS", "R" : "GLONASS", "C" : "BeiDou", "E" : "Galileo"}
        self.sys_names = [self.system_code_mapper[sys_code] for sys_code in self.available_systems]
        self.total_sats = sum(len(self.prn_overview[sys])for sys in self.available_systems)


    def read_a_list_of_nav_files(self, rinex_nav_file, data_rate):
        """
        Reads in a list of RINEX navigation files and merge the data
        into one data array.
        """
        nav_files = [nav_file for nav_file in rinex_nav_file if nav_file != ""]
        first_nav_data = Rinex_v3_Reader().read_rinex_nav(nav_files[0], data_rate=data_rate)
        data = np.concatenate([first_nav_data['ephemerides']] + [Rinex_v3_Reader().read_rinex_nav(nav_file, data_rate=data_rate)['ephemerides'] for nav_file in nav_files[1:]], axis=0)
        glonass_fcn = first_nav_data.get('glonass_fcn', None)
        data = np.unique(data, axis=0) # ensure no duplicates
        return data, glonass_fcn


    def filter_array_on_PRN(self, PRN:str) -> ndarray:
        """
        Filtering the array with ephemerides based on system code "G", "E" "C" or "R"
        based on PRN number. PRN needs to be in the form "G04".
        """
        mask = np.char.startswith(self.ephemerides[:,0], PRN)
        return self.ephemerides[mask]


    def filter_array_on_system(self, sys_code:str) -> ndarray:
        """ Filtering the array with ephemerides based on system code "G", "E" "C" or "R" """
        mask = np.char.startswith(self.ephemerides[:,0], sys_code)
        return self.ephemerides[mask]

    def get_availible_satellites_for_a_system(self) -> dict:
        """Creates a dict of list of available satellites for each availible system """
        aviliable_sats = {sys:None for sys in self.available_systems}
        for sys in self.available_systems:
            curr_sys = self.filter_array_on_system(sys)
            aviliable_sats[sys] = list(np.unique(curr_sys[:, 0]))
        return aviliable_sats

    def find_available_systems_in_eph_data(self):
        """Returns a list of all avaiable systems in navigation file"""
        satcodes = self.ephemerides[:, 0]
        all_sys_codes = [code[0] for code in satcodes if code and code[0].isalpha()]
        return list(set(all_sys_codes))

    def get_closest_ephemerides_for_PRN_at_time(self, PRN:str, desired_tow:ndarray) -> ndarray:
        """
        Collects the ephemerides for the specified PRN number,
        and retrieves those that are closest to the specified time.

        INPUT:
        -----
        ephemerids : array containing ephemerides for all available systems and satellites
        PRN        : str that specifies which satellite ehemerides shall be extracted for
        desired_tow: array the time you are looking for a match for

        OUTPUT:
        ------
        eph_closest_in_time: array containing ephemerides (closest in time) for the specific PRN

        """
        ephemerids_filtered = self.filter_array_on_PRN(PRN)
        if ephemerids_filtered.size == 0:
            return None
        epochs_dates = ephemerids_filtered[:, 1:7].astype(int)
        week_sat, tow_sat = date2gpstime_vectorized(epochs_dates)
        diff = np.abs(tow_sat[:, np.newaxis] - desired_tow)
        closest_indices  = np.argmin(diff, axis=0)
        closest_indices_repeated = closest_indices[:, np.newaxis].repeat(desired_tow.shape[0], axis=1)[:,0]
        modified_eph_array = ephemerids_filtered[closest_indices_repeated] # repeting array that contains all ephemerides needed in correct order wrt to observation epochs
        return modified_eph_array



    def compute_satellite_azimut_and_elevation_angle(self, drop_below_horizon=False):
        """
        Computes the satellites azimute and elevation angle based on satellitte and
        receiver ECEF-coordinates. Utilizes vectorization (no for loops) for
        better performance.

        drop_below_horizon : Sets values to np.nan for satellites below horizon
        """
        # Check if satellite coordinates are computed
        if not self.sat_coord_computed:
            raise ValueError('Satellite coordinates are not computed. Please compute coordinates by calling the method "get_sat_ecef_coordinates" before performing this operation.')


        ## -- WGS 84 ellipsoid:
        a   =  6378137.0         # semi-major ax
        b   =  6356752.314245    # semi minor ax

        # Compute latitude and longitude for the receiver
        lat,lon,h = ECEF2geodb(a, b, self.x_rec, self.y_rec, self.z_rec)

        bar_format = '{desc}:{percentage:3.0f}%|{bar}|({n_fmt}/{total_fmt} satellites)'
        desc = ', '.join(self.sys_names[:-1]) + (' and ' + self.sys_names[-1] if len(self.sys_names) > 1 else self.sys_names[0])
        with tqdm(total = self.total_sats, desc =f"Satellite azimuth and elevation angles are being computed for {desc}" , position=0, leave=True, bar_format=bar_format) as pbar:
            for sys in self.available_systems:
                az_array = np.full((self.nepochs, self.max_sat_per_sys[sys]+1), np.nan)
                el_array = np.full((self.nepochs, self.max_sat_per_sys[sys]+1), np.nan)
                for sys_code in self.prn_overview[sys]:
                    PRN = int(sys_code[1::])
                    if self.sat_coord[sys]['position'][str(PRN)] is None:
                        pbar.update(1) # Update the progress bar for each satellite processed
                        continue
                    X,Y,Z = self.sat_coord[sys]['position'][str(PRN)].T
                    # Find coordinate difference between satellite and receiver
                    dX = (X - self.x_rec)
                    dY = (Y - self.y_rec)
                    dZ = (Z - self.z_rec)

                    # Convert from ECEF to ENU (east,north, up)
                    east, north, up = np.vectorize(ECEF2enu)(lat,lon,dX,dY,dZ)

                    # Calculate azimuth angle and correct for quadrants
                    azimuth = np.rad2deg(np.arctan(east/north))
                    azimuth = np.where((east > 0) & (north < 0) | ((east < 0) & (north < 0)), azimuth + 180, azimuth)
                    azimuth = np.where((east < 0) & (north > 0), azimuth + 360, azimuth)

                    # Calculate elevation angle
                    elevation = np.rad2deg(np.arctan(up / np.sqrt(east**2 + north**2)))
                    if drop_below_horizon:
                        elevation = np.where((elevation <= 0) | (elevation >= 90), np.nan, elevation) # Set elevation angle to NaN if not in the range (0, 90)
                        azimuth = np.where((elevation <= 0) | (elevation >= 90), np.nan, azimuth) # Set elevation angle to NaN if not in the range (0, 90)
                    # Store azimuth and elevation in the numpy arrays
                    az_array[:,PRN] = azimuth
                    el_array[:,PRN] = elevation
                    pbar.update(1) # Update the progress bar for each satellite processed
                self.sat_coord[sys]['azimuth'] = az_array
                self.sat_coord[sys]['elevation'] = el_array

        return self.sat_coord


    @staticmethod
    def compute_azimuth_and_elevation_static(
        X: np.ndarray,
        Y: np.ndarray,
        Z: np.ndarray,
        x_rec: float,
        y_rec: float,
        z_rec: float
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Computes the azimuth and elevation angles of satellites relative to the receiver.

        Parameters:
        -----------
        X, Y, Z     : ndarray
                      Satellite ECEF coordinates (arrays of the same size)
        x_rec, y_rec, z_rec : float
                              Receiver ECEF coordinates (single values)

        Returns:
        --------
        azimuth, elevation : tuple of ndarray
                             Azimuth and elevation angles (in degrees) for each satellite
        """
        # Calculate differences
        dX = X - x_rec
        dY = Y - y_rec
        dZ = Z - z_rec

        # Compute the receiver's geodetic latitude and longitude using WGS-84 ellipsoid parameters
        a = 6378137.0  # Semi-major axis (meters)
        b = 6356752.314245  # Semi-minor axis (meters)
        lat_rec, lon_rec, _ = ECEF2geodb(a, b, x_rec, y_rec, z_rec)

        # Convert differences to local ENU (East, North, Up) coordinates
        east, north, up = np.vectorize(ECEF2enu)(lat_rec, lon_rec, dX, dY, dZ)

        # Calculate azimuth angle and correct for quadrants
        azimuth = np.rad2deg(np.arctan(east / north))
        azimuth = np.where((east > 0) & (north < 0) | ((east < 0) & (north < 0)), azimuth + 180, azimuth)
        azimuth = np.where((east < 0) & (north > 0), azimuth + 360, azimuth)

        # Calculate elevation angle
        elevation = np.rad2deg(np.arctan(up / np.sqrt(east**2 + north**2)))

        return azimuth, elevation


    def get_sat_ecef_coordinates(self, desired_time:ndarray, time_fmt:Literal["TOW", "GREGORIAN"] = 'TOW', PRN: Optional[str] = None) -> ndarray:
        """
        Main method for converting from Kepler to ECEF and interpolate to the desired time. Input time can be both "Time of Week" (seconds) or
        gregorian time. This is set by time_format. Set to "TOW" by default.". If no PRN is set, the method will compute the satellite coordiantes for
        all available systems and satellites.

        INPUT:
        -----
        desired_time: The desired time for interpolation. Given in "Time-of-Week" or gregorian time.

                      Example on desired_time input fmt:
                      ---------------------------------
                      TOW: array([480134, 480135, 480136, 480137, 480138])

                      GREOGORIAN: array([[2020,   10,   30,   13,   22,   14],
                                         [2020,   10,   30,   13,   22,   15],
                                         [2020,   10,   30,   13,   22,   16]])

        time_fmt : Either Time-Of-Week (TOW) or gregorian calender (year, month, day, hh, min, sec)
        PRN : Optional. The PRN number of the desired satellite included system code in front (ex: "G04" or "E19").

        OUTPUT:
        ------
        sat_coord : dict containing X,Y,Z -coordinate for choosen system and satellites

        """

        if time_fmt == "GREGORIAN":
            desired_time = np.atleast_2d(date2gpstime_vectorized(desired_time)).T
        if PRN and 'R' not in PRN:
            filtered_eph_data = self.get_closest_ephemerides_for_PRN_at_time(PRN, desired_time)
            xs, ys, zs, dTrel = Kepler2ECEF(self.x_rec, self.y_rec, self.z_rec).kepler2ecef(filtered_eph_data,desired_time)
            return xs, ys, zs, dTrel
        elif PRN and 'R' in PRN:
            filtered_eph_data = self.get_closest_ephemerides_for_PRN_at_time(PRN, desired_time)
            sat_coords, sat_velocity, dT, dT_rate = GLOStateVec2ECEF().interpolate_glonass_coord_runge_kutta(filtered_eph_data, desired_time)
            xs, ys, zs = sat_coords.T
            return xs, ys, zs, dT
        else:
            bar_format = '{desc}:{percentage:3.0f}%|{bar}|({n_fmt}/{total_fmt} satellites)'
            desc = ', '.join(self.sys_names[:-1]) + (' and ' + self.sys_names[-1] if len(self.sys_names) > 1 else self.sys_names[0])
            with tqdm(total = self.total_sats, desc =f"Satellite coordinates are being computed for {desc}" , position=0, leave=True, bar_format=bar_format) as pbar:
                for sys in self.available_systems:
                    for sys_code in self.prn_overview[sys]:
                        PRN = int(sys_code[1::])
                        filtered_eph_data = self.get_closest_ephemerides_for_PRN_at_time(sys_code, desired_time)
                        if filtered_eph_data is None:
                            pbar.update(1)
                            continue
                        if 'R' in sys_code:
                            current_sat_coord,_,_,_ = GLOStateVec2ECEF().interpolate_glonass_coord_runge_kutta(filtered_eph_data, desired_time)
                        else:
                            xs, ys, zs, dTrel = Kepler2ECEF(self.x_rec, self.y_rec, self.z_rec).kepler2ecef(filtered_eph_data, desired_time)
                            current_sat_coord = np.array([xs, ys, zs]).T

                        self.sat_coord[sys]['position'][str(PRN)] = current_sat_coord
                        pbar.update(1) # Update the progress bar for each satellite processed
        self.sat_coord_computed = True
        self.nepochs = len(desired_time)
        return self.sat_coord



if __name__ == "__main__":
    # nav1 = r"C:\Users\perhe\Desktop\TEST VEC\test_5\OPEC00NOR_S_20220010000_01D_GN.rnx"
    # nav2 = r"C:\Users\perhe\Desktop\TEST VEC\test_5\OPEC00NOR_S_20220010000_01D_RN.rnx"
    # nav3 = r"C:\Users\perhe\Desktop\TEST VEC\test_5\OPEC00NOR_S_20220010000_01D_CN.rnx"
    # nav4 = r"C:\Users\perhe\Desktop\TEST VEC\test_5\OPEC00NOR_S_20220010000_01D_EN.rnx"

    # rinnav = [nav1,nav2,nav3,nav4]
    # x_rec = 3149785.9652
    # y_rec = 598260.8822
    # z_rec = 5495348.4927
    # time_epochs = np.array([[2199, 45345]])
    # time_epochs = np.array([[2020,   10,   30,   13,   22,   14],
    #                         [2020,   10,   30,   13,   22,   15],
    #                         [2020,   10,   30,   13,   22,   16]])
    # time_epochs = np.load(r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\Admin\gammel kode ifm vectorisering\time_epochs.npy")
    # CONVERTER = SatelliteEphemerisToECEF(rinnav, x_rec, y_rec, z_rec, data_rate=120)
    # CONVERTER.get_sat_ecef_coordinates(time_epochs[:,1])
    # sat_pos = CONVERTER.compute_satellite_azimut_and_elevation_angle()
    # sat_pos = SatelliteEphemerisToECEF(rinnav, x_rec, y_rec, z_rec).get_sat_ecef_coordinates(time_epochs[:,1], PRN = 'R10')
    # sat_pos = SatelliteEphemerisToECEF(rinnav, x_rec, y_rec, z_rec).get_sat_ecef_coordinates(time_epochs, time_fmt='GREGORIAN' ,PRN = 'R10')


    pass