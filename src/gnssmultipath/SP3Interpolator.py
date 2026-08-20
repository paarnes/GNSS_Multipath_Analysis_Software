import pandas as pd
import numpy as np
from datetime import datetime
from typing import Literal, Tuple
from gnssmultipath.utils.PickleHandler import PickleHandler
from gnssmultipath.readers.readRinexObs import readRinexObs
from gnssmultipath.readers.SP3Reader import SP3Reader
from gnssmultipath.SatelliteEphemerisToECEF import SatelliteEphemerisToECEF, Kepler2ECEF
from gnssmultipath.Geodetic_functions import date2gpstime, date2gpstime_vectorized, gpstime2date_arrays, gpstime2date_arrays_with_microsec
from gnssmultipath.constants import SPEED_OF_LIGHT as c
from tqdm import tqdm


class SP3Interpolator:
    """
    SP3 files already provide satellite positions in the Earth-Centered Earth-Fixed (ECEF) frame.
    These positions are valid for the provided epoch timestamp and do not require
    an additional correction for Earth's rotation when interpolated directly.

    Example:
    -------
    .. code-block:: python
        interpolator = SP3Interpolator(sp3_df, epochInterval_sp3)
        interpolated_positions = interpolator.interpolate_sat_coordinates(time_epochs, gnss_systems)
    """

    def __init__(self, sp3_dataframe, epoch_interval):
        """
        Initializes the SP3 Interpolator with the provided SP3 DataFrame.

        Parameter:
        ----------
        - sp3_dataframe: Pandas DataFrame containing SP3 data (columns: ['Epoch', 'Satellite', 'X', 'Y', 'Z', 'Clock Bias']).
        - epoch_interval: Interval between each epoch in seconds.
        """
        self.sp3_dataframe = sp3_dataframe
        self.epoch_interval = epoch_interval
        # Lazy cache: PRN -> (sat_seconds, sat_xyz, sat_clock) sorted by epoch.
        # Built once on first access via _get_satellite_arrays()/_build_satellite_arrays_cache().
        self._sat_arrays_cache = None
        # Origin for the interpolation time axis, see _seconds_from_reference().
        self._reference_ns = None

    _GPS_EPOCH_NS = np.datetime64('1980-01-06T00:00:00', 'ns').astype('int64')

    def _seconds_from_reference(self, epoch_ns: np.ndarray) -> np.ndarray:
        """
        Convert integer nanoseconds to float seconds relative to the first SP3 epoch.

        Only differences between epochs matter for the interpolation, and seconds counted
        from an absolute epoch land on the float64 resolution limit: one ULP is 119 ns at
        "seconds since 2000", which is enough to move a GPS satellite half a millimetre and
        makes the result depend on the platform. Counting from the first SP3 epoch keeps the
        magnitude near 1e5 s, where an ULP is ~1e-11 s.
        """
        return (epoch_ns - self._reference_ns) / 1e9

    @classmethod
    def _gpstime_to_ns(cls, weeks: np.ndarray, tows: np.ndarray) -> np.ndarray:
        """GPS week and time-of-week to integer nanoseconds since the GPS epoch."""
        weeks = np.asarray(weeks, dtype=float)
        tows = np.asarray(tows, dtype=float)
        whole_seconds = np.floor(tows)
        # The seconds and the sub-second part are kept apart, so the 0.1 us resolution of the
        # RINEX epoch field survives instead of being rounded away by the week offset.
        return ((weeks * 604800.0 + whole_seconds).astype('int64') * 1_000_000_000
                + np.round((tows - whole_seconds) * 1e9).astype('int64'))

    @staticmethod
    def epoch_to_seconds(epoch):
        """
        Convert an epoch (datetime object) into seconds since 2000-01-01.

        Parameter:
        ----------
        - epoch: Datetime object representing the epoch.

        Return:
        ------
        Total seconds since the reference epoch (January 1st, 2000).
        """
        base_time = datetime(2000, 1, 1)
        delta = epoch - base_time
        return delta.total_seconds()

    @staticmethod
    def _neville_vectorized(x, y):
        """
        Batched Neville's polynomial interpolation.

        Parameter:
        ----------
        - x: Array of shape (n_obs, n) with x-values per observation epoch.
        - y: Array of shape (n_obs, n, k) with y-values for k channels per epoch.

        Return:
        ------
        - Array of shape (n_obs, k) with the interpolated values at x = 0.
        """
        y = y.copy()
        n = x.shape[1]
        for j in range(1, n):
            for i in range(n - j):
                denom = x[:, i + j] - x[:, i]                       # (n_obs,)
                num = (x[:, i + j:i + j + 1] * y[:, i, :]
                       - x[:, i:i + 1] * y[:, i + 1, :])             # (n_obs, k)
                y[:, i, :] = num / denom[:, None]
        return y[:, 0, :]

    def _build_satellite_arrays_cache(self):
        """Build and cache per-PRN sorted (epoch_seconds, xyz, clock) NumPy arrays."""
        df = self.sp3_dataframe
        epochs = df['Epoch'].to_numpy()
        epoch_ns_all = pd.to_datetime(epochs).to_numpy().astype('datetime64[ns]').astype('int64') - self._GPS_EPOCH_NS
        self._reference_ns = int(epoch_ns_all.min()) if epoch_ns_all.size else 0
        epoch_seconds_all = self._seconds_from_reference(epoch_ns_all)

        sats = df['Satellite'].to_numpy()
        xyz_all = df[['X', 'Y', 'Z']].to_numpy()
        clk_all = df['Clock Bias'].to_numpy()

        cache = {}
        # Group indices by PRN (single pass over the array)
        unique_sats, inverse = np.unique(sats, return_inverse=True)
        for k, prn in enumerate(unique_sats):
            mask = inverse == k
            t = epoch_seconds_all[mask]
            order = np.argsort(t, kind='stable')
            cache[str(prn)] = (
                t[order],
                xyz_all[mask][order],
                clk_all[mask][order],
            )
        self._sat_arrays_cache = cache

    def _get_satellite_arrays(self, prn):
        """Return (sat_seconds, sat_xyz, sat_clock) for *prn*, building cache if needed."""
        if self._sat_arrays_cache is None:
            self._build_satellite_arrays_cache()
        if prn not in self._sat_arrays_cache:
            raise ValueError(f"No data found for satellite {prn}.")
        return self._sat_arrays_cache[prn]

    def _observation_seconds_from_time_epochs(self, time_epochs):
        """Convert observation time_epochs (GPS week + TOW) to the interpolation time axis."""
        time_epochs = np.asarray(time_epochs)
        if time_epochs.ndim > 1 and time_epochs.shape[-1] == 2 and time_epochs.shape[0] != 2:
            weeks = time_epochs[:, 0]
            tows  = time_epochs[:, 1]
        elif time_epochs.ndim > 1:
            weeks = time_epochs[0]
            tows  = time_epochs[1]
        else:
            weeks = np.atleast_1d(time_epochs[0])
            tows  = np.atleast_1d(time_epochs[1])

        if self._sat_arrays_cache is None:
            self._build_satellite_arrays_cache()
        return self._seconds_from_reference(self._gpstime_to_ns(weeks, tows))

    def _interpolate_satellite_vec(self, sat_seconds, sat_xyz, sat_clock,
                                   obs_seconds, n_interpol_points=7):
        """
        Vectorized interpolation of (X, Y, Z, clock_bias) for one satellite at *obs_seconds*.

        Parameter:
        ----------
        - sat_seconds: (n_sp3,) sorted SP3 epoch seconds for the satellite.
        - sat_xyz:     (n_sp3, 3) SP3 positions for the satellite.
        - sat_clock:   (n_sp3,)  SP3 clock biases for the satellite.
        - obs_seconds: (n_obs,)  observation epoch seconds.
        - n_interpol_points: number of nearest SP3 points to use (default 7).

        Return:
        ------
        - (positions, clock_biases) with shapes (n_obs, 3) and (n_obs,).
        """
        n_sp3 = sat_seconds.shape[0]
        n = min(n_interpol_points, n_sp3)
        n_obs = obs_seconds.shape[0]

        # Find nearest n SP3 indices per observation (O(n_sp3 * n_obs), but cheap in NumPy)
        abs_diffs = np.abs(sat_seconds[:, None] - obs_seconds[None, :])  # (n_sp3, n_obs)
        if n_sp3 > n:
            nearest = np.argpartition(abs_diffs, n - 1, axis=0)[:n, :]
        else:
            nearest = np.broadcast_to(np.arange(n_sp3)[:, None], (n_sp3, n_obs)).copy()
        # Sort along window axis so x is monotonic (better numerical behaviour)
        nearest = np.sort(nearest, axis=0)                                # (n, n_obs)

        # Build x of shape (n_obs, n) and y of shape (n_obs, n, 4)
        x = (sat_seconds[nearest] - obs_seconds[None, :]).T               # (n_obs, n)
        gathered_xyz = sat_xyz[nearest]                                   # (n, n_obs, 3)
        gathered_clk = sat_clock[nearest]                                 # (n, n_obs)
        y = np.empty((n_obs, n, 4), dtype=np.float64)
        y[:, :, 0:3] = np.transpose(gathered_xyz, (1, 0, 2))
        y[:, :, 3]   = gathered_clk.T

        result = self._neville_vectorized(x, y)                           # (n_obs, 4)
        return result[:, 0:3], result[:, 3]

    def compute_relativistic_correction_single_sat(self, prn, time_epochs):
        """
        Compute the relativistic clock correction for a single satellite.

        Parameter:
        ----------
        - prn: PRN of the satellite (e.g., 'G12').
        - time_epochs: Array of observation times in GPS time format (week, TOW).

        Return:
        ------
        Relativistic clock correction values as a NumPy array.
        """
        delta_t = self.epoch_interval  # Time step for velocity approximation (seconds)

        # Times for velocity computation
        t_delta_plus = time_epochs + np.array([[0, delta_t]]).T
        t_delta_minus = time_epochs - np.array([[0, delta_t]]).T

        # Interpolate satellite positions at t+delta_t, t, and t-delta_t
        pos_plus, _ = self.interpolate_single_satellite(prn, t_delta_plus)
        pos_minus, _ = self.interpolate_single_satellite(prn, t_delta_minus)
        pos_now, _ = self.interpolate_single_satellite(prn, time_epochs)

        # Compute satellite velocities
        velocities = (pos_plus - pos_minus) / (2 * delta_t)

        # Compute relativistic corrections
        corrections = -2 * np.sum(pos_now * velocities, axis=1) / c**2

        return corrections

    def interpolate_single_satellite(self, prn, time_epochs, n_interpol_points=7):
        """
        Interpolates satellite positions and clock biases for a single satellite specified by its PRN.

        Parameter:
        ----------
        - prn: PRN of the satellite to interpolate (e.g., 'G12').
        - time_epochs: Array of observation times in GPS time format (week, TOW).
        - n_interpol_points: Number of nearest points for interpolation.

        Return:
        ------
        - Interpolated positions and clock biases as a dictionary.
        """
        time_epochs = np.asarray(time_epochs)
        observation_seconds = self._observation_seconds_from_time_epochs(time_epochs)
        sat_seconds, sat_xyz, sat_clock = self._get_satellite_arrays(prn)
        return self._interpolate_satellite_vec(
            sat_seconds, sat_xyz, sat_clock, observation_seconds, n_interpol_points
        )


    def interpolate_sat_coordinates(self, time_epochs, gnss_systems, n_interpol_points=7, output_format: Literal["pd.DataFrame", "dict"] = "pd.DataFrame"):
        """
        Interpolates satellite positions and clock biases for all systems and satellites for given time epochs.

        Parameter:
        ----------
        - time_epochs: Array of observation times in GPS time format (week, TOW).
        - gnss_systems: List of GNSS systems to include (e.g., ['G', 'R', 'E']).
        - n_interpol_points: Number of nearest points for interpolation.
        - output_format: Desired output format. Options are 'dict' or 'dataframe'.


        Return:
        -------
        - Interpolated positions and clock biases in the specified output format.
        """

        # Convert GPS time to datetime objects (kept for the optional DataFrame output below)
        if len(time_epochs) > 2:
            observation_times = gpstime2date_arrays_with_microsec(time_epochs[:, 0], time_epochs[:, 1])
        else:
            observation_times = gpstime2date_arrays_with_microsec(time_epochs[0], time_epochs[1])

        # Convert observation times to seconds since the reference epoch (vectorized)
        observation_seconds = self._observation_seconds_from_time_epochs(np.asarray(time_epochs))

        # Build per-satellite NumPy cache once (cheap if already built)
        if self._sat_arrays_cache is None:
            self._build_satellite_arrays_cache()

        # Create a dictionary to store results for all GNSS systems
        interpolated_positions = {}

        # Pre-group cached PRNs by GNSS system letter
        prns_by_system = {g: sorted([p for p in self._sat_arrays_cache if p[:1] == g])
                          for g in gnss_systems}
        total_satellites = sum(len(v) for v in prns_by_system.values())
        bar_format = '{desc}: {percentage:3.0f}%|{bar}| ({n_fmt}/{total_fmt})'
        pbar = tqdm(total=total_satellites, desc="Interpolating satellite coordinates",
                    position=0, leave=True, bar_format=bar_format)

        for gnss in gnss_systems:
            interpolated_positions[gnss] = {}
            for satellite in prns_by_system[gnss]:
                sat_seconds, sat_xyz, sat_clock = self._sat_arrays_cache[satellite]
                positions, clock_bias = self._interpolate_satellite_vec(
                    sat_seconds, sat_xyz, sat_clock,
                    observation_seconds, n_interpol_points,
                )
                interpolated_positions[gnss][satellite] = {
                    "positions": positions,
                    "clock_bias": clock_bias,
                }
                pbar.update(1)

        pbar.close()

        # If output_format is 'dataframe', convert the dictionary to a DataFrame
        if output_format == "pd.DataFrame":
            all_epochs = []
            all_sats = []
            all_x = []
            all_y = []
            all_z = []
            all_clk = []
            epoch_datetimes = [datetime(*t) for t in observation_times]
            for gnss, satellites in interpolated_positions.items():
                for satellite, data in satellites.items():
                    positions = data["positions"]
                    clock_biases = data["clock_bias"]
                    n = len(epoch_datetimes)
                    all_epochs.append(epoch_datetimes)
                    all_sats.append([satellite] * n)
                    all_x.append(positions[:, 0])
                    all_y.append(positions[:, 1])
                    all_z.append(positions[:, 2])
                    all_clk.append(clock_biases)

            return pd.DataFrame({
                "Epoch": np.concatenate(all_epochs) if all_epochs else [],
                "Satellite": np.concatenate(all_sats) if all_sats else [],
                "X": np.concatenate(all_x) if all_x else [],
                "Y": np.concatenate(all_y) if all_y else [],
                "Z": np.concatenate(all_z) if all_z else [],
                "Clock Bias": np.concatenate(all_clk) if all_clk else [],
            })

        return interpolated_positions


    @staticmethod
    def filter_by_prn(interpolated_data, prn_list:list):
        """
        Filters the interpolated data for specific PRN numbers.

        Parameter:
        ----------
        - interpolated_data: Interpolated data as a dictionary or DataFrame.
        - prn_list: List of PRN numbers to include (e.g., ['G01', 'E02', 'R03']).

        Return:
        ------
        - Filtered data as the same type as input (dict or DataFrame).
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

        Parameter:
        ----------
        - interpolated_data: Interpolated data as a dictionary or DataFrame.
        - gnss_systems: List of GNSS systems to include (e.g., ['G', 'R', 'E']).

        Return:
        ------
        - Filtered data as the same type as input (dict or DataFrame).
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

    obs_data = readRinexObs(rinObs)
    GNSS_obs = obs_data.GNSS_obs
    time_epochs = obs_data.time_epochs
    nepochs = obs_data.nepochs
    GNSSsystems = obs_data.GNSSsystems
    obsCodes = obs_data.obsCodes
    tInterval = obs_data.tInterval

    observation_times = gpstime2date_arrays(time_epochs[:,0],time_epochs[:,1])

    # sat_coord_sp3, epoch_dates_sp3, navGNSSsystems_sp3, nEpochs_sp3, epochInterval_sp3, _= readSP3Nav(sp3)
    sp3_reader = SP3Reader(sp3, coords_in_meter=True, desiredGNSSsystems=gnss_systems)
    sp3_df = sp3_reader.read()

    # Print metadata
    metadata = sp3_reader.get_metadata()
    nEpochs_sp3 = metadata["n_epochs"]
    epochInterval_sp3 = metadata["epoch_interval_sec"]
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


    # interpolated_positions = interpolator.interpolate_sat_coordinates(time_epochs, gnss_systems, output_format="dict")
    # diff = interpolated_positions["G"]["G01"]["positions"] - results_rnav["Sat_position"]["G"]["position"]["1"][0:440]
    # diff_E = interpolated_positions["E"]["E01"]["positions"] - results_rnav["Sat_position"]["E"]["position"]["1"][0:440]
    # diff_R = interpolated_positions["R"]["R01"]["positions"] - results_rnav["Sat_position"]["R"]["position"]["1"][0:440]
    # print(diff_R)





