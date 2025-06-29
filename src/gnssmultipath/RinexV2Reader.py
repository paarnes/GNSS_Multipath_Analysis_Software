import re
from typing import Tuple
from datetime import datetime
from collections import Counter, defaultdict
import pandas as pd
import numpy as np
from gnssmultipath.Geodetic_functions import date2gpstime_vectorized
from tqdm.auto import tqdm
from itertools import chain


class RinexUtils:

    @staticmethod
    def get_gpstime_array(time_epochs: list) -> np.ndarray:
        week, tow = date2gpstime_vectorized(np.array(time_epochs))
        return np.stack((week, tow), axis=1) 
    @staticmethod
    def save_dataframe_as_csv(df, fname="rinex_data.csv", sep=";"):
        df.to_csv(fname, sep=sep, index=False, encoding='utf-8')
    
    @staticmethod
    def filter_by_satellite(df: pd.DataFrame, satellite: str) -> pd.DataFrame:
        """Filter the DataFrame to include only rows from a specific satellite (e.g. 'G10')."""
        return df[df['sat'] == satellite].copy()
    
    @staticmethod
    def filter_by_epoch(df: pd.DataFrame, epoch_str: str) -> pd.DataFrame:
        """
        Filter the DataFrame to a specific epoch.
        Example format: '2024:12:03:02:23:50.000000'
        """
        return df[df['epoch'] == epoch_str].copy()
    
    @staticmethod
    def filter_by_signal_type(df: pd.DataFrame, signal_type: str) -> pd.DataFrame:
        """
        Filter the DataFrame to include only columns related to a specific signal type
        (e.g., 'L1', 'C1'). Keeps epoch and sat columns.
        """
        base_cols = ['epoch', 'sat']
        matching_cols = [col for col in df.columns if col.endswith(f"_{signal_type}") or col == signal_type]
        return df[base_cols + matching_cols].copy()
    
    @staticmethod
    def filter_by_system(df: pd.DataFrame, system_letter: str) -> pd.DataFrame:
        """
        Filter the DataFrame to include only rows for a specific GNSS system.
    
        Parameters:
            df (pd.DataFrame): The observation DataFrame.
            system_letter (str): A single character identifying the system ('G', 'E', 'R', etc.)
        """
        return df[df['sat'].str.startswith(system_letter)].copy()

    @staticmethod
    def create_gnss_obs_structure(df: pd.DataFrame, obs_types: list) -> dict:
        """
        Create a nested dictionary GNSS_obs[system][epoch_index][PRN] = np.array([observations]).
        
        Parameters:
            df (pd.DataFrame): DataFrame with 'epoch', 'sat', and observation code columns.
            obs_types (List[str]): List of observation types (e.g. ['C1', 'L1', ...]).

        Returns:
            dict: Nested GNSS observation structure.
        """
        from collections import defaultdict
        import numpy as np
        import pandas as pd

        GNSS_obs = defaultdict(lambda: defaultdict(dict))
        unique_epochs = sorted(df['epoch'].unique())
        epoch_index_map = {epoch: idx for idx, epoch in enumerate(unique_epochs)}

        for _, row in df.iterrows():
            sat = row['sat']
            system = sat[0]  # e.g., 'G'
            prn = int(sat[1:])  # e.g., '10' from 'G10'
            epoch_idx = epoch_index_map[row['epoch']]
            
            obs_values = np.array([
                row.get(code, 0.0) if pd.notna(row.get(code)) else 0.0
                for code in obs_types
            ])

            GNSS_obs[system][epoch_idx][prn] = obs_values

        return GNSS_obs
    



class RinexV2Reader:

    EPOCH_LINE_PATTERN = re.compile(r'^\s*\d{2,4}\s+\d{2}\s+\d{2}\s+\d{2}\s+\d{2}\s+\d{1,2}\.\d+')

    def __init__(self, filepath):
        self.filepath = filepath
        self.header = {}

    def infer_interval_from_epochs(self, time_epochs):
        if len(time_epochs) < 2:
            return None  # not enough data

        def to_datetime(t):
            return datetime(t[0], t[1], t[2], t[3], t[4], int(t[5])) \
                .timestamp() + (t[5] % 1)

        times = [to_datetime(t) for t in time_epochs[:10]]
        diffs = [round(times[i+1] - times[i], 3) for i in range(len(times) - 1)]
        interval, _ = Counter(diffs).most_common(1)[0]
        return interval
    
    def get_nepochs(self, lines):
        return sum(1 for line in lines if self.EPOCH_LINE_PATTERN.match(line))

    def parse(
        self,
        readSS=False,
        readLLI=True,
        includeAllGNSSsystems=True,
        includeAllObsCodes=True,
        desiredGNSSsystems=None,
        desiredObsCodes=None,
        desiredObsBands=None
    ) -> Tuple[dict, pd.DataFrame]:
        
        """
        Parse the RINEX 2.x observation file and return the header and the 
        observation data as a pandas DataFrame.

        Parameter:
        -----------
        - readSS (bool, optional): If True, include Signal Strength (SS) columns (1=weak, 9=very strong).
        - readLLI (bool, optional): If True, include Loss of Lock Indicator (LLI) columns.
        - includeAllGNSSsystems (bool, optional): If True, include all GNSS systems. If False, filter by desiredGNSSsystems.
        - includeAllObsCodes (bool, optional): If True, include all observation codes. If False, filter by desiredObsCodes or desiredObsBands.
        - desiredGNSSsystems (list, optional): List of GNSS system letters to include (e.g., ['G', 'E']).
        - desiredObsCodes (list, optional): List of observation codes to include (e.g., ['C1', 'L1']).
        - desiredObsBands (list, optional): List of band identifiers to include (e.g., ['1', '5']).

        Returns:
        --------
        - header: Parsed RINEX header as a dictionary.
        - df: DataFrame with columns for epoch, satellite, and selected observation codes.
            If readLLI/readSS is True, includes LLI_*/SS_* columns.

        """
        
        with open(self.filepath, 'r') as f:
            lines = f.readlines()
        header_end_idx = self._parse_header(lines)

        full_obs_types = self.header.get('OBS_TYPES', [])
        tInterval = self.header.get('INTERVAL', None)
    

        # GNSS system filtering (e.g., ['G', 'E'])
        system_filter = None
        if not includeAllGNSSsystems and desiredGNSSsystems:
            system_filter = desiredGNSSsystems

        # Observation code filtering (e.g., ['C1', 'L1'])
        obs_types = full_obs_types
        if not includeAllObsCodes and desiredObsCodes:
            obs_types = [code for code in full_obs_types if code[0] in desiredObsCodes]
        elif not includeAllObsCodes and desiredObsBands:
            # Filter based on band prefixes like '1', '5', etc.
            obs_types = [code for code in full_obs_types if len(code) > 1 and code[1] in desiredObsBands]


        obs_idx = header_end_idx + 1
        records = []
        time_epochs = []
        GNSS_obs = defaultdict(lambda: defaultdict(dict))
        GNSS_LLI = defaultdict(lambda: defaultdict(dict))
        GNSS_SS = defaultdict(lambda: defaultdict(dict))
        GNSS_SVs = defaultdict(list)
        epoch_index_map = {}
        nepochs = self.get_nepochs(lines[header_end_idx::])
        self.header["NEPOCHS"] = nepochs

        def format_epoch(t):
            return f"{t[0]:04}:{t[1]:02}:{t[2]:02}:{t[3]:02}:{t[4]:02}:{t[5]:09.6f}"

        bar_format = '{desc}: {percentage:3.0f}%|{bar}| ({n_fmt}/{total_fmt})'
        pbar = tqdm(total=nepochs, desc="Parsing epochs",position=0, leave=True,  bar_format=bar_format)

        while obs_idx < len(lines):
            line = lines[obs_idx]
            if len(line) < 32:
                obs_idx += 1
                continue

            dt = self._parse_epoch_time(line[:26])
            epoch_str = format_epoch(dt)
            nsats = int(line[29:32])
            sats = []

            # First line may have some sats after column 32
            sats_line = line[32:].rstrip()
            sats.extend([sats_line[i:i+3] for i in range(0, len(sats_line), 3)])

            # Read extra lines if needed
            while len(sats) < nsats:
                obs_idx += 1
                extra_line = lines[obs_idx].rstrip()
                sats.extend(re.findall(r'[A-Z]\d{2}', extra_line))

            obs_idx += 1  # Move to start of observation data

            time_epochs.append(dt)
            epoch_idx = len(epoch_index_map)
            epoch_index_map[epoch_str] = epoch_idx
            GNSS_SVs[epoch_idx] = sats  # Store satellites per epoch

            pbar.update(1)  # Progress bar increment

            for sv in sats:
                n_obs_lines = (len(full_obs_types) + 4) // 5
                data_lines = lines[obs_idx:obs_idx + n_obs_lines]
                obs_idx += n_obs_lines  # Always skip based on full obs types

                if system_filter and sv[0] not in system_filter:
                    continue  # skip storing

                data = []
                for line in data_lines:
                    data.extend([line[i:i + 16] for i in range(0, 80, 16) if i < len(line)])

                obs_vals_all, lli_vals_all, ss_vals_all = [], [], []
                for d in data[:len(full_obs_types)]:
                    val = d[0:14].strip()
                    lli = d[14:15].strip()
                    ss = d[15:16].strip()
                    obs_vals_all.append(float(val) if val else np.nan)
                    lli_vals_all.append(int(lli) if lli else np.nan)
                    ss_vals_all.append(int(ss) if ss else np.nan)

                row = {'epoch': epoch_str, 'sat': sv}
                for i, code in enumerate(full_obs_types):
                    if code not in obs_types:
                        continue  # only store selected codes
                    row[code] = obs_vals_all[i]
                    if code.startswith('L') and lli_vals_all[i] is not None and readLLI:
                        row[f'LLI_{code}'] = lli_vals_all[i]
                    if ss_vals_all[i] is not None and readSS:
                        row[f'SS_{code}'] = ss_vals_all[i]

                # Build the GNSS_obs structure
                if epoch_str not in epoch_index_map:
                    epoch_index_map[epoch_str] = len(epoch_index_map)

                epoch_idx = epoch_index_map[epoch_str]
                system = sv[0]
                prn = int(sv[1:])

                obs_array = np.array([
                    obs_vals_all[i] if code in obs_types and pd.notna(obs_vals_all[i]) else np.nan
                    for i, code in enumerate(full_obs_types)
                    if code in obs_types
                ])

                # Prepare LLI and SS arrays 
                if readLLI:
                    lli_array = np.array([
                        lli_vals_all[i] if code in obs_types else np.nan
                        for i, code in enumerate(full_obs_types)
                        if code in obs_types
                    ], dtype=np.float64)

                if readSS:
                    ss_array = np.array([
                        ss_vals_all[i] if code in obs_types else np.nan
                        for i, code in enumerate(full_obs_types)
                        if code in obs_types
                    ], dtype=np.float64)

                # Populate nested dicts
                system = sv[0]
                try:
                    prn = int(sv[1:])
                except ValueError:
                    continue  # skip invalid PRNs

                GNSS_obs[system][epoch_idx][prn] = obs_array

                if readLLI:
                    GNSS_LLI[system][epoch_idx][prn] = lli_array
                if readSS:
                    GNSS_SS[system][epoch_idx][prn] = ss_array

                records.append(row)
        pbar.close()

        # Determine the interval if None
        if tInterval is None:
            tInterval = self.infer_interval_from_epochs(time_epochs)

        self.header['TIME_EPOCHS'] = time_epochs
        df = pd.DataFrame.from_records(records)

        return self.header, df, GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs

    def _parse_header(self, lines):
        obs_types = []
        self.header['GNSS_SYSTEMS'] = self.get_gnss_systems_in_file(lines)
        obs_types_dict = {}
        sys_idx = 1
        for i, line in enumerate(lines):
            label = line[60:].strip()
            if label == 'RINEX VERSION / TYPE':
                self.header['RINEX_VERSION'] = float(line[:9])
                self.header['GNSS_TYPE'] = line[40:41]
            elif label == 'PGM / RUN BY / DATE':
                self.header['RINEX_PROGR'] = line[:20].strip()
                self.header['RINEX_DATE'] = line[40:60].strip()
            elif label == 'MARKER NAME':
                self.header['MARKER_NAME'] = line[:60].strip()
            elif label == 'APPROX POSITION XYZ':
                self.header['APPROX_POSITION'] = tuple(map(float, line[:60].split()))
            elif label == 'ANTENNA: DELTA H/E/N':
                self.header['ANTENNA_DELTA'] = tuple(map(float, line[:60].split()))
            elif label == '# / TYPES OF OBSERV':
                obs_types += line[6:60].split()
                for sys in self.header['GNSS_SYSTEMS']:
                    obs_types_dict[sys_idx] = {sys: obs_types}
                    sys_idx += 1
            elif label == 'INTERVAL':
                self.header['INTERVAL'] = float(line[:10])
            elif label == 'TIME OF FIRST OBS':
                self.header['TIME_FIRST_OBS'] = self._parse_time_fields(line[:43])
                self.header['TIME_SYSTEM'] = line[48:51].strip()
            elif label == 'TIME OF LAST OBS':
                self.header['TIME_LAST_OBS'] = self._parse_time_fields(line[:43])
            elif label == 'RCV CLOCK OFFS APPL':
                self.header['RCV_CLOCK_OFFS_APPL'] = int(line[:6])
            elif label == 'LEAP SECONDS':
                self.header['LEAP_SECONDS'] = int(line[:6])
            elif label == 'REC # / TYPE / VERS':
                self.header['REC_TYPE'] = line[20:40].strip()
            elif label == 'END OF HEADER':
                self.header['OBS_TYPES'] = obs_types
                self.header['GNSS_OBS_CODES'] = obs_types_dict
                return i
        return i

    def _parse_time_fields(self, line):
        fields = line.strip().split()[:6]
        year, month, day, hour, minute = map(int, fields[:5])
        second = float(fields[5])
        return (year, month, day, hour, minute, second)

    def _parse_epoch_time(self, line):
        fields = line.strip().split()[:6]
        year, month, day, hour, minute = map(int, fields[:5])
        second = float(fields[5])
        if year < 80:  # handle Y2K case
            year += 2000
        elif year < 100:
            year += 1900
        return (year, month, day, hour, minute, second)


    # def get_gnss_systems_in_file(self, GNSS_SVs: dict) -> list:
    #     """
    #     Optimized version: extract unique GNSS system identifiers from GNSS_SVs dict.
    #     """
    #     systems = {sv[0] for sv in chain.from_iterable(GNSS_SVs.values()) if isinstance(sv, str) and sv}
    #     systems_list = sorted(systems)
    #     self.header["GNSS_SYSTEMS"] = systems_list
    #     return systems_list

    def get_gnss_systems_in_file(self, lines: list) -> list:
        """
        Fast detection of unique GNSS system identifiers (e.g., 'G', 'E', 'R') 
        using a single regex pass on all RINEX lines.
        """
        joined = ''.join(lines)
        systems = sorted(set(re.findall(r'([A-Z])\d{2}', joined)))
        self.header["GNSS_SYSTEMS"] = systems
        return systems


if __name__=="__main__":

    ## Rinex observation file
    path_to_testdata = "TestData/rin_211/"
    rin_obs = path_to_testdata + 'OPEC00NOR_S_20220010000_01D_30S_MO_v211.obs'
    # rin_obs = path_to_testdata + 'p0803430.24o'

    reader = RinexV2Reader(rin_obs)
    # system_filter = ['E']  # Example filter for GPS, Galileo, and GLONASS
    # obs_filter = ['C1', 'L1']  # Example filter for specific observation codes
    # header, df = reader.parse(system_filter=system_filter, obs_filter=obs_filter)
    header, df, GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs = reader.parse(
        # readSS=False,
        # readLLI=True,
        # includeAllGNSSsystems=False,
        # includeAllObsCodes=False,
        # desiredGNSSsystems=['G', 'E'],
        # desiredObsCodes=None,
        # desiredObsBands=['1', '5']  # Example for L1 and L5 bands
    )
    print(df)
    
