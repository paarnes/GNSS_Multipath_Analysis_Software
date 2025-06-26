import re
from collections import defaultdict
from datetime import datetime
from collections import Counter
import pandas as pd
import numpy as np
from Geodetic_functions import date2gpstime_vectorized

class Rinex2Reader:
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
    
    def get_gpstime_array(self, time_epochs: list) -> np.ndarray:
        week, tow = date2gpstime_vectorized(np.array(time_epochs))
        return np.stack((week, tow), axis=1) 
    
    def save_dataframe_as_csv(self, df, fname="rinex_data.csv", sep=";"):
        df.to_csv(fname, sep=sep, index=False, encoding='utf-8')
    
    def filter_by_satellite(self, df: pd.DataFrame, satellite: str) -> pd.DataFrame:
        """Filter the DataFrame to include only rows from a specific satellite (e.g. 'G10')."""
        return df[df['sat'] == satellite].copy()
    
    def filter_by_epoch(self, df: pd.DataFrame, epoch_str: str) -> pd.DataFrame:
        """
        Filter the DataFrame to a specific epoch.
        Example format: '2024:12:03:02:23:50.000000'
        """
        return df[df['epoch'] == epoch_str].copy()

    def filter_by_signal_type(self, df: pd.DataFrame, signal_type: str) -> pd.DataFrame:
        """
        Filter the DataFrame to include only columns related to a specific signal type
        (e.g., 'L1', 'C1'). Keeps epoch and sat columns.
        """
        base_cols = ['epoch', 'sat']
        matching_cols = [col for col in df.columns if col.endswith(f"_{signal_type}") or col == signal_type]
        return df[base_cols + matching_cols].copy()
    
    def filter_by_system(self, df: pd.DataFrame, system_letter: str) -> pd.DataFrame:
        """
        Filter the DataFrame to include only rows for a specific GNSS system.
    
        Parameters:
            df (pd.DataFrame): The observation DataFrame.
            system_letter (str): A single character identifying the system ('G', 'E', 'R', etc.)
        """
        return df[df['sat'].str.startswith(system_letter)].copy()

    def parse(self, system_filter=None):
        with open(self.filepath, 'r') as f:
            lines = f.readlines()

        header_end_idx = self._parse_header(lines)

        obs_types = self.header.get('OBS_TYPES', [])
        tInterval = self.header.get('INTERVAL', None)

        obs_idx = header_end_idx + 1
        records = []
        time_epochs = []

        def format_epoch(t):
            return f"{t[0]:04}:{t[1]:02}:{t[2]:02}:{t[3]:02}:{t[4]:02}:{t[5]:09.6f}"

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

            for sv in sats:
                n_obs_lines = (len(obs_types) + 4) // 5
                data = []
                for _ in range(n_obs_lines):
                    data.extend(lines[obs_idx][i:i+16] for i in range(0, 80, 16) if i < len(lines[obs_idx]))
                    obs_idx += 1

                if system_filter and sv[0] not in system_filter:
                    continue  # skip storing, but lines already read

                obs_vals, lli_vals, ss_vals = [], [], []
                for d in data[:len(obs_types)]:
                    val = d[0:14].strip()
                    lli = d[14:15].strip()
                    ss = d[15:16].strip()
                    obs_vals.append(float(val) if val else None)
                    lli_vals.append(int(lli) if lli else 0)
                    ss_vals.append(int(ss) if ss else None)

                row = {
                    'epoch': epoch_str,
                    'sat': sv
                }
                for i, code in enumerate(obs_types):
                    row[code] = obs_vals[i]
                    if code.startswith('L') and lli_vals[i] is not None:
                        row[f'LLI_{code}'] = lli_vals[i]
                    if ss_vals[i] is not None:
                        row[f'SS_{code}'] = ss_vals[i]
                records.append(row)

        # Determine the interval if None
        if tInterval is None:
            tInterval = self.infer_interval_from_epochs(time_epochs)

        df = pd.DataFrame.from_records(records)

        return self.header, df

    def _parse_header(self, lines):
        obs_types = []
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



if __name__=="__main__":

    ## Rinex observation file
    path_to_testdata = "TestData/rin_211/"
    rin_obs = path_to_testdata + 'OPEC00NOR_S_20220010000_01D_30S_MO_v211.obs'

    reader = Rinex2Reader(rin_obs)
    header, df = reader.parse(system_filter=None)
    
    print(len(df['epoch'].unique()))
    

    # from readRinexObs import readRinexObs304

    # rin_obs = r"TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx"
    # GNSS_obs_304, GNSS_LLI_304, GNSS_SS_304, GNSS_SVs_304, time_epochs_304, nepochs_304, GNSSsystems_304,\
    #     obsCodes_304, approxPosition_304, max_sat_304, tInterval_304, markerName_304, rinexVersion_304, recType_304, timeSystem_304, leapSec_304, gnssType_304,\
    #     rinexProgr_304, rinexDate_304, antDelta_304, tFirstObs_304, tLastObs_304, clockOffsetsON_304, GLO_Slot2ChannelMap_304, success_304 = readRinexObs304(rin_obs)

    # print(rin_obs)