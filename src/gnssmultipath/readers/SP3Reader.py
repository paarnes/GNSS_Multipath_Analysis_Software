import pandas as pd
import numpy as np
import re
from typing import Union, List
from datetime import datetime


def _parse_sp3_float(field: str) -> float:
    """Parse a fixed-width SP3 numeric field, returning NaN when it is blank."""
    text = field.strip()
    if not text:
        return np.nan
    try:
        return float(text)
    except ValueError:
        return np.nan


class SP3Reader:
    """
    SP3Reader Class

    A utility class to read and parse SP3 files containing satellite ephemeris data. The class supports
    multiple SP3 files, GNSS system filtering, and coordinate conversion.

    Attributes:
    -----------
    - filepaths: List of SP3 file paths to process.
    - coords_in_meter: Boolean indicating whether coordinates should be converted to meters.
    - clock_bias_in_sec: Boolean indicating whether clock bias should be converted to seconds.
    - num_epochs: Integer indicating the total number of epochs across all SP3 files.
    - epoch_interval: Float representing the interval between epochs (in seconds).
    - gnss_systems: Set containing GNSS systems present in the SP3 files.
    - desiredGNSSsystems: List of GNSS systems to filter and include in the output.

    Example:
    --------
    .. code-block:: python

        file_paths = ["file1.sp3", "file2.sp3"]
        sp3_reader = SP3Reader(file_paths, coords_in_meter=True, desiredGNSSsystems=["G", "R"])
        sp3_df = sp3_reader.read()
        metadata = sp3_reader.get_metadata()
    """

    def __init__(self, filepaths: Union[str, List]=None, coords_in_meter: bool = True, clock_bias_in_sec: bool = True, desiredGNSSsystems: list = ["G", "R", "E", "C"]):
        """
        Initialize the SP3Reader class with optional file paths, coordinate scaling, and GNSS filtering.

        Parameters:
        ----------
        - filepaths: str or list of str, optional. Single SP3 file path or a list of SP3 file paths.
        - coords_in_meter: bool, optional. Scale coordinates to meters if True (default is True).
        - clock_bias_in_sec: bool, optional. Convert clock bias to seconds if True (default is True).
        - desiredGNSSsystems: list of str, optional. GNSS systems to include (default is ["G", "R", "E", "C"]).
        """
        if isinstance(filepaths, str):
            self.filepaths = [filepaths]
        elif isinstance(filepaths, list):
            self.filepaths = filepaths
        else:
            self.filepaths = []

        self.coords_in_meter = coords_in_meter
        self.clock_bias_in_sec = clock_bias_in_sec
        self.num_epochs = 0
        self.epoch_interval = None
        self.gnss_systems = set()
        self.desiredGNSSsystems = desiredGNSSsystems


    def _read_file(self, filepath):
        """
        Read a single SP3 file and extract its satellite data into a DataFrame.

        Parameters:
        ----------
        filepath: str. Path to the SP3 file.

        Returns:
        -------
        pandas.DataFrame: DataFrame containing satellite ephemeris data with columns ['Epoch', 'Satellite', 'X', 'Y', 'Z', 'Clock Bias'].
        """
        with open(filepath, 'r') as file:
            lines = file.readlines()

        satellite_data = []
        current_epoch = None

        for line in lines:
            # Read metadata from the header
            if line.startswith('##'):
                self._parse_header(line)
            elif line.startswith('+'):
                self._parse_gnss_systems(line)
            # Detect epoch line
            elif line.startswith('*'):
                current_epoch = self._parse_epoch(line.strip().split()[1:])
            # Detect satellite data line (e.g., PG01, etc.)
            elif line.startswith('P'):
                satellite = line[1:4]  # Satellite ID (e.g., G01)
                gnss_system = satellite[0]  # GNSS system type (e.g., 'G' for GPS)

                # Skip satellite if not in desired GNSS systems
                if gnss_system not in self.desiredGNSSsystems:
                    continue

                x = _parse_sp3_float(line[4:18])   # X coordinate
                y = _parse_sp3_float(line[18:32])  # Y coordinate
                z = _parse_sp3_float(line[32:46])  # Z coordinate
                clk = _parse_sp3_float(line[46:60])  # Clock bias

                # SP3 marks bad/missing satellite positions with exactly
                # 0.000000 in all three coordinate columns. Without this
                # check, downstream interpolation would treat the geocenter
                # as a real satellite position and corrupt the orbit.
                if x == 0.0 and y == 0.0 and z == 0.0:
                    x = y = z = np.nan

                # Convert coordinates from kilometers to meters if required
                if self.coords_in_meter:
                    x *= 1000
                    y *= 1000
                    z *= 1000

                # Convert clock bias to seconds if required and not a placeholder value.
                # SP3 uses 999999.999999 as the bad-clock sentinel; tolerate
                # small floating-point drift around that value.
                if abs(clk - 999999.999999) < 1e-6:
                    clk = np.nan
                elif self.clock_bias_in_sec:
                    clk *= 1e-6

                satellite_data.append([current_epoch, satellite, x, y, z, clk])

        # Create a DataFrame
        return pd.DataFrame(
            satellite_data,
            columns=['Epoch', 'Satellite', 'X', 'Y', 'Z', 'Clock Bias']
        )

    def read(self):
        """
        Combine satellite data from multiple SP3 files into a single DataFrame.

        Returns:
        -------
        pandas.DataFrame: Combined DataFrame with data from all provided SP3 files.
        """
        if not self.filepaths:
            raise ValueError("No SP3 files provided for combination.")

        combined_data = []
        for filepath in self.filepaths:
            file_data = self._read_file(filepath)
            combined_data.append(file_data)

        # Concatenate all DataFrames and sort by epoch
        combined_df = pd.concat(combined_data, ignore_index=True)
        combined_df.sort_values(by='Epoch', inplace=True)

        return combined_df

    def _parse_header(self, line):
        """
        Parse the header line (##) to extract the number of epochs and the epoch interval.
        Parameter:
        ---------
        line: Header line starting with '##'.
        """
        parts = line.split()
        self.num_epochs += int(parts[1])
        self.epoch_interval = float(parts[3])

    def _parse_gnss_systems(self, line):
        """
        Parse the GNSS systems from the satellite list lines (+).
        :param line: Line starting with '+' containing satellite IDs.
        """
        for satellite in re.findall(r'[A-Z]', line[3:]):
            self.gnss_systems.add(satellite)

    @staticmethod
    def _parse_epoch(epoch_parts):
        """
        Convert SP3 epoch components into a datetime object.

        Parameters:
        ----------
        epoch_parts: list of str. Components of the epoch in the format [YYYY, MM, DD, HH, MM, SS].

        Returns:
        -------
        datetime: Datetime object representing the epoch.
        """
        year, month, day, hour, minute = map(int, epoch_parts[:5])
        second = float(epoch_parts[5])
        microsecond = int((second % 1) * 1e6)
        return datetime(year, month, day, hour, minute, int(second), microsecond)

    def get_metadata(self):
        """
        Retrieve metadata about the parsed SP3 files.

        Returns:
        -------
        dict: Dictionary containing metadata with keys ['n_epochs', 'epoch_interval_sec', 'gnss_systems_in_file', 'desired_gnss_systems'].
        """
        return {
            "n_epochs": self.num_epochs,
            "epoch_interval_sec": self.epoch_interval,
            "gnss_systems_in_file": sorted(self.gnss_systems),
            "desired_gnss_systems": self.desiredGNSSsystems
        }




if __name__ == "__main__":

    # Example usage
    file_path = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\SP3\Testfile_20220101.eph"
    file_path2 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\SP3\com21910.eph"
    file_paths = [file_path, file_path2]
    desiredGNSSsystems= ["G", "R", "E", "C"]
    sp3_reader = SP3Reader(file_paths, coords_in_meter=True, desiredGNSSsystems=desiredGNSSsystems)
    sp3_df = sp3_reader.read()

    # Print metadata
    metadata = sp3_reader.get_metadata()
    print("SP3 Metadata:", metadata)






