import pandas as pd
import re
from datetime import datetime

class SP3Reader:
    def __init__(self, filepaths=None, coords_in_meter: bool = True, clock_bias_in_sec: bool = True, desiredGNSSsystems: list = ["G", "R", "E", "C"]):
        """
        Initialize the SP3 reader with file paths and various options.

        :param filepaths: Single SP3 file path (string) or a list of SP3 file paths.
        :param coords_in_meter: Boolean to decide whether to scale coordinates to meters.
        :param clock_bias_in_sec: Boolean to convert clock bias to seconds.
        :param desiredGNSSsystems: List of GNSS systems to include (e.g., ["G", "R", "E", "C"]).
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
        Read a single SP3 file and extract satellite data into a pandas DataFrame.

        :param filepath: Path to the SP3 file.
        :return: pandas DataFrame containing satellite ephemeris data.
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

                x = float(line[4:18].strip())  # X coordinate
                y = float(line[18:32].strip())  # Y coordinate
                z = float(line[32:46].strip())  # Z coordinate
                clk = float(line[46:60].strip())  # Clock bias

                # Convert coordinates from kilometers to meters if required
                if self.coords_in_meter:
                    x *= 1000
                    y *= 1000
                    z *= 1000

                # Convert clock bias to seconds if required and not a placeholder value
                if self.clock_bias_in_sec and clk != 999999.999999:
                    clk *= 1e-6

                satellite_data.append([current_epoch, satellite, x, y, z, clk])

        # Create a DataFrame
        return pd.DataFrame(
            satellite_data,
            columns=['Epoch', 'Satellite', 'X', 'Y', 'Z', 'Clock Bias']
        )

    def read(self):
        """
        Combines multiple SP3 files into a single DataFrame.

        :return: DataFrame containing combined SP3 data from all files.
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
        :param line: Header line starting with '##'.
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
        Parse the epoch components into a datetime object.
        :param epoch_parts: List of epoch components [YYYY, MM, DD, HH, MM, SS].
        :return: datetime object representing the epoch.
        """
        year, month, day, hour, minute = map(int, epoch_parts[:5])
        second = float(epoch_parts[5])
        microsecond = int((second % 1) * 1e6)
        return datetime(year, month, day, hour, minute, int(second), microsecond)

    def get_metadata(self):
        """
        Get the parsed metadata: number of epochs, epoch interval, and GNSS systems.
        :return: Dictionary containing metadata.
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






