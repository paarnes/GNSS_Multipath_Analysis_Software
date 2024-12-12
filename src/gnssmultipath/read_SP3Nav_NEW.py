import pandas as pd
import re
from datetime import datetime

class SP3Reader:
    def __init__(self, filepath, coords_in_meter: bool = True, clock_bias_in_sec: bool = True, desiredGNSSsystems: list = ["G", "R", "E", "C"]):
        """
        Initialize the SP3 reader with the file path and various options.
        :param filepath: Path to the SP3 file.
        :param coords_in_meter: Boolean to decide whether to scale coordinates to meters.
        :param clock_bias_in_sec: Boolean to convert clock bias to seconds.
        :param desiredGNSSsystems: List of GNSS systems to include (e.g., ["G", "R", "E", "C"]).
        """
        self.filepath = filepath
        self.coords_in_meter = coords_in_meter
        self.clock_bias_in_sec = clock_bias_in_sec
        self.num_epochs = None
        self.epoch_interval = None
        self.gnss_systems = set()
        self.desiredGNSSsystems = desiredGNSSsystems

    def read_file(self):
        """
        Read the SP3 file, parse metadata, and extract satellite data into a pandas DataFrame.
        :return: pandas DataFrame containing satellite ephemeris data.
        """
        with open(self.filepath, 'r') as file:
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
        df = pd.DataFrame(
            satellite_data, 
            columns=['Epoch', 'Satellite', 'X', 'Y', 'Z', 'Clock Bias']
        )

        return df

    def _parse_header(self, line):
        """
        Parse the header line (##) to extract the number of epochs and the epoch interval.
        :param line: Header line starting with '##'.
        """
        parts = line.split()
        self.num_epochs = int(parts[1])
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
            "Number of Epochs": self.num_epochs,
            "Epoch Interval (s)": self.epoch_interval,
            "GNSS Systems": sorted(self.gnss_systems)
        }

# Example usage
file_path = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\SP3\Testfile_20220101.eph"
sp3_reader = SP3Reader(file_path, coords_in_meter=True, desiredGNSSsystems=["G"])
sp3_df = sp3_reader.read_file()

# Print metadata
metadata = sp3_reader.get_metadata()
print("SP3 Metadata:", metadata)






