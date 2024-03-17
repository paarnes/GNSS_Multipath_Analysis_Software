import pandas as pd
import numpy as np
import os
import sys



class createCSVfile:
    def __init__(self, analysisResults, output_dir, column_delimiter=';'):
        self.analysisResults = analysisResults
        self.output_dir = output_dir
        self.column_delimiter = column_delimiter
        self.results_dict = {}
        self.format_rules = {
            "Azimuth": "{:.2f}",
            "Elevation": "{:.2f}",
            "MP_": "{:.4f}",
            "SNR_": "{:.1f}"
        }
        self.time_stamps = self.analysisResults["ExtraOutputInfo"]["time_epochs_utc_time"]
        self.mp_data_lst = ["PRN; Time_UTC; Elevation; Azimuth"]
        self.GNSSsystemCode2Fullname = {'G': 'GPS', 'R': 'GLONASS', 'E': 'Galileo', 'C': 'BeiDou'}
        self.GNSS_Name2Code = {v: k for k, v in self.GNSSsystemCode2Fullname.items()}
        self.results_dict = self.build_results_dict()


    def flatten_result_array(self, arr):
        """
        Flatten a numpy array to 1D
        """
        flatten_array = arr[:, 1:].T.ravel().tolist()
        return flatten_array

    def set_float_fmt_dataframe(self, df):
        for column in df.columns:
            for prefix, fmt in self.format_rules.items():
                if column.startswith(prefix):
                    df[column] = df[column].map(lambda x: fmt.format(x))
                    break

    def extract_multipath_and_put_in_result_dict(self, sys_name):
        """
        Extract the multipath values from "analysisResults" dictionary
        and gather them in the results_dict.
        """

        for band in self.analysisResults[sys_name]["Bands"]:
            for code in self.analysisResults[sys_name][band]["Codes"]:
                curr_code_data = self.analysisResults[sys_name][band].get(code, None) if not isinstance(code, list) else None
                if curr_code_data is not None:
                    rms_multipath_avg = np.round(curr_code_data["multipath_range1"], 4)
                else:
                    continue
                mp_header = f"MP_{code}"
                self.mp_data_lst.append(mp_header)
                self.results_dict[sys_name][mp_header] = rms_multipath_avg


    def extract_SNR_and_put_in_result_dict(self, sys_name):
        """
        Extract the SNR values from "analysisResults" dictionary
        and gather them in the results_dict.
        """
        SNR_dict = self.analysisResults[sys_name].get("SNR", None)
        if SNR_dict is not None:
            for signal_code, snr_array in SNR_dict.items():
                snr_array[snr_array == 0] = np.nan # convert null to np.nan
                if not np.all(np.isnan(snr_array)):
                    signal_header = f"SNR_{signal_code}"
                    self.mp_data_lst.append(signal_header)
                    # Append SNR data to results_dict if not all elements are np.nan
                    self.results_dict[sys_name][signal_header] = snr_array



    def build_results_dict(self):
        # results_dict = {}
        GNSS_systems = list(self.analysisResults["Sat_position"].keys())

        for gnss_sys in GNSS_systems:
            sys_name = self.GNSSsystemCode2Fullname[gnss_sys]
            self.results_dict[sys_name] = {"Elevation": [], "Azimuth": []}

        # Build up result dict
        for gnss_sys in GNSS_systems:
            curr_sys = self.analysisResults["Sat_position"][gnss_sys]
            sys_name = self.GNSSsystemCode2Fullname[gnss_sys]
            self.results_dict[sys_name]["Azimuth"] = curr_sys["azimuth"]
            self.results_dict[sys_name]["Elevation"] = curr_sys["elevation"]
            
            # Extract multipath and SNR values and store in results_dict
            self.extract_multipath_and_put_in_result_dict(sys_name)
            self.extract_SNR_and_put_in_result_dict(sys_name)


        return self.results_dict

    def write_results_to_csv(self):
        for sys_name, sys_data in self.results_dict.items():
            # Extract data into arrays
            sys_code = self.GNSS_Name2Code[sys_name]
            timestamps = self.time_stamps * sys_data["Elevation"][:, 1:].shape[1]

            prns = list(range(1, sys_data["Azimuth"].shape[1]))
            prns = [f"{sys_code}{prn:02d}" for prn in prns]
            prn_repeated = list(np.repeat(prns, sys_data["Azimuth"].shape[0]))

            # Flatten numpy array to 1D
            az = self.flatten_result_array(np.round(sys_data["Azimuth"], 2))
            el = self.flatten_result_array(np.round(sys_data["Elevation"], 2))

            # Create a DataFrame for the current system
            df = pd.DataFrame({
                "PRN": prn_repeated,
                "Time_UTC": timestamps,
                "Azimuth": az,
                "Elevation": el
            }, dtype=object)

            
            # Add SNR columns to the DataFrame
            if any(key.startswith("MP") for key in sys_data.keys()):
                mp_headers = [header for header in sys_data.keys() if header.startswith("MP_")]
                for i, header in enumerate(mp_headers):
                    df[header] = self.flatten_result_array(sys_data[header])

            # Add SNR columns to the DataFrame
            if any(key.startswith("SNR") for key in sys_data.keys()):
                snr_headers = [header for header in sys_data.keys() if header.startswith("SNR_")]
                for i, header in enumerate(snr_headers):
                    df[header] = self.flatten_result_array(sys_data[header])
                    
            # Set specified float formatting on the columns
            self.set_float_fmt_dataframe(df)

            # Remove rows where all values except PRN and time are np.nan
            df[df.columns[2:]] = df[df.columns[2:]].apply(pd.to_numeric, errors='coerce') # Convert selected columns to numeric
            df.dropna(subset=['Elevation'], inplace=True) # drop rows where the satellite is below the horizon
            output_file = os.path.join(self.output_dir, f"{sys_name}_results.csv")
            print(f'INFO: The result CSV file {output_file} has been written')
            df.to_csv(output_file, index=False, sep=self.column_delimiter)


if __name__ =="__main__":
    sys.path.append(r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\src")
    from gnssmultipath import PickleHandler
    analysisResults = PickleHandler.read_zstd_pickle(r"C:\Users\perhe\Desktop\CSV_export\analysisResults.pkl")
    outputDir = r"C:\Users\perhe\Desktop\CSV_export\TEST"
    createCSV = createCSVfile(analysisResults, outputDir)
    createCSV.write_results_to_csv()
