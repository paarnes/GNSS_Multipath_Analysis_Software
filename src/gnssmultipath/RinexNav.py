"""
This module contains a parent class called RinexNav with two subclasess for reading RINEX
navigation files.

Example on how to use it:
from gnssmultipath import Rinex_v3_Reader
navdata = Rinex_v3_Reader().read_rinex_nav(nav_file)



Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""

import re
import warnings
from datetime import datetime, timedelta
import numpy as np
from pandas import DataFrame
from tqdm import tqdm

warnings.filterwarnings("ignore")

class RinexNav:

    def __init__(self):
        self.dataframe = None
        self.block_len = {'G' : 7, 'R': 3, 'E' : 7, 'C': 7}
        self.glo_fcn = None


    def filter_data_rinex_nav(self, filename, desired_GNSS, data_rate=30):
        """
        This function is creating a list of lists that contain ephemeride data
        for the desired GNSS systems only. The rate of data can be set by the user.
        Default value is 30 min.
        """
        pattern = r'^[' + ''.join(desired_GNSS) + ']\d{2}'
        with open(filename, "r") as f:
            lines = f.readlines()
        desired_lines = [lines[idx:idx + self.block_len[line[0]] + 1] for idx, line in enumerate(lines) if re.match(pattern, line) is not None]
        desired_lines = self.filter_ephemeris_data_on_time(desired_lines, time_difference_minutes=data_rate)  # remove epochs with time difference less than specified time
        return desired_lines


    def filter_ephemeris_data_on_time(self, ephemeris_list, time_difference_minutes=30):
        """
        Filter rinex nav data based on time. The desired data rate can be set
        by time_difference_minutes. Default value is 30 min.
        """
        filtered_data = []
        time_difference = timedelta(minutes=time_difference_minutes)
        previous_epoch_str = ephemeris_list[0][0].split()[1:6]
        previous_epoch = datetime(*map(int, previous_epoch_str))
        prev_sys = ephemeris_list[0][0].split()[0][0]
        prev_sat = ephemeris_list[0][0].split()[0]
        for ephemeris_sublist in ephemeris_list:
            current_sys = ephemeris_sublist[0].split()[0][0]
            current_sat = ephemeris_sublist[0].split()[0]
            epoch_str = ephemeris_sublist[0].split()[1:6]
            epoch = datetime(*map(int, epoch_str))
            if prev_sys != current_sys: # check if new system
                filtered_data.append(ephemeris_sublist)
                previous_epoch = epoch
                prev_sys = current_sys
            elif prev_sys == current_sys and prev_sat != current_sat: # check if new satellite within same sys
                filtered_data.append(ephemeris_sublist)
                previous_epoch = epoch
                prev_sat = current_sat

            # Calculate the time difference between the current epoch and the previous one
            # Appends the data if the time diff is greater than the set limit, or if it is the first epoch (filtered_data is empty)
            time_diff = epoch - previous_epoch
            if time_diff >= time_difference or not filtered_data:
                filtered_data.append(ephemeris_sublist)
                previous_epoch = epoch

        return filtered_data


    def extract_glonass_fcn_from_rinex_nav(self, data_array):
        """
        Extract the GLONASS frequency channels numbers (FCN) from
        RINEX navigation file.

        Input:
        ------
        data_array: numpy array containing ephemerides for all avaible systems

        Output:
        ------
        fcn_dict: dictionary with PRN as keys and FCN as values {'R01': 1,'R02': -4, 'R03': 5,...}
        """
        # Get the PRN and FCN columns
        prn_column = data_array[:, 0]  # PRN numbers have index 0
        glo_rows = np.char.startswith(prn_column, 'R') # Create array with GLONASS data only
        glo_data = data_array[glo_rows]
        unique_prns_name, unique_glo_prns_idx = np.unique(glo_data[:, 0], return_index=True)
        PRN_data = glo_data[unique_glo_prns_idx] # ephemeride array with glonass only
        # Create dictionary with PRN as keys and FCN as values
        fcn_dict = {int(prn[1::]): int(fcn) for prn, fcn in zip(PRN_data[:,0], PRN_data[:,17].astype(float).astype(int))}
        return fcn_dict



    def read_header_lines(self, filename):
        """
        Read content of the header of a RINEX navigation file.
        """
        header = []
        try:
            with open(filename, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    line = line.rstrip()
                    header.append(line)
                    if 'END OF HEADER' in line:
                        break
        except OSError as e:
            print(f"Could not open/read file: {filename}\nError: {e}")
            return
        return header






class Rinex_v2_Reader(RinexNav):
    """
    Class for reading RINEX v2 navigation files.

    Example on how to use it:
    from gnssmultipath import Rinex_v2_Reader
    navdata = Rinex_v2_Reader().read_rinex_nav(nav_file)

    """
    def __init__(self):
        super().__init__()

    def read_rinex_nav(self, filename, dataframe = None):
        """
        Reads the navigation message from GPS broadcast efemerids in RINEX v.2 format. (GPS only)

        Reads one navigation message at a time until the end of the row. Accumulate in
        in a common matrix, "data", where there is a line for each message.
        Note that each message forms a block in the file and that the same satellite
        can have several messages, usually with an hour's difference in reference time.


        Parameters
        ----------
        filename : Filename of the RINEX navigation file
        dataframe : Set to 'yes' or 'YES' to get the data output as a pandas DataFrame (array as default)

        Returns
        -------
        data : Matrix with data for all epochs
        header: List with header content
        n_eph: Number of epochs


        """


        try:
            print('Reading broadcast ephemeris from RINEX-navigation file.....')
            filnr = open(filename, 'r')
        except OSError:
            print("Could not open/read file: %s", filename)


        line = filnr.readline().rstrip()
        header = []
        while 'END OF HEADER' not in line:
            line = filnr.readline().rstrip()
            header.append(line)

        data  = np.zeros((1,36))

        while line != '':
            block_arr = np.array([])

            ## -- Read first line of navigation message
            line = filnr.readline().rstrip()

            # Replacing 'D' with 'E' ('D' is fortran syntax for exponentiall form)
            line = line.replace('D','E')
            ## -- Have to add space between datacolums where theres no whitespace
            for idx, val in enumerate(line):
                if line[0:2] != ' ' and line[22] != ' ':
                    line = line[:22] + " " + line[22:]
                if line[idx] == 'e' or line[idx] == 'E' and idx !=0:
                    line = line[:idx+4] + " " + line[idx+4:]

            fl = [el for el in line.split(" ") if el != ""]
            block_arr =np.append(block_arr,np.array([fl]))
            block_arr = block_arr.reshape(1,len(block_arr))

            ## Looping throug the next 7-lines for current message (satellitte)
            for i in np.arange(0,7):
                line = filnr.readline().rstrip()
                ## -Replacing 'D' with 'E'
                line = line.replace('D','E')

                ## -- Have to add space between datacolums where theres no whitespace
                for idx, val in enumerate(line):
                    if line[idx] == 'E':
                        line = line[:idx+4] + " " + line[idx+4:]

                ## --Reads the line vector nl from the text string line and adds navigation
                # message for the relevant satellite n_sat. It becomes a long line vector
                # for the relevant message and satellite.
                nl = [el for el in line.split(" ") if el != ""]
                block_arr = np.append(block_arr,np.array([nl]))
                block_arr = block_arr.reshape(1,len(block_arr))

            ## -- Collecting all data into common variable
            if block_arr.shape[1] > 36:
                block_arr = block_arr[:,0:36]
            if np.size(block_arr) != 0:
                data  = np.concatenate([data , block_arr], axis=0)
            else:
                data  = np.delete(data , (0), axis=0)
                print('File %s is read successfully!' % (filename))


        filnr.close()
        n_eph = len(data)
        data = data.astype(float)
        if dataframe == 'yes' or dataframe == 'YES':
            data = DataFrame(data)

        # Create a dictinary for the data
        nav_data = {'ephemerides':data,
                    'header':header,
                    'nepohs': n_eph
                    }

        return nav_data



class Rinex_v3_Reader(RinexNav):
    """
    Class for reading RINEX v3 navigation files.

    Example on how to use it:
    from gnssmultipath import Rinex_v3_Reader
    navdata = Rinex_v3_Reader().read_rinex_nav(nav_file)
    """

    def __init__(self):
        super().__init__()
        self.valid_systems = {'G', 'R', 'E', 'C'}


    def read_rinex_nav(self, filename, desired_GNSS: list = ['G','R','E','C'], dataframe = False, data_rate = 30):
        """
        Reads the navigation message from broadcast efemerids in RINEX v.3 format.
        Support all global systems: GPS, GLONASS, Galileo and BeiDou

        Reads one navigation message at a time until the end of the row. Accumulate in
        in a common matrix, "data", where there is a line for each message.
        Note that each message forms a block in the file and that the same satellite
        can have several messages, usually with an hour's difference in reference time.


        Parameters
        ----------
        filename : Filename of the RINEX navigation file
        desired_GNSS: List of desired systems. Ex desired_GNSS = ['G','R','E']
        dataframe : Set to True to get the data output as a pandas DataFrame (array as default)
        data_rate: The desired data rate of ephemerides given in minutes. Default is 30 min.

        Returns
        -------
        data : Matrix with data for all epochs (or dataframe)
        header: List with header content
        n_eph: Number of epochs
        glo_fcn: Dictionary containing Glonass FCN is GLONASS included in rinex nav


        """

        for sys in desired_GNSS:
            if not all(letter in self.valid_systems for letter in sys):
                raise ValueError("Invalid GNSS system in desired_GNSS list. Must be one these ['G','R','E','C'].")

        header = self.read_header_lines(filename)
        nav_lines = self.filter_data_rinex_nav(filename, desired_GNSS, data_rate=data_rate)
        current_epoch = 0
        n_update_break = max(1, len(nav_lines) // 10)
        n_ep = 100 if len(nav_lines)>10 else len(nav_lines)
        bar_format = '{desc}: {percentage:3.0f}%|{bar}| ({n_fmt}/{total_fmt})'
        with tqdm(total= n_ep, desc ="Rinex navigation file is being read" , position=0, leave=True, bar_format=bar_format) as pbar:
            data  = np.zeros((1,36))
            for lines in nav_lines:
                line = lines[0].rstrip()
                block_arr = np.array([])
                sys_PRN = line[0:3]
                line = line.replace('D','E') # Replacing 'D' with 'E' ('D' is fortran syntax for exponentiall form)
                ## -- Have to add space between datacolums where theres no whitespace
                for idx, val in enumerate(line):
                    if line[0:2] != ' ' and line[23] != ' ':
                        line = line[:23] + " " + line[23:]
                    if line[idx] == 'e' or line[idx] == 'E' and idx !=0:
                        line = line[:idx+4] + " " + line[idx+4:]

                fl = [el for el in line.split(" ") if el != ""]
                block_arr = np.append(block_arr,np.array([fl]))
                block_arr = block_arr.reshape(1,len(block_arr))

                ## Looping throug the next 3-lines for current message (satellitte) (GLONASS)
                if 'R' in sys_PRN:
                    for i in np.arange(1,len(lines)):
                        line = lines[i].rstrip()
                        line = line.replace('D','E') # Replacing 'D' with 'E' ('D' is fortran syntax for exponentiall form)
                        ## -- Have to add space between datacolums where theres no whitespace
                        for idx, val in enumerate(line):
                            if line[idx] == 'e' or line[idx] == 'E':
                                line = line[:idx+4] + " " + line[idx+4:]

                        ## --Reads the line vector nl from the text string line and adds navigation
                        # message for the relevant satellite n_sat. It becomes a long line vector
                        # for the relevant message and satellite.
                        nl = [el for el in line.split(" ") if el != ""]
                        block_arr = np.append(block_arr,np.array([nl]))
                        block_arr = block_arr.reshape(1,len(block_arr))
                else:
                    ## Looping throug the next 7-lines for current message (satellitte) (GPS,Galileo,BeiDou)
                    for i in np.arange(1,len(lines)):
                        line = lines[i].rstrip()
                        line = line.replace('D','E') # Replacing 'D' with 'E' ('D' is fortran syntax for exponentiall form)
                        if line == '':
                            continue
                        else:
                            ## -- Have to add space between datacolums where theres no whitespace
                            for idx, val in enumerate(line):
                                if line[idx] == 'e' or line[idx] == 'E':
                                    line = line[:idx+4] + " " + line[idx+4:]

                            ## Reads the line vector nl from the text string line and adds navigation
                            # message for the relevant satellite n_sat. It becomes a long line vector
                            # for the relevant message and satellite.
                            ## Runs through line to see if each line contains 4 objects. If not, adds nan.
                            if i < 7 and line.lower().count('e') < 4:
                                if line[10:20].strip() == '':
                                    line = line[:10] +  'nan' + line[10:]
                                if line[30:40].strip() == '':
                                    line = line[:30] +  'nan' + line[30:]
                                if line[50:60].strip() == '':
                                    line = line[:50] +  'nan' + line[50:]
                                if line[70:80].strip() == '':
                                    line = line[:70] +  'nan' + line[70:]


                            if i == 7 and line.lower().count('e') < 2 and 'E' not in sys_PRN:
                                if line[10:20].strip() == '':
                                    line = line[:10] +  'nan' + line[10:]
                                if line[30:40].strip() == '':
                                    line = line[:30] +  'nan' + line[30:]

                            if i == 7 and line.lower().count('e') < 1 and 'E' in sys_PRN: #only one object in last line for Galileo
                                if line[10:20].strip() == '':
                                    line = line[:10] +  'nan' + line[10:]

                            if i == 7 and 'E' in sys_PRN: #only one object in last line for Galileo, but add nan to get 36 in total
                                line = line + 'nan'


                            nl = [el for el in line.split(" ") if el != ""]
                            block_arr = np.append(block_arr,np.array([nl]))
                            block_arr = block_arr.reshape(1,len(block_arr))

                ## -- Collecting all data into common variable
                if block_arr.shape[1] > 36:
                    block_arr = block_arr[:,0:36]
                try:
                    if np.size(block_arr) != 0 and 'R' not in sys_PRN:
                        data  = np.concatenate([data , block_arr], axis=0)

                    elif np.size(block_arr) != 0 and 'R' in sys_PRN:
                        GLO_dum = np.zeros([1,np.size(data,axis=1) - np.size(block_arr,axis=1)])
                        block_arr = np.append(block_arr,GLO_dum) # adding emtpy columns to match size of other systems
                        block_arr = block_arr.reshape(1,len(block_arr))
                        data  = np.concatenate([data , block_arr], axis=0)
                    else:
                        data  = np.delete(data , (0), axis=0)
                except:
                    print("ERROR! Sure this is a RINEX v.3 navigation file?")
                    break

                current_epoch += 1
                if len(nav_lines) >=10 and np.mod(current_epoch, n_update_break) == 0:  # Update progress bar every n_update_break epochs
                    pbar.update(10)
                elif len(nav_lines) < 10 and np.mod(current_epoch, n_update_break) == 0:
                    pbar.update(1)

        # Remove first row if contains only zeros
        if np.all(data[0,:] == '0.0'):
            data  = np.delete(data , (0), axis=0)

        if np.any(np.char.startswith(data[:, 0], 'R')):
            self.glo_fcn = self.extract_glonass_fcn_from_rinex_nav(data)

        n_eph = len(data)
        if dataframe:
            data = DataFrame(data)
            data_columns = list(data.columns)
            data_columns.pop(0) # Removing index that contains PRN nr ex 'G01'
            data[data_columns] = data[data_columns].astype(float) # Change the other values to float


        # Create a dictinary for the data
        nav_data = {'ephemerides':data,
                    'header':header,
                    'nepohs': n_eph,
                    'glonass_fcn': self.glo_fcn}

        return nav_data





if __name__=="__main__":
    # brod1 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\BRDC00IGS_R_20220010000_01D_MN.rnx"
    # data = Rinex_v3_Reader().read_rinex_nav(brod1, dataframe=True, data_rate=60)
    pass
