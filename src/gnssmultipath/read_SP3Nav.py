"""
Module for reading SP3 files.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""
import copy
import numpy as np


def readSP3Nav(filename, desiredGNSSsystems=None):
    """
    Function that reads the GNSS satellite position data from a SP3 position file.

    INPUTS:
    -------

    filename:             path and filename of sp3 position file, string

    desiredGNSSsystems:   List of strings. Each string is a code for a
                          GNSS system that should have its position data stored
                          in sat_pos. Must be one of: "G", "R", "E",
                          "C". If left undefined, it is automatically set to
                          ["G", "R", "E", "C"]

    OUTPUTS:
    -------

    sat_pos:          dict. Each elements contains position data for a
                      specific GNSS system. Order is defined by order of
                      navGNSSsystems. Each key element is another dict that
                      stores position data of a specific epoch of that
                      GNSS system. Each of these dict  is an array
                      with [X, Y, Z] position for each satellite. Each satellite
                      has their PRN number as a key.

                      sat_pos[GNSSsystem][epoch][PNR] = [X, Y, Z]
                      Ex:
                          sat_pos['G'][100][24] = [X, Y, Z]
                      This command will extract the coordinates for GPS at epoch
                      100 for satellite PRN 24.

    epoch_dates:      matrix. Each row contains date of one of the epochs.
                      [nEpochs x 6]

    navGNSSsystems:   list of strings. Each string is a code for a
                      GNSS system with position data stored in sat_pos.
                      Must be one of: "G", "R", "E", "C"

    nEpochs:          number of position epochs, integer

    epochInterval:    interval of position epochs, seconds

    success:          boolean, 1 if no error occurs, 0 otherwise

   """

    #Initialize variables
    epochInterval = None
    success = 1

    ## --- Open nav file
    try:
        fid = open(filename, 'r', encoding='utf-8')
    except Exception as exc:
        success = 0
        raise ValueError('No file selected!') from exc


    if desiredGNSSsystems is None:
        desiredGNSSsystems = ["G", "R", "E", "C"]

    navGNSSsystems = ["G", "R", "E", "C"] # GNSS system order
    GNSSsystem_map = dict(zip(navGNSSsystems,[1, 2, 3, 4])) # Mapping GNSS system code to GNSS system index
    sat_pos = {}
    # Read header
    headerLine = 0
    line = fid.readline().rstrip()
    # All header lines begin with '*'
    while '*' not in line[0]:
        headerLine = headerLine + 1
        if headerLine == 1:
            sp3Version = line[0:2]
            # Control sp3 version
            if '#c' not in sp3Version and '#d' not in sp3Version:
                print(f'ERROR(readSP3Nav): SP3 Navigation file is version {sp3Version}, must be version c or d!')
                success = 0
                return success
            # Control that sp3 file is a position file and not a velocity file
            Pos_Vel_Flag = line[2]

            if 'P' not in Pos_Vel_Flag:
                print('ERROR(readSP3Nav): SP3 Navigation file is has velocity flag, should have position flag!')
                success = 0
                return success

            #Store coordinate system and amount of epochs
            CoordSys = line[46:51]
            nEpochs = int(line[32:39])

        if headerLine == 2:
            # Store GPS-week, "time-of-week" and epoch interval[seconds]
            GPS_Week = int(line[3:7])
            tow      = float(line[8:23])
            epochInterval = float(line[24:38])


        if headerLine == 3:
            #initialize list for storing indices of satellites to be excluded
            RemovedSatIndex = []
            if '#c' in sp3Version:
                nSat = int(line[4:6])
            else:
                nSat = int(line[3:6])

            line = line[9:60] # Remove beginning of line
            # Initialize array for storing the order of satellites in the SP3 file (ie. what PRN and GNSS system index)
            GNSSsystemIndexOrder = []
            PRNOrder = []
            # Keep reading lines until all satellite IDs have been read
            for k in range(0,nSat):
                # Control that current satellite is among desired systems
                if np.in1d(line[0], desiredGNSSsystems):
                    ## -- Get GNSSsystemIndex from map container
                    GNSSsystemIndex = GNSSsystem_map[line[0]]
                    #Get PRN number/slot number
                    PRN = int(line[1:3])
                    #remove satellite that has been read from line
                    line = line[3::]
                    #Store GNSSsystemIndex and PRN in satellite order arrays
                    GNSSsystemIndexOrder.append(GNSSsystemIndex)
                    PRNOrder.append(PRN)
                    #if current satellite ID was last of a line, read next line
                    #and increment number of headerlines
                    if np.mod(k+1,17)==0 and k != 0:
                        line = fid.readline().rstrip()
                        line = line[9:60]
                        headerLine = headerLine + 1
                #If current satellite ID is not among desired GNSS systems,
                #append its index to the array of undesired satellites
                else:
                    RemovedSatIndex.append(k)
                    GNSSsystemIndexOrder.append(np.nan)
                    PRNOrder.append(np.nan)
                    #if current satellite ID was last of a line, read next line
                    #and increment number of headerlines
                    if np.mod(k+1,17)==0 and k != 0:
                        line = fid.readline().rstrip()
                        line = line[9:60]
                        headerLine = headerLine + 1
        # Read next line
        line = fid.readline().rstrip()

    # Initialize matrix for epoch dates
    epoch_dates = []
    sys_dict = {}
    PRN_dict_GPS = {}
    PRN_dict_Glonass = {}
    PRN_dict_Galileo = {}
    PRN_dict_BeiDou = {}

    # Read satellite positions of every epoch
    ini_sys = list(GNSSsystem_map.keys())[0]
    for k in range(0,nEpochs):
        #Store date of current epoch
        epochs = line[3:31].split(" ")
        epochs = [x for x in epochs if x != "" ] # removing ''
        ## -- Make a check if there's a new line. (if the header is not giving the correct nepochs)
        if epochs == []:
            print(f'The number of epochs given in the headers is not correct!\nInstead of {str(nEpochs)} epochs, the file contains {str(k+1)} epochs.\nSP3-file {filename} has been read successfully')
            return sat_pos, epoch_dates, navGNSSsystems, nEpochs, epochInterval,success
        epoch_dates.append(epochs)

        # Store positions of all satellites for the current epoch
        obs_dict_GPS = {}
        obs_dict_Glonass = {}
        obs_dict_Galileo = {}
        obs_dict_BeiDou = {}
        for i in range(0,nSat):
            line = fid.readline().rstrip()
            #if current satellite is among desired systems, store positions
            if np.in1d(i, RemovedSatIndex,invert = True):
                #Get PRN and GNSSsystemIndex of current satellite for previously stored order
                PRN = PRNOrder[i]
                GNSSsystemIndex = GNSSsystemIndexOrder[i]
                # Store position of current satellite in the correct location in
                sys_keys = list(GNSSsystem_map.keys())
                sys_values = list(GNSSsystem_map.values())
                sys_inx = sys_values.index(GNSSsystemIndex)
                sys = sys_keys[sys_inx]
                obs = line[5:46].split(" ")
                obs = [float(x)*1000 for x in obs if x != "" ] # multiplying with 1000 to get meters
                if sys != ini_sys:
                    ini_sys = sys
                if sys == 'G':
                    obs_G = [x for x in obs if x != "" ]
                    obs_dict_GPS[PRN]  = np.array([obs_G])
                    PRN_dict_GPS[k] = obs_dict_GPS
                elif sys =='R':
                    obs_R = [x for x in obs if x != "" ]
                    obs_dict_Glonass[PRN]  = np.array([obs_R])
                    PRN_dict_Glonass[k] = obs_dict_Glonass
                elif sys =='E':
                    obs_E = [x for x in obs if x != "" ]
                    obs_dict_Galileo[PRN]  = np.array([obs_E])
                    PRN_dict_Galileo[k] = obs_dict_Galileo
                elif sys =='C':
                    obs_C = [x for x in obs if x != "" ]
                    obs_dict_BeiDou[PRN]  = np.array([obs_C])
                    PRN_dict_BeiDou[k] = obs_dict_BeiDou

            sys_dict['G'] = PRN_dict_GPS
            sys_dict['R'] = PRN_dict_Glonass
            sys_dict['E'] = PRN_dict_Galileo
            sys_dict['C'] = PRN_dict_BeiDou
            sat_pos[sys] = sys_dict[sys]

        #Get the next line
        line = fid.readline().rstrip()

    # The next line should be eof. If not, raise a warning
    try:
        line = fid.readline().rstrip()
    except EOFError:
        print('ERROR(readSP3Nav): End of file was not reached when expected!')
        success = 0
        return success

    # Remove NaN values
    GNSSsystemIndexOrder = [x for x in GNSSsystemIndexOrder if x != 'nan']
    PRNOrder = [x for x in GNSSsystemIndexOrder if x != 'nan']
    epoch_dates = np.array(epoch_dates)
    print(f'SP3 Navigation file "{filename}" has been read successfully.')

    return sat_pos, epoch_dates, navGNSSsystems, nEpochs, epochInterval, success












def combineSP3Nav(three_sp3_files, sat_positions_1, epoch_dates_1, navGNSSsystems_1, nEpochs_1, epochInterval_1,
                  sat_positions_2, epoch_dates_2, navGNSSsystems_2, nEpochs_2, epochInterval_2,
                  sat_positions_3, epoch_dates_3, navGNSSsystems_3, nEpochs_3, epochInterval_3, GNSSsystems):

    """
    Function that combines the precise orbital data of two or three SP3
    files. Note that the SP3 files should first be read by the function
    readSP3Nav.

    INPUTS:
    -------

    three_sp3_files:      boolean. 1 if there are three SP3 files to be
                          combined, 0 otherwise

    sat_positions_1:      dict. Each element contains position data for a
                          specific GNSS system. Order is defined by the order of
                          navGNSSsystems. Each key element is another dict that
                          stores position data of a specific epoch of that
                          GNSS system. Each of these dict is an array
                          with [X, Y, Z] position for each satellite. Each satellite
                          has its PRN number as a key.

                          sat_positions_1[GNSSsystem][epoch][PNR] = [X, Y, Z]

    epoch_dates_1:        matrix. Each row contains the date of one of the epochs
                          from the first SP3 file
                          [nEpochs_1 x 6]

    navGNSSsystems_1:     array. Contains strings. Each string is a code for a
                          GNSS system with position data stored in sat_positions.
                          Must be one of: "G", "R", "E", "C"

    nEpochs_1:            number of position epochs in the first SP3 file, integer

    epochInterval_1:      interval of position epochs in the first SP3 file, seconds

    sat_positions_2:      dict. Each element contains position data for a
                          specific GNSS system. Order is defined by the order of
                          navGNSSsystems. Each key element is another dict that
                          stores position data of a specific epoch of that
                          GNSS system. Each of these dict is an array
                          with [X, Y, Z] position for each satellite. Each satellite
                          has its PRN number as a key.

                          sat_positions_2[GNSSsystem][epoch][PNR] = [X, Y, Z]

    epoch_dates_2:        matrix. Each row contains the date of one of the epochs
                          from the second SP3 file
                          [nEpochs_1 x 6]

    navGNSSsystems_2:     array. Contains strings. Each string is a code for a
                          GNSS system with position data stored in sat_positions.
                          Must be one of: "G", "R", "E", "C"

    nEpochs_2:            number of position epochs in the first SP3 file, integer

    epochInterval_2:      interval of position epochs in the second SP3 file, seconds

    sat_positions_3:      dict. Each element contains position data for a
                          specific GNSS system. Order is defined by the order of
                          navGNSSsystems. Each key element is another dict that
                          stores position data of a specific epoch of that
                          GNSS system. Each of these dict is an array
                          with [X, Y, Z] position for each satellite. Each satellite
                          has its PRN number as a key.

                          sat_positions_3[GNSSsystem][epoch][PNR] = [X, Y, Z]

    epoch_dates_3:        matrix. Each row contains the date of one of the epochs
                          from the third SP3 file
                          [nEpochs_1 x 6]

    navGNSSsystems_3:     array. Contains strings. Each string is a code for a
                          GNSS system with position data stored in sat_positions.
                          Must be one of: "G", "R", "E", "C"

    nEpochs_3:            number of position epochs in the third SP3 file, integer

    epochInterval_3:      interval of position epochs in the third SP3 file, seconds



    OUTPUTS:
    --------

    sat_positions:        dict. Contains data from all two/three SP3 files. Each cell
                          element contains position data for a specific GNSS
                          system. Order is defined by the order of navGNSSsystems_1.
                          Each cell element is another cell that stores
                          position data of specific satellites of that
                          GNSS system. Each of these cell elements is a matrix
                          with [X, Y, Z] position of an epoch in each row.

    epoch_dates:          matrix. Each row contains the date of one of the epochs
                          from all two/three SP3 files
                          [nEpochs_1 x 6]

    navGNSSsystems:       array. Contains strings. Each string is a code for a
                          GNSS system with position data stored in sat_positions.
                          Must be one of: "G", "R", "E", "C"

    nEpochs:              number of position epochs in all two/three SP3 files, integer

    epochInterval:        interval of position epochs in all SP3 files, seconds

    success:              boolean, 1 if no error occurs, 0 otherwise
    """

    success = 1  # Setting success to 1
    ## -- Check that first two SP3 files have the same interval
    if epochInterval_1 != epochInterval_2:
        print('ERROR(combineSP3Nav): The first and second SP3 files do not have the same epoch interval')
        sat_positions, epoch_dates, navGNSSsystems, nEpochs, epochInterval = np.nan, np.nan, np.nan, np.nan, np.nan
        success = 0
        return success

    ## -- Check that first two SP3 files have the same GNSS systems
    if list(set(navGNSSsystems_1) - set(navGNSSsystems_2)):
        print('ERROR(CombineSP3Nav): SP3 file 1 and 2 do not contain the same GNSS systems')
        sat_positions, epoch_dates, navGNSSsystems, nEpochs, epochInterval = np.nan, np.nan, np.nan, np.nan, np.nan
        success = 0
        return success

    if three_sp3_files:
        ## -- Check that the last SP3 file has the same GNSS systems as the others
        if list(set(navGNSSsystems_2) - set(navGNSSsystems_3)):
            print('ERROR(CombineSP3Nav): SP3 file 3 does not contain the same GNSS systems as SP3 file 1 and 2')
            sat_positions, epoch_dates, navGNSSsystems, nEpochs, epochInterval = np.nan, np.nan, np.nan, np.nan, np.nan
            success = 0
            return success

    navGNSSsystems = navGNSSsystems_1
    epochInterval = epochInterval_1

    epoch_dates = np.vstack([epoch_dates_1, epoch_dates_2]) # Combine epoch dates from the first and second SP3 file
    nEpochs = nEpochs_1 + nEpochs_2  # Compute the total number of epochs
    sat_positions = copy.deepcopy(sat_positions_1)  # Initialize dictionary for storing combined satellite positions

    # Combine satellite positions from the first and second SP3 files
    for k in range(0, len(GNSSsystems)):
        curr_sys = GNSSsystems[k + 1]
        len_sat = len(sat_positions_1[curr_sys])
        for ep in range(0, len(sat_positions_2[curr_sys].keys())):
            sat_positions[curr_sys].update({ep + len_sat: sat_positions_2[curr_sys][ep]})

    # If three SP3 files are present
    if three_sp3_files:
        epoch_dates = np.vstack([epoch_dates, epoch_dates_3]) # Add epoch dates from the third SP3 file to the first two
        nEpochs = nEpochs + nEpochs_3  # Compute the total number of epochs
        # Check that the last SP3 file has the same interval as the others
        if epochInterval_2 != epochInterval_3:
            print('ERROR(combineSP3Nav): The second and third SP3 files do not have the same epoch interval')
            sat_positions, epoch_dates, navGNSSsystems, nEpochs, epochInterval = np.nan, np.nan, np.nan, np.nan, np.nan
            success = 0
            return success

        # Check that the last SP3 file has the same GNSS systems as the others
        if list(set(navGNSSsystems_2) - set(navGNSSsystems_3)):
            print('Warning (CombineSP3Nav): SP3 file 2 and 3 do not contain the same amount of GNSS systems')

        # Combine satellite positions from the first, second, and third SP3 files
        sat_positions_dum = copy.deepcopy(sat_positions)
        for k in range(0, len(GNSSsystems)):
            curr_sys = GNSSsystems[k + 1]
            len_sat = len(sat_positions_dum[curr_sys])
            for ep in range(0, len(sat_positions_3[curr_sys].keys())):
                sat_positions[curr_sys].update({ep + len_sat: sat_positions_3[curr_sys][ep]})

    return sat_positions, epoch_dates, navGNSSsystems, nEpochs, epochInterval, success
