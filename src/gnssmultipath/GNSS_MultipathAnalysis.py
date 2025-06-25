"""
This is the main module for running the software GNSS_Multipath_Analysis.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""

import os
import warnings
import logging
import time
from typing import Union, List
import numpy as np
from tqdm import tqdm
from gnssmultipath.readRinexObs import readRinexObs
from gnssmultipath.Geodetic_functions import gpstime_to_utc_datefmt, gpstime2date
from gnssmultipath.computeSatElevAzimuth_fromNav import computeSatElevAzimuth_fromNav
from gnssmultipath.signalAnalysis import signalAnalysis
from gnssmultipath.detectClockJumps import detectClockJumps
from gnssmultipath.writeOutputFile import writeOutputFile
from gnssmultipath.createCSVfile import createCSVfile
from gnssmultipath.make_polarplot import make_polarplot,make_skyplot, make_polarplot_SNR, plot_SNR_wrt_elev
from gnssmultipath.make_polarplot_dont_use_TEX import make_polarplot_dont_use_TEX, make_skyplot_dont_use_TEX, make_polarplot_SNR_dont_use_TEX, plot_SNR_wrt_elev_dont_use_TEX
from gnssmultipath.plotResults import plotResults, plotResults_dont_use_TEX, make_barplot, make_barplot_dont_use_TEX
from gnssmultipath.PickleHandler import PickleHandler
from gnssmultipath.PreciseSatCoords import PreciseSatCoords
from gnssmultipath.SP3PositionEstimator import SP3PositionEstimator

warnings.filterwarnings("ignore")



def GNSS_MultipathAnalysis(rinObsFilename: str,
                          broadcastNav1: Union[str, None] = None,
                          broadcastNav2: Union[str, None] = None,
                          broadcastNav3: Union[str, None] = None,
                          broadcastNav4: Union[str, None] = None,
                          sp3NavFilename_1: Union[str, None] = None,
                          sp3NavFilename_2: Union[str, None] = None,
                          sp3NavFilename_3: Union[str, None] = None,
                          desiredGNSSsystems: Union[List[str], None] = None,
                          phaseCodeLimit: Union[float, int, None] = None,
                          ionLimit: Union[float, None] = None,
                          cutoff_elevation_angle: Union[int, None] = None,
                          outputDir: Union[str, None] = None,
                          plotEstimates: bool = True,
                          plot_polarplot: bool = True,
                          include_SNR: bool = True,
                          save_results_as_pickle: bool = True,
                          save_results_as_compressed_pickle: bool = False,
                          write_results_to_csv: bool = True,
                          output_csv_delimiter: str = ';',
                          nav_data_rate: int = 60,
                          includeResultSummary: Union[bool, None] = None,
                          includeCompactSummary: Union[bool, None] = None,
                          includeObservationOverview: Union[bool, None] = None,
                          includeLLIOverview: Union[bool, None] = None,
                          use_LaTex: bool = True
                          ) -> dict:

    """
    GNSS Multipath Analysis
    ------------------------
    Made by: Per Helge Aarnes
    E-mail: per.helge.aarnes@gmail.com

    GNSS_MultipathAnalysis is a software for analyzing the multipath effect on Global Navigation Satellite Systems (GNSS) and
    is based on the MATLAB software "GNSS_Receiver_QC_2020" made by Bjørn-Eirik Roald.
    This is the main function of the software that, through the help of other functions:

      - reads RINEX observation files
      - reads SP3 satellite navigation files (if a SP3 file is fed in)
      - reads rinex navigation files (if a rinex navigation file is fed in)
      - computes elevation angles of satellites for every epoch
      - makes estimates of multipath, ionospheric delay, and cycle slips
          for all signals in RINEX observation file
      - plots estimates if user choose to do it
      - computes and stores statistics on estimates
      - writes an output files containing results

    This function calls on a range of functions. These in turn call on
    further more functions. The function called upon directly in
    GNSS_MultipathAnalysis are:

      - readRinexObs.py
      - computeSatElevations.py
      - computeSatElevAzimuth_fromNav.py
      - signalAnalysis.py
      - plotEstimates.py
      - make_barplot.py
      - make_polarplot.py
      - writeOutputFile.py

    --------------------------------------------------------------------------------------------------------------------------

    INPUTS:
    ------

    rinObsFilename:           string. Path to RINEX 3 observation file

    sp3NavFilename_1:         string. Path to first SP3 navigation file

    sp3NavFilename_2:         string. Path to second SP3 navigation file (optional).

    sp3NavFilename_3:         string. Path to third SP3 navigation file (optional).

    desiredGNSSsystems:       List with the desired GNSS system. Ex ['G','R'] if you want to
                              only run the analysis on GPS and GLONASS. Default: All systems. (if set to None) (optional)

    phaseCodeLimit:           critical limit that indicates cycle slip for
                              phase-code combination. Unit: m/s. If set to 0,
                              default value of 6.667 m/s will be used (optional)

    ionLimit:                 critical limit that indicates cycle slip for
                              the rate of change of the ionopheric delay.
                              Unit: m/s. If set to 0, default value of 0.0667 m/s will be used (optional)

    cutoff_elevation_angle    Critical cutoff angle for satellite elevation angles, degrees
                              Estimates where satellite elevation angle
                              is lower than cutoff are removed, so are
                              estimated slip periods (optional)

    outputDir:                string. Path to directory where output file
                              should be generated. If user does not wish to
                              specify output directory, this variable should
                              be empty string, "". In this case the output file
                              will be generated in sub-directory inside same
                              directory as GNSS_Receiver_QC_2020.m (optional)

    plotEstimates:            boolean. False if user desires estimates not to be
                              ploted. True by default. (optional)

    plot_polarplot:           boolean. True or False. If not defined polarplots will be made (optional)


    include_SNR:              boolean. If not defined, SNR from Rinex obs file will NOT be used (optional)

    save_results_as_pickle:   boolean. If True, the results will be stored as dictionary in form of a binary pickle file. Default set to True.


    save_results_as_compressed_pickle : boolean. If True, the results will be stored as dictionary in form of a binary compressed pickle file (zstd compression). Default set to False.

    write_results_to_csv: boolean. If True, a subset of the results will be exported as a CSV file. Default is True.

    output_csv_delimiter:     str. Set the delimiter of the CSV file. Default is semi colon (;).


    nav_data_rate:            integer. The desired data rate of ephemerides given in minutes. Default is 60 min. The purpose with this
                              parameter is to speed up processing time. Both reading the RINEX navigation file and looping through the
                              ephemerides aftwerward will be significatnly faster by increasing this value. Note: A too high value will
                              degrade the accuracy of the interploated satellite coordinates.

    includeResultSummary:     boolean. 1 if user desires output file to
                              include more detailed overview of statistics,
                              including for each individual satellites.
                              0 otherwise (optional)

    includeCompactSummary:    boolean. 1 if user desired output file to
                              include more compact overview og statistics. (optional)

    includeObservationOverview: boolean. 1 if user desires output file to
                                include overview of obseration types observed
                                by each satellite. 0 otherwise (optional)

    use_LaTex:                 boolean. Will use TeX as an interpreter in plots. Default set to true. "Requires TeX installed on computer".

    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS:

    analysisResults:          A dictionary that contains alls results of all analysises, for all GNSS systems.


    The software is also returning results file. A report provided as a text file, and a CSV file with the estimated values.
    --------------------------------------------------------------------------------------------------------------------------
    """
    start_time = time.time()


    if broadcastNav1 is None and sp3NavFilename_1 is None:
        raise RuntimeError("No SP3 or navigation file is defined! This is \
                           mandatory for this software, so please add one of them.")

    if broadcastNav1 is not None and sp3NavFilename_1 is not None:
        raise RuntimeError("You defined both a navigation file and a SP3 file. Please\
                           choose between using broadcast ephemerides or precise.")

    if broadcastNav1 is None:
        broadcastNav1 = ""
    if broadcastNav2 is None:
        broadcastNav2 = ""
    if broadcastNav3 is None:
        broadcastNav3 = ""
    if broadcastNav4 is None:
        broadcastNav4 = ""
    if sp3NavFilename_1 is None:
        sp3NavFilename_1 = ""
    if sp3NavFilename_2 is None:
        sp3NavFilename_2 = ""
    if sp3NavFilename_3 is None:
        sp3NavFilename_3 = ""
    if phaseCodeLimit is None:
        phaseCodeLimit =  4/60 *100

    if ionLimit is None:
        ionLimit = 4/60

    if cutoff_elevation_angle is None:
        cutoff_elevation_angle = 0

    #  Create output file
    if outputDir is None:
        outputDir = 'Output_Files'

    if not os.path.isdir(outputDir):
        os.mkdir(outputDir)

    # Unless graph directory already exists, create directory
    graphDir = outputDir + '/Graphs'
    if not os.path.isdir(graphDir):
        os.mkdir(graphDir)


    if includeResultSummary is None:
        includeResultSummary = 1

    if includeCompactSummary is None:
        includeCompactSummary = 1

    if includeObservationOverview is None:
        includeObservationOverview = 1

    if includeLLIOverview is None:
        includeLLIOverview = 1

    if desiredGNSSsystems is None:
        includeAllGNSSsystems   = 1
        desiredGNSSsystems = ["G", "R", "E", "C"]  # All GNSS systems.
    else:
        includeAllGNSSsystems   = 0

    estimated_position, stats = None, None


    ## ---  Control of the user input arguments
    if not isinstance(sp3NavFilename_1, str):
        print('ERROR(GNSS_MultipathAnalysis): The input variable sp3NavFilename_1 must be a string\n' \
            'Argument is now of type %s\n' %  (type(sp3NavFilename_1)))
        analysisResults = np.nan
        return

    if not isinstance(sp3NavFilename_2, str):
        print('ERROR(GNSS_MultipathAnalysis): The input variable sp3NavFilename_2 must be a string\n' \
            'Argument is now of type %s\n' %  (type(sp3NavFilename_2)))
        analysisResults = np.nan
        return


    if not isinstance(sp3NavFilename_3, str):
        print('ERROR(GNSS_MultipathAnalysis): The input variable sp3NavFilename_3must be a string\n' \
            'Argument is now of type %s\n' %  (type(sp3NavFilename_3)))
        analysisResults = np.nan
        return

    if not isinstance(rinObsFilename, str):
        print('ERROR(GNSS_MultipathAnalysis): The input variable rinObsFilename must be a string\n' \
            'Argument is now of type %s\n' %  (type(rinObsFilename)))
        analysisResults = np.nan
        return


    if not os.path.isfile(rinObsFilename):
        print('ERROR(GNSS_MultipathAnalysis): RINEX observation file can not be found. Please check that the path is correct.\n')
        analysisResults = np.nan
        return


    if not os.path.isfile(sp3NavFilename_1) and len(sp3NavFilename_1) != 0:
        print('WARNING: Second SP3 Navigation file can not be found.\n')

    if not os.path.isfile(sp3NavFilename_2) and len(sp3NavFilename_2) != 0:
        print('WARNING: Second SP3 Navigation file can not be found.\n')


    if not os.path.isfile(sp3NavFilename_3) and len(sp3NavFilename_3) != 0:
        print('WARNING: Third SP3 Navigation file can not be found.\n')


    # Check for conflicting save options
    if save_results_as_pickle and save_results_as_compressed_pickle:
        save_results_as_pickle = False


    latex_installed = True
    glo_fcn = None
    # A frequency overview for the different systems.
    frequencyOverview_temp2 = {
        'G': np.array([[1.57542e+09],[1.22760e+09],[np.nan],[np.nan],[1.17645e+09],[np.nan],[np.nan],[np.nan],[np.nan]]),
        'R': np.array([[1.602000e+09, 5.625000e+05],[1.246000e+09, 4.375000e+05],[1.202025e+09, 0.000000e+00],[1.600995e+09, 0.000000e+00],
                        [np.nan, 0.000000e+00],[1.248060e+09, 0.000000e+00],[np.nan, 0.000000e+00],[np.nan, 0.000000e+00],[np.nan, 0.000000e+00]]),
        'E': np.array([[1.575420e+09],[np.nan],[np.nan],[np.nan],[1.176450e+09],[1.278750e+09],[1.207140e+09],[1.191795e+09],[np.nan]]),
        'C': np.array([[1.575420e+09],[1.561098e+09],[np.nan],[np.nan],[1.176450e+09],[1.268520e+09],[1.207140e+09],[1.191795e+09],[np.nan]])
        }

    ## -- Create a logger instance (logging.INFO, which will include INFO, WARNING, ERROR, and CRITICAL (not DEBUG))
    path_logfile = os.path.join(outputDir,'Logfile.log')
    logging.basicConfig(filename=path_logfile, filemode='w', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)


    ## Make dictionaries
    GNSSsystemCode2Fullname = dict(zip(['G','R','E','C'],['GPS','GLONASS','Galileo','BeiDou']))
    GNSSsystem2BandsMap = dict(zip(['GPS', 'GLONASS', 'Galileo', 'BeiDou'],\
        [[1, 2, 5], [1, 2, 3, 4, 6], [1, 5, 6, 7, 8], [1, 2, 5, 6, 7, 8]]))


    ## --- Read observation file
    includeAllObsCodes  = 0

    if include_SNR:
        desiredObsCodes = ["C", "L", "S"]
    else:
        desiredObsCodes = ["C", "L"] # only code and phase observations


    desiredObsBands = list(np.arange(1,10)) # all carrier bands. Tot 9, but arange stops at 8 -> 10
    readSS = 1
    readLLI = 1

    ## --- Read RINEX 3.0x observation file
    [GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
        obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
        rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success] = \
        readRinexObs(rinObsFilename, readSS=readSS, readLLI=readLLI, includeAllGNSSsystems=includeAllGNSSsystems,includeAllObsCodes=includeAllObsCodes, desiredGNSSsystems=desiredGNSSsystems,\
        desiredObsCodes=desiredObsCodes, desiredObsBands=desiredObsBands)


    sat_pos = {}
    if sp3NavFilename_1 != '':
        ## -- Compute satellite elevation angles from SP3 files
        # sat_elevation_angles, sat_azimut_angles, sat_coordinates = computeSatElevations(GNSS_SVs, GNSSsystems, approxPosition,\
        #     nepochs, time_epochs, max_sat, sp3NavFilename_1, sp3NavFilename_2, sp3NavFilename_3)
        # sat_elevation_angles1, sat_azimut_angles1, sat_coordinates1 = computeSatElevations(GNSS_SVs, GNSSsystems, approxPosition,\
        #     nepochs, time_epochs, max_sat, sp3NavFilename_1, sp3NavFilename_2, sp3NavFilename_3)
        x_rec_approx, y_rec_approx, z_rec_approx = np.atleast_2d(approxPosition).flatten()

        # Filter out empty strings using a list comprehension
        sp3_files = [sp3NavFilename_1, sp3NavFilename_2, sp3NavFilename_3]
        sp3_files = [file for file in sp3_files if file]

        sat_obj = PreciseSatCoords(sp3_files, time_epochs=time_epochs, GNSSsystems= GNSSsystems)
        df_sat_coordinates = sat_obj.satcoords

        if all(coord == 0 for coord in [x_rec_approx, y_rec_approx, z_rec_approx]):
                desired_time = np.array(gpstime2date(time_epochs[0,0], round(time_epochs[0,1],6)))
                position_estimator = SP3PositionEstimator(df_sat_coordinates, desired_time= desired_time,
                                                          GNSS_obs = GNSS_obs, time_epochs = time_epochs,
                                                          GNSSsystems = GNSSsystems, obsCodes = obsCodes,
                                                          sp3_metadata_dict = sat_obj.sp3_metadata_dict)
                estimated_position, stats = position_estimator.estimate_position()
                x_rec_approx, y_rec_approx, z_rec_approx,_ = estimated_position.flatten()

        df_az_el = sat_obj.compute_azimuth_and_elevation(receiver_position=(x_rec_approx, y_rec_approx, z_rec_approx), drop_below_horizon=True)
        sat_dict = sat_obj.create_satellite_data_dict(df_sat_coordinates, df_az_el)

        # Create dicts for each data type
        sat_coordinates = {}
        sat_elevation_angles = {}
        sat_azimut_angles = {}

        # Loop through sat_dict, enumerate to get an index for each key
        for idx, (system, data) in enumerate(sat_dict.items()):
            sat_coordinates[system] = data.get('coordinates', {})
            sat_elevation_angles[idx] = data.get('elevation', None)
            sat_azimut_angles[idx] = data.get('azimuth', None)

    else:
        nav_files = [broadcastNav1,broadcastNav2,broadcastNav3,broadcastNav4]
        sat_pos, glo_fcn, estimated_position, stats = computeSatElevAzimuth_fromNav(nav_files, approxPosition, GNSS_SVs, time_epochs, nav_data_rate, GNSS_obs, GNSSsystems, obsCodes)

        ## -- Build same struture for satellit elevation angles if broadcast nav defined
        sat_elevation_angles = {}
        sat_pos_dummy = sat_pos.copy()

        # Intersect the keys in the sat_pos_dummy dict and the GNSSsystems dict
        GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, GNSSsystems = filter_common_gnss_keys(sat_pos_dummy, GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, GNSSsystems)

        for sys in np.arange(0,len(GNSSsystems)):
            currentGNSSsystem = GNSSsystems[sys+1]

            # Raise an error if currentGNSSsystem is not in sat_pos_dummy keys
            if currentGNSSsystem not in sat_pos_dummy.keys():
                sys_name = GNSSsystemCode2Fullname[currentGNSSsystem]
                sys_in_rin_nav = [GNSSsystemCode2Fullname[sys] for sys in sat_pos_dummy.keys()]
                raise KeyError(f'GNSS system "{sys_name}" is not present in the RINEX navigation file: {sys_in_rin_nav}')

            if currentGNSSsystem != 'C':
                sat_elevation_angles[sys] = sat_pos_dummy[currentGNSSsystem]['elevation'][:,0:37]
            else:
                sat_elevation_angles[sys] = sat_pos_dummy[currentGNSSsystem]['elevation']

        ## - Check for missing systems in navigation file, and remove if found
        missing_sys = []
        dummy_GNSSsystems = GNSSsystems.copy()
        for key,sys in dummy_GNSSsystems.items():
            if len(sat_pos[sys]['position']) == 0:
                del sat_pos[sys]
                del GNSSsystems[key]
                missing_sys.append(sys)
        if len(missing_sys) != 0:
            print('\n\nSystems %s does not exist in navigation file! \nMultipath analysis for these systems is therefore not possible. \nConsider using another navigation file.\n\n' % (missing_sys))


    ## Define carrier frequencies for every GNSS system. Note: Carrier band numbers follow RINEX 3 convention
    nGNSSsystems = len(GNSSsystems)
    max_GLO_ID = 36

    ## -- Initialize cell for storing carrier bands for all systems
    frequencyOverview = {}
    ## -- Read frequenxies from overview file
    for i in np.arange(0,nGNSSsystems):
        curr_sys = GNSSsystems[i+1]
        frequencyOverview[i+1] = frequencyOverview_temp2[curr_sys]

    ## -- Change cell describing GLONASS carrier frequencies so that each satellite frequency is described.
    ## NOTE: This uses the GLONASS frequency channel information from RINEX 3

    ## -- Observation header
    if "R" in list(GNSSsystems.values()):
        GNSSsystemIndex = [k for k in GNSSsystems if GNSSsystems[k] == 'R'][0]
        try:
            GLOSatID = list(GLO_Slot2ChannelMap.keys())
        except:
            if glo_fcn:
                GLO_Slot2ChannelMap = glo_fcn
                GLOSatID = list(GLO_Slot2ChannelMap.keys())
            else:
                raise ValueError("ERROR! GLONASS k-numbers do not exist. This is mandatory to be able to run analysis for GLONASS. Please add GLONASS SLOT / FRQ  to RINEX header.\
                                or use a rinex navigation file instead of SP3.")
        frequencyOverviewGLO = np.full([9,max_GLO_ID+1], np.nan)
        for k in np.arange(0,9):
            for j in np.arange(0,max_GLO_ID):
                if j in GLOSatID:
                    frequencyOverviewGLO[k, j] = frequencyOverview[GNSSsystemIndex][k,0] + \
                    GLO_Slot2ChannelMap[j] * frequencyOverview[GNSSsystemIndex][k,1] # added +1 and remove axis 1

        # store GLONASS carrier frequencies in their new dicture
        frequencyOverview[GNSSsystemIndex] = frequencyOverviewGLO


    ## Create observation type overview
    ## create overview for each GNSS system. Each overview gives the observation
    ## types for each of the RINEX 3 convention carrier bands. As no system has
    ## observations on all 9 RINEX convention bands many of these "bands" will
    ## be empty.

    # Remove obscode that exist in RINEX file header but dont contain any data
    # obsCodes = remove_obscodes_missing_data_in_rinexobs(GNSS_obs, GNSSsystems, obsCodes)

    obsCodeOverview = {}
    for i in np.arange(0,nGNSSsystems):
        GNSSsystemIndex = list(GNSSsystems.keys())[i]
        curr_sys = GNSSsystems[GNSSsystemIndex]
        obsCodeOverview[GNSSsystemIndex] = {}
        CODES = [x for x in obsCodes[GNSSsystemIndex][curr_sys] if 'C' in x[0] or 'P' in x[0]] #P is used in RINEX v2
        band_list = [band[1] for band in CODES]
        for j in np.arange(1,10):
            obsCodeOverview[GNSSsystemIndex][str(j)] = [] # preallocating slots for band (make 9 slots anyway)
        for band in band_list:
            obs_dum = [obs for obs in CODES if obs[1] == band]
            obsCodeOverview[GNSSsystemIndex][band] = obs_dum


    ### --- Build the dicture of the dict used for storing results ----

    ## --Initialize variable storing total number of observation codes processed
    nCodes_Total = 0
    ## -- Initialize results dictionary
    analysisResults = {}
    analysisResults['nGNSSsystem'] = nGNSSsystems
    analysisResults['GNSSsystems'] = list(GNSSsystems.values())
    # for sys in range(0,nGNSSsystems):
    for sys in np.arange(0, nGNSSsystems):
        GNSSsystemName = GNSSsystemCode2Fullname[GNSSsystems[sys+1]] # Full name of current GNSS system
        analysisResults[GNSSsystemName] = {}  # Include full name of current GNSS system
        current_sys_dict = {} # Initialize dict for current GNSS system
        current_sys_dict['observationOverview'] = {} # Initialize observationOverview field as a dict
        ## -- Extract the possible bands of current GNSS system, example GPS: 1,2,5
        GNSSsystemPossibleBands = GNSSsystem2BandsMap[GNSSsystemName]
        nPossibleBands = len(GNSSsystemPossibleBands)

        for i in np.arange(0,int(max_sat[sys][0])):
            i = i + 1 # dont want sat_0 but sat_1
            ## create field for current sat. Field is dict
            current_sys_dict['observationOverview']['Sat_'+ str(i)] = {}
            current_sys_dict['observationOverview']['Sat_'+ str(i)]['Bands'] = []
            Bands_list = []
            for j in np.arange(0,nPossibleBands):
                current_sys_dict['observationOverview']['Sat_'+ str(i)]['n_possible_bands'] = nPossibleBands
                Bands_list.append("Band_" + str(GNSSsystemPossibleBands[j]))
                current_sys_dict['observationOverview']['Sat_'+ str(i)]['Band_' + str(GNSSsystemPossibleBands[j])] = ""
                current_sys_dict['observationOverview']['Sat_'+ str(i)]['Bands']  = Bands_list

        ## -- Initilize fields for current system dict
        current_sys_dict['nBands'] = 0
        current_sys_dict['Bands'] = {}
        Bands_list = []
        for bandNumInd in np.arange(0,9):
            ## See if current system has any observations in in carrier band(bandnum)
            bandNumInd = str(bandNumInd+1) # because python nullindexed
            nCodes_currentBand = len(obsCodeOverview[sys+1][bandNumInd])
            if nCodes_currentBand > 0:
                current_sys_dict['nBands'] += 1 # Increment number of bands for current system dict
                Bands_list.append("Band_" + str(bandNumInd)) # Append current band to list of band for this system dict
                current_sys_dict['Bands'] =Bands_list
                current_band_dict = {} # Create field for this band, field dict
                current_band_dict['nCodes'] = nCodes_currentBand  # Store number of codes in current band
                nCodes_Total = nCodes_Total + nCodes_currentBand  # Increment total number of codes processed
                curr_band = current_sys_dict['Bands'][current_sys_dict['nBands']-1]
                current_band_dict['Codes'] = obsCodeOverview[sys+1][bandNumInd] # Store codes for this band
                current_sys_dict[curr_band] = current_band_dict # Store current band dict as field in current system dict
        # Store current systems dict as field in results dict
        analysisResults[GNSSsystemName] = current_sys_dict

    #### -------- Execute analysis of current data, produce results and store results in results dictionary ------
    codeNum = 0
    ## -- Defining frrmat of progressbar
    bar_format = '{desc}: {percentage:3.0f}%|{bar}| ({n_fmt}/{total_fmt})'
    for sys in np.arange(0,nGNSSsystems):    # replaced "range" with np.arange for speed
        currentGNSSsystem = GNSSsystems[sys+1]  # Get current GNSS system code, example GPS: G
        GNSSsystemName = GNSSsystemCode2Fullname[GNSSsystems[sys+1]]
        current_sys_dict = analysisResults[GNSSsystemName]
        nBands = current_sys_dict['nBands'] # Get number of carrier bands in current system dict
        ## -- Itterate through Bands in system dict.
        ## NOTE variable "bandNumInd" is NOT the carrier band number, but the index of that band in this system dict
        n_signals= sum(current_sys_dict[current_sys_dict['Bands'][bandNumInd]]['nCodes'] for bandNumInd  in range(0,nBands))
        pbar = tqdm(total=n_signals, desc='Currently processing all available signals for %s' % (GNSSsystemName), position=0, leave=True, bar_format=bar_format)
        # for bandNumInd in trange(0,nBands,initial=0, desc='Currently processing all available bands for %s' % (GNSSsystemName), leave=False,bar_format=bar_format,position=0):
        for bandNumInd in np.arange(0,nBands): #,initial=0, desc='Currently processing all available bands for %s' % (GNSSsystemName), leave=False,bar_format=bar_format,position=0):
            current_band_dict = current_sys_dict[current_sys_dict['Bands'][bandNumInd]] # Make HARD copy of current band dict
            nCodes = current_band_dict['nCodes'] # Get number of codes for current band dict
            currentBandName = current_sys_dict['Bands'][bandNumInd] # Get current band full name
            ## For each code pseudorange observation in current band dict,
            ## execute analysis once with every other signal in othe band to
            ## create linear combination. The analysis with the most estimates
            ## is the analysis that is stored.
            for i in np.arange(0,nCodes): # replaced "range" with np.arange for speed
                ## -- Get code(range) and phase obervation codes
                range1_Code = current_band_dict['Codes'][i]
                phase1_Code = "L" + range1_Code[1::]
                ## --Increment code counter and update waitbar
                codeNum = codeNum + 1
                if phase1_Code in obsCodes[sys+1][currentGNSSsystem]:
                    ## Initialize variable storing the best number for estimates
                    ## for the different analysis on current code
                    best_nEstimates = 0
                    best_currentStats = np.nan
                    # Itterate through the codes in the other bands to execute analysis with them and current range1 code
                    for secondBandnum in np.arange(0,nBands):  # replaced "range" with np.arange for speed
                        # Disregard observation code in same carrier band as current range1 observation
                        if secondBandnum != bandNumInd:
                            other_band_dict = current_sys_dict[current_sys_dict['Bands'][secondBandnum]] # Make HARD copy of the other band dict
                            nCodesOtherBand = other_band_dict['nCodes']  # Get number of codes in other band dict
                            # Itterate through codes in other band
                            for k in np.arange(0,nCodesOtherBand):
                                ## Get code(range) and phase obsertion codes from the other band
                                range2_Code = other_band_dict['Codes'][k]
                                if range2_Code == []:
                                    continue
                                phase2_Code = "L" + range2_Code[1::]
                                ## Check if phase2 observation was read from RINEX 3 observtaion file
                                if phase2_Code in obsCodes[sys+1][currentGNSSsystem]:
                                    # Test if some signals cotains only zeros
                                    range1_Code_idx = obsCodes[sys+1][currentGNSSsystem].index(range1_Code)
                                    phase1_Code_idx = obsCodes[sys+1][currentGNSSsystem].index(phase1_Code)
                                    phase2_Code_idx = obsCodes[sys+1][currentGNSSsystem].index(phase2_Code)
                                    obs_values = np.stack(list(GNSS_obs[currentGNSSsystem].values()))
                                    obs_range1 = obs_values[:, :, range1_Code_idx]
                                    obs_phase1 = obs_values[:, :, phase1_Code_idx]
                                    obs_phase2 = obs_values[:, :, phase2_Code_idx]
                                    if np.all(obs_range1 == 0) or np.all(obs_phase1 == 0) or np.all(obs_phase2 == 0):
                                        logger.warning(f"INFO(GNSS_MultipathAnalysis): One or more of the following observation codes {range1_Code},{phase1_Code} and {phase2_Code} ({GNSSsystemName}),"\
                                                       " lack data for the entire observation period. Therefore, this linear combination cannot be utilized.")
                                        continue

                                    ## Execute the analysis of current combination of observations. Return statistics on analysis
                                    current_max_sat = int(max_sat[sys].item()) if hasattr(max_sat[sys], "item") else int(max_sat[sys])
                                    currentStats, success = signalAnalysis(
                                        currentGNSSsystem, range1_Code, range2_Code, GNSSsystems, frequencyOverview, nepochs, tInterval,
                                        current_max_sat, GNSS_SVs[currentGNSSsystem], obsCodes[sys+1], GNSS_obs[currentGNSSsystem], GNSS_LLI[currentGNSSsystem],
                                        sat_elevation_angles[sys], phaseCodeLimit, ionLimit, cutoff_elevation_angle
                                    )

                                    if not success:
                                        return success

                                    ##  -- Get number of estimates produced from analysis
                                    current_nEstimates = currentStats['nEstimates']
                                    if current_nEstimates == 0:
                                        logger.warning(f'INFO(GNSS_MultipathAnalysis): Estimates for signal combination {range1_Code}-{phase1_Code}-{phase2_Code} were not possible.'\
                                                       ' The reason could be a lack of simultaneous observations from the three signals.')

                                    ## -- Check if current analysis has more estimate than previous
                                    if current_nEstimates > best_nEstimates:
                                        ## store current analysis results as "best so far"
                                        best_nEstimates = current_nEstimates
                                        best_currentStats = currentStats

                                ## If phase2 observation is not read from RINEX 3 observation file
                                else:
                                    pbar.update(1)
                                    logger.warning(f"INFO(GNSS_MultipathAnalysis): {range2_Code} code exists in RINEX observation file, but not {phase2_Code}. Linear combinations using this signal are not used.")
                                    other_band_dict['Codes'][ismember(other_band_dict['Codes'], range2_Code)] = [] # Remove range1 observation dict from other band dict, as it can not be used later
                                    other_band_dict["nCodes"] -= 1  #  Deincrement number of codes in other band dict
                                    current_sys_dict[current_sys_dict['Bands'][secondBandnum]] = other_band_dict # replace the, now altered, hard copy of other band dict in its original place in system dict

                    ## -- Store best analysis result dict in current band dict
                    current_code_dict = best_currentStats
                    ## For every satellite that had an observation of range1, it
                    ## is stored in an overview. Hence the user can get overview
                    ## of which satellites have transmitted which observations
                    ## number of satellites for current system, observation or no
                    if type(current_code_dict) == dict:
                        nSat = len(current_code_dict['range1_slip_distribution_per_sat'])
                    else:
                        pbar.update(1)
                        continue

                    for sat in np.arange(0,nSat):
                        ## -- If current satellite had observation of range1 code
                        if current_code_dict['n_range1_obs_per_sat'][0,sat+1] > 0:
                            ## Name of satellite
                            satCode = 'Sat_' + str(sat+1)
                            ## -- Check that code has not been added to list by fault
                            if current_sys_dict['observationOverview'][satCode][currentBandName] !=  current_code_dict['range1_Code']:
                                ## -- Add current range1 code to string of codes for  current satellite, sorted into bands
                                if current_sys_dict['observationOverview'][satCode][currentBandName] == "":
                                    current_sys_dict['observationOverview'][satCode][currentBandName] = current_sys_dict['observationOverview'][satCode][currentBandName] + current_code_dict['range1_Code']
                                else:
                                    current_sys_dict['observationOverview'][satCode][currentBandName] =  current_sys_dict['observationOverview'][satCode][currentBandName] + ', ' + current_code_dict['range1_Code']


                    if plotEstimates:
                        ## -- Plot and save graphs
                        if use_LaTex: #check if use TEX. Adding try/except to handle if user dont have TEX installed
                            try:
                                plotResults(current_code_dict['ion_delay_phase1'], current_code_dict['multipath_range1'], \
                                    current_code_dict['sat_elevation_angles'], tInterval, currentGNSSsystem, \
                                    current_code_dict['range1_Code'], current_code_dict['range2_Code'], \
                                    current_code_dict['phase1_Code'], current_code_dict['phase2_Code'], graphDir)
                            except:
                                latex_installed = False
                                plotResults_dont_use_TEX(current_code_dict['ion_delay_phase1'], current_code_dict['multipath_range1'], \
                                    current_code_dict['sat_elevation_angles'], tInterval, currentGNSSsystem, \
                                    current_code_dict['range1_Code'], current_code_dict['range2_Code'], \
                                    current_code_dict['phase1_Code'], current_code_dict['phase2_Code'], graphDir)
                        else:
                            plotResults_dont_use_TEX(current_code_dict['ion_delay_phase1'], current_code_dict['multipath_range1'], \
                                  current_code_dict['sat_elevation_angles'], tInterval, currentGNSSsystem, \
                                  current_code_dict['range1_Code'], current_code_dict['range2_Code'], \
                                  current_code_dict['phase1_Code'], current_code_dict['phase2_Code'], graphDir)

                      ## -- Place the current code dict in its original place in current band dict
                    current_band_dict[range1_Code] = current_code_dict
                    pbar.update(1)
                else:
                    ## If phase1 observation is not read from RINEX observation file
                    pbar.update(1)
                    logger.warning(f"INFO(GNSS_MultipathAnalysis): {range1_Code} code exists in RINEX observation file, but not {phase1_Code}\n'\
                                   'Linear combination using this signal is not used.")

                    current_band_dict['Codes'][ismember(current_band_dict['Codes'], range1_Code)] = []
                    current_band_dict['nCodes'] = current_band_dict['nCodes'] - 1


            ## -- Replace the, now altered, hard copy of current band dict in its original place in system dict
            current_sys_dict[current_sys_dict['Bands'][bandNumInd]] = current_band_dict

        # Store the satellite position,azimuths and elevation angles if SP3 files is used
        if sp3NavFilename_1 != '':
            try:
                sat_pos[currentGNSSsystem] = {}
                sat_pos[currentGNSSsystem]['position']  = sat_coordinates[currentGNSSsystem]
                sat_pos[currentGNSSsystem]['azimuth']    = sat_azimut_angles[sys]
                sat_pos[currentGNSSsystem]['elevation'] = sat_elevation_angles[sys]
            except:
                pass


        ## -- Replace the, now altered, hard copy of current system dict in its original place in results dict
        analysisResults[GNSSsystemName] = current_sys_dict
        ## -- Store information needed for output file in result dict
        rinex_obs_filename = rinObsFilename.split('/')
        rinex_obs_filename = rinex_obs_filename[-1]
        analysisResults['ExtraOutputInfo']  = {}
        analysisResults['ExtraOutputInfo']['rinex_obs_filename']       = rinex_obs_filename
        analysisResults['ExtraOutputInfo']['markerName']               = markerName
        analysisResults['ExtraOutputInfo']['rinexVersion']             = rinexVersion
        analysisResults['ExtraOutputInfo']['rinexProgr']               = rinexProgr
        analysisResults['ExtraOutputInfo']['recType']                  = recType
        analysisResults['ExtraOutputInfo']['Rinex_Receiver_Approx_Pos'] = np.atleast_2d(approxPosition).flatten().tolist()
        analysisResults['ExtraOutputInfo']['tFirstObs']                = tFirstObs
        analysisResults['ExtraOutputInfo']['tLastObs']                 = tLastObs
        analysisResults['ExtraOutputInfo']['tInterval']                = tInterval
        analysisResults['ExtraOutputInfo']['time_epochs_gps_time']     = time_epochs
        analysisResults['ExtraOutputInfo']['time_epochs_utc_time']     = gpstime_to_utc_datefmt(time_epochs)
        analysisResults['ExtraOutputInfo']['GLO_Slot2ChannelMap']      = GLO_Slot2ChannelMap
        analysisResults['ExtraOutputInfo']['nEpochs']                  = nepochs
        analysisResults['ExtraOutputInfo']['elevation_cutoff']         = cutoff_elevation_angle

        if estimated_position is not None:
            analysisResults['ExtraOutputInfo']['Estimated_Receiver_Approx_Pos'] = np.round(estimated_position.flatten(),4).tolist()[:-1]
            analysisResults['ExtraOutputInfo']['Estimated_Receiver_Approx_Pos_stats'] = stats


        ## -- Store default limits or user set limits in dict
        if phaseCodeLimit == 0:
            analysisResults['ExtraOutputInfo']['phaseCodeLimit']  = 4/60*100
        else:
            analysisResults['ExtraOutputInfo']['phaseCodeLimit']  = phaseCodeLimit

        if ionLimit == 0:
            analysisResults['ExtraOutputInfo']['ionLimit']  = 4/60
        else:
            analysisResults['ExtraOutputInfo']['ionLimit']  = ionLimit

        if sp3NavFilename_1 != "":
            sp3_list = [sp3NavFilename_1,sp3NavFilename_2,sp3NavFilename_3]
            analysisResults['ExtraOutputInfo']['SP3_filename'] = [os.path.basename(sp3) for sp3 in sp3_list if sp3 !=""]

        if broadcastNav1 != "":
            nav_list = [broadcastNav1,broadcastNav2,broadcastNav3,broadcastNav4]
            analysisResults['ExtraOutputInfo']['rinex_nav_filename'] = [os.path.basename(nav) for nav in nav_list if nav !=""] #added 19.02.2023

        ## -- Compute number of receiver clock jumps and store
        nClockJumps, meanClockJumpInterval, stdClockJumpInterval = detectClockJumps(GNSS_obs, nGNSSsystems, obsCodes, time_epochs, tInterval,GNSSsystems)
        analysisResults['ExtraOutputInfo']['nClockJumps'] = nClockJumps
        analysisResults['ExtraOutputInfo']['meanClockJumpInterval'] = meanClockJumpInterval
        analysisResults['ExtraOutputInfo']['stdClockJumpInterval']  = stdClockJumpInterval

        pbar.close()
    if 'sat_pos' in locals(): # add satellite position,azimut, elevation to analysResults
        analysisResults['Sat_position'] = sat_pos

    print('\n\nINFO: Analysis complete!\n')
    if latex_installed == False:
        logger.warning("INFO(GNSS_MultipathAnalysis): Use of TEX was enabled, but not installed on your computer! Install that to get prettier text formatting in plots.")

    baseFileName = os.path.basename(rinObsFilename)
    outputFilename = baseFileName.split('.')[0] +   '_Report.txt'
    writeOutputFile(outputFilename, outputDir, analysisResults, includeResultSummary, includeCompactSummary, includeObservationOverview, includeLLIOverview)
    print('INFO: The output file %s has been written.\n' % (outputFilename))


    ## -- Make barplot if plotEstimates is True
    if plotEstimates:
        print('INFO: Making bar plot. Please wait...\n')
        if use_LaTex:
            try:
                make_barplot(analysisResults,graphDir)
            except:
                make_barplot_dont_use_TEX(analysisResults,graphDir)
        else:
            make_barplot_dont_use_TEX(analysisResults,graphDir)


        ## -- Make skyplot of all systems
        GNSS_Name2Code =  dict(zip(['GPS', 'GLONASS', 'Galileo', 'BeiDou'], ['G', 'R', 'E', 'C'])) #mmaoung from name to code
        for sys in analysisResults['GNSSsystems']:
            curr_sys = GNSS_Name2Code[sys]
            try:
                azimut_currentSys = analysisResults['Sat_position'][curr_sys]['azimuth']
                elevation_currentSys = analysisResults['Sat_position'][curr_sys]['elevation']
                print('INFO: Making a regular polar plot for showing azimut and elevation angle for each satellite. Please wait...')
                if use_LaTex:
                    try:
                        make_skyplot(azimut_currentSys,elevation_currentSys,sys,graphDir)
                    except:
                        make_skyplot_dont_use_TEX(azimut_currentSys,elevation_currentSys,sys,graphDir)
                else:
                    make_skyplot_dont_use_TEX(azimut_currentSys,elevation_currentSys,sys,graphDir)

            except:
                print('Skyplot is not possible for %s! Missing data.' % (sys))
                pass

        if plot_polarplot:
            print('INFO: Making a polar plot of the multipath effect. Please wait ...')
            if use_LaTex:
                try:
                    make_polarplot(analysisResults, graphDir)
                except:
                    make_polarplot_dont_use_TEX(analysisResults, graphDir)
            else:
                make_polarplot_dont_use_TEX(analysisResults, graphDir)


    if include_SNR:
        # Seaching for SNR codes
        for sys in GNSS_obs.keys():
            GNSSsystemIndex = [k for k in GNSSsystems if GNSSsystems[k] == sys][0]
            SNR_codes = [SNR_code for SNR_code in obsCodes[GNSSsystemIndex][sys] if 'S' in SNR_code[0]]
        if plotEstimates and len(SNR_codes) != 0:
            print('INFO: Making a plot of the Signal To Noise Ration (SNR). Please wait ...')
            if use_LaTex:
                try:
                    make_polarplot_SNR(analysisResults,GNSS_obs,GNSSsystems,obsCodes,graphDir)
                    plot_SNR_wrt_elev(analysisResults,GNSS_obs,GNSSsystems,obsCodes,graphDir,tInterval)
                except:
                    make_polarplot_SNR_dont_use_TEX(analysisResults,GNSS_obs,GNSSsystems,obsCodes,graphDir)
                    plot_SNR_wrt_elev_dont_use_TEX(analysisResults,GNSS_obs,GNSSsystems,obsCodes,graphDir,tInterval)
            else:
                make_polarplot_SNR_dont_use_TEX(analysisResults,GNSS_obs,GNSSsystems,obsCodes,graphDir)
                plot_SNR_wrt_elev_dont_use_TEX(analysisResults,GNSS_obs,GNSSsystems,obsCodes,graphDir,tInterval)
        else:
            logger.warning("INFO(GNSS_MultipathAnalysis): There is no SNR codes available in the RINEX files. Plot of the Signal To Noise Ration is not possible.")

        for syst in GNSSsystems.keys():
            curr_syst =GNSSsystemCode2Fullname[GNSSsystems[syst]]
            analysisResults[curr_syst]['SNR']  = {}
            curr_obscodes = obsCodes[syst][GNSSsystems[syst]]
            snr_codes_idx =  [n for n, l in enumerate(curr_obscodes) if l.startswith('S')]
            for code_idx in snr_codes_idx:
                signal = curr_obscodes[code_idx]
                curr_ban = [element for element in analysisResults[curr_syst].keys() if element.endswith(signal[1])][0]
                curr_signal = np.stack(list(GNSS_obs[GNSSsystems[syst]].values()))[:, :, code_idx]
                curr_signal = np.squeeze(curr_signal)
                analysisResults[curr_syst]['SNR'][signal] = curr_signal


    if write_results_to_csv:
        result_dir = os.path.join(outputDir, "Result_files_CSV")
        if not os.path.exists(result_dir):
            os.makedirs(result_dir)
        createCSV = createCSVfile(analysisResults, result_dir, output_csv_delimiter)
        createCSV.write_results_to_csv()


    ## -- Saving the workspace as a binary pickle file ---
    if save_results_as_pickle:
        pickle_filename = 'analysisResults.pkl'
        print(f'\nINFO: The analysis results are being written to the file {pickle_filename}. Please wait..')
        results_name = os.path.join(outputDir, pickle_filename)
        PickleHandler.write_pickle(analysisResults, results_name)
        print(f'INFO: The analysis results has been written to the file {pickle_filename}.\n')
    elif save_results_as_compressed_pickle:
        pickle_filename = 'analysisResults.pkl.zst'
        print(f'\nINFO: The analysis results are being written to the file {pickle_filename}. Please wait..')
        results_name = os.path.join(outputDir, pickle_filename)
        PickleHandler.write_zstd_pickle(analysisResults, results_name)
        print(f'INFO: The analysis results has been written to the file {pickle_filename}.\n')

    end_time = time.time()
    compute_processing_time(start_time, end_time)
    logging.shutdown()

    return analysisResults


def compute_processing_time(start_time, end_time):
    """ Computes the processing time"""
    total_time_seconds = end_time - start_time
    hours = str(int(total_time_seconds // 3600)).zfill(2)
    minutes = str(int((total_time_seconds % 3600) // 60)).zfill(2)
    seconds = str(int(total_time_seconds % 60)).zfill(2)
    print(f"INFO: Finished! Processing time: {hours}:{minutes}:{seconds}")
    return


def ismember(list_,code):
    """
    The function takes in a string and a list, and finds the index of
    """
    indx = [idx for idx, val in enumerate(list_) if val == code]
    if indx != []:
        indx = indx[0]
    return indx


def filter_common_gnss_keys(sat_pos, GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, GNSSsystems):
    """
    Filters the GNSS-related dictionaries to include only common systems present in `sat_pos_dummy` keys.

    Parameters:
    ----------
    - sat_pos (dict): Dictionary with satellite position data.
    - GNSS_obs (dict): Dictionary of GNSS observations.
    - GNSS_LLI (dict): Dictionary of GNSS Loss of Lock Indicators.
    - GNSS_SS (dict): Dictionary of GNSS Signal Strengths.
    - GNSS_SVs (dict): Dictionary of GNSS Space Vehicles.
    - GNSSsystems (dict): Dictionary of GNSS systems.

    Returns:
    --------
    - Tuple[dict, dict, dict, dict, dict]: Filtered GNSS dictionaries with common systems.
    """
    # Intersect the keys in sat_pos_dummy and GNSSsystems
    common_systems = set(sat_pos.keys()).intersection(set(GNSSsystems.values()))

    # Filter each dictionary to only include common systems
    GNSSsystems_filtered = {k: v for k, v in GNSSsystems.items() if v in common_systems}
    GNSS_obs_filtered = {k: v for k, v in GNSS_obs.items() if k in common_systems}
    GNSS_LLI_filtered = {k: v for k, v in GNSS_LLI.items() if k in common_systems}
    GNSS_SS_filtered = {k: v for k, v in GNSS_SS.items() if k in common_systems}
    GNSS_SVs_filtered = {k: v for k, v in GNSS_SVs.items() if k in common_systems}

    return GNSS_obs_filtered, GNSS_LLI_filtered, GNSS_SS_filtered, GNSS_SVs_filtered, GNSSsystems_filtered




if __name__== "__main__":
    pass
