import os, sys,numpy as np,pickle

from readRinexObs304 import *
from Geodetic_functions import *

from computeSatElevations import computeSatElevations
from computeSatElevAimut_fromNav import computeSatElevAimut_fromNav
from readFrequencyOverview import readFrequencyOverview
from signalAnalysis import signalAnalysis
from plotResults import plotResults
from detectClockJumps import detectClockJumps
from tqdm import tqdm, trange
from writeOutputFile import writeOutputFile
from make_polarplot import make_polarplot
from plotResults import *

def GNSS_MultipathAnalysis(rinObsFilename,
                          broadcastNav1=None,
                          broadcastNav2=None,
                          broadcastNav3=None,
                          broadcastNav4=None,
                          sp3NavFilename_1=None,
                          sp3NavFilename_2=None, 
                          sp3NavFilename_3=None,
                          desiredGNSSsystems=None,
                          phaseCodeLimit= None, 
                          ionLimit=None, 
                          cutoff_elevation_angle=None,
                          outputDir=None,
                          plotEstimates= None,
                          plot_polarplot=None,
                          includeResultSummary= None,
                          includeCompactSummary=None,
                          includeObservationOverview=None,
                          includeLLIOverview= None
                          ):
    
    """
    Function that, through the help of other functions: 
      
      - reads RINEX 3 observation files
      - reads SP3 satellite navigation files
      - computes elevation angles of satellites for every epoch
      - makes estimates of multipath, ionospheric delay, and ambiguity slips
          for all signals in RINEX 3 observation file
      - plots estimates in graphs, if user desired it to
      - computes and stores statistics on estimates
      - writes an output files containing results
    
    This function calls on a range of functions. These in turn call on
    further more functions. The function called upon directly in 
    GNSS_MultipathAnalysis are:
    
      - readRinexObs304.m
      - computeSatElevations.m
      - readFrequencyOverview.m
      - signalAnalysis.m
      - writeOutputFile.m
    
    --------------------------------------------------------------------------------------------------------------------------
    INPUTS:
    
    rinObsFilename:           string. Path to RINEX 3 observation file
    
    sp3NavFilename_1:         string. Path to first SP3 navigation file
    
    sp3NavFilename_2:         string. Path to second SP3 navigation file. If
                              no second SP3 file is to be used, this variable 
                              should be empty string, ""  
    
    sp3NavFilename_3:         string. Path to third SP3 navigation file. If
                              no third SP3 file is to be used, this variable 
                              should be empty string, ""
                              
    desiredGNSSsystems:       List with the desired GNSS system. Ex ['G','R'] if you want to 
                              only run the analysis on GPS and GLONASS. Default: All systems. (if set to None) 
    
    phaseCodeLimit:           critical limit that indicates cycle slip for
                              phase-code combination. Unit: m/s. If set to 0,
                              default value will be used
    
    ionLimit:                 critical limit that indicates cycle slip for
                              the rate of change of the ionopheric delay. 
                              Unit: m/s. If set to 0, default value will be used
    
    cutoff_elevation_angle    Critical cutoff angle for satellite elevation angles, degrees
                              Estimates where satellite elevation angle
                              is lower than cutoff are removed, so are
                              estimated slip periods
    
    outputDir:                string. Path to directory where output file 
                              should be generated. If user does not wish to 
                              specify output directory, this variable should 
                              be empty string, "". In this case the output file 
                              will be generated in sub-directory inside same 
                              directory as GNSS_Receiver_QC_2020.m
    
    plotEstimates:            boolean. 1 if user desires estimates to be
                              ploted in figures. 0 otherwise
    

    
    includeResultSummary:     boolean. 1 if user desires output file to
                              include more detailed overview of statistics, 
                              including for each individual satellites. 
                              0 otherwise
    
    includeCompactSummary:    booleab. 1 if user desired output file to
                              include more compact overview og statistics.
    
    includeObservationOverview:     boolean. 1 if user desires output file to
                                      include overview of obseration types observed
                                      by each satellite. 0 otherwise
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS:
    
    analysisResults:          A dictionary that contains alls results of all analysises, for all GNSS systems.
    --------------------------------------------------------------------------------------------------------------------------
    """
    if broadcastNav1 == None and sp3NavFilename_1 == None:
        raise RuntimeError("No SP3 or navigation file is defined! This is \
                           mandatory for this software, so please add one of them.")

    if broadcastNav1 !=None and sp3NavFilename_1 != None:
        raise RuntimeError("You defined both a navigation file and a SP3 file. Please\
                           choose between using broadcast ephemerides or precise.")
        
    
    if broadcastNav1 == None:    
        broadcastNav1 = ""  
    # if broadcastNav2 == None:    
    #     broadcastNav2 = 
    # if broadcastNav3 == None:    
    #     broadcastNav3 = ""
    
    if sp3NavFilename_1 == None:    
        sp3NavFilename_1 = ""  
    if sp3NavFilename_2 == None:    
        sp3NavFilename_2 = ""
    if sp3NavFilename_3 == None:    
        sp3NavFilename_3 = ""

    if phaseCodeLimit == None:
        phaseCodeLimit = 0

    if ionLimit == None:
        ionLimit = 0
    
    if cutoff_elevation_angle == None:
        cutoff_elevation_angle = 0
    
    if outputDir == None:
        outputDir = "" 
        
    if plotEstimates == None:
        plotEstimates = 1 
    if plot_polarplot == None:
        plot_polarplot = 1

    if includeResultSummary == None:
        includeResultSummary = 1 
        
    if includeCompactSummary == None:
        includeCompactSummary = 1 
        
    if includeObservationOverview == None:
        includeObservationOverview = 1 

    if includeLLIOverview == None:
        includeLLIOverview = 1

    if desiredGNSSsystems == None:
        includeAllGNSSsystems   = 1
        desiredGNSSsystems = ["G", "R", "E", "C"];  # All GNSS systems.
    else:
        includeAllGNSSsystems   = 0
    
    
    ## ---  Control of the user input arguments 
    if type(sp3NavFilename_1) != str:
        print('ERROR(GNSS_Receiver_QC_2020): The input variable sp3NavFilename_1 must be a string\n' \
            'Argument is now of type %s\n' %  (type(sp3NavFilename_1)))
        analysisResults = np.nan
        return
    
    if type(sp3NavFilename_2) != str:
        print('ERROR(GNSS_Receiver_QC_2020): The input variable sp3NavFilename_2 must be a string\n' \
            'Argument is now of type %s\n' %  (type(sp3NavFilename_2)))
        analysisResults = np.nan
        return
    
    
    if type(sp3NavFilename_3) != str:
        print('ERROR(GNSS_Receiver_QC_2020): The input variable sp3NavFilename_3must be a string\n' \
            'Argument is now of type %s\n' %  (type(sp3NavFilename_3)))
        analysisResults = np.nan
        return
    
    if type(rinObsFilename) != str:
        print('ERROR(GNSS_Receiver_QC_2020): The input variable rinObsFilename must be a string\n' \
            'Argument is now of type %s\n' %  (type(rinObsFilename)))
        analysisResults = np.nan
        return
 
    
    if not os.path.isfile(rinObsFilename):
        print('ERROR(GNSS_Receiver_QC_2020): RINEX observation file can not be found. Please check that the path is correct.\n') 
        analysisResults = np.nan
        return 
    
    
    if not os.path.isfile(sp3NavFilename_1) and len(sp3NavFilename_1) != 0:
        print('WARNING: Second SP3 Navigation file can not be found.\n')
    
    if not os.path.isfile(sp3NavFilename_2) and len(sp3NavFilename_2) != 0:
        print('WARNING: Second SP3 Navigation file can not be found.\n')
    
    
    if not os.path.isfile(sp3NavFilename_3) and len(sp3NavFilename_3) != 0:
        print('WARNING: Third SP3 Navigation file can not be found.\n')
    
    
    
    
    ## Make dictionaries
    GNSSsystemCode2Fullname = dict(zip(['G','R','E','C'],['GPS','GLONASS','Galileo','BeiDou']))
    GNSSsystem2BandsMap = dict(zip(['GPS', 'GLONASS', 'Galileo', 'BeiDou'],\
        [[1, 2, 5], [1, 2, 3, 4, 6], [1, 5, 6, 7, 8], [1, 2, 5, 6, 7, 8]]))
    
    
    ## --- Read observation file
    
    # includeAllGNSSsystems   = 0
    includeAllObsCodes      = 0
    desiredObsCodes = ["C", "L"];               # code and phase observations
    desiredObsBands = list(np.arange(1,10))     # all carrier bands. Tot 9, but arange stops at 8 -> 10
    # desiredGNSSsystems = ["G", "R", "E", "C"];  # All GNSS systems. 
    
    readSS = 1
    readLLI = 1
    
    ## --- Read RINEX 3.0x observation file
    [GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
        obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
        rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success] = \
        readRinexObs304(rinObsFilename, readSS, readLLI, includeAllGNSSsystems,includeAllObsCodes, desiredGNSSsystems,\
        desiredObsCodes, desiredObsBands)
    
    if sp3NavFilename_1 != '':
        # ## -- Compute satellite elevation angles from SP3 files
        # sat_elevation_angles = computeSatElevations(GNSS_SVs, GNSSsystems, approxPosition,\
        #     nepochs, time_epochs, max_sat, sp3NavFilename_1, sp3NavFilename_2, sp3NavFilename_3)
        ## -- Compute satellite elevation angles from SP3 files
        sat_elevation_angles,sat_azimut_angles = computeSatElevations(GNSS_SVs, GNSSsystems, approxPosition,\
            nepochs, time_epochs, max_sat, sp3NavFilename_1, sp3NavFilename_2, sp3NavFilename_3)
    else:
        nav_files = [broadcastNav1,broadcastNav2,broadcastNav3,broadcastNav4]
        sat_pos = computeSatElevAimut_fromNav(nav_files,approxPosition,GNSS_SVs,GNSS_obs,time_epochs)
        # sat_pos = computeSatElevAimut_fromNav(broadcastNav1,approxPosition,GNSS_SVs,GNSS_obs,time_epochs)
        
        ## -- Build same struture for satellit elevation angles if broadcast nav defined
        sat_elevation_angles = {}
        sat_pos_dummy = sat_pos.copy()
        for sys in range(0,len(GNSSsystems)):
            currentGNSSsystem = GNSSsystems[sys+1]
            if currentGNSSsystem != 'C':
                sat_elevation_angles[sys] = sat_pos_dummy[currentGNSSsystem]['Elevation'][:,0:37]
            else:
                sat_elevation_angles[sys] = sat_pos_dummy[currentGNSSsystem]['Elevation']
            
        ## - Check for missing systems in navigation file, and remove if found
        missing_sys = []
        dummy_GNSSsystems = GNSSsystems.copy()
        for key,sys in dummy_GNSSsystems.items():
            if len(sat_pos[sys]['Position']) == 0:
                del sat_pos[sys]
                del GNSSsystems[key]
                missing_sys.append(sys)
        if len(missing_sys) != 0:
            print('\n\nSystems %s does not exist in navigation file! \nMultipath analysis for these systems is therefore not possible. \nConsider using another navigation file.\n\n' % (missing_sys))

                
                

        
    ## Define carrier frequencies for every GNSS system. Note: Carrier band numbers follow RINEX 3 convention
    nGNSSsystems = len(GNSSsystems)
    max_GLO_ID = 36
    
    ## -- Read frequency overview file
    frequencyOverviewFilename = 'Rinex_Frequency_Overview.txt'
    frequencyOverview_temp, frequencyGNSSsystemOrder, _ = readFrequencyOverview(frequencyOverviewFilename)
    
    ## -- Making a new dummy frequencyOverview to rename key values in dict (make it able to run only one system)
    systems_in_overview = ['G','R','E','C'] # KAN FJERNE DETTE ETTER TESTING AV GLONASS ER FERDI. ELLER PRØV Å LA PROGRAMMET SLIK AT ANALYSEN KAN GJENNOMFØRS KUN PÅ ET FORHÅNDSBESTEMT SYSTEM
    frequencyOverview_temp2 =  dict(zip(systems_in_overview, list(frequencyOverview_temp.values())))
    
    ## -- Initialize cell for storing carrier bands for all systems
    frequencyOverview = {}
    ## -- Read frequenxies from overview file
    # for i in range(0,nGNSSsystems):
    #     GNSSsystemIndex = [k for k in GNSSsystems if GNSSsystems[k]==GNSSsystems[i+1][0]][0]
    #     frequencyOverview[i+1] = frequencyOverview_temp[GNSSsystemIndex]
    for i in range(0,nGNSSsystems):
        curr_sys = GNSSsystems[i+1]
        frequencyOverview[i+1] = frequencyOverview_temp2[curr_sys]
    
    ## -- Change cell describing GLONASS carrier frequencies so that each satellite frequency is described. 
    ## NOTE: This uses the GLONASS frequency channel information from RINEX 3
    
    ## -- Observation header
    # if ismember("R", string(GNSSsystems))
    if "R" in list(GNSSsystems.values()):
       # GNSSsystemIndex = find([GNSSsystems{:}] == "R"); 
       GNSSsystemIndex = [k for k in GNSSsystems if GNSSsystems[k] == 'R'][0]
       
       # GLOSatID = cell2mat(keys(GLO_Slot2ChannelMap));
       GLOSatID = list(GLO_Slot2ChannelMap.keys())
       frequencyOverviewGLO = np.full([9,max_GLO_ID+1], np.nan)
       for k in range(0,9):
           for j in range(0,max_GLO_ID):
               if j in GLOSatID:
                   frequencyOverviewGLO[k, j] = frequencyOverview[GNSSsystemIndex][k,0] + \
                   GLO_Slot2ChannelMap[j] * frequencyOverview[GNSSsystemIndex][k,1] # added +1 and remove axis 1
      
      
        
       # store GLONASS carrier frequencies in their new structure
       # frequencyOverview{GNSSsystemIndex} = frequencyOverviewGLO;
       frequencyOverview[GNSSsystemIndex] = frequencyOverviewGLO
    
    
    ## Create observation type overview
    ## create overview for each GNSS system. Each overview gives the observation
    ## types for each of the RINEX 3 convention carrier bands. As no system has 
    ## observations on all 9 RINEX convention bands many of these "bands" will 
    ## be empty.
    
    
    def ismember(list_,code):
        """
        The function takes in a string and a list, and finds the index of 
        """
        indx = [idx for idx, val in enumerate(list_) if val == code]
        if indx != []:
            indx = indx[0]
        return indx
    
    ## MIDLERTIDIG FJERN ETTERPÅ!!!
    
    
    
    obsCodeOverview = {}
    for i in range(0,nGNSSsystems):
        # obsCodeOverview{i} = cell(1,9);
        # GNSSsystemIndex = [k for k in GNSSsystems if GNSSsystems[k]==GNSSsystems[i+1]][0]
        GNSSsystemIndex = list(GNSSsystems.keys())[i]
        curr_sys = GNSSsystems[GNSSsystemIndex]
        obsCodeOverview[GNSSsystemIndex] = {}
        # CODES = [x for x in obsCodes[1]['G'] if 'C' in x[0]]
        CODES = [x for x in obsCodes[GNSSsystemIndex][curr_sys] if 'C' in x[0]]
        band_list = [band[1] for band in CODES]
        for j in range(1,10):
            obsCodeOverview[GNSSsystemIndex][str(j)] = [] # preallocating slots for band (make 9 slots anyway) 
    
        for band in band_list:
            obs_dum = [obs for obs in CODES if obs[1] == band]
            obsCodeOverview[GNSSsystemIndex][band] = obs_dum
    
    
    
    ### --- Build the structure of the struct used for storing results ----
    
    ## --initialize variable storing total number of observation codes processed
    nCodes_Total = 0;
    
    ## -- initialize results struct
    analysisResults = {}
    analysisResults['nGNSSsystem'] = nGNSSsystems
    analysisResults['GNSSsystems'] = list(GNSSsystems.values())
    # for sys in range(0,nGNSSsystems):
    for sys in np.arange(0,nGNSSsystems):
        ## -- Full name of current GNSS system
        GNSSsystemName = GNSSsystemCode2Fullname[GNSSsystems[sys+1]]
        ## -- Include full name of current GNSS system 
        analysisResults[GNSSsystemName] =  {}
        ## -- Initialize struct for current GNSS system
        current_sys_struct = {}
        ## -- Initialize observationOverview field as a struct
        current_sys_struct['observationOverview'] = {}
        ## -- Extract the possible bands of current GNSS system, example GPS: 1,2,5
        GNSSsystemPossibleBands = GNSSsystem2BandsMap[GNSSsystemName]
        nPossibleBands = len(GNSSsystemPossibleBands)
        
        # for i in range(0,int(max_sat[sys][0])):
        for i in np.arange(0,int(max_sat[sys][0])):
            i = i + 1 # dont want sat_0 but sat_1
            ## create field for current sat. Field is struct
            current_sys_struct['observationOverview']['Sat_'+ str(i)] = {}
            current_sys_struct['observationOverview']['Sat_'+ str(i)]['Bands'] = []
            Bands_list = []
            for j in range(0,nPossibleBands):
                current_sys_struct['observationOverview']['Sat_'+ str(i)]['n_possible_bands'] = nPossibleBands
                ## in sat. struct create 1 field for every possible band for this
                ## system as empty string
                Bands_list.append("Band_" + str(GNSSsystemPossibleBands[j]))
                current_sys_struct['observationOverview']['Sat_'+ str(i)]['Band_' + str(GNSSsystemPossibleBands[j])] = ""
                current_sys_struct['observationOverview']['Sat_'+ str(i)]['Bands']  = Bands_list
        ## -- Initilize fields for current system struct
        current_sys_struct['nBands'] = 0
        current_sys_struct['Bands'] = {}
    
        Bands_list = []
        for bandNumInd in range(0,9):
            ## See if current system has any observations in in carrier band(bandnum)
            # nCodes_currentBand = length(obsCodeOverview{sys}{bandNumInd});
            bandNumInd = str(bandNumInd+1) # because python nullindexed
            nCodes_currentBand = len(obsCodeOverview[sys+1][bandNumInd])
    
            if nCodes_currentBand > 0:
                ## --Increment number of bands for current system struct
                # current_sys_struct.nBands = current_sys_struct.nBands + 1;
                current_sys_struct['nBands'] = current_sys_struct['nBands'] + 1
                ## -- Append current band to list of band for thi system struct
                # current_sys_struct.Bands = [current_sys_struct.Bands, strcat("Band_", string(bandNumInd))];
                Bands_list.append("Band_" + str(bandNumInd))
                # current_sys_struct.Bands = [current_sys_struct.Bands, strcat("Band_", string(bandNumInd))];
                current_sys_struct['Bands'] =Bands_list
    
                ## -- Create field for this band, field struct
                # current_band_struct = struct;
                current_band_struct = {}
               
                ## -- Store number of codes in current band
                # current_band_struct.nCodes = nCodes_currentBand;
                current_band_struct['nCodes'] = nCodes_currentBand
    
                ## -- Increment total number of codes processed
                nCodes_Total = nCodes_Total + nCodes_currentBand
               
                ## -- Store codes for this band
                # current_band_struct.Codes = cellstr(obsCodeOverview{sys}{bandNumInd});
                curr_band = current_sys_struct['Bands'][current_sys_struct['nBands']-1]
                current_band_struct['Codes'] = obsCodeOverview[sys+1][bandNumInd]
    
                ## -- Store current band struct as field in current system struct 
                # current_sys_struct.(current_sys_struct.Bands{current_sys_struct.nBands}) = current_band_struct;
                # current_sys_struct['Bands'] = current_band_struct
                # current_sys_struct['Bands'][current_sys_struct['nBands']-1] = current_band_struct # subtrackting 1 to get index
                current_sys_struct[curr_band] = current_band_struct # subtrackting 1 to get index
    
            
        
        # Store current systems struct as field in results struct
        # analysisResults.(analysisResults.GNSSsystems{sys}) = current_sys_struct;
        # analysisResults['analysisResults']['GNSSsystems']['sys'] = current_sys_struct
        analysisResults[GNSSsystemName] = current_sys_struct # Her henter koden til bjørn Eirik verdien fra analysResult.Vurder å legge til hele systemnavnene i analyses result. 
                                                              #Ikke bare første bosktav
    
    
    #### -------- Execute analysis of current data, produce results and store results in results dictionary ------
    
    ## -- Intialize wait bar
    # wbar = waitbar(0, 'INFO(GNSS\\_Receiver\\_QC\\_2020): Data analysis is being executed. Please wait.');
    # waitbar = tqdm(0, 'INFO(GNSS\\_Receiver\\_QC\\_2020): Data analysis is being executed. Please wait.')
    ## --Initialize counter of number of codes processed so far. This is used  mainly for waitbar
    codeNum = 0
    
    ## -- Defining frrmat of progressbar
    bar_format = '{desc}: {percentage:3.0f}%|{bar}| ({n_fmt}/{total_fmt})'
    # for sys in range(0,nGNSSsystems):
    for sys in np.arange(0,nGNSSsystems):    # replaced "range" with np.arange for speed      
    # for sys in trange(0,nGNSSsystems,initial=0,desc='INFO: Data analysis is being executed...',bar_format=bar_format,position=0):
        ## --Get current GNSS system code, example GPS: G
        currentGNSSsystem = GNSSsystems[sys+1]
        GNSSsystemName = GNSSsystemCode2Fullname[GNSSsystems[sys+1]] # I USE THIS INSTEAD
                
        # print('INFO: Data analysis is being executed for %s' % (GNSSsystemName))
        ## --Make HARD copy of current system struct
        # current_sys_struct = analysisResults.(analysisResults.GNSSsystems{sys});
        current_sys_struct = analysisResults[GNSSsystemName]
        ## -- Get number of carrier bands in current system struct
        # nBands = current_sys_struct.nBands
        nBands = current_sys_struct['nBands']
        ## -- Itterate through Bands in system struct. 
        ## NOTE variable "bandNumInd" is NOT the carrier band number, but the index of that band in this system struct  
    
    
        # for bandNumInd in range(0,nBands): 
        for bandNumInd in trange(0,nBands,initial=0, desc='Currently processing all available bands for %s' % (GNSSsystemName), leave=False,bar_format=bar_format,position=0): 
            ## Make HARD copy of current band struct
            # current_band_struct = current_sys_struct.(current_sys_struct.Bands{bandNumInd});
            current_band_struct = current_sys_struct[current_sys_struct['Bands'][bandNumInd]]
            
            ## Get number of codes for current band struct
            # nCodes = current_band_struct.nCodes;
            nCodes = current_band_struct['nCodes']
            # Get current band full name
            # currentBandName = current_sys_struct.Bands{bandNumInd};
            currentBandName = current_sys_struct['Bands'][bandNumInd]
            
            ## For each code pseudorange observation in current band struct,
            ## execute analysis once with every other signal in othe band to
            ## create linear combination. The analysis with the most estimates
            ## is the analysis that is stored.
            
            # for i in range(0,nCodes):
            for i in np.arange(0,nCodes): # replaced "range" with np.arange for speed
            # for i in trange(0,nCodes,initial=1):
            # for i in tqdm(range(0,nCodes),desc='INFO(GNSS\\_Receiver\\_QC\\_2020): Data analysis is being executed.',position=1,leave=False):
                ## -- Get code(range) and phase obervation codes
                # range1_Code = current_band_struct.Codes{i};
                range1_Code = current_band_struct['Codes'][i]
                # phase1_Code = "L"+extractAfter(range1_Code, 1);
                phase1_Code = "L" + range1_Code[1::]
               
                ## --Increment code counter and update waitbar
                codeNum = codeNum + 1
                # msg = sprintf(['INFO(GNSS\\_Receiver\\_QC\\_2020): Data analysis is being executed. Please wait.\nCurrently processing ',...
                #     '%s %s signal'], GNSSsystemCode2Fullname(currentGNSSsystem), range1_Code);
                # waitbar(codeNum/nCodes_Total, wbar, msg);
                
                # waitbar.update(codeNum/nCodes_Total)
                # waitbar.set_description('INFO(GNSS\\_Receiver\\_QC\\_2020): Data analysis is being executed. Please wait.\nCurrently processing '\
                #     '%s %s signal' % (GNSSsystemCode2Fullname[currentGNSSsystem], range1_Code))
                ## --Check if phase observation matching code observation has been
                ## -- read from RINEX 3 observation file
                # if ismember(phase1_Code, obsCodes{sys})
                if phase1_Code in obsCodes[sys+1][currentGNSSsystem]:
                    ## Initialize variable storing the best number for estimates
                    ## for the different analysis on current code
                    best_nEstimates = 0
                    best_currentStats = np.nan
                   
                    ## Itterate through the codes in the other bands to execute
                    ## analysis with them and current range1 code
                    # for secondBandnum = 1:nBands
                    # for secondBandnum in range(0,nBands):
                    for secondBandnum in np.arange(0,nBands):  # replaced "range" with np.arange for speed
                        # Disregard observation code in same carrier band as
                        # current range1 observation
                        if secondBandnum != bandNumInd:
                            ## Make HARD copy of the other band struct
                            # other_band_struct = current_sys_struct.(current_sys_struct.Bands{secondBandnum});
                            other_band_struct = current_sys_struct[current_sys_struct['Bands'][secondBandnum]]
                            ## Get number of codes in other band struct
                            # nCodesOtherBand = other_band_struct.nCodes; 
                            nCodesOtherBand = other_band_struct['nCodes'] 
                         
                            # Itterate through codes in other band
                            # for k in range(0,nCodesOtherBand):
                            for k in np.arange(0,nCodesOtherBand):                                
                                  ## Get code(range) and phase obsertion codes from the other band
                                  # range2_Code = other_band_struct.Codes{k};
                                  # phase2_Code = "L"+extractAfter(range2_Code, 1);
                                  range2_Code = other_band_struct['Codes'][k]
                                  phase2_Code = "L" + range2_Code[1::]
                                  ## Check if phase2 observation was read from RINEX 3 observtaion file
    
                                  if phase2_Code in obsCodes[sys+1][currentGNSSsystem]:
                                       ## Execute the analysis of current combination of
                                       ## observations. Return statistics on analysis
                                       # [currentStats, success] = signalAnalysis(currentGNSSsystem, range1_Code, range2_Code, GNSSsystems, frequencyOverview, nepochs, tInterval, ...
                                       # max_sat(sys), GNSS_SVs{sys}, obsCodes{sys}, GNSS_obs{sys}, GNSS_LLI{sys}, sat_elevation_angles{sys}, phaseCodeLimit, ionLimit, cutoff_elevation_angle);
                                       # currentStats, success = signalAnalysis(currentGNSSsystem, range1_Code, range2_Code, GNSSsystems, frequencyOverview, nepochs, tInterval, \
                                       # int(max_sat[sys]), GNSS_SVs[currentGNSSsystem], obsCodes[sys+1], GNSS_obs[currentGNSSsystem], GNSS_LLI[currentGNSSsystem],\
                                       #     sat_elevation_angles[sys], phaseCodeLimit, ionLimit, cutoff_elevation_angle)
                                       
                                       currentStats, success = signalAnalysis(currentGNSSsystem, range1_Code, range2_Code, GNSSsystems, frequencyOverview, nepochs, tInterval, \
                                       int(max_sat[sys]), GNSS_SVs[currentGNSSsystem], obsCodes[sys+1], GNSS_obs[currentGNSSsystem], GNSS_LLI[currentGNSSsystem],\
                                           sat_elevation_angles[sys], phaseCodeLimit, ionLimit, cutoff_elevation_angle)
                
                                       if not success:
                                           pass ##   REMOVE AND ADD RETUNR!!!!!!
                                            # return success
                                       
                                    
                                    ##  -- Get number of estimates produced from analysis
                                    # current_nEstimates = currentStats.nEstimates;
                                       current_nEstimates = currentStats['nEstimates']
                                
                                       ## -- Check if current analysis has more estimate than previous
                                       if current_nEstimates > best_nEstimates:
                                            ## store current analysis results as "best so far"
                                            best_nEstimates = current_nEstimates
                                            best_range1 = range1_Code
                                            best_range2 = range2_Code
                                            best_currentStats = currentStats
                                    
                                    ## If phase2 observation is not read from RINEX 3 observation file
                                  else:
                                        print('\nINFO(GNSS_MultipathAnalysis): %s code exists in RINEX observation file, but not %s\n'\
                                            'Linear combinations using this signal is not used.\n\n' % (range2_Code, phase2_Code))
                                        ## -- Remove range1 observation struct from other band struct, as it can not be used later
                                        # other_band_struct.Codes(ismember(other_band_struct.Codes, range2_Code)) = [];
                                        other_band_struct['Codes'][ismember(other_band_struct['Codes'], range2_Code)] = []
                                        ## -- Deincrement numbe rof codes in otehr band struct
                                        # other_band_struct.nCodes = other_band_struct.nCodes - 1; 
                                        other_band_struct['nCodes'] = other_band_struct['nCodes'] - 1 
                                        ## -- replace the, now altered, hard copy of other band struct in its original place in system struct
                                        # current_sys_struct.(current_sys_struct.Bands{secondBandnum}) = other_band_struct;
                                        current_sys_struct[current_sys_struct['Bands'][secondBandnum]] = other_band_struct
    
                          
                    ## -- Display which analysis is the best(REMOVE LATER)
                    # fprintf('%s %s Winner: %s\n\n', currentGNSSsystem, best_range1, best_range2)
                  
                    ## -- Store best analysis result struct in current band struct
                    current_code_struct = best_currentStats
                  
                  
                    ## For every satellite that had an observation of range1, it
                    ## is stored in an overview. Hence the user can get overview
                    ## of which satellites have transmitted which observations
                    ## number of satellites for current system, observation or no
                    # nSat = length(current_code_struct.range1_slip_distribution_per_sat);
                    # if not np.isnan(current_code_struct):
                    #     nSat = len(current_code_struct['range1_slip_distribution_per_sat'])
                    # else:
                    #     # nSat = 0
                    #     break
                    if type(current_code_struct) == dict:
                        nSat = len(current_code_struct['range1_slip_distribution_per_sat'])
                    else:
                        print('\n\nWARNING! No estimates for band: "%s" and system: "%s". Probably because of only one obscode available for the current band. Observations from two different frequency is needed.' % (currentBandName,GNSSsystemName))
                        break
     
                    # for sat = 1:nSat
                    for sat in range(0,nSat):
                        ## -- If current satellite had observation of range1 code
                        # if current_code_struct.nRange1Obs_Per_Sat(sat)>0
                        if current_code_struct['n_range1_obs_per_sat'][0,sat+1] > 0:
                            ## Name of satellite struct
                            # satCode = strcat('Sat_', string(sat));
                            # satCode = strcat('Sat_', string(sat));
                            satCode = 'Sat_' + str(sat+1) 
                    
                            ## -- Check that code has not been added to list by fault
                            # if ~contains(current_sys_struct.observationOverview.(satCode).(currentBandName), current_code_struct.range1_Code)
                            if current_sys_struct['observationOverview'][satCode][currentBandName] !=  current_code_struct['range1_Code']:
                                ## -- Add current range1 code to string of codes for  current satellite, sorted into bands
                                # if strcmp(current_sys_struct.observationOverview.(satCode).(currentBandName), ""):
                                if current_sys_struct['observationOverview'][satCode][currentBandName] == "":
                                    # current_sys_struct.observationOverview.(satCode).(currentBandName) = strcat(current_sys_struct.observationOverview.(satCode).(currentBandName), current_code_struct.range1_Code);
                                    current_sys_struct['observationOverview'][satCode][currentBandName] = current_sys_struct['observationOverview'][satCode][currentBandName] + current_code_struct['range1_Code']
    
                                else:
                                  # current_sys_struct.observationOverview.(satCode).(currentBandName) = strcat(current_sys_struct.observationOverview.(satCode).(currentBandName), ', ', current_code_struct.range1_Code);
                                  current_sys_struct['observationOverview'][satCode][currentBandName] =  current_sys_struct['observationOverview'][satCode][currentBandName] + ', ' + current_code_struct['range1_Code']
    
                  
                   ## -- If plotEstimates boolean is 1, plot estimates from best analysis and store figures
                    if plotEstimates:
                      ## If user has not specified an output directory, set output directory to "Output_Files" 
                      # if isempty(outputDir)
                      #   outputDir = 'Output_Files';
                      # end
                      if outputDir == "":
                          outputDir = 'Output_Files'
                         
                      ## -- Unless output directory already exists, create output directory
                      if not os.path.isdir(outputDir):
                          os.mkdir(outputDir)
                      
                     
                      ## --Unless graph directory already exists, create directory
                      graphDir = outputDir + '/Graphs' 
                      if not os.path.isdir(graphDir):
                          os.mkdir(graphDir)
                      
                     
                      ## -- Plot and save graphs
                      plotResults(current_code_struct['ion_delay_phase1'], current_code_struct['multipath_range1'], \
                          current_code_struct['sat_elevation_angles'], tInterval, currentGNSSsystem, \
                          current_code_struct['range1_Code'], current_code_struct['range2_Code'], \
                          current_code_struct['phase1_Code'], current_code_struct['phase2_Code'], graphDir)
                      
    
                 
                      ## -- Place the current code struct in its original place in current band struct
                      current_band_struct[range1_Code] = current_code_struct
                  
                    else:
                        ## If phase1 observation is not read from RINEX observation file
                        print('\nINFO(GNSS_MultipathAnalysis): %s code exists in RINEX observation file, but not %s\n',\
                                        'Linear combination using this signal is not used.\n\n' % (range1_Code, phase1_Code))

                        current_band_struct['Codes'][ismember(current_band_struct['Codes'], range1_Code)] = []
                        current_band_struct['nCodes'] = current_band_struct['nCodes'] - 1 
                        
                    if plot_polarplot:
                        try:
                            azimut_currentSys = sat_pos[currentGNSSsystem]['Azimut']
                        except:
                            azimut_currentSys = sat_azimut_angles[sys]
                        make_polarplot(current_code_struct['multipath_range1'],current_code_struct['sat_elevation_angles'],\
                                       azimut_currentSys,GNSSsystemName,range1_Code ,graphDir)
                        
    
               
            ## -- Replace the, now altered, hard copy of current band struct in its original place in system struct
            # current_sys_struct.(current_sys_struct.Bands{bandNumInd}) = current_band_struct;
            current_sys_struct[current_sys_struct['Bands'][bandNumInd]] = current_band_struct
           
        
        
        ## -- Replace the, now altered, hard copy of current system struct in its original place in results struct
        # analysisResults.(analysisResults.GNSSsystems{sys}) = current_sys_struct;
        analysisResults[GNSSsystemName] = current_sys_struct
    
    
        ## -- Store information needed for output file in result struct
        # rinex_obs_filename = strsplit(rinObsFilename, '/');
        # rinex_obs_filename = rinex_obs_filename{end};
        rinex_obs_filename = rinObsFilename.split('/')
        rinex_obs_filename = rinex_obs_filename[-1]
        
        analysisResults['ExtraOutputInfo']  = {}
        analysisResults['ExtraOutputInfo']['rinex_obs_filename']  = rinex_obs_filename
        analysisResults['ExtraOutputInfo']['markerName']          = markerName
        analysisResults['ExtraOutputInfo']['rinexVersion']        = rinexVersion
        analysisResults['ExtraOutputInfo']['rinexProgr']          = rinexProgr
        analysisResults['ExtraOutputInfo']['recType']             = recType
        analysisResults['ExtraOutputInfo']['tFirstObs']           = tFirstObs
        analysisResults['ExtraOutputInfo']['tLastObs']            = tLastObs
        analysisResults['ExtraOutputInfo']['tInterval']           = tInterval
        analysisResults['ExtraOutputInfo']['GLO_Slot2ChannelMap'] = GLO_Slot2ChannelMap
        analysisResults['ExtraOutputInfo']['nEpochs']             = nepochs
        
        ## -- Store default limits or user set limits in struct
        if phaseCodeLimit == 0:
            analysisResults['ExtraOutputInfo']['phaseCodeLimit']  = 4/60*100
        else:
            analysisResults['ExtraOutputInfo']['phaseCodeLimit']  = phaseCodeLimit
        
        
        if ionLimit == 0:
            analysisResults['ExtraOutputInfo']['ionLimit']  = 4/60
        else:
            analysisResults['ExtraOutputInfo']['ionLimit']  = ionLimit
        
        
        
        ## -- Compute number of receiver clock jumps and store
        # [nClockJumps, meanClockJumpInterval, stdClockJumpInterval] = detectClockJumps(GNSS_obs, nGNSSsystems, obsCodes, time_epochs, tInterval);
        # analysisResults.ExtraOutputInfo.nClockJumps = nClockJumps;
        
        # analysisResults.ExtraOutputInfo.meanClockJumpInterval = meanClockJumpInterval;
        # analysisResults.ExtraOutputInfo.stdClockJumpInterval = stdClockJumpInterval;
        
        nClockJumps, meanClockJumpInterval, stdClockJumpInterval = detectClockJumps(GNSS_obs, nGNSSsystems, obsCodes, time_epochs, tInterval,GNSSsystems)
        # analysisResults.ExtraOutputInfo.nClockJumps = nClockJumps;
        analysisResults['ExtraOutputInfo']['nClockJumps'] = nClockJumps
        analysisResults['ExtraOutputInfo']['meanClockJumpInterval'] = meanClockJumpInterval
        analysisResults['ExtraOutputInfo']['stdClockJumpInterval']  = stdClockJumpInterval
        
    if 'sat_pos' in locals(): # add satellite position,azimut, elevation to analysResults
        analysisResults['Sat_position'] = sat_pos
    

        
    print('\n\nINFO: Analysis complete!')
    ## -- Create output file    
    baseFileName = os.path.basename(rinObsFilename)
    outputFilename = baseFileName.split('.')[0] +   '_Report.txt'
    writeOutputFile(outputFilename, outputDir, analysisResults, includeResultSummary, includeCompactSummary, includeObservationOverview, includeLLIOverview)
    print('INFO: The output file %s has been written.' % (outputFilename))
    ## -- Make barplot if plotEstimates is True 
    if plotEstimates:
        make_barplot(analysisResults,graphDir)
    ## -- Saving the workspace as a binary pickle file ---
    results_name = os.path.join(outputDir, 'analysisResults.pkl')
    f = open(results_name,"wb")
    ## -- write the python object (dict) to pickle file
    pickle.dump(analysisResults,f)
    print('INFO: The analysis results has been written to the file %s.' % (results_name))
    f.close()


    return analysisResults



def ismember(list_,code):
    """
    The function takes in a string and a list, and finds the index of 
    """
    indx = [idx for idx, val in enumerate(list_) if val == code]
    if indx != []:
        indx = indx[0]
    return indx





