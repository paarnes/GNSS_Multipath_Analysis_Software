import os, numpy as np, pickle
from readRinexObs import *
from Geodetic_functions import *
from computeSatElevations import computeSatElevations
from computeSatElevAzimuth_fromNav import computeSatElevAzimuth_fromNav
from readFrequencyOverview import readFrequencyOverview
from signalAnalysis import signalAnalysis
from detectClockJumps import detectClockJumps
from tqdm import tqdm, trange
from writeOutputFile import writeOutputFile
from make_polarplot import make_polarplot,make_skyplot, make_polarplot_SNR, plot_SNR_wrt_elev
from make_polarplot_dont_use_TEX import make_polarplot_dont_use_TEX, make_skyplot_dont_use_TEX, make_polarplot_SNR_dont_use_TEX, plot_SNR_wrt_elev_dont_use_TEX
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
                          include_SNR=None,
                          tLim_R   = None,
                          tLim_GEC = None,
                          includeResultSummary= None,
                          includeCompactSummary=None,
                          includeObservationOverview=None,
                          includeLLIOverview= None,
                          use_LaTex =None
                          ):
    
    """
    GNSS Multipath Analysis
    ------------------------
    Made by: Per Helge Aarnes
    E-mail: per.helge.aarnes@gmail.com
    
    Based on GNSS_Reciever_QC_2020 software made by Bjørn-Eirik Roald.
    
    The main function of the software that, through the help of other functions: 
        
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
    
      - readRinexObs304.py
      - computeSatElevations.py
      - readFrequencyOverview.py
      - signalAnalysis.py
      - plotEstimates.py
      - make_barplot.py
      - make_polarplot.py
      - writeOutputFile.py
    
    --------------------------------------------------------------------------------------------------------------------------
    INPUTS:
    
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
    
    
    includeResultSummary:     boolean. 1 if user desires output file to
                              include more detailed overview of statistics, 
                              including for each individual satellites. 
                              0 otherwise (optional)
    
    includeCompactSummary:    boolean. 1 if user desired output file to
                              include more compact overview og statistics. (optional)
    
    includeObservationOverview:     boolean. 1 if user desires output file to
                                      include overview of obseration types observed
                                      by each satellite. 0 otherwise (optional)
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
    if broadcastNav2 == None:    
        broadcastNav2 = "" 
    if broadcastNav3 == None:    
        broadcastNav3 = "" 
    if broadcastNav4 == None:    
        broadcastNav4 = ""         
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
        
    if tLim_R == None:
        tLim_R = 1800 # 30 min
    if tLim_GEC == None:
        tLim_GEC = 7200 # 2 hours 

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
    
    if use_LaTex == None:
        use_LaTex = True
    
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
    includeAllObsCodes  = 0
    
    if include_SNR == None:
        desiredObsCodes = ["C", "L"] # only code and phase observations
    elif include_SNR == True:
        desiredObsCodes = ["C", "L", "S"]
    # desiredObsCodes = ["C", "L"] # only code and phase observations
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
        sat_elevation_angles, sat_azimut_angles, sat_coordinates = computeSatElevations(GNSS_SVs, GNSSsystems, approxPosition,\
            nepochs, time_epochs, max_sat, sp3NavFilename_1, sp3NavFilename_2, sp3NavFilename_3)
    else:
        nav_files = [broadcastNav1,broadcastNav2,broadcastNav3,broadcastNav4]
        sat_pos = computeSatElevAzimuth_fromNav(nav_files,approxPosition,GNSS_SVs,GNSS_obs,time_epochs,tLim_GEC,tLim_R)
        
        ## -- Build same struture for satellit elevation angles if broadcast nav defined
        sat_elevation_angles = {}
        sat_pos_dummy = sat_pos.copy()
        for sys in np.arange(0,len(GNSSsystems)):
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
           raise ValueError("ERROR! GLONASS k-numbers do not exist. This is mandatory to be able to run analysis for GLONASS. Please add GLONASS SLOT / FRQ  to RINEX header.")
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
    
    ## --initialize variable storing total number of observation codes processed
    nCodes_Total = 0
    
    ## -- initialize results dictionary
    analysisResults = {}
    analysisResults['nGNSSsystem'] = nGNSSsystems
    analysisResults['GNSSsystems'] = list(GNSSsystems.values())
    # for sys in range(0,nGNSSsystems):
    for sys in np.arange(0,nGNSSsystems):
        ## -- Full name of current GNSS system
        GNSSsystemName = GNSSsystemCode2Fullname[GNSSsystems[sys+1]]
        ## -- Include full name of current GNSS system 
        analysisResults[GNSSsystemName] =  {}
        ## -- Initialize dict for current GNSS system
        current_sys_dict = {}
        ## -- Initialize observationOverview field as a dict
        current_sys_dict['observationOverview'] = {}
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
                ## --Increment number of bands for current system dict
                current_sys_dict['nBands'] = current_sys_dict['nBands'] + 1
                ## -- Append current band to list of band for this system dict
                Bands_list.append("Band_" + str(bandNumInd))
                current_sys_dict['Bands'] =Bands_list
    
                ## -- Create field for this band, field dict
                current_band_dict = {}
               
                ## -- Store number of codes in current band
                current_band_dict['nCodes'] = nCodes_currentBand
    
                ## -- Increment total number of codes processed
                nCodes_Total = nCodes_Total + nCodes_currentBand
               
                ## -- Store codes for this band
                # current_band_dict.Codes = cellstr(obsCodeOverview{sys}{bandNumInd});
                curr_band = current_sys_dict['Bands'][current_sys_dict['nBands']-1]
                current_band_dict['Codes'] = obsCodeOverview[sys+1][bandNumInd]
    
                ## -- Store current band dict as field in current system dict 
                current_sys_dict[curr_band] = current_band_dict # subtrackting 1 to get index
    
            
        
        # Store current systems dict as field in results dict
        analysisResults[GNSSsystemName] = current_sys_dict 
    
    #### -------- Execute analysis of current data, produce results and store results in results dictionary ------
    
    ## -- Intialize wait bar
    # waitbar = tqdm(0, 'INFO(GNSS\\_Receiver\\_QC\\_2020): Data analysis is being executed. Please wait.')
    ## --Initialize counter of number of codes processed so far. This is used  mainly for waitbar
    codeNum = 0
    
    ## -- Defining frrmat of progressbar
    bar_format = '{desc}: {percentage:3.0f}%|{bar}| ({n_fmt}/{total_fmt})'
    for sys in np.arange(0,nGNSSsystems):    # replaced "range" with np.arange for speed      
        ## --Get current GNSS system code, example GPS: G
        currentGNSSsystem = GNSSsystems[sys+1]
        GNSSsystemName = GNSSsystemCode2Fullname[GNSSsystems[sys+1]] # I USE THIS INSTEAD                
        current_sys_dict = analysisResults[GNSSsystemName]
        ## -- Get number of carrier bands in current system dict
        nBands = current_sys_dict['nBands']
        ## -- Itterate through Bands in system dict. 
        ## NOTE variable "bandNumInd" is NOT the carrier band number, but the index of that band in this system dict  
    
        for bandNumInd in trange(0,nBands,initial=0, desc='Currently processing all available bands for %s' % (GNSSsystemName), leave=False,bar_format=bar_format,position=0): 
            ## Make HARD copy of current band dict
            current_band_dict = current_sys_dict[current_sys_dict['Bands'][bandNumInd]]
            
            ## Get number of codes for current band dict
            nCodes = current_band_dict['nCodes']
            # Get current band full name
            currentBandName = current_sys_dict['Bands'][bandNumInd]
            
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
                   
                    ## Itterate through the codes in the other bands to execute
                    ## analysis with them and current range1 code
                    for secondBandnum in np.arange(0,nBands):  # replaced "range" with np.arange for speed
                        # Disregard observation code in same carrier band as current range1 observation
                        if secondBandnum != bandNumInd:
                            ## Make HARD copy of the other band dict
                            other_band_dict = current_sys_dict[current_sys_dict['Bands'][secondBandnum]]
                            ## Get number of codes in other band dict
                            nCodesOtherBand = other_band_dict['nCodes'] 
                         
                            # Itterate through codes in other band
                            for k in np.arange(0,nCodesOtherBand):                                
                                  ## Get code(range) and phase obsertion codes from the other band
                                  range2_Code = other_band_dict['Codes'][k]
                                  phase2_Code = "L" + range2_Code[1::]
                                  ## Check if phase2 observation was read from RINEX 3 observtaion file
                                  if phase2_Code in obsCodes[sys+1][currentGNSSsystem]:
                                       ## Execute the analysis of current combination of observations. Return statistics on analysis
                                       currentStats, success = signalAnalysis(currentGNSSsystem, range1_Code, range2_Code, GNSSsystems, frequencyOverview, nepochs, tInterval, \
                                       int(max_sat[sys]), GNSS_SVs[currentGNSSsystem], obsCodes[sys+1], GNSS_obs[currentGNSSsystem], GNSS_LLI[currentGNSSsystem],\
                                           sat_elevation_angles[sys], phaseCodeLimit, ionLimit, cutoff_elevation_angle)
                
                                       if not success:
                                            return success
                                       
                                       ##  -- Get number of estimates produced from analysis
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
                                        ## -- Remove range1 observation dict from other band dict, as it can not be used later
                                        other_band_dict['Codes'][ismember(other_band_dict['Codes'], range2_Code)] = []
                                        ## -- Deincrement numbe rof codes in otehr band dict
                                        other_band_dict['nCodes'] = other_band_dict['nCodes'] - 1 
                                        ## -- replace the, now altered, hard copy of other band dict in its original place in system dict
                                        current_sys_dict[current_sys_dict['Bands'][secondBandnum]] = other_band_dict
        
                    ## -- Store best analysis result dict in current band dict
                    current_code_dict = best_currentStats                  
                    ## For every satellite that had an observation of range1, it
                    ## is stored in an overview. Hence the user can get overview
                    ## of which satellites have transmitted which observations
                    ## number of satellites for current system, observation or no
                    if type(current_code_dict) == dict:
                        nSat = len(current_code_dict['range1_slip_distribution_per_sat'])
                    else:
                        print('\n\nWARNING! No estimates for band: "%s" and system: "%s". Probably because of only one obscode available for the current band. Observations from two different frequency is needed.' % (currentBandName,GNSSsystemName))
                        break
     
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
    
                  
                   ## -- If plotEstimates boolean is 1, plot estimates from best analysis and store figures
                    if plotEstimates:
                      ## If user has not specified an output directory, set output directory to "Output_Files" 
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
                      if use_LaTex: #check if use TEX. Adding try/except to handle if user dont have TEX installed
                          try:
                              plotResults(current_code_dict['ion_delay_phase1'], current_code_dict['multipath_range1'], \
                                  current_code_dict['sat_elevation_angles'], tInterval, currentGNSSsystem, \
                                  current_code_dict['range1_Code'], current_code_dict['range2_Code'], \
                                  current_code_dict['phase1_Code'], current_code_dict['phase2_Code'], graphDir)
                          except:
                              print("TEX not installed on your computer! Install that to get prettier text formatting in plots.")
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
                  
                else:
                    ## If phase1 observation is not read from RINEX observation file
                    print('\nINFO(GNSS_MultipathAnalysis): %s code exists in RINEX observation file, but not %s\n'\
                                    'Linear combination using this signal is not used.\n\n' % (range1_Code, phase1_Code))


                    current_band_dict['Codes'][ismember(current_band_dict['Codes'], range1_Code)] = []
                    current_band_dict['nCodes'] = current_band_dict['nCodes'] - 1 
                        
                         
            ## -- Replace the, now altered, hard copy of current band dict in its original place in system dict
            current_sys_dict[current_sys_dict['Bands'][bandNumInd]] = current_band_dict
            if sp3NavFilename_1 != '':
                try:
                    sat_pos[currentGNSSsystem] = {}
                    sat_pos[currentGNSSsystem]['Position']  = sat_coordinates[currentGNSSsystem]                    
                    sat_pos[currentGNSSsystem]['Azimut']    = sat_azimut_angles[sys]
                    sat_pos[currentGNSSsystem]['Elevation'] = current_code_dict['sat_elevation_angles']
                except:
                    pass
        
        
        ## -- Replace the, now altered, hard copy of current system dict in its original place in results dict
        analysisResults[GNSSsystemName] = current_sys_dict
    
        ## -- Store information needed for output file in result dict
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
        analysisResults['ExtraOutputInfo']['elevation_cutoff']    = cutoff_elevation_angle
        
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
        
    if 'sat_pos' in locals(): # add satellite position,azimut, elevation to analysResults
        analysisResults['Sat_position'] = sat_pos
    

        
    print('\n\nINFO: Analysis complete!\n')
    ## -- Create output file    
    if outputDir == "":
        outputDir = 'Output_Files'
       
    ## -- Unless output directory already exists, create output directory
    if not os.path.isdir(outputDir):
        os.mkdir(outputDir)
    
   
    ## --Unless graph directory already exists, create directory
    graphDir = outputDir + '/Graphs' 
    if not os.path.isdir(graphDir):
        os.mkdir(graphDir)
        
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
                azimut_currentSys = analysisResults['Sat_position'][curr_sys]['Azimut']
                elevation_currentSys = analysisResults['Sat_position'][curr_sys]['Elevation'] 
                print('INFO: Making a regular polar plot for showing azimut and elevation angle for each satellite. Please wait...\n')
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
            print('INFO: Making a polar plot of the multipath effect. Please wait ...\n')
            if use_LaTex:
                try:
                    make_polarplot(analysisResults, graphDir)
                except:
                    make_polarplot_dont_use_TEX(analysisResults, graphDir)  
            else:
                make_polarplot_dont_use_TEX(analysisResults, graphDir)  
            
        if include_SNR:
            print('INFO: Making a plot of the Signal To Noise Ration (SNR). Please wait ...\n')
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
                
                    

            
            
            
    ## -- Saving the workspace as a binary pickle file ---
    results_name = os.path.join(outputDir, 'analysisResults.pkl')
    f = open(results_name,"wb")
    ## -- write the python object (dict) to pickle file
    pickle.dump(analysisResults,f)
    print('INFO: The analysis results has been written to the file %s.\n' % ('analysisResults.pkl'))
    print('INFO: Finished!')
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





