import os, pickle,numpy as np
from GNSS_MultipathAnalysis import GNSS_MultipathAnalysis
from readRinexObs import readRinexObs

abs_path = os.path.abspath("GNSS_MultipathAnalysis.py") #find path to the main file
base_path = os.path.normpath(os.getcwd() + os.sep + os.pardir)
## Path to TestData
relpath_to_testdata = 'TestData'
full_path_testdata = os.path.join(base_path, relpath_to_testdata) 

#Path to outputdir 
relpath_to_outputdir = 'Results'
full_path_ouputdir = os.path.join(base_path, relpath_to_outputdir) 

## -----  Defining input data --------
## Rinex observation file
rinObsFilename1 = full_path_testdata  + '/ObservationFiles/' + 'OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx'
rinObsFilename2 = full_path_testdata  + '/ObservationFiles/' + 'NMBUS_SAMSUNG_S20.20o'

## SP3 files
sp3NavFilename_1      = full_path_testdata  + '/SP3/' + 'NMBUS_2020 10 30.SP3'
sp3NavFilename_1_opec = full_path_testdata  + '/SP3/' + 'Testfile_20220101.eph'


## Broadcast ephemerides
broadcastNav1 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_GN.rnx'
broadcastNav2 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_RN.rnx'
broadcastNav3 =  full_path_testdata  + '/NavigationFiles/' + 'BCEmerge_30_10_2020.20p'
broadcastNav4 =  full_path_testdata  + '/NavigationFiles/' + 'BRDC00IGS_R_20220010000_01D_MN.rnx'



#%% #####-- Simple examples for running analysis (few userdefined settings) 
# Using SP3  OPEC

analysisResults = GNSS_MultipathAnalysis(rinObsFilename1,
                                         desiredGNSSsystems=["G"],
                                         sp3NavFilename_1 = sp3NavFilename_1_opec)


#%% Using broadcasted eph OPEC
analysisResults = GNSS_MultipathAnalysis(rinObsFilename1,
                                         desiredGNSSsystems=['G'], # Run analysis on GPS only ['G','R','E'] for GPS,GLONASS and Galileo 
                                         outputDir = full_path_ouputdir,
                                         broadcastNav1=broadcastNav4
                                         )

#%% Using broadcasted eph NMBUS
analysisResults = GNSS_MultipathAnalysis(rinObsFilename2,
                                          desiredGNSSsystems=['E'],
                                          outputDir = full_path_ouputdir,
                                          broadcastNav1=broadcastNav3
                                          )
#%% Using SP3  NMBUS
analysisResults = GNSS_MultipathAnalysis(rinObsFilename2,
                                         desiredGNSSsystems=['G'],
                                         outputDir = full_path_ouputdir,
                                         sp3NavFilename_1=sp3NavFilename_1,
                                         plotEstimates=True
                                         )

#%% Using SP3  NMBUS without plotting
analysisResults = GNSS_MultipathAnalysis(rinObsFilename2,
                                         desiredGNSSsystems=['G'],
                                         outputDir = full_path_ouputdir,
                                         broadcastNav1=broadcastNav3,
                                         plotEstimates=False,
                                         plot_polarplot=False
                                         )

#%% -- Advanced example (more user defined settings)

## Parameters
GNSSsystems                 = ["R"] # run analysis in GLONASS only
phaseCodeLimit              = 6.667
ionLimit                    = 0.0667
cutoff_elevation_angle      = 10 # 10 degree elevation cutoff
outputDir                   = full_path_ouputdir
plotEstimates               = True
plot_polarplot              = True
includeResultSummary        = True
includeCompactSummary       = True
includeObservationOverview  = True
includeLLIOverview          = True



analysisResults = GNSS_MultipathAnalysis(rinObsFilename1,
                       desiredGNSSsystems=GNSSsystems,
                       sp3NavFilename_1 = sp3NavFilename_1_opec,
                       phaseCodeLimit = phaseCodeLimit,
                       ionLimit = ionLimit,
                       cutoff_elevation_angle = cutoff_elevation_angle,
                       outputDir = outputDir, 
                       plotEstimates = plotEstimates,
                       plot_polarplot=plot_polarplot,
                       includeResultSummary = includeResultSummary, 
                       includeCompactSummary = includeCompactSummary, 
                       includeObservationOverview = includeObservationOverview,
                       includeLLIOverview = includeLLIOverview
                       )
 
  

#%% ------------ How to read in the result file from a analysis
path_to_resFile = os.path.join(full_path_ouputdir, 'analysisResults.pkl')
file_to_read = open(path_to_resFile, "rb")
loaded_dictionary = pickle.load(file_to_read)


#%% ----------- Examplecode for using readrinex to read in GPS and Galielo observation
rinObsFilename = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx'
desiredGNSSsystems = ["G","E"]          # GPS and Galileo
desiredObsCodes = ["C", "L"]       # code and phase observations
desiredObsBands = list(np.arange(1,10)) # all carrier bands. 

## --- Read RINEX 3.0x observation file
[GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
    obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
    rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success] = \
    readRinexObs(rinObsFilename,  
                 desiredGNSSsystems = desiredGNSSsystems,
                 desiredObsCodes = desiredObsCodes,
                 desiredObsBands = desiredObsBands)
