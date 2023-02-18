import os, sys, pickle,numpy as np
from readRinexObs import readRinexObs304
from GNSS_MultipathAnalysis import GNSS_MultipathAnalysis

abs_path = os.path.abspath("GNSS_MultipathAnalysis.py")
base_path = os.path.normpath(os.getcwd() + os.sep + os.pardir)
## Path to TestData
relpath_to_testdata = 'TestData'
full_path_testdata = os.path.join(base_path, relpath_to_testdata) 

#Path to outputdir 
relpath_to_outputdir = 'Results_OPEC_eph'
full_path_ouputdir = os.path.join(base_path, relpath_to_outputdir) 

## -----  Defining input data --------
## Rinex observation file
rinObsFilename1 = full_path_testdata  + '/ObservationFiles/' + 'OPEC00NOR_S_20220010000_01D_30S_MO.rnx' # with SNR
rinObsFilename2 = full_path_testdata  + '/ObservationFiles/' + 'OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx'
rinObsFilename3 = full_path_testdata  + '/ObservationFiles/' + 'NMBUS_SAMSUNG_S20.20o'

## SP3 files
sp3NavFilename_1      = full_path_testdata  + '/SP3/' + 'NMBUS_2020 10 30.SP3'
sp3NavFilename_1_opec = full_path_testdata  + '/SP3/' + 'test1.eph'
sp3NavFilename_2_opec = full_path_testdata  + '/SP3/' + 'test2.SP3'
sp3NavFilename_3_opec = full_path_testdata  + '/SP3/' + 'test3.SP3'

## Broadcast ephemerides
broadcastNav1 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_GN.rnx'
broadcastNav2 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_RN.rnx'
broadcastNav3 =  full_path_testdata  + '/NavigationFiles/' + 'ENDRET_BCEmerge_30_10_2020.20p'
broadcastNav4 =  full_path_testdata  + '/NavigationFiles/' + 'BCEmerge_30_10_2020.20p'
broadcastNav5 =  full_path_testdata  + '/NavigationFiles/' + 'BRDC00IGS_R_20220010000_01D_MN.rnx'
#%% -- Simple example for running analysis (no userdefined settings) using SP3 files
analysisResults = GNSS_MultipathAnalysis(rinObsFilename1,
                                          desiredGNSSsystems=['G'],
                                          outputDir = full_path_ouputdir,
                                          sp3NavFilename_1=sp3NavFilename_2_opec
                                          )
#%% Using broadcasted eph NMBUS
analysisResults = GNSS_MultipathAnalysis(rinObsFilename3,
                                           desiredGNSSsystems=['G'],
                                          outputDir = full_path_ouputdir,
                                          broadcastNav1=broadcastNav4
                                          )
#%% Using SP3  NMBUS
analysisResults = GNSS_MultipathAnalysis(rinObsFilename3,
                                           desiredGNSSsystems=['G'],
                                          outputDir = full_path_ouputdir,
                                          sp3NavFilename_1=sp3NavFilename_1
                                          )
#%% Using SP3  OPEC
analysisResults = GNSS_MultipathAnalysis(rinObsFilename2,
                                          # desiredGNSSsystems=['G'],
                                          outputDir = full_path_ouputdir,
                                          sp3NavFilename_1=sp3NavFilename_1_opec
                                          )
#%% Using broadcasted eph OPEC
analysisResults = GNSS_MultipathAnalysis(rinObsFilename2,
                                            desiredGNSSsystems=['G'],
                                          outputDir = full_path_ouputdir,
                                          broadcastNav1=broadcastNav5
                                          )

#%%
sp3 =  r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData'  + '/SP3/' + 'COD0R03FIN_20100010000_01D_05M_ORB.SP3'
rinObsFilename = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\opec0010_3.10o'
analysisResults = GNSS_MultipathAnalysis(rinObsFilename,
                                         desiredGNSSsystems=['G'],
                                          outputDir = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo/test211',
                                          sp3NavFilename_1=sp3
                                          )
#%% -- Advanced example (more user defined settings)

## Parameters
GNSSsystems                 = ["R"] # run analysis in GLONASS only
phaseCodeLimit              = 0
ionLimit                    = 0
cutoff_elevation_angle      = 0
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
                       sp3NavFilename_2 = sp3NavFilename_2_opec,
                       sp3NavFilename_3 = sp3NavFilename_3_opec,
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
 
  

#%% Read in the result file from a analysis
path_to_resFile = os.path.join(full_path_ouputdir, 'analysisResults.pkl')
file_to_read = open(path_to_resFile, "rb")
loaded_dictionary = pickle.load(file_to_read)


#%% Examplecode for using readrinex
import os, sys, pickle,numpy as np
from readRinexObs import readRinexObs

rinObsFilename = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO.rnx'
rinObsFilename = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\opec0010_2.10o'



phaseCodeLimit              = 0
ionLimit                    = 0
cutoff_elevation_angle      = 0 
outputDir                   = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\Results'
plotEstimates               = 1
saveWorkspace               = 1
includeResultSummary        = 1
includeCompactSummary       = 1
includeObservationOverview  = 1
includeLLIOverview          = 1
readLLI=1
readSS = 1
includeAllGNSSsystems   = 0
includeAllObsCodes      = 0
desiredGNSSsystems = ["G"]
desiredObsCodes = ["C", "L", "S"];               # code and phase observations
desiredObsBands = list(np.arange(1,10))     # all carrier bands. Tot 9, but arange stops at 8 -> 10

## --- Read RINEX 3.0x observation file
[GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
    obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
    rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success] = \
    readRinexObs(rinObsFilename, readSS, readLLI, includeAllGNSSsystems,includeAllObsCodes, desiredGNSSsystems,\
    desiredObsCodes, desiredObsBands)
