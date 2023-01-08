import os, sys, pickle,numpy as np
from readRinexObs304 import readRinexObs304
from GNSS_MultipathAnalysis import GNSS_MultipathAnalysis

abs_path = os.path.abspath("GNSS_Receiver_QC_2020.py")
base_path = os.path.normpath(os.getcwd() + os.sep + os.pardir)
## Path to TestData
relpath_to_testdata = 'TestData'
full_path_testdata = os.path.join(base_path, relpath_to_testdata) 

#Path to outputdir 
relpath_to_outputdir = 'Results'
full_path_ouputdir = os.path.join(base_path, relpath_to_outputdir) 
## -----  Defining input data --------

## Rinex observation file
rinObsFilename = full_path_testdata  + '/ObservationFiles/' + 'OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx'
# rinObsFilename = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Multipath_analysis/NMBUS_SAMSUNG_S20.20o'


## SP3 files
sp3NavFilename_1 = full_path_testdata  + '/SP3/' + 'test1.eph'
sp3NavFilename_2 = full_path_testdata  + '/SP3/' + 'test2.SP3'
sp3NavFilename_3 = full_path_testdata  + '/SP3/' + 'test3.SP3'


broadcastNav1 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_GN.rnx'
broadcastNav2 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_RN.rnx'
broadcastNav3 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_EN.rnx'
broadcastNav4 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_CN.rnx'

broadcastNav1 =  full_path_testdata  + '/NavigationFiles/' + 'BRDC00IGS_R_20220010000_01D_MN.rnx'

## Broadcast ephemerides
# broadcastNav1 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_GN.rnx'
# broadcastNav2 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_RN.rnx'
# broadcastNav3 =  full_path_testdata  + '/NavigationFiles/' + 'ENDRET_BCEmerge_30_10_2020.20p'
# broadcastNav3 =  full_path_testdata  + '/NavigationFiles/' + 'BCEmerge_30_10_2020.20p'


## -- Simple example for running analysis (no userdefined settings)
# analysisResults = GNSS_MultipathAnalysis(rinObsFilename, 
#                                          broadcastNav1=broadcastNav1,
#                                          broadcastNav2=broadcastNav3
#                                          )


analysisResults = GNSS_MultipathAnalysis(rinObsFilename, 
                                          broadcastNav1=broadcastNav3,
                                          desiredGNSSsystems = ['E']
                                          )

# analysisResults = GNSS_MultipathAnalysis(rinObsFilename, 
#                                          broadcastNav1=broadcastNav1,
#                                          broadcastNav2=broadcastNav2,
#                                          broadcastNav3=broadcastNav3,
#                                          broadcastNav4=broadcastNav4
#                                          )


# from readRinexNav import *
# data, header, n_eph = read_rinex3_nav(broadcastNav1)
#%% -- Advanced example (more user defined settings)

## Parameters
phaseCodeLimit              = 0
ionLimit                    = 0
cutoff_elevation_angle      = 0
outputDir                   = full_path_ouputdir
plotEstimates               = True
includeResultSummary        = True
includeCompactSummary       = True
includeObservationOverview  = True
includeLLIOverview          = True


GNSSsystems = ["R"] # run analysis in GLONASS only
analysisResults = GNSS_MultipathAnalysis(rinObsFilename,
                       desiredGNSSsystems=GNSSsystems,
                       sp3NavFilename_1 = sp3NavFilename_1,
                       sp3NavFilename_2 = sp3NavFilename_2,
                       sp3NavFilename_3 = sp3NavFilename_3,
                       phaseCodeLimit = phaseCodeLimit,
                       ionLimit = ionLimit,
                       cutoff_elevation_angle = cutoff_elevation_angle,
                       outputDir = outputDir, 
                       plotEstimates = plotEstimates , 
                       includeResultSummary = includeResultSummary , 
                       includeCompactSummary = includeCompactSummary , 
                       includeObservationOverview = includeObservationOverview ,
                       includeLLIOverview = includeLLIOverview
                       )
 
  
## --  Read in the result file from a analysis
# result_path = os.path.join(full_path_ouputdir, "analysisResults.pkl")
file_to_read = open(r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo/analysisResults.pkl', "rb")
loaded_dictionary = pickle.load(file_to_read)


#%% Read multiple nav
import os, sys, pickle,numpy as np
abs_path = os.path.abspath("GNSS_Receiver_QC_2020.py")
base_path = os.path.normpath(os.getcwd() + os.sep + os.pardir)
## Path to TestData
relpath_to_testdata = 'TestData'
full_path_testdata = os.path.join(base_path, relpath_to_testdata) 

from readRinexNav import *

broadcastNav1 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_GN.rnx'
broadcastNav2 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_RN.rnx'
broadcastNav3 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_EN.rnx'
broadcastNav4 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_CN.rnx'


data, header, n_eph = read_rinex3_nav(broadcastNav2)
nav_list = [broadcastNav1,broadcastNav2,None,None]
nav_list = [i for i in nav_list if i is not None]
sat_pos2 = {}
for file in nav_list:
    data, header, n_eph =read_rinex3_nav(file)
    curr_sys = data[0][0][0]
    sat_pos2[curr_sys] = data 



# phaseCodeLimit              = 0;
# ionLimit                    = 0;
# cutoff_elevation_angle      = 0; 
# # outputDIR                   = r'C:\Users\perhe\Desktop\Data\Test2/'
# outputDir                   = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Multipath_analysis/Test2'
# plotEstimates               = 1;
# saveWorkspace               = 1;
# includeResultSummary        = 1;
# includeCompactSummary       = 1;
# includeObservationOverview  = 1;
# includeLLIOverview          = 1;
# readLLI=1
# readSS = 1
# includeAllGNSSsystems   = 0
# includeAllObsCodes      = 0
# desiredGNSSsystems = ["G"]
# desiredObsCodes = ["C", "L"];               # code and phase observations
# desiredObsBands = list(np.arange(1,10))     # all carrier bands. Tot 9, but arange stops at 8 -> 10

# ## --- Read RINEX 3.0x observation file
# [GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
#     obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
#     rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success] = \
#     readRinexObs304(rinObsFilename, readSS, readLLI, includeAllGNSSsystems,includeAllObsCodes, desiredGNSSsystems,\
#     desiredObsCodes, desiredObsBands)
