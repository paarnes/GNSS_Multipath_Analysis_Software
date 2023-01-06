
import os, sys, pickle,numpy as np

sys.path.append('..Multipath_analysis/')
sys.path.append('..Read_SP3/')
sys.path.append('..\Read_RINEX_OBS/')
sys.path.append('..\Read_RINEX_nav/')
sys.path.append('..\TestData/')
sys.path.append('..\Kepler2ECEF/')
from Read_RINEX_OBS.readRinexObs304 import readRinexObs304
from Multipath_analysis.GNSS_Receiver_QC_2020 import GNSS_Receiver_QC_2020

rinObsFilename = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Multipath_analysis/OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx'
# rinObsFilename = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Multipath_analysis/OPEC00NOR_S_20180010000_01D_01S_MO.rnx'
# rinObsFilename = 'NMBUS_SAMSUNG_S20.20o'
sp3NavFilename_1 = 'test1.eph'
# sp3NavFilename_1 = 'OPEC_2018.SP3'
# sp3NavFilename_1 = 'samsung2.SP3'

# broadcastNav = 'NMBU_SAMSUNG_S20_ALL.20p'
# broadcastNav = 'OPEC00NOR_S_20220010000_01D_GN.rnx'

# sp3NavFilename_2 = 'test2.SP3'
# sp3NavFilename_3 = 'test3.SP3'
sp3NavFilename_2 = ""
sp3NavFilename_3 = ""

## Parameters
# phaseCodeLimit              = 0
# ionLimit                    = 0
# cutoff_elevation_angle      = 0
outputDir                   = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Multipath_analysis/Test1'
# plotEstimates               = True
# includeResultSummary        = True
# includeCompactSummary       = True
# includeObservationOverview  = True
# includeLLIOverview          = True



GNSSsystems = ["G","R","E","C"]

# GNSS_Receiver_QC_2020(rinObsFilename, sp3NavFilename_1, sp3NavFilename_2, sp3NavFilename_3, phaseCodeLimit, ionLimit, cutoff_elevation_angle,\
                          # outputDir, plotEstimates, includeResultSummary, includeCompactSummary,  includeObservationOverview, includeLLIOverview)
    

# analysisResults = GNSS_Receiver_QC_2020(rinObsFilename, broadcastNav1=broadcastNav)
analysisResults = GNSS_Receiver_QC_2020(rinObsFilename, sp3NavFilename_1=sp3NavFilename_1)




# create a binary pickle file 
f = open("file.pkl","wb")
# write the python object (dict) to pickle file
pickle.dump(analysisResults,f)
f.close()

# read in the pickle file
file_to_read = open("file.pkl", "rb")
loaded_dictionary = pickle.load(file_to_read)




phaseCodeLimit              = 0;
ionLimit                    = 0;
cutoff_elevation_angle      = 0; 
# outputDIR                   = r'C:\Users\perhe\Desktop\Data\Test2/'
outputDir                   = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Multipath_analysis/Test2'
plotEstimates               = 1;
saveWorkspace               = 1;
includeResultSummary        = 1;
includeCompactSummary       = 1;
includeObservationOverview  = 1;
includeLLIOverview          = 1;
readLLI=1
readSS = 1
includeAllGNSSsystems   = 0
includeAllObsCodes      = 0
desiredGNSSsystems = ["G"]
desiredObsCodes = ["C", "L"];               # code and phase observations
desiredObsBands = list(np.arange(1,10))     # all carrier bands. Tot 9, but arange stops at 8 -> 10

## --- Read RINEX 3.0x observation file
[GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
    obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
    rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success] = \
    readRinexObs304(rinObsFilename, readSS, readLLI, includeAllGNSSsystems,includeAllObsCodes, desiredGNSSsystems,\
    desiredObsCodes, desiredObsBands)
