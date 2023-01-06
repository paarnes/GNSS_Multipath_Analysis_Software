import numpy as np
from rinexFindNEpochs304 import rinexFindNEpochs304
from rinexReadObsFileHeader304 import rinexReadObsFileHeader304
from tqdm import tqdm
from rinexReadObsBlock304 import rinexReadObsBlock304
from kepler2ecef import date2gpstime
from rinexReadObsBlockHead304 import rinexReadObsBlockHead304
from readRinexObs304 import readRinexObs304

readLLI = 1
readSS = 1
includeAllGNSSsystems = 0
includeAllObsTypes = 0
includeAllObsCodes = 0
desiredGNSSsystems = ["G","R","E","C"]
desiredObsCodes = ["C","L"]
desiredObsBands = [1,2, 5]
# desiredObsTypes = ["C", "L", "S", "D"]


# filename = 'opec0020_3.04_kort.10o'

# filename = 'OPEC00NOR_S_20220010000_01D_30S_MO.rnx'
filename = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Multipath_analysis/OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx'

GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
    obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
    rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success =\
    readRinexObs304(filename, readSS, readLLI, includeAllGNSSsystems, \
                includeAllObsCodes,desiredGNSSsystems, desiredObsCodes, desiredObsBands)



C1C = []
systems = list(GNSS_obs.keys())
for sys in systems:
    for ep in range(1,len(GNSS_obs[sys])+1):
        C1C.append(GNSS_obs[sys][ep][:,0].reshape(1,len(GNSS_obs[sys][1][:,0])))
    
