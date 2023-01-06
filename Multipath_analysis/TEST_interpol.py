import os, sys,numpy as np
sys.path.append(r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Read_RINEX_OBS')
sys.path.append(r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Read_RINEX_nav')
sys.path.append(r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Kepler2ECEF')
from read_rinex3_nav import read_rinex3_nav
from read_SP3Nav import readSP3Nav
from kepler2ecef import Satkoord2
from kepler2ecef import date2gpstime
from kepler2ecef import extract_nav_message
from readRinexObs304 import readRinexObs304
from compute_GLO_coord_from_nav import compute_GLO_coord_from_nav
from scipy.interpolate import lagrange
import matplotlib.pyplot as plt
from tqdm import trange,tqdm
from compute_azimut_elev import compute_azimut_elev
import pandas as pd

# # sp3NavFilename_1 = 'test1.eph'
# # sp3NavFilename_1 = 'OPEC_2018.SP3'
# sp3NavFilename = 'samsung2.SP3'
# # desiredGNSSsystems = ["G", "R", "E", "C"];  # All GNSS systems. 
# sat_positions, epoch_dates, navGNSSsystems, nEpochs, epochInterval, success = readSP3Nav(sp3NavFilename)

## Read RINEX

# rinObsFilename = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Multipath_analysis/OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx'
rinObsFilename = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Multipath_analysis/NMBUS_SAMSUNG_S20.20o'


phaseCodeLimit              = 0;
ionLimit                    = 0;
cutoff_elevation_angle      = 0; 
outputDir                   = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Multipath_analysis/Test2'
plotEstimates               = 1;
saveWorkspace               = 1;
includeResultSummary        = 1;
includeCompactSummary       = 1;
includeObservationOverview  = 1;
includeLLIOverview          = 1;
readSS = 1
readLLI = 1
includeAllGNSSsystems=0
includeAllObsCodes = 0

# desiredGNSSsystems = ["R"]
desiredGNSSsystems = ["G","R","E","C"]
# desiredGNSSsystems = ["E","C"]
desiredObsCodes = ["C", "L"]         
desiredObsBands = [1,2, 5]

[GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
    obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
    rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success] = \
    readRinexObs304(rinObsFilename, readSS, readLLI, includeAllGNSSsystems,includeAllObsCodes, desiredGNSSsystems,\
    desiredObsCodes, desiredObsBands)


#%%
##--- Read nav
# sp3 = 'opec_3.18n'
# sp3 = 'NMBUS_SAMSUNG_S20.20n'
sp3 = 'NMBU_SAMSUNG_S20_ALL.20p'
data,header, n_eph = read_rinex3_nav(sp3,dataframe='no')


#%%
# test_data = [line for indx, line in enumerate(data) if 'G' in data[indx,0]]
# test = np.array([])
# test = []
# for indx, line in enumerate(data):
#     if 'G' in data[indx,0]:
#         # test = np.append(test,data[indx,:].reshape(1,len(data[indx,:])))
#         test.append(data[indx,:].reshape(1,len(data[indx,:]))[0])


# data2 = data[]
#%%




### ------- Compute satellite coordiantes in ECEF and compute azimut and elevation angels

## .-- Extracting approx postion from RINEX obs-file
x = float(approxPosition[0])
y = float(approxPosition[1])
z = float(approxPosition[2])

GNSS_FullName = dict(list(zip(['G','R','E','C'],['GPS','GLONASS','Galileo','BeiDou'])))

## -- Find availible satellittes for the whole RINEX file
sat_pos = {}            # Dict for storing all data
# aktuelle_sat_list = []  # Dummy list for availible satellites
for curr_sys in GNSS_SVs:
    sat_pos[curr_sys] = {}
    aktuelle_sat_list = []  # Dummy list for availible satellites
    for epoch in GNSS_obs[curr_sys]:
        aktuelle_sat = GNSS_SVs[curr_sys][epoch-1].astype(int)
        aktuelle_sat_list.append(aktuelle_sat[aktuelle_sat != 0])
    aktuelle_sat_list_dum = [list(x) for x in set(tuple(x) for x in aktuelle_sat_list)] #get list of sets
    sat_pos[curr_sys]['Available_Sat'] = list(sorted(set(x for l in aktuelle_sat_list_dum for x in l)))
    # sat_pos[curr_sys]['Available_Sat'] =  list([x for x in set(tuple(x) for x in aktuelle_sat_list)][0])     


## -- Compute satellite coordinates for all availible satellites   
t = time_epochs[:,1]      # Extracting time for RINEX obs-file
for sys in list(sat_pos.keys()):
    aktuelle_sat_list = sat_pos[sys]['Available_Sat']
    df_data = pd.DataFrame(data) # making dataframe of data
    curr_data = df_data[df_data.iloc[:,0].str.contains(sys)].to_numpy() # extaction data for current system only
    counter = 0
    curr_pos = {}
    curr_azimut = {}
    curr_elevation = {}
    # t = time_epochs[:,1]      # Extracting time for RINEX obs-file
    X = np.zeros([len(t),61]) # Array for storing posistions and angles
    Y = np.zeros([len(t),61])
    Z = np.zeros([len(t),61])
    azimut    = np.zeros([len(t),61]) # Satellites azimut
    elevation = np.zeros([len(t),61]) # Satellites elevation 
    for PRN in aktuelle_sat_list:
        counter = counter + 1
        print("\rCurrently computing coordinates for the %s system. Progress: %.1f%%" %(GNSS_FullName[curr_sys],counter/len(aktuelle_sat_list)*100), end='\r',flush=True)  # \r makes the line get overwritten
        try:
            # ephemerides = extract_nav_message(data,PRN,t[epoch-1])
            ephemerides = extract_nav_message(curr_data,PRN,t[epoch-1]) # passing curr_data instead to get correct system
            ephemerides[0] = ephemerides[0][1::] # Removing system letter from number. Ex G10 -> 10
            ephemerides = ephemerides.astype(float)
        except:
            ephemerides = np.nan
            continue
        ## Beregner satellittkoordinatene

        for i in range(0,len(t)):
            curr_time = t[i]
            if sys != 'R':
                X[i,PRN], Y[i,PRN], Z[i,PRN],_ = Satkoord2(ephemerides, curr_time, x, y, z)
            else:
                # If current system is GLONASS
                curr_time = time_epochs[i]
                pos, _, _, _= compute_GLO_coord_from_nav(ephemerides, curr_time)
                X[i,PRN] = pos[0]
                Y[i,PRN] = pos[1]
                Z[i,PRN] = pos[2]
            ## - Compute azimut and elevation angle
            azimut[i,PRN],elevation[i,PRN] = compute_azimut_elev(X[i,PRN], Y[i,PRN], Z[i,PRN], x, y, z)
            
        ## -- Assign the computed variable to temporarly dicts for storing results   
        # curr_pos[str(PRN)] = np.array([X[:,PRN],Y[:,PRN],Z[:,PRN]]).reshape(len(X[:,PRN]),3)
        curr_pos[str(PRN)] = np.array([X[:,PRN],Y[:,PRN],Z[:,PRN]]).T

    ## -- Update dictionary with coordinates    
    sat_pos[sys]['Position']  = curr_pos
    sat_pos[sys]['Azimut']    = azimut
    sat_pos[sys]['Elevation'] = elevation
    
            
