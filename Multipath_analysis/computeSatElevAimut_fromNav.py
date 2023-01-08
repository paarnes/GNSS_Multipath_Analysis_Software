import sys,numpy as np
from readRinexNav import read_rinex3_nav
from read_SP3Nav import readSP3Nav
from Geodetic_functions import *
from readRinexObs304 import *
import pandas as pd

def computeSatElevAimut_fromNav(navigationFile,approxPosition,GNSS_SVs,GNSS_obs,time_epochs):
    """
    A function for computing satellite elevations and azimut angles based on
    broadcasted ephemerides. Support all global navigation systems.
    
    Input:
        navigationFile: list of navigation files
    """
        
    ##--- Read rinex nav
    # data_,header, n_eph = read_rinex3_nav(navigationFile,dataframe='no')
    nav_list = [i for i in navigationFile if i is not None] #remove "NONE" if exist in list

    for idx, nav_file in enumerate(nav_list):
        data_,header, n_eph = read_rinex3_nav(nav_file,dataframe='no')
        if idx == 0:
            data = data_
        else:
            data = np.append(data,data_,axis=0)
    ### ------- Compute satellite coordiantes in ECEF and compute azimut and elevation angels
    
    ## -- Extracting approx postion from RINEX obs-file
    x = float(approxPosition[0])
    y = float(approxPosition[1])
    z = float(approxPosition[2])
    
    GNSS_FullName = dict(list(zip(['G','R','E','C'],['GPS','GLONASS','Galileo','BeiDou'])))
    
    ## -- Find availible satellittes for the whole RINEX file
    sat_pos = {}  # Dict for storing all data
    for curr_sys in GNSS_SVs:
        sat_pos[curr_sys] = {}
        aktuelle_sat_list = []  # Dummy list for availible satellites
        for epoch in GNSS_obs[curr_sys]:
            aktuelle_sat = GNSS_SVs[curr_sys][epoch-1].astype(int)
            aktuelle_sat_list.append(aktuelle_sat[aktuelle_sat != 0])
        aktuelle_sat_list_dum = [list(x) for x in set(tuple(x) for x in aktuelle_sat_list)] #get list of sets
        sat_pos[curr_sys]['Available_Sat'] = list(sorted(set(x for l in aktuelle_sat_list_dum for x in l)))
        
    ## -- Compute satellite coordinates for all availible satellites   
    t = time_epochs[:,1]  # Extracting time for RINEX obs-file
    for sys in list(sat_pos.keys()):
        aktuelle_sat_list = sat_pos[sys]['Available_Sat']
        df_data = pd.DataFrame(data) # making dataframe of data
        curr_data = df_data[df_data.iloc[:,0].str.contains(sys)].to_numpy() # extaction data for current system only
        counter = 0
        curr_pos = {}       # Cells for storing data
        curr_azimut = {}
        curr_elevation = {}
        X = np.zeros([len(t),61]) # Array for storing posistions and angles
        Y = np.zeros([len(t),61])
        Z = np.zeros([len(t),61])
        azimut    = np.zeros([len(t),61]) # Satellites azimut
        elevation = np.zeros([len(t),61]) # Satellites elevation 
        for PRN in aktuelle_sat_list:
            counter = counter + 1
            print("\rCurrently computing coordinates for the %s system. Progress: %.1f%%" %(GNSS_FullName[sys],counter/len(aktuelle_sat_list)*100), end='\r',flush=True)  # \r makes the line get overwritten

            # print("\rCurrently computing coordinates for the %s system. Progress: %.1f%%" %(GNSS_FullName[sys],counter/len(aktuelle_sat_list)*100), end='\r',flush=True)  # \r makes the line get overwritten
            # try:
            #     ephemerides = extract_nav_message(curr_data,PRN,t[epoch-1]) # passing curr_data instead to get correct system
            #     ephemerides[0] = ephemerides[0][1::] # Removing system letter from number. Ex G10 -> 10
            #     ephemerides = ephemerides.astype(float)
            # except:
            #     ephemerides = np.nan
            #     continue
            ## Computing satellite coordinates
            # for i in range(0,len(t)):
            for i in np.arange(0,len(t)):
                curr_time = t[i]
                try:
                    ephemerides = extract_nav_message(curr_data,PRN,curr_time) # passing curr_data instead to get correct system
                    ephemerides[0] = ephemerides[0][1::] # Removing system letter from number. Ex G10 -> 10
                    ephemerides = ephemerides.astype(float)
                except:
                    ephemerides = np.nan
                    continue
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
            curr_pos[str(PRN)] = np.array([X[:,PRN],Y[:,PRN],Z[:,PRN]]).T
    
        ## -- Update dictionary with coordinates    
        sat_pos[sys]['Position']  = curr_pos
        sat_pos[sys]['Azimut']    = azimut
        sat_pos[sys]['Elevation'] = elevation
        
    return sat_pos
    
            
