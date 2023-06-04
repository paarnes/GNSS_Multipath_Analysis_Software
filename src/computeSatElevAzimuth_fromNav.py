import numpy as np
from readRinexNav import read_rinex3_nav
from Geodetic_functions import *
import pandas as pd, re, os

def computeSatElevAzimuth_fromNav(navigationFile,approxPosition,GNSS_SVs,GNSS_obs,time_epochs, tLim_GEC=None,tLim_R=None):
    """
    A function for computing satellite elevations and azimut angles based on
    broadcasted ephemerides. Support all global navigation systems (GPS,GLONASS,Galileo & BeiDou).
    
    Input:
    -----
    navigationFile: list of navigation files
    """
        
    ##--- Read rinex navigation file/files. If several files defined, the will
    ## be added to the same data array for further analysis
    
    nav_list = [i for i in navigationFile if i !=""] #remove "NONE" if exist in list
    for idx, nav_file in enumerate(nav_list):
        shorten_navigation_file(nav_file) # removes system 'S','J' and 'I' and ephermerids of same TOC to increase speed
        base_path = os.path.split(nav_file)[0]
        org_name = os.path.basename(nav_file)
        new_name = org_name.split('.')[0] + '_temp.' + org_name.split('.')[-1]
        full_path = os.path.join(base_path,new_name)
        data_,header, n_eph = read_rinex3_nav(full_path,dataframe='no')
        os.remove(full_path) # removes the temp broadcasted file 
        if idx == 0:
            data = data_
        else:
            data = np.append(data,data_,axis=0)
            
    ### ------- Compute satellite coordiantes in ECEF and compute azimut and elevation angels
    ## -- Extracting approx postion from RINEX obs-file
    x = float(approxPosition[0])
    y = float(approxPosition[1])
    z = float(approxPosition[2])
    
    ## -- Running the computations for each system
    GNSS_FullName = dict(list(zip(['G','R','E','C'],['GPS','GLONASS','Galileo','BeiDou']))) # dict for mapping from system code to full name 
    sat_pos = {}                           # Dict for storing all data
    for sys in GNSS_SVs:
        sat_pos[sys] = {}                  # new dict for each system
        curr_pos = {}                      # Cells for storing data
        nepochs = len(GNSS_obs[sys])       # total nr of epochs 
        X = np.zeros([nepochs,61])         # Array for storing X-coordinate
        Y = np.zeros([nepochs,61])         # Array for storing X-coordinate
        Z = np.zeros([nepochs,61])         # Array for storing X-coordinate
        azimut    = np.zeros([nepochs,61]) # Array for storing satellites azimut angle
        elevation = np.zeros([nepochs,61]) # Array for storing Satellites elevation angle
        aktuelle_sat_list_all = []         # Dummy list for availible satellites
        for ep in np.arange(0,nepochs):
            aktuelle_sat = [PRN for PRN in list(GNSS_SVs[sys][ep][1::].astype(int)) if PRN !=0] ## added [1::] because first ele is nr of satellites
            aktuelle_sat_list_all.append(aktuelle_sat)
        aktuelle_sat_list_all = sorted(list({x for l in aktuelle_sat_list_all for x in l})) # finding uniqe value in list in list
            
        ## -- Run through all epochs
        null_epoch = time_epochs[:,1][0]
        for epoch in np.arange(0,nepochs):
            aktuelle_sat_list = [PRN for PRN in list(GNSS_SVs[sys][epoch][1::].astype(int)) if PRN !=0]   ## added [1::] because first ele is nr of satellites 30.01.2023
            print("\rCurrently computing coordinates for the %s system. Progress: %.1f%%" %(GNSS_FullName[sys],epoch/len(GNSS_obs[sys])*100), end='\r',flush=True)  # \r makes the line get overwritten
            ## -- Compute satellite coordinates for all availible satellites   
            curr_time = time_epochs[:,1][epoch] # Extracting time for RINEX obs-file
            df_data = pd.DataFrame(data) # making dataframe of data
            curr_data = df_data[df_data.iloc[:,0].str.contains(sys)].to_numpy() # extaction data for current system only
            ## -- Extract satellites ephemerids
            time_diff = curr_time - null_epoch
            ### Add if test here to speed up the process. Ephemerides get extracted only if time_diff is greater
            ### than the pre defined limit tLim_GEC and tLim_R
            if sys != 'R' and time_diff > tLim_GEC or sys == 'R' and time_diff > tLim_R or time_diff == 0:
                eph_data = np.empty((36,1))        # Preallocation array for storing ephemerids
                null_epoch = curr_time # resetting null_epoch
                for PRN in aktuelle_sat_list_all:
                    try: ## need try/except for satellites that in rinex obs, but not in nav file. consider adding functionality for his on the top script
                        eph = extract_nav_message(curr_data,PRN,curr_time)
                        eph[0] = eph[0][1::] # Removing system letter from number. Ex G10 -> 10
                    except:
                        eph = np.nan
                        continue
                    eph_data = np.column_stack((eph_data, eph))
                eph_data = pd.DataFrame(eph_data)
                eph_data.columns = eph_data.iloc[0].astype(int) # set first row as header
            ## -- Updating aktuelle_sat_list based on availible sat from both obs file and nav file 
            aktuelle_sat_list = list(set(aktuelle_sat_list) & set(list(eph_data.columns))) # (prevent error message if in obs but not in nav)
            for PRN in aktuelle_sat_list:
                ephemerides = eph_data[PRN].to_numpy().reshape(36,1) # passing curr_data instead to get correct system
                ephemerides = ephemerides.astype(float)
                if sys != 'R':
                    ## -- Computing satellite coordinates for GPS,Galileo and BeiDou if in Rinex file
                    X[epoch,PRN], Y[epoch,PRN], Z[epoch,PRN],_ = Satkoord2(ephemerides, curr_time, x, y, z)
                else:
                    ## -- Computing satellite coordinates for GLONASS if in Rinex file
                    curr_time = time_epochs[epoch]
                    pos, _, _, _= compute_GLO_coord_from_nav(ephemerides, curr_time)
                    X[epoch,PRN] = pos[0]
                    Y[epoch,PRN] = pos[1]
                    Z[epoch,PRN] = pos[2]
                ## - Compute azimut and elevation angle
                azimut[epoch,PRN],elevation[epoch,PRN] = compute_azimut_elev(X[epoch,PRN], Y[epoch,PRN], Z[epoch,PRN], x, y, z)
                ## -- Assign the computed variable to temporarly dicts for storing results   
                curr_pos[str(PRN)] = np.array([X[:,PRN],Y[:,PRN],Z[:,PRN]]).T
        
            ## -- Update dictionary with coordinates    
            sat_pos[sys]['Position']  = curr_pos
            sat_pos[sys]['Azimut']    = azimut
            sat_pos[sys]['Elevation'] = elevation
        
    return sat_pos
    
            

def shorten_navigation_file(navigationFile):
    """
    Function that is parsing a navigation file and removes system that
    not GPS,GLONASS, GALILEO and BeiDou. In addition it removes epochs that
    have same TOC (Time of clock) as the previous epoch. This function is only used
    to shorted down the reading time of nav file. In addition it will shorten the time
    it takes to compute satellite coordinates. (Shorter extraction time)
    
    Parameters:
    ----------
    tLim_GEC: Time limit for GPS,Galileo and BeiDou. If the time difference between 
              two epochs in the navigation file is less than this value, it will be removed
              to speed up the processing time.
              
    tLim_R  : Time limit for GLONASS. If the time difference between 
              two epochs in the navigation file is less than this value, it will be removed
              to speed up the processing time. It recommended to have a smalle limit
              for GLONASS since state vector are less suitet for interpolartion over
              larger time windows. 
    
    """
    fd = open(navigationFile,'r')
    lines = fd.readlines()
    del_start_indx = []
    del_end_indx = []
    pass_header = False

    for idx in np.arange(0,len(lines)):
        ## -- Checking if file contains 'S' code
        if re.match(r'S\d+',lines[idx]):
            del_start_indx.append(idx)
            IDX = idx + 1
            while 'S' not in lines[IDX]:
                IDX = IDX + 1
                if IDX == len(lines) or re.match(r'[J,S,G,R,I,E,C]\d+',lines[IDX]):
                    break
            del_end_indx.append(IDX)
                
        ## -- Checking if file contains 'I' code
        if re.match(r'I\d+',lines[idx]):
            del_start_indx.append(idx)
            IDX = idx + 1
            while 'I' not in lines[IDX]:
                IDX = IDX + 1
                if IDX == len(lines) or re.match(r'[J,S,G,R,I,E,C]\d+',lines[IDX]):
                    break
            del_end_indx.append(IDX)
        ## -- Checking if file contains 'J' code
        if re.match(r'J\d+',lines[idx]):
            del_start_indx.append(idx)
            IDX = idx + 1
            while 'J' not in lines[IDX]:
                IDX = IDX + 1
                if IDX == len(lines) or re.match(r'[J,S,G,R,I,E,C]\d+',lines[IDX]):
                    break
            del_end_indx.append(IDX)
            
    ## -- Deleting the indencies detected above
    del_start_indx.sort(reverse=True)
    del_end_indx.sort(reverse=True)
    for index,val in enumerate(del_start_indx):
        idx_start = del_start_indx[index]
        idx_end = del_end_indx[index]
        del lines[idx_start:idx_end]
        
    ## -- Some RINEX navigation files have several lines with same TOC
    ## This makes the reading and extraction time of ephemerides in connention
    ## with satelitte coordinates computation much slower. The script under
    ## is removing epochs if it has the same TOC as the previous epoch. 
    del_start_indx = []
    del_end_indx = []
    for idx in np.arange(0,len(lines)):
        ## -- Test if passed header
        if 'END OF HEADER' in lines[idx]:
            pass_header = True
        ## -- Checking for GPS, Galileo and BeiDou
        if pass_header==True and idx+8 < len(lines):
            if re.match(r'G\d+',lines[idx]) or re.match(r'E\d+',lines[idx]) or re.match(r'C\d+',lines[idx]):
                # line_curr = lines[idx][0:24] 
                line_curr = lines[idx][0:23] # change 29.01.2023 to remove "-" in end of seconds
                line_curr = [el for el in line_curr.split(" ") if el != ""]
                # line_next = lines[idx+8][0:24] 
                line_next = lines[idx+8][0:23] # change 29.01.2023 to remove "-" in end of seconds
                line_next = [el for el in line_next.split(" ") if el != ""]
                if line_curr == line_next:
                    IDX_start = idx + 8
                    IDX_end = IDX_start + 8
                    del_start_indx.append(IDX_start)
                    del_end_indx.append(IDX_end)
        ## -- Checking for GLONASS (less lines in file)     
        if pass_header==True and idx+4 < len(lines):
            if re.match(r'R\d+',lines[idx]):
                # line_curr = lines[idx][0:24] 
                line_curr = lines[idx][0:23] # change 29.01.2023  to remove "-" in end of seconds
                line_curr = [el for el in line_curr.split(" ") if el != ""]
                # line_next = lines[idx+4][0:24] 
                line_next = lines[idx+4][0:23] # change 29.01.2023  to remove "-" in end of seconds
                line_next = [el for el in line_next.split(" ") if el != ""]
                if line_curr == line_next:
                    IDX_start = idx + 4
                    IDX_end = IDX_start + 4
                    del_start_indx.append(IDX_start)
                    del_end_indx.append(IDX_end)
            
    ## -- Deleting the indencies that have equal TOC       
    del_start_indx.sort(reverse=True)
    del_end_indx.sort(reverse=True)
    for index,val in enumerate(del_start_indx):
        idx_start = del_start_indx[index]
        idx_end = del_end_indx[index]
        del lines[idx_start:idx_end]
        
        
    # ## -- Removing epoch that are within the time limit. This step is only for increasing speed.
    # del_start_indx = []
    # del_end_indx = []
    # dum_del = []
    # for idx in np.arange(0,len(lines)):
    #     ## -- Test if passed header
    #     if 'END OF HEADER' in lines[idx]:
    #         pass_header = True
    #     ## -- Checking for GPS, Galileo and BeiDou
    #     if pass_header==True and idx+8 < len(lines):
    #         if re.match(r'G\d+',lines[idx]) or re.match(r'E\d+',lines[idx]) or re.match(r'C\d+',lines[idx]):
    #             # line_curr = lines[idx][0:24] 
    #             line_curr = lines[idx][0:23] # change 29.01.2023 to remove "-" in end of seconds
    #             line_curr = [el for el in line_curr.split(" ") if el != ""]
    #             # line_next = lines[idx+8][0:24] 
    #             line_next = lines[idx+8][0:23] # change 29.01.2023 to remove "-" in end of seconds
    #             line_next = [el for el in line_next.split(" ") if el != ""]
    #             if abs(time_difference(line_curr, line_next)) < tLim_GEC and line_curr[0] == line_next[0]: # check if time differnce is less than 1 hour
    #                 dum = ['']    
    #                 for PRN_line in dum_del:
    #                     dum.append(PRN_line)
    #                 dum_del.append(line_next)
    #                 PRN_line = dum[-1]   
    #                 if line_next[0] in PRN_line:
    #                     if abs(time_difference(PRN_line, line_next)) < tLim_GEC:
    #                         # dum = ['']  
    #                         # dum_del = []
    #                         # dum_del.append(line_next)
    #                         IDX_start = idx + 8
    #                         IDX_end = IDX_start + 8
    #                         del_start_indx.append(IDX_start)
    #                         del_end_indx.append(IDX_end)
    #                 else:                              
    #                     IDX_start = idx + 8
    #                     IDX_end = IDX_start + 8
    #                     del_start_indx.append(IDX_start)
    #                     del_end_indx.append(IDX_end)
    
    #     ## -- Checking for GLONASS (less lines in file)     
    #     if pass_header==True and idx+4 < len(lines):
    #         if re.match(r'R\d+',lines[idx]):
    #             # line_curr = lines[idx][0:24] 
    #             line_curr = lines[idx][0:23] # change 29.01.2023  to remove "-" in end of seconds
    #             line_curr = [el for el in line_curr.split(" ") if el != ""]
    #             # line_next = lines[idx+4][0:24] 
    #             line_next = lines[idx+4][0:23] # change 29.01.2023  to remove "-" in end of seconds
    #             line_next = [el for el in line_next.split(" ") if el != ""]
    #             if time_difference(line_curr, line_next) < tLim_R and line_curr[0] == line_next[0]: # check if time differnce is less than 1 hour
    #                 IDX_start = idx + 4
    #                 IDX_end = IDX_start + 4
    #                 del_start_indx.append(IDX_start)
    #                 del_end_indx.append(IDX_end)
              
    # ## -- Deleting the indencies that have equal TOC       
    # del_start_indx.sort(reverse=True)
    # del_end_indx.sort(reverse=True)
    # for index,val in enumerate(del_start_indx):
    #     idx_start = del_start_indx[index]
    #     idx_end = del_end_indx[index]
    #     del lines[idx_start:idx_end]
    

    ## -- Write to file 
    base_path = os.path.split(navigationFile)[0]
    org_fileName = os.path.basename(navigationFile)
    temp_fileName = org_fileName.split('.')[0] + '_' + 'temp.'+ org_fileName.split('.')[-1] 
    full_path = os.path.join(base_path,temp_fileName)
    with open(full_path,'w') as fid:
        for line in lines:
            fid.write(line)
    fid.close()
    return



# def time_difference(list1, list2):
#     """
#     This function is compute the time difference between two list on the format
#     [year, month, day, hour, minute, and seconds]. 
#     """
#     list1 = [int(i) for i in list1[1::]]
#     list2 = [int(i) for i in list2[1::]]
#     date1 = datetime(list1[0], list1[1], list1[2], list1[3], list1[4], list1[5])
#     date2 = datetime(list2[0], list2[1], list2[2], list2[3], list2[4], list2[5])
#     difference = date2 - date1
#     return difference.total_seconds()

# time1 = [2022, 1, 1, 0, 0, 0]
# time2 = [2022, 1, 1, 0, 3, 0]
# print(time_difference(line_curr, line_next)) # should print 1.0
# eph_data = np.empty((36, 0))
# for sat in aktuelle_sat_list:
#     eph = extract_nav_message(curr_data,sat,t[i])
#     eph_data = np.column_stack((eph_data, eph))
#     eph_data = pd.DataFrame(eph_data)
#     df.columns = df.iloc[0] # set first row as header