import numpy as np
import re
import os.path

file = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\BRDC00IGS_R_20220010000_01D_MN.rnx'

fd = open(file,'r')
lines = fd.readlines()
# lines = list(reversed(lines_org))
del_start_indx = []
del_end_indx = []
pass_header = False



for idx in np.arange(0,len(lines)):
    # if 'END OF HEADER' in lines[idx]:
    #     pass_header = True
    ## -- Append index that contains 'S' code
    if re.match(r'S\d+',lines[idx]):
        del_start_indx.append(idx)
        # del_end_indx.append(idx + 3)
        IDX = idx + 1
        while 'S' not in lines[IDX]:
            IDX = IDX + 1
            if IDX == len(lines):
                break
        del_end_indx.append(IDX)
            
            
    ## -- Append index that contains 'I' code
    if re.match(r'I\d+',lines[idx]):
        del_start_indx.append(idx)
        IDX = idx + 1
        while 'I' not in lines[IDX]:
            IDX = IDX + 1
            if IDX == len(lines):
                break
        del_end_indx.append(IDX)
    ## -- Append index that contains 'J' code
    if re.match(r'J\d+',lines[idx]):
        del_start_indx.append(idx)
        IDX = idx + 1
        while 'J' not in lines[IDX]:
            IDX = IDX + 1
            if IDX == len(lines):
                break
        del_end_indx.append(IDX)
        
    


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
for idx in range(0,len(lines)):
    ## -- Test if passed header
    if 'END OF HEADER' in lines[idx]:
        pass_header = True
    ## -- Checking for GPS, Galileo and BeiDou
    if pass_header==True and idx+8 < len(lines):
        if re.match(r'G\d+',lines[idx]) or re.match(r'E\d+',lines[idx]) or re.match(r'C\d+',lines[idx]):
            line_curr = lines[idx][0:24] 
            line_curr = [el for el in line_curr.split(" ") if el != ""]
            line_next = lines[idx+8][0:24] 
            line_next = [el for el in line_next.split(" ") if el != ""]
            if line_curr == line_next:
                IDX_start = idx + 8
                IDX_end = IDX_start + 8
                del_start_indx.append(IDX_start)
                del_end_indx.append(IDX_end)
    ## -- Checking for GLONASS (less lines in file)     
    if pass_header==True and idx+4 < len(lines):
        if re.match(r'R\d+',lines[idx]):
            line_curr = lines[idx][0:24] 
            line_curr = [el for el in line_curr.split(" ") if el != ""]
            line_next = lines[idx+4][0:24] 
            line_next = [el for el in line_next.split(" ") if el != ""]
            if line_curr == line_next:
                IDX_start = idx + 4
                IDX_end = IDX_start + 4
                del_start_indx.append(IDX_start)
                del_end_indx.append(IDX_end)
        
        
del_start_indx.sort(reverse=True)
del_end_indx.sort(reverse=True)
for index,val in enumerate(del_start_indx):
    idx_start = del_start_indx[index]
    idx_end = del_end_indx[index]
    del lines[idx_start:idx_end]   
    
# for index,val in enumerate(sorted(del_start_indx, reverse=True)):
#     start_indx = del_start_indx[index]
#     end_indx = del_end_indx[index]
#     dum.append(lines[start_indx:end_indx])
#     del lines[start_indx:end_indx]
    
with open(r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\TEST.rnx','w') as fid:
    for line in lines:
        fid.write(line)
fid.close()

#%%
import numpy as np
import re
import os
file = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\BRDC00IGS_R_20220010000_01D_MN.rnx'

def shorten_navigation_file(navigationFile):
    """
    Function that is parsing a navigation file and removes system that
    not GPS,GLONASS, GALILEO and BeiDou. In addition it removes epochs that
    have same TOC (Time of clock) as the previous epoch. This function is only used
    to shorted down the reading time of nav file. In addition it will shorten the time
    it takes to compute satellite coordinates. (Shorter extraction time)
    """
    fd = open(navigationFile,'r')
    lines = fd.readlines()
    # lines = list(reversed(lines_org))
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
                if IDX == len(lines):
                    break
            del_end_indx.append(IDX)
                
        ## -- Checking if file contains 'I' code
        if re.match(r'I\d+',lines[idx]):
            del_start_indx.append(idx)
            IDX = idx + 1
            while 'I' not in lines[IDX]:
                IDX = IDX + 1
                if IDX == len(lines):
                    break
            del_end_indx.append(IDX)
        ## -- Checking if file contains 'J' code
        if re.match(r'J\d+',lines[idx]):
            del_start_indx.append(idx)
            IDX = idx + 1
            while 'J' not in lines[IDX]:
                IDX = IDX + 1
                if IDX == len(lines):
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
    for idx in range(0,len(lines)):
        ## -- Test if passed header
        if 'END OF HEADER' in lines[idx]:
            pass_header = True
        ## -- Checking for GPS, Galileo and BeiDou
        if pass_header==True and idx+8 < len(lines):
            if re.match(r'G\d+',lines[idx]) or re.match(r'E\d+',lines[idx]) or re.match(r'C\d+',lines[idx]):
                line_curr = lines[idx][0:24] 
                line_curr = [el for el in line_curr.split(" ") if el != ""]
                line_next = lines[idx+8][0:24] 
                line_next = [el for el in line_next.split(" ") if el != ""]
                if line_curr == line_next:
                    IDX_start = idx + 8
                    IDX_end = IDX_start + 8
                    del_start_indx.append(IDX_start)
                    del_end_indx.append(IDX_end)
        ## -- Checking for GLONASS (less lines in file)     
        if pass_header==True and idx+4 < len(lines):
            if re.match(r'R\d+',lines[idx]):
                line_curr = lines[idx][0:24] 
                line_curr = [el for el in line_curr.split(" ") if el != ""]
                line_next = lines[idx+4][0:24] 
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
    base_path = os.path.split(navigationFile)[0]
    org_fileName = os.path.basename(navigationFile)
    temp_fileName = org_fileName.split('.')[0] + '_' + 'temp.'+ org_fileName.split('.')[-1] 
    full_path = os.path.join(base_path,temp_fileName)
    with open(full_path,'w') as fid:
        for line in lines:
            fid.write(line)
    fid.close()
    
    return

file = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles/OPEC00NOR_S_20220010000_01D_GN.rnx'
shorten_navigation_file(file)
#%%
from readRinexNav import *
data,_,_ = read_rinex3_nav(r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\TEST.rnx')

#%%
import pickle
file = r'C:\Users\perhe\OneDrive\Documents\Python_skript\Output_Files_beforatan2/analysisResults.pkl'
file = open(file, "rb")
loaded_dictionary = pickle.load(file)




