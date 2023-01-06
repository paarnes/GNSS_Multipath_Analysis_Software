# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 19:14:31 2022

@author: perhe
"""

import os, sys
os.chdir(r"C:\Users\perhe\OneDrive\Documents\Python skripts\GNSS\Read SP3\24")

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import os


# %%
def readMatlabResults(headerFilename, filename1, filename2):
    
    # initialize dictionary to store all results
    results = {}
    
    #read lines from tct files
    with open(headerFilename, 'r') as hfid:
        header = hfid.readlines()
    
    # load csv files
    data1 = pd.read_csv(filename1, sep=',', header= None).to_numpy()
    data2 = pd.read_csv(filename2, sep=',', header= None).to_numpy()
    
    
    # extract first line of header
    linenum = 0
    line = header[linenum]
    
    # move cursor to METADATA section
    while not('METADATA' in line): 
        # next line
        linenum += 1
        line = header[linenum]
    
    
    # store metadata    
    metadata_dict = {}
    
    linenum += 1
    line = header[linenum]
    
    # read all of metadata section
    while not('METADATA ENDS HERE' in line): 
        
        if not('*'  in line[0]):
            print('ERROR: Line in metadata does not have an expected asteric indicating variable format.')
            return
        
        infoName = line.split(':')[0]
        info =     ':'.join(line.split(':')[1:])
        
        if '****' in infoName:  # GLONASS satellite channel map
            infoName = infoName.split('**** ')[1]
            info = info.split(',')
            
            GLONASS_channel_map = {}
            for satChannelPair in info:
                satID = int(satChannelPair.split('/')[0])
                channel = int(satChannelPair.split('/')[1])
                GLONASS_channel_map[satID] = channel
            
            metadata_dict[infoName] = GLONASS_channel_map            
        elif '***' in infoName: # float values
            infoName = infoName.split('*** ')[1]
            info = float(info)
            metadata_dict[infoName] = info
        elif '**' in infoName: # date values
            infoName = infoName.split('** ')[1]
            info = info.split('/')
            
            year = int(info[0])
            month = int(info[1])
            rest = info[2]
            
            day = int(rest.split(' ')[0])
            
            rest = rest.split(' ')[1]
            rest = rest.split(':')
            
            hour = int(rest[0])
            minute = int(rest[1])
            sec = float(rest[2])
            
            fullDate = [year, month, day, hour, minute, sec]
            
            metadata_dict[infoName] = fullDate
        else: # string values
            infoName = infoName.split('* ')[1]
            metadata_dict[infoName] = info.rstrip('\n')
        
        #next line
        linenum += 1
        line = header[linenum]
        
    # store metadata dictionarie
    results['metadata'] = metadata_dict
        
    # %%
    # move cursor to GNSS systems list
    while not('GNSS systems:' in line): 
        # next line
        linenum += 1
        line = header[linenum]
    
    # Store GNSS systems
    line = line.split(':')[1]
    line = line.rstrip('\n')
    line = line.split(',')
    
    GNSS_systems = []
    for sys in line:
        GNSS_systems.append(sys)
    
    results['GNSS_systems'] = GNSS_systems
    results['nGNSS_systems'] = len(GNSS_systems)
    
    # iterate through GNSS systems
    for sys in GNSS_systems:
        
        # system dictionary
        sys_dict = {}
        
        # move cursor to system start
        while not(sys in line): 
            # next line
            linenum += 1
            line = header[linenum]
        
         # move cursor to Band list 
        while not('Bands:' in line): 
            # next line
            linenum += 1
            line = header[linenum]
        
        # Store Bands
        line = line.split(':')[1]
        line = line.rstrip('\n')
        line = line.split(',')
    
        bands = []
        for band in line:
            bands.append(band)
            
        sys_dict['Bands'] = bands
        sys_dict['nBands'] = len(bands)
        
        # iterate through bands
        for band in bands:
            
            band_dict = {}
            
            # next line
            linenum += 1
            line = header[linenum]
            
            # move cursor to band
            while not(band in line): 
                # next line
                linenum += 1
                line = header[linenum]
                
            # next line
            linenum += 1
            line = header[linenum]
            
            # Store codes
            line = line.split(':')[1]
            line = line.rstrip('\n')
            line = line.split(',')
        
            codes = []
            for code in line:
                codes.append(code)
                
            band_dict['Codes'] = codes
            band_dict['nCodes'] = len(codes)
            
            # next line
            linenum += 1
            line = header[linenum]
            
            # iterate through codes
            for code in codes:
                
                code_dict = {}
                 # move cursor to code
                while not(code in line): 
                    # next line
                    linenum += 1
                    line = header[linenum]
                    
                # next line
                linenum += 1
                line = header[linenum]
                
                # store data from data1 file
                while not(code + ' ENDS HERE' in line):
                    
                    infoName = line.split(':')[0]
                    columns =     ':'.join(line.split(':')[1:])
                    start_column = int(columns.split(',')[0]) - 1 # -1 since Matlab indexes from 1
                    end_column = int(columns.split(',')[1]) - 1 # -1 since Matlab indexes from 1
                    
                    code_dict[infoName] = data1[:, start_column:end_column+1]
                    
                    
                    # next line
                    linenum += 1
                    line = header[linenum]
                    
                band_dict[code] = code_dict
            sys_dict[band] = band_dict
        results[sys] = sys_dict
                
    #%%
    
    linenum += 1
    line = header[linenum]    
    
     # move cursor to observation overview
    while not('DATA FILE 2 DESCRIPTION' in line): 
        # next line
        linenum += 1
        line = header[linenum]
        
        #iterate through systems
    for sys in GNSS_systems:
        
        #systems dictionaries 
        sys_dict = results[sys]
        
        # move cursor to system start
        while not(sys in line): 
            # next line
            linenum += 1
            line = header[linenum]
            
        # next line
        linenum += 1
        line = header[linenum]
            
        bands = results[sys]['Bands']
        
        # iterate through bands
        for band in bands:
            
            band_dict = sys_dict[band]
            
            # next line
            linenum += 1
            line = header[linenum]
            
            # move cursor to band
            while not(band in line): 
                # next line
                linenum += 1
                line = header[linenum]
                
            # next line
            linenum += 1
            line = header[linenum]
            # next line
            linenum += 1
            line = header[linenum]
            
            codes = results[sys][band]['Codes']
            
            # iterate through codes
            for code in codes:
                
                code_dict = band_dict[code]
                
                 # move cursor to code
                while not(code in line): 
                    # next line
                    linenum += 1
                    line = header[linenum]
                    
                # next line
                linenum += 1
                line = header[linenum]
                
                # load data from data2 file
                while not(code + ' ENDS HERE' in line):
                    
                    infoName = line.split(':')[0]
                    info =     line.split(':')[1]
                    
                    # if info is range/phase codes, convert from ascii values to string
                    if infoName == 'range1/phase1/range2/phase2 Codes':
                        nRows = int(info.split(',')[0])
                        start_column = int(info.split(',')[1]) - 1 # -1 since Matlab indexes from 1
                        end_column = int(info.split(',')[2]) - 1 # -1 since Matlab indexes from 1
                        
                        obsCodes = data2[:nRows, start_column:end_column+1]
                        
                        code_dict['range1_Code'] =  chr(int(obsCodes[0, 0])).upper() + str(int(obsCodes[1, 0])) + chr(int(obsCodes[2, 0])).upper()  
                        code_dict['phase1_Code'] =  chr(int(obsCodes[0, 1])).upper() + str(int(obsCodes[1, 1])) + chr(int(obsCodes[2, 1])).upper() 
                        code_dict['range2_Code'] =  chr(int(obsCodes[0, 2])).upper() + str(int(obsCodes[1, 2])) + chr(int(obsCodes[2, 2])).upper() 
                        code_dict['phase2_Code'] =  chr(int(obsCodes[0, 3])).upper() + str(int(obsCodes[1, 3])) + chr(int(obsCodes[2, 3])).upper() 
                    else:
                    
                        nRows = int(info.split(',')[0])
                        start_column = int(info.split(',')[1]) - 1 # -1 since Matlab indexes from 1
                        end_column = int(info.split(',')[2]) - 1 # -1 since Matlab indexes from 1
                        
                        code_dict[infoName] = data2[:nRows, start_column:end_column+1]
                    
                    
                    # next line
                    linenum += 1
                    line = header[linenum]
                    
                band_dict[code] = code_dict
            sys_dict[band] = band_dict
        results[sys] = sys_dict
                    
                    
                    
          #%%      
    linenum += 1
    line = header[linenum]    
    
     # move cursor to observation overview
    while not('OBSERVATION OVERVIEW' in line): 
        # next line
        linenum += 1
        line = header[linenum]
        
    
            
            
    for sys in GNSS_systems:
        
        
        obsOverview_dict = {}
            
        # next line
        linenum += 1
        line = header[linenum]   
        
        # move cursor to observation overview
        while not(sys in line): 
            # next line
            linenum += 1
            line = header[linenum]
            
        # next line
        linenum += 1
        line = header[linenum]   
        
        nSat = int(line.split(':')[1])
        
        # move cursor to first sat
        while not('Sat_1' in line): 
            # next line
            linenum += 1
            line = header[linenum]
           
        
        
        for satI in range(nSat):
            linenum += 1
            line = header[linenum]
    
            sat = 'Sat_'+str(satI)
            sat_dict = {}
            line = line.split('Bands:')[1].rstrip('\n')
            bands = line.split(',')
            sat_dict['Bands'] = bands
            
            for i, band in enumerate(bands):
                linenum += 1
                line = header[linenum]
                codes = line.split(':')[1]
                sat_dict[band] = codes
            
            linenum += 1
            line = header[linenum]
            
            obsOverview_dict[sat] = sat_dict
        results[sys]['observationOverview'] = obsOverview_dict
        
    return results
   
    
    
if __name__ == '__main__':
    
    headerFilename = "OPEC00NOR_S_20220240000_01D_30S_MO_3.04_Data_Output_Header.txt"
    filename1 = 'OPEC00NOR_S_20220240000_01D_30S_MO_3.04_Data_Output_File_1.csv'
    filename2 = 'OPEC00NOR_S_20220240000_01D_30S_MO_3.04_Data_Output_File_2.csv'
    
    results = readMatlabResults(headerFilename, filename1, filename2)