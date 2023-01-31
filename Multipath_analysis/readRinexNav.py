def read_rinex2_nav(filename, dataframe = None):    
    """
    Reads the navigation message from GPS broadcast efemerids in RINEX v.2 format. (GPS only)
  
    Reads one navigation message at a time until the end of the row. Accumulate in
    in a common matrix, "data", where there is a line for each message. 
    Note that each message forms a block in the file and that the same satellite
    can have several messages, usually with an hour's difference in reference time.
    
    
    Parameters
    ----------
    filename : Filename of the RINEX navigation file
    dataframe : Set to 'yes' or 'YES' to get the data output as a pandas DataFrame (array as default)
    
    Returns
    -------
    data : Matrix with data for all epochs
    header: List with header content
    n_eph: Number of epochs


    """
  
    import numpy as np
    from pandas import DataFrame
    
    try:
        print('Reading broadcast ephemeris from RINEX-navigation file.....')
        filnr = open(filename, 'r')
    except OSError:
        print("Could not open/read file: %s", filename)
        

    line = filnr.readline().rstrip()
    header = []
    while 'END OF HEADER' not in line:
        line = filnr.readline().rstrip()
        header.append(line)
        
    data  = np.zeros((1,36))
   
    while line != '':
        block_arr = np.array([])
        
        ## -- Read first line of navigation message
        line = filnr.readline().rstrip()
        
        # Replacing 'D' with 'E' ('D' is fortran syntax for exponentiall form)
        line = line.replace('D','E')
        ## -- Have to add space between datacolums where theres no whitespace
        for idx, val in enumerate(line):
            if line[0:2] != ' ' and line[22] != ' ':
                line = line[:22] + " " + line[22:]
            if line[idx] == 'E':
                line = line[:idx+4] + " " + line[idx+4:]
        
        fl = [el for el in line.split(" ") if el != ""]
        block_arr =np.append(block_arr,np.array([fl]))
        block_arr = block_arr.reshape(1,len(block_arr))
        
        ## Looping throug the next 7-lines for current message (satellitte)
        for i in np.arange(0,7):
            line = filnr.readline().rstrip()
            ## -Replacing 'D' with 'E'
            line = line.replace('D','E')
            
            ## -- Have to add space between datacolums where theres no whitespace
            for idx, val in enumerate(line):
                if line[idx] == 'E':
                    line = line[:idx+4] + " " + line[idx+4:]
            
            ## --Reads the line vector nl from the text string line and adds navigation 
            # message for the relevant satellite n_sat. It becomes a long line vector 
            # for the relevant message and satellite.
            nl = [el for el in line.split(" ") if el != ""]
            block_arr = np.append(block_arr,np.array([nl]))
            block_arr = block_arr.reshape(1,len(block_arr))
    
        ## -- Collecting all data into common variable    
        if np.size(block_arr) != 0:
            data  = np.concatenate([data , block_arr], axis=0)
        else:
            data  = np.delete(data , (0), axis=0)
            print('File %s is read successfully!' % (filename))
            
    filnr.close()
    n_eph = len(data)
    data = data.astype(float)
    if dataframe == 'yes' or dataframe == 'YES':
        data = DataFrame(data)
        

    return data, header, n_eph

# data, header, n_eph = read_rinex2_nav('testfile.20n')

def read_rinex3_nav(filename, dataframe = None):    
    """
    Reads the navigation message from GPS broadcast efemerids in RINEX v.3 format. (NOT GLONASS!)
  
    Reads one navigation message at a time until the end of the row. Accumulate in
    in a common matrix, "data", where there is a line for each message. 
    Note that each message forms a block in the file and that the same satellite
    can have several messages, usually with an hour's difference in reference time.
    
    
    Parameters
    ----------
    filename : Filename of the RINEX navigation file
    dataframe : Set to 'yes' or 'YES' to get the data output as a pandas DataFrame (array as default) 
    
    Returns
    -------
    data : Matrix with data for all epochs (or dataframe)
    header: List with header content
    n_eph: Number of epochs


    """

    import numpy as np
    from pandas import DataFrame

    
    try:
        filnr = open(filename, 'r')
        print('Reading broadcast ephemeris from RINEX-navigation file.....')
    except OSError:
        print("Could not open/read file: %s", filename)
        return
        
    
    line = filnr.readline().rstrip()
    header = []
    while 'END OF HEADER' not in line:
        line = filnr.readline().rstrip()
        header.append(line)
        
    data  = np.zeros((1,36))
    while line != '':
        block_arr = np.array([])
        ## -- Read first line of navigation message
        line = filnr.readline().rstrip()
        # ## --- Replace "E" with "e" i academic notation
        # line.replace('E+', 'e+')
        # line.replace('E-', 'e-')
        
        sys_PRN = line[0:3]            
        while 'S' in sys_PRN or 'I' in sys_PRN or 'J' in sys_PRN: ## trying to add GLONASS #PH
            line = filnr.readline().rstrip()
            sys_PRN = line[0:3]
            while sys_PRN == '   ':
                line = filnr.readline().rstrip()
                sys_PRN = line[0:3]
            
        
        ## -- Have to add space between datacolums where theres no whitespace
        for idx, val in enumerate(line):
            if line[0:2] != ' ' and line[23] != ' ':
                line = line[:23] + " " + line[23:]
            if line[idx] == 'e' or line[idx] == 'E' and idx !=0:
                line = line[:idx+4] + " " + line[idx+4:]
        
        fl = [el for el in line.split(" ") if el != ""]
        block_arr = np.append(block_arr,np.array([fl]))
        block_arr = block_arr.reshape(1,len(block_arr))
        
        
        ## Looping throug the next 3-lines for current message (satellitte) (GLONASS)
        if 'R' in sys_PRN:  #PH added if test to add GLONASS 12.12.2022
            for i in np.arange(0,3):
                line = filnr.readline().rstrip()
                ## -- Have to add space between datacolums where theres no whitespace
                for idx, val in enumerate(line):
                    if line[idx] == 'e' or line[idx] == 'E':
                        line = line[:idx+4] + " " + line[idx+4:]
    
                ## --Reads the line vector nl from the text string line and adds navigation 
                # message for the relevant satellite n_sat. It becomes a long line vector 
                # for the relevant message and satellite.
            
                nl = [el for el in line.split(" ") if el != ""]
        
                block_arr = np.append(block_arr,np.array([nl]))
                block_arr = block_arr.reshape(1,len(block_arr))
        else:   
            ## Looping throug the next 7-lines for current message (satellitte) (GPS,Galileo,BeiDou)
            for i in np.arange(0,7):
                line = filnr.readline().rstrip()
                if line == '': #PH 07.01.2023
                    continue
                else:
                    ## -- Have to add space between datacolums where theres no whitespace
                    for idx, val in enumerate(line):
                        if line[idx] == 'e' or line[idx] == 'E':
                            line = line[:idx+4] + " " + line[idx+4:]
                    
                    ## --Reads the line vector nl from the text string line and adds navigation 
                    # message for the relevant satellite n_sat. It becomes a long line vector 
                    # for the relevant message and satellite.
                    
                    ## Runs through line to see if each line contains 4 objects. If not, adds nan.
                    if i < 6 and line.lower().count('e') < 4:
                        if line[10:20].strip() == '':
                            line = line[:10] +  'nan' + line[10:]
                        if line[30:40].strip() == '':
                            line = line[:30] +  'nan' + line[30:]
                        if line[50:60].strip() == '':
                            line = line[:50] +  'nan' + line[50:]
                        if line[70:80].strip() == '':
                            line = line[:70] +  'nan' + line[70:]
                        
                    
                    if i == 6 and line.lower().count('e') < 2 and 'E' not in sys_PRN:
                        if line[10:20].strip() == '':
                            line = line[:10] +  'nan' + line[10:]
                        if line[30:40].strip() == '':
                            line = line[:30] +  'nan' + line[30:]
                            
                    if i == 6 and line.lower().count('e') < 1 and 'E' in sys_PRN: #only one object in last line for Galileo
                        if line[10:20].strip() == '':
                            line = line[:10] +  'nan' + line[10:]
                    
                    if i == 6 and 'E' in sys_PRN: #only one object in last line for Galileo, but add nan to get 36 in total
                        line = line + 'nan'
                        
                    
                    # if i == 6 and line.lower().count('e') == 1:
                    #     line = line + 'nan'
                        
                    # if i == 4 and 'E' in sys_PRN and line.lower().count('e') < 4:
                    #     line = line + 'nan'
            
                    nl = [el for el in line.split(" ") if el != ""]
            
                    block_arr = np.append(block_arr,np.array([nl]))
                    block_arr = block_arr.reshape(1,len(block_arr))
        
    
    
        ## -- Collecting all data into common variable
        if block_arr.shape[1] > 36:
            block_arr = block_arr[:,0:36]
        try:
            if np.size(block_arr) != 0 and 'R' not in sys_PRN:
                data  = np.concatenate([data , block_arr], axis=0)
                
            elif np.size(block_arr) != 0 and 'R' in sys_PRN:
                GLO_dum = np.zeros([1,np.size(data,axis=1) - np.size(block_arr,axis=1)])
                block_arr = np.append(block_arr,GLO_dum) # adding emtpy columns to match size of other systems
                block_arr = block_arr.reshape(1,len(block_arr))
                data  = np.concatenate([data , block_arr], axis=0)
            else:
                data  = np.delete(data , (0), axis=0)
                print('File %s is read successfully!' % (filename))
        except:
            print("ERROR! Sure this is a RINEX v.3 navigation file?")
            break
            

    filnr.close()
    n_eph = len(data)
    if dataframe == 'yes' or dataframe == 'YES':
        data = DataFrame(data)
        data_columns = list(data.columns)
        data_columns.pop(0) # Removing index that contains PRN nr ex 'G01'
        data[data_columns] = data[data_columns].astype(float) # Change the other values to float
        

    return data, header, n_eph




# data, header, n_eph = read_rinex3_nav('opec_3.18n',dataframe='yes')

# data, header, n_eph = read_rinex3_nav('OPEC00NOR_S_20220010000_01D_GN.rnx',dataframe='no')



