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
        block = []
        block_arr = np.array([])
        ## -- Read first line of navigation message
        line = filnr.readline().rstrip()
        
        sys_PRN = line[0:3]
        # -- Passing S, I, R and J systems
        # while 'S' in sys_PRN or 'I' in sys_PRN or 'R' in sys_PRN or 'J' in sys_PRN:
        #     line = filnr.readline().rstrip()
        #     sys_PRN = line[0:3]
        #     while sys_PRN == '   ':
        #         line = filnr.readline().rstrip()
        #         sys_PRN = line[0:3]
            
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
            if line[idx] == 'e' or line[idx] == 'E':
                line = line[:idx+4] + " " + line[idx+4:]
        
        fl = [el for el in line.split(" ") if el != ""]
        block_arr =np.append(block_arr,np.array([fl]))
        block_arr = block_arr.reshape(1,len(block_arr))
        
        
        ## Looping throug the next 3-lines for current message (satellitte) (GLONASS)
        if 'R' in sys_PRN:  #PH added if test to add GLONASS 12.12.2022
            for i in range(0,3):
                line = filnr.readline().rstrip()
                ## -- Have to add space between datacolums where theres no whitespace
                for idx, val in enumerate(line):
                    if line[idx] == 'e':
                        line = line[:idx+4] + " " + line[idx+4:]
    
                ## --Reads the line vector nl from the text string line and adds navigation 
                # message for the relevant satellite n_sat. It becomes a long line vector 
                # for the relevant message and satellite.
            
                nl = [el for el in line.split(" ") if el != ""]
        
                block_arr = np.append(block_arr,np.array([nl]))
                block_arr = block_arr.reshape(1,len(block_arr))
        else:   
            ## Looping throug the next 7-lines for current message (satellitte) (GPS,Galileo,BeiDou)
            for i in range(0,7):
                line = filnr.readline().rstrip()
                ## -- Have to add space between datacolums where theres no whitespace
                for idx, val in enumerate(line):
                    if line[idx] == 'e' or line[idx] == 'E':
                        line = line[:idx+4] + " " + line[idx+4:]
                
                ## --Reads the line vector nl from the text string line and adds navigation 
                # message for the relevant satellite n_sat. It becomes a long line vector 
                # for the relevant message and satellite.
                
                if i == 6 and line.count('e') == 1:
                    line = line + 'nan'
                    
                if i == 4 and 'E' in sys_PRN and line.count('e') < 4:
                    line = line + 'nan'
        
                nl = [el for el in line.split(" ") if el != ""]
        
                block_arr = np.append(block_arr,np.array([nl]))
                block_arr = block_arr.reshape(1,len(block_arr))
        
    
    
        ## -- Collecting all data into common variable
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

