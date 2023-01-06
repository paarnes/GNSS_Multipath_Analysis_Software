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
        block = []
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
        for i in range(0,7):
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



