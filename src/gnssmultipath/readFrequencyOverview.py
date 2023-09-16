def readFrequencyOverview(filename):
    """
     Function that reads current frequencies of GNSS satellite carrier bands,
     from the text file "Rinex_Frequency_Overview.txt". The band numbers follow the conventions of RINEX 3.
    --------------------------------------------------------------------------------------------------------------------------
     INPUTS:
    
     filename:             string, path and filename of frequency overview file             
    --------------------------------------------------------------------------------------------------------------------------
     OUTPUTS:
    
     frequencyOverview:    dict. Each key conatins frequencies of one
                           GNSS system. The order is decided by GNSSsystems.
                           If the cell is for any system but GLONASS, the cell
                           cotains an array with the frequency of band i at
                           index i.
    
                           frequencyOverview[GNSSsystemIndex][bandnumber]
                           if not GLONASS
    
                           For GLONASS, the value for the key is a array with two columns.
                           First column gives core frequency of that band.
                           Second column gives increment frequency of that
                           band. 
    
     GNSSsystems:          dict. Contains codes for each GNSS system. Must be
                           either 'G', 'R', 'E', 'C'
    
     success:              boolean, 0 if error is thrown, 1 otherwise
    --------------------------------------------------------------------------------------------------------------------------
    """
    import numpy as np

    success = 1
    fid = open(filename, 'r')
    
    ## -- Initialize variables
    nGNSSsystems = 0
    GNSSsystems = {}
    frequencyOverview = {}
    
    line = fid.readline().rstrip()    
    ## Continue until end of file marker reached
    while 'eof' not in line:
        ## read past header. first line after header has >
        while '>' not in line:    
            line = fid.readline().rstrip()
            if 'eof' in line:
                print('ERROR(readFrequencyOverview): End of file reached unexpectedly!')
                success = 0
                return success
    
        # increment number of GNSS systems
        nGNSSsystems = nGNSSsystems + 1;
        # Store current GNSS system
        GNSSsystems[nGNSSsystems] = str(line[1])
        ## -- Current GNSS system 
        currentGNSSsystem = GNSSsystems[nGNSSsystems]
       
        if currentGNSSsystem in ["G", "R", "E", "C"]:
            if currentGNSSsystem != 'R':
                frequencyOverview[nGNSSsystems] = np.zeros([9,1])
                for k in range(0,9):
                    line = fid.readline().rstrip()
                    line = line[2::]
                    line_ = [el for el in line.split(" ") if el != ""]
                    frequencyOverview[nGNSSsystems][k] = float(line_.pop(0))
                
            else: # Reading in GLONASS frequency
                frequencyOverview[nGNSSsystems] = np.zeros([9,2]);
    
                for k in range(0,9):
                    line = fid.readline().rstrip()
                    line = line[2::]
                    line_ = [el for el in line.split(" ") if el != ""]
                    freq = float(line_.pop(0))
                    if line_ != []:
                        dFreq = float(line_.pop(0))
                    else:
                        dFreq = np.nan
                        
                    frequencyOverview[nGNSSsystems][k,:] = [freq, dFreq]
    
        ## -- Next line
        line = fid.readline().rstrip()
    
    
    print('INFO(readFrequencyOverview): Frequency overview file has been read')

    return frequencyOverview, GNSSsystems, success

