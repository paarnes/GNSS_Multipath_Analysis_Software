def combineSP3Nav(three_sp3_files,sat_positions_1, epoch_dates_1, navGNSSsystems_1, \
                  nEpochs_1, epochInterval_1,sat_positions_2, epoch_dates_2, navGNSSsystems_2,\
                  nEpochs_2, epochInterval_2, sat_positions_3, epoch_dates_3, navGNSSsystems_3,\
                  nEpochs_3, epochInterval_3,GNSSsystems):
    
    """
    # Function that combines the precise orbital data of two or three SP3
    # files. Note that the SP3 files should first be read by the function
    # readSP3Nav.m. 
    #--------------------------------------------------------------------------------------------------------------------------
    # INPUTS:
    
    # three_sp3_files:      boolean. 1 if there are three SP3 files to be
    #                       combined, 0 otherwise
    
    # sat_positions_1:      cell. Conatains data from first SP3 file. Each cell 
    #                       elements contains position data for a specific GNSS 
    #                       system. Order is defined by order of navGNSSsystems_1. 
    #                       Each cell element is another cell that stores 
    #                       position data of specific satellites of that 
    #                       GNSS system. Each of these cell elements is a matrix 
    #                       with [X, Y, Z] position of a epoch in each row.
    
    #                       sat_positions_1{GNSSsystemIndex}{PRN}(epoch, :) = [X, Y, Z]
    
    # epoch_dates_1:        matrix. Each row contains date of one of the epochs 
    #                       from the first SP3 file
    #                       [nEpochs_1 x 6]
    
    # navGNSSsystems_1:     array. Contains string. Each string is a code for a
    #                       GNSS system with position data stored in sat_positions.
    #                       Must be one of: "G", "R", "E", "C"
    
    # nEpochs_1:            number of position epochs in first SP3 file, integer
    
    # epochInterval_1:      interval of position epochs in first SP3 file, seconds
    
    # sat_positions_2:      cell. Conatains data from second SP3 file. Each cell 
    #                       elements contains position data for a specific GNSS 
    #                       system. Order is defined by order of navGNSSsystems_2. 
    #                       Each cell element is another cell that stores 
    #                       position data of specific satellites of that 
    #                       GNSS system. Each of these cell elements is a matrix 
    #                       with [X, Y, Z] position of a epoch in each row.
    
    #                       sat_positions_2{GNSSsystemIndex}{PRN}(epoch, :) = [X, Y, Z]
    
    # epoch_dates_2:        matrix. Each row contains date of one of the epochs 
    #                       from the second SP3 file
    #                       [nEpochs_1 x 6]
    
    # navGNSSsystems_2:     array. Contains string. Each string is a code for a
    #                       GNSS system with position data stored in sat_positions.
    #                       Must be one of: "G", "R", "E", "C"
    
    # nEpochs_2:            number of position epochs in first SP3 file, integer
    
    # epochInterval_2:      interval of position epochs in second SP3 file, seconds
    
    # sat_positions_3:      cell. Conatains data from third SP3 file. Each cell 
    #                       elements contains position data for a specific GNSS 
    #                       system. Order is defined by order of navGNSSsystems_3. 
    #                       Each cell element is another cell that stores 
    #                       position data of specific satellites of that 
    #                       GNSS system. Each of these cell elements is a matrix 
    #                       with [X, Y, Z] position of a epoch in each row.
    
    #                       sat_positions_3{GNSSsystemIndex}{PRN}(epoch, :) = [X, Y, Z]
    
    # epoch_dates_3:        matrix. Each row contains date of one of the epochs 
    #                       from the third SP3 file
    #                       [nEpochs_1 x 6]
    
    # navGNSSsystems_3 :    array. Contains string. Each string is a code for a
    #                       GNSS system with position data stored in sat_positions.
    #                       Must be one of: "G", "R", "E", "C"
    
    # nEpochs_3:            number of position epochs in third SP3 file, integer
    
    # epochInterval_3:      interval of position epochs in third SP3 file, seconds
    #--------------------------------------------------------------------------------------------------------------------------
    # OUTPUTS:
    
    # sat_positions:        cell. Conatains data from all two/three SP3 file. Each cell 
    #                       elements contains position data for a specific GNSS 
    #                       system. Order is defined by order of navGNSSsystems_1. 
    #                       Each cell element is another cell that stores 
    #                       position data of specific satellites of that 
    #                       GNSS system. Each of these cell elements is a matrix 
    #                       with [X, Y, Z] position of a epoch in each row.
    
    # epoch_dates:          matrix. Each row contains date of one of the epochs 
    #                       from all two/three SP3 file
    #                       [nEpochs_1 x 6]
    
    # navGNSSsystems:       array. Contains string. Each string is a code for a
    #                       GNSS system with position data stored in sat_positions.
    #                       Must be one of: "G", "R", "E", "C"
    
    # nEpochs:              number of position epochs in all two/three SP3 file, integer
    
    # epochInterval:        interval of position epochs in all SP3 file, seconds
    
    
    # success:              boolean, 1 if no error occurs, 0 otherwise
    #--------------------------------------------------------------------------------------------------------------------------
    
    """
    
    import numpy as np
    import copy 
    
    success = 1 # Setting success to 1 
    
    ## -- Check that first two SP3 files have same interval
    if epochInterval_1 != epochInterval_2:
        print('ERROR(combineSP3Nav):The first and second SP3 files do not have the same epoch interval') 
        sat_positions, epoch_dates, navGNSSsystems, nEpochs, epochInterval = np.nan,np.nan,np.nan,np.nan,np.nan
        success = 0
        return success
    
    
    ## -- Check that first two SP3 files have same GNSS systems
    # if ~isempty(setdiff(navGNSSsystems_1, navGNSSsystems_2)):
    if list(set(navGNSSsystems_1) - set(navGNSSsystems_2)):
        print('ERROR(CombineSP3Nav): SP3 file 1 and 2 do not contain the same GNSS systems')
        sat_positions, epoch_dates, navGNSSsystems, nEpochs, epochInterval = np.nan,np.nan,np.nan,np.nan,np.nan
        success = 0
        return success
    
    
    if three_sp3_files:
        ## -- Check that last SP3 file has same GNSS systems as the others
        if list(set(navGNSSsystems_2) - set(navGNSSsystems_3)):
            print('ERROR(CombineSP3Nav): SP3 file 3 does not contain the same GNSS systems as SP3 file 1 and 2')
            sat_positions, epoch_dates, navGNSSsystems, nEpochs, epochInterval = np.nan,np.nan,np.nan,np.nan,np.nan
            success = 0
            return success

    
    
    
    nGNSSsystems = len(navGNSSsystems_1)
    max_GPS_PRN     = 36 # Max number of GPS PRN in constellation
    max_GLONASS_PRN = 36 # Max number of GLONASS PRN in constellation
    max_Galileo_PRN = 36 # Max number of Galileo PRN in constellation
    max_Beidou_PRN  = 60 # Max number of BeiDou PRN in constellation
    max_sat = np.array([max_GPS_PRN, max_GLONASS_PRN, max_Galileo_PRN, max_Beidou_PRN])
    
    navGNSSsystems = navGNSSsystems_1
    epochInterval = epochInterval_1
    
    ## -- Combine epoch dates from first and second SP3 file
    epoch_dates = np.vstack([epoch_dates_1,epoch_dates_2])  
    
    ## -- Compute total amount of epochs
    nEpochs = nEpochs_1 + nEpochs_2
    
    ## -- Initialize cell structure for storing combined satellite positions
    sat_positions = copy.deepcopy(sat_positions_1)
    
    # Combine satellite positions from first and second SP3 file
    for k in range(0,nGNSSsystems):
       curr_sys = GNSSsystems[k+1]
       len_sat = len(sat_positions_1[curr_sys])
       for ep in range(0,len(sat_positions_2[curr_sys].keys())):
           sat_positions[curr_sys].update({ep+len_sat: sat_positions_2[curr_sys][ep]})
       
    
    # If three SP3 files are present
    if three_sp3_files:
        # Add epoch dates from thirs SP3 file to the first two 
        epoch_dates = np.vstack([epoch_dates, epoch_dates_3])
        # Compute total amount of epochs
        nEpochs = nEpochs + nEpochs_3
        
        # check that last SP3 file has same interval as the others
        if epochInterval_2 != epochInterval_3:
            print('ERROR(combineSP3Nav):The second and third SP3 files do not have the same epoch interval') 
            sat_positions, epoch_dates, navGNSSsystems, nEpochs, epochInterval = np.nan,np.nan,np.nan,np.nan,np.nan
            success = 0
            return success
        
        
        # check that last SP3 file has same GNSS systems as the others
        if list(set(navGNSSsystems_2) - set(navGNSSsystems_3)):
        # if ~isempty(setdiff(navGNSSsystems_2, navGNSSsystems_3))
            print('Warning (CombineSP3Nav): SP3 file 2 and 3 do not contain the same amount of GNSS systems') 
        
        
        # Combine satellite positions from first, second and third SP3 files           
        sat_positions_dum = copy.deepcopy(sat_positions)
        for k in range(0,nGNSSsystems):
           curr_sys = GNSSsystems[k+1]
           len_sat = len(sat_positions_dum[curr_sys])
           for ep in range(0,len(sat_positions_3[curr_sys].keys())):
               sat_positions[curr_sys].update({ep+len_sat: sat_positions_3[curr_sys][ep]})

    return sat_positions, epoch_dates, navGNSSsystems, nEpochs, epochInterval, success