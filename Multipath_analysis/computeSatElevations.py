def computeSatElevations(GNSS_SVs, GNSSsystems, approxPosition,\
    nepochs, time_epochs, max_sat, sp3_nav_filename_1, sp3_nav_filename_2, sp3_nav_filename_3):
    """
    # Function that computes the satellite elevation angles of all satellites
    # of all GNSS systems at each epoch of the observations period.
    #--------------------------------------------------------------------------------------------------------------------------
    # INPUTS
    
    # GNSS_SVs:                 Cell containing number of satellites with 
    #                           observations and PRN for each satellite, 
    #                           for each GNSS system
    
    #                           GNSS_SVs{GNSSsystemIndex}(epoch,j)     
    #                                       j=1: number of observed GPS satellites
    #                                       j>1: PRN of observed satellites
    
    # GNSSsystems:              cell array containing different GNSS 
    #                           systems included in RINEX file. Elements are strings. 
    #                           Must be either "G","R","E" or "C".  
    
    # approxPosition:           array containing approximate position from rinex
    #                           observation file header. [X, Y, Z]
    
    # nepochs:                  number of epochs with observations in 
    #                           rinex observation file.
    
    # time_epochs:              matrix conatining gps-week and "time of week" 
    #                           for each epoch
    #                           time_epochs(epoch,i),   i=1: week
    #                                                   i=2: time-of-week in seconds (tow)
    
    # max_sat:                  array that stores max satellite PRN number for 
    #                           each of the GNSS systems. Follows same order as GNSSsystems
    
    # nav_filename:             string, path and filename of rinex3.xx navigation file.
    
    # almanac_nav_filename:     string, path and filename of sen almanac
    #                           navigation filename for GLONASS
    #--------------------------------------------------------------------------------------------------------------------------
    # OUTPUTS
    
    # sat_elevation_angles:     Cell contaning satellite elevation angles at each
    #                           epoch, for each GNSS system. 
    #                           sat_elevation_angles{GNSSsystemIndex}(epoch, PRN)
    #--------------------------------------------------------------------------------------------------------------------------
    """
    
    # import sys
    # sys.path.append(r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Read SP3')
    from read_SP3Nav import readSP3Nav
    from combineSP3Nav import combineSP3Nav
    from get_elevation_angle import get_elevation_angle
    from gpstime2date import gpstime2date
    from tqdm import tqdm
    import numpy as np
    
    nGNSSsystems = len(GNSSsystems);
    # sat_elevation_angles = cell(nGNSSsystems,1);
    sat_elevation_angles = {}
    sat_azimut_angles    = {}
    
    if all(approxPosition == 0):
        print('ERROR(computeSatElevations): The approximate receiver position is missing.\n',\
            'Please chack that APPROX POSITION XYZ is in header of Rinex file.\n',\
            'Elevation angles will not be computed\n\n.')
        # return
    
    ## -- Testing whether two and tree sp3 files are defined
    two_sp3_files = 0
    three_sp3_files = 0
    
    if sp3_nav_filename_2 != "":
        two_sp3_files = 1
        if sp3_nav_filename_3 != "":
          three_sp3_files = 1
       
        
    
    ## ---  Read first SP3 file
    sat_positions_1, epoch_dates_1, navGNSSsystems_1, nEpochs_1, epochInterval_1, success = readSP3Nav(sp3_nav_filename_1)
    # if two SP3 files inputted by user, read second SP3 file
    if two_sp3_files:
        sat_positions_2, epoch_dates_2, navGNSSsystems_2, nEpochs_2, epochInterval_2, success = readSP3Nav(sp3_nav_filename_2)
    else:
        # sat_positions_2, epoch_dates_2, navGNSSsystems_2, nEpochs_2, epochInterval_2 = deal(NaN)
        sat_positions_2, epoch_dates_2, navGNSSsystems_2, nEpochs_2, epochInterval_2 = np.nan,np.nan,np.nan,np.nan,np.nan
    
    # if three SP3 files inputted by user, read third SP3 file
    if three_sp3_files:
        sat_positions_3, epoch_dates_3, navGNSSsystems_3, nEpochs_3, epochInterval_3, success = readSP3Nav(sp3_nav_filename_3)
    else:
        # sat_positions_3, epoch_dates_3, navGNSSsystems_3, nEpochs_3, epochInterval_3 = deal(NaN);
        sat_positions_3, epoch_dates_3, navGNSSsystems_3, nEpochs_3, epochInterval_3 = np.nan,np.nan,np.nan,np.nan,np.nan
    
    ## Combine data from different SP3 files 
    if two_sp3_files:
        sat_positions, epoch_dates, navGNSSsystems, nEpochs, epochInterval, success = combineSP3Nav(three_sp3_files,\
            sat_positions_1, epoch_dates_1, navGNSSsystems_1, nEpochs_1, epochInterval_1,\
            sat_positions_2, epoch_dates_2, navGNSSsystems_2, nEpochs_2, epochInterval_2,\
            sat_positions_3, epoch_dates_3, navGNSSsystems_3, nEpochs_3, epochInterval_3,GNSSsystems)
    else:
        sat_positions = sat_positions_1;
        epoch_dates = epoch_dates_1;
        navGNSSsystems = navGNSSsystems_1;
        nEpochs = nEpochs_1;
        epochInterval = epochInterval_1;
    
    
    # sat_positions, epoch_dates_, navGNSSsystems, nEpochs, epochInterval, success = readSP3Nav(sp3_nav_filename_1)
    """
     creating a modified version of GNSS_SVs. Instead of only giving
     satellites with observations for each epoch, FirstLastObsEpochOverview
     gives every satellite who's first and last observations are before and
     after current epoch. Hence if a temporary period of time where a
     satellite does not have observations occurs the elevation angle of that
     satellite will still be computed.
    """
    # FirstLastObsEpochOverview = cell(1, nGNSSsystems);
    FirstLastObsEpochOverview = {}
    
    
    
    for i in range(0,nGNSSsystems):
       curr_sys = GNSSsystems[i+1]
       # nepochs_current_sys, max_sat_current_sys = size(GNSS_SVs{i})
       nepochs_current_sys, max_sat_current_sys = GNSS_SVs[curr_sys].shape
       # max_sat_current_sys = max_sat_current_sys - 1;
       # FirstLastObsEpochOverview{i} = zeros(nepochs_current_sys, max_sat_current_sys);
       FirstLastObsEpochOverview[i] = np.zeros([nepochs_current_sys, max_sat_current_sys])
       
       for PRN_ in range(0,max_sat_current_sys-1):
          PRN = PRN_ + 1 # add 1 because of python null-indexed
          # logical, 1 for every epoch with observation for this PRN
          # dummy = sum(GNSS_SVs{i}(:, 2:end)==PRN,2);
          dumm = GNSS_SVs[curr_sys][:, 1::] ==PRN
          dummy = np.sum(dumm.astype(np.int8),axis=1).reshape(len(dumm.astype(np.int8)),1)
          # find first and last epoch with observation for this PRN
          # firstObsEpoch_current_sys = find(dummy==1, 1, 'first');
          # lastObsEpoch_current_sys = find(dummy==1, 1, 'last');
          if not all(dummy ==0):
              firstObsEpoch_current_sys = np.where(dummy == 1)[0][0]
              lastObsEpoch_current_sys = np.where(dummy == 1)[0][-1]
              FirstLastObsEpochOverview[i][firstObsEpoch_current_sys:lastObsEpoch_current_sys, PRN] = PRN
          else:
              firstObsEpoch_current_sys = np.nan
              lastObsEpoch_current_sys = np.nan
              FirstLastObsEpochOverview[i][:, PRN] = 0
          
          # FirstLastObsEpochOverview{i}(firstObsEpoch_current_sys:lastObsEpoch_current_sys, PRN) = PRN;
          # FirstLastObsEpochOverview[i][firstObsEpoch_current_sys:lastObsEpoch_current_sys, PRN] = PRN
    
    # Initialize progress bar
    # wbar = waitbar(0, 'INFO(computeSatElevations): Satellite elevation angles are being calculated. Please wait');
    # with tqdm(total=1000,desc ='INFO(computeSatElevations): Satellite elevation angles are being calculated. Please wait' , position=0, leave=True) as wbar:
    # overview of satellites without navigation data in sp3 file
    # satMissingData = {};
    satMissingData = []
    # for k = 1:nGNSSsystems:
    for k in tqdm(range(0,nGNSSsystems),desc='Looping thorugh the systems',position=0,leave=True):
       # Initialize data matrix for current GNSSsystem
       sat_elevation_angles[k] = np.zeros([int(nepochs), int(max_sat[k])+1]) # added +1 20.11.2022
       sat_azimut_angles[k]    = np.zeros([int(nepochs), int(max_sat[k])+1]) # added +1 20.11.2022
       sys = GNSSsystems[k+1]
       
       if sys in navGNSSsystems: # lat til tqdm
           for epoch in tqdm(range(0,nepochs),desc='Satellite elevation angles are being calculated for system %s of %s' %(k+1, nGNSSsystems),position=0, leave=False):
              # Update progress bar 
              # if mod(epoch, 500) == 0
              #      waitbar((epoch+nepochs*(k-1))/(nepochs*nGNSSsystems), wbar, ...
              #          'INFO(computeSatElevations): Satellite elevation angles are being calculated. Please wait');
               # if np.mod(epoch, 500) == 0:
               #      # waitbar((epoch+nepochs*(k-1))/(nepochs*nGNSSsystems), wbar,'INFO(computeSatElevations): Satellite elevation angles are being calculated. Please wait');
               #      wbar.update(10)
                    # wbar.update((epoch+nepochs*(k-1))/(nepochs*nGNSSsystems))
             
              # GPS Week and time of week of current epoch
               # week = time_epochs(epoch,1);
               # tow  = time_epochs(epoch,2);
               week = time_epochs[epoch,0]
               tow  = time_epochs[epoch,1]
               # Satellites in current epoch that should have elevation computed
               # SVs = nonzeros(FirstLastObsEpochOverview{k}(epoch, :));
               SVs = np.nonzero(FirstLastObsEpochOverview[k][epoch, :])[0]
               n_sat = len(SVs)
    
              # for sat = 1:n_sat
               for sat in range(0,n_sat):
                   # Get satellite elevation angle of current sat at current epoch
                   # elevation_angle, missing_nav_data = get_elevation_angle(GNSSsystems{k}, SVs(sat), week, tow, sat_positions, nEpochs,\
                   #    epoch_dates, epochInterval, navGNSSsystems, approxPosition);
                   # elevation_angle, missing_nav_data = get_elevation_angle(GNSSsystems[k+1], SVs[sat], week, tow, sat_positions, nEpochs,\
                   #    epoch_dates, epochInterval, navGNSSsystems, approxPosition);
                   elevation_angle, azimut_angle, missing_nav_data = get_elevation_angle(GNSSsystems[k+1], SVs[sat], week, tow, sat_positions, nEpochs,\
                      epoch_dates, epochInterval, navGNSSsystems, approxPosition);
    
                  # sat_elevation_angles{k}(epoch,SVs(sat)) = elevation_angle;
                   sat_elevation_angles[k][epoch,SVs[sat]] = elevation_angle
                   sat_azimut_angles[k][epoch,SVs[sat]] = azimut_angle
                   if missing_nav_data:
                       # combine PRN number and GNSS system char
                       # SVN = strcat(GNSSsystems{k},num2str(SVs(sat)));
                       SVN = str(GNSSsystems[k]) + str(SVs(sat))
                       # if ~any(strcmp(SVN, satMissingData))
                       if not any(SVN in satMissingData):
                           # satMissingData = [satMissingData, SVN]; 
                           satMissingData.append(SVN) 
    
    
    # if ~isempty(satMissingData)!:
    if satMissingData:
       print('INFO(computeSatElevations): The following satellites had missing orbit data in SP3 file.',\
       'Their elevation angles were set to 0 for these epochs') 
       print(satMissingData)
    
        
        # close(wbar)
    print('INFO(computeSatElevations): Satellite elevation angles have been computed')

    return sat_elevation_angles,sat_azimut_angles