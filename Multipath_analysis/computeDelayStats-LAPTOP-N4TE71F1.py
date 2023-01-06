def computeDelayStats(ion_delay_phase1, multipath_range1, current_sat_elevation_angles, range1_slip_periods, LLI_slip_periods, range1_observations, tInterval):
    """
    Function that computes statistical values on estimates of multipath delay, ionospheric delay and satellite elevation angles
    
    #--------------------------------------------------------------------------------------------------------------------------
    # INPUTS 
    
    #ion_delay_phase1:      matrix containing estimates of ionospheric delays 
    #                       on the first phase signal for each PRN, at each epoch.
    
    #                       ion_delay_phase1(epoch, PRN)
    
    # multipath_range1:     matrix containing estimates of multipath delays 
    #                       on the first range signal for each PRN, at each epoch.
    
    #                       multipath_range1(epoch, PRN)
    
    # current_sat_elevation_angles: Array contaning satellite elevation angles at each
    #                               epoch, for current GNSS system. 
    
    #                               sat_elevation_angles(epoch, PRN)
    
    # range1_slip_periods_per_sat:  cell, each cell element contains one matrix for
    #                               every PRN. Each matrix contains epochs of range1
    #                               slip starts in first column and slip ends in second
    #                               column.
    
    # LLI_slip_periods_per_sat:     cell, each cell element contains one matrix for
    #                               every PRN. Each matrix contains epochs of LLI
    #                               slip starts in first column and LLI slip ends in second
    #                               column. These may be often be empty depending on if
    #                               RINEX observation files include LLI indicators or
    #                               not.
    
    
    # range1_observations:  matrix. Contains all range1 observations for all
    #                       epochs and for all SVs. 
    
    #                       range1_observations(epoch, PRN)
    
    # tInterval:            observation interval in seconds
    #--------------------------------------------------------------------------------------------------------------------------
    # OUTPUTS
    
    # mean_multipath_range1:        array. contains mean values of estimates of 
    #                               multipath delay of the first range signal,
    #                               for each satellite
    #                               
    #                               mean_multipath_range1(PRN)                 
                     
    
    # overall_mean_multipath_range1: overall mean of all estimates of 
    #                                multipath delay of the first range signal        
        
    # rms_multipath_range1:         rms values of estimates of multipath delay 
    #                               of the first range signal, for each 
    #                               satellite.
    
    #                               rms_multipath_range1(PRN)
    
    # average_rms_multipath_range1: average rms of estimates of multipath delay 
    #                               of the first range signal 
    
    # elevation_weighted_rms_multipath_range1:  rms_multipath_range1, but
    #                                           weighted based on elevation
    #                                           angles
    
    # elevation_weighted_average_rms_multipath_range1:  average_rms_multipath_range1, but
    #                                                   weighted based on elevation
    #                                                   angles
        
    # mean_ion_delay_phase1:        array. contains mean values of estimates of 
    #                               ionospheric delay on the first phase signal, 
    #                               for each satellite
    #                               
    #                               mean_ion_delay_phase1(PRN)  
    
    # overall_mean_ion_delay_phase1: overall mean of all estimates of 
    #                                ionospheric delay on the first phase signal  
    
    # mean_sat_elevation_angles:    array. contains mean elevation angle, for
    #                               each satellite.
    
    #                               mean_sat_elevation_angles(PRN)
    
    # nEstimates:                   total amount of epochs with multipath estimates, double.  
    
    # nEstimates_per_sat:           amount of epochs with multipath estimates
    #                               per sat.
    
    # nRange1Obs_Per_Sat:           array. each elements says how many range1
    #                               observations for that SV
    
    #                               nRange1Obs_Per_Sat(PRN)
    
    # nRange1Obs:                   total number of range1 observations, all SVs
    
    # range1_slip_distribution_per_sat:     cell. each element contains a structure for
    #                               each satellite. each structure has
    #                               information on number of range1 slips
    #                               distributed over groups depending on the
    #                               elevation angle of the satellite at the
    #                               time of the slip.
    
    # range1_slip_distribution:     structure. Stores information on number of 
    #                               range1 slips distributed over groups depending 
    #                               on the elevation angle of the satellite at the
    #                               time of the slip.
    
    # LLI_slip_distribution_per_sat:     cell. each element contains a structure for
    #                               each satellite. each structure has
    #                               information on number of LLI detected slips
    #                               distributed over groups depending on the
    #                               elevation angle of the satellite at the
    #                               time of the slip.
    
    # LLI_slip_distribution:        structure. Stores information on number of 
    #                               LLI detected slips distributed over groups depending 
    #                               on the elevation angle of the satellite at the
    #                               time of the slip.
    
    # combined_slip_distribution_per_sat:     cell. each element contains a structure for
    #                               each satellite. each structure has
    #                               information on number of slips detected
    #                               both by LLI and this siftware analysis,
    #                               distributed over groups depending on the
    #                               elevation angle of the satellite at the
    #                               time of the slip.
    
    # combined_slip_distribution:   structure. Stores information on number of 
    #                               slips detected both ny LLI and this software 
    #                               analysis, distributed over groups depending 
    #                               on the elevation angle of the satellite at the
    #                               time of the slip.
    #--------------------------------------------------------------------------------------------------------------------------
    """
    from numpy import mean, sqrt, sin, nan,pi,isnan, sum, concatenate, zeros, where, intersect1d, nanmean
    # nSat = length(range1_slip_periods);
    nSat = len(range1_slip_periods)

    ## set all 0 values to NaN so they are excluded from stats calculation
    # ion_delay_phase1(ion_delay_phase1==0) = NaN;
    # multipath_range1(multipath_range1==0) = NaN;
    # current_sat_elevation_angles(current_sat_elevation_angles==0) = NaN;
    ion_delay_phase1[ion_delay_phase1==0] = nan
    multipath_range1[multipath_range1==0] = nan
    current_sat_elevation_angles[current_sat_elevation_angles==0] = nan
    
    ## -- Compute mean
    # mean multipath for each satellite, excluding NaN values
    # mean_multipath_range1 = mean(multipath_range1, 'omitnan');
    mean_multipath_range1 = mean(multipath_range1)
    
    ## overall mean multipath, excluding NaN values
    # overall_mean_multipath_range1 = mean(mean_multipath_range1, 'omitnan');
    overall_mean_multipath_range1 = mean(mean_multipath_range1)
    
    ## RMS multipath of each satellite, excluding NaN
    rms_multipath_range1 = sqrt(mean(multipath_range1*multipath_range1)) # matlab uses dot notation here
    
    ## Average RMS multipath, excluding NaN
    average_rms_multipath_range1 = sqrt(mean(multipath_range1*multipath_range1)) # matlab uses dot notation here
    
    
    ##
    ## Weighted RMS multipath 
    # weights = current_sat_elevation_angles;
    # crit_weight = 4*sin(30*pi/180)^2; # 1
    # weights = 4*sin(weights*pi/180).^2; 
    # weights(weights > crit_weight) = 1;
    # elevation_weighted_multipath_range1 = multipath_range1.*weights;
    weights = current_sat_elevation_angles
    crit_weight = 4*sin(30*pi/180)**2 # 1
    weights = 4*sin(weights*pi/180)**2 
    weights[weights > crit_weight] = 1
    elevation_weighted_multipath_range1 = multipath_range1*weights
    
    ## RMS multipath of each satellite, excluding NaN
    # elevation_weighted_rms_multipath_range1 = sqrt(mean((elevation_weighted_multipath_range1.*elevation_weighted_multipath_range1), 'omitnan'));
    elevation_weighted_rms_multipath_range1 = sqrt(nanmean(elevation_weighted_multipath_range1*elevation_weighted_multipath_range1))
    
    ## Average RMS multipath, excluding NaN
    # elevation_weighted_average_rms_multipath_range1 = sqrt(mean((elevation_weighted_multipath_range1.*elevation_weighted_multipath_range1), 'all', 'omitnan'));
    elevation_weighted_average_rms_multipath_range1 = sqrt(nanmean(elevation_weighted_multipath_range1*elevation_weighted_multipath_range1))
    
    ## --Ionosphere 
    ## mean ionospheric delay for each satellite, excluding NaN
    # mean_ion_delay_phase1 = mean(ion_delay_phase1, 'omitnan');
    mean_ion_delay_phase1 = mean(ion_delay_phase1)
    
    ## Overall mean ionospheric delay, excluding NaN
    # overall_mean_ion_delay_phase1 = mean(mean_ion_delay_phase1, 'omitnan');
    overall_mean_ion_delay_phase1 = mean(mean_ion_delay_phase1)
    
    ##
    ## Average elevation angle for each satellite, excluding NaN
    # dummy = (range1_observations~= 0 & ~isnan(range1_observations));
    # obs_elevations = current_sat_elevation_angles.*dummy;
    # obs_elevations(obs_elevations==0) = NaN;
    # mean_sat_elevation_angles = mean(obs_elevations,'omitnan');
    dumm1 = (range1_observations!= 0)*1 # multiplying with 1 to get from True/False -> 1/0
    dumm2 = (~isnan(range1_observations))*1
    dummy = (dumm1 & dumm2)
    obs_elevations = current_sat_elevation_angles*dummy
    obs_elevations[obs_elevations==0] =nan
    mean_sat_elevation_angles = mean(obs_elevations)
    
    
    ##
    ## Amount of epochs with estimates
    # nEstimates = sum((multipath_range1 ~= 0 & ~isnan(multipath_range1) ),'all');
    # nEstimates_per_sat = sum((multipath_range1 ~= 0 & ~isnan(multipath_range1) ));
    nEstimates = sum(multipath_range1 != 0 & ~isnan(multipath_range1))
    nEstimates_per_sat = sum(multipath_range1 != 0 & ~isnan(multipath_range1))
        
    ## -- Slip distribution
    # distribution of slips over elevation angles AND
    # distribution of LLI slips over elevation angles
    
    # range1_slip_distribution_per_sat = cell(1, nSat);
    # range1_slip_distribution         = struct;
    range1_slip_distribution_per_sat = {}
    range1_slip_distribution         = {}
    
    range1_slip_distribution['n_slips_0_10']   = 0
    range1_slip_distribution['n_slips_10_20']  = 0
    range1_slip_distribution['n_slips_20_30']  = 0
    range1_slip_distribution['n_slips_30_40']  = 0
    range1_slip_distribution['n_slips_40_50']  = 0
    range1_slip_distribution['n_slips_over50'] = 0
    range1_slip_distribution['n_slips_NaN']    = 0
    range1_slip_distribution['n_slips_Tot']    = 0
    
    
    
    # LLI_slip_distribution_per_sat = cell(1, nSat);
    # LLI_slip_distribution         = struct;
    LLI_slip_distribution_per_sat = {}
    LLI_slip_distribution         = {}
    
    LLI_slip_distribution['n_slips_0_10']   = 0
    LLI_slip_distribution['n_slips_10_20']  = 0
    LLI_slip_distribution['n_slips_20_30']  = 0
    LLI_slip_distribution['n_slips_30_40']  = 0
    LLI_slip_distribution['n_slips_40_50']  = 0
    LLI_slip_distribution['n_slips_over50'] = 0
    LLI_slip_distribution['n_slips_NaN']    = 0
    LLI_slip_distribution['n_slips_Tot']    = 0
    
    

    
    # combined_slip_distribution_per_sat = cell(1, nSat);
    # combined_slip_distribution         = struct;
    combined_slip_distribution_per_sat = {}
    combined_slip_distribution         = {}
    
    combined_slip_distribution['n_slips_0_10']   = 0
    combined_slip_distribution['n_slips_10_20']  = 0
    combined_slip_distribution['n_slips_20_30']  = 0
    combined_slip_distribution['n_slips_30_40']  = 0
    combined_slip_distribution['n_slips_40_50']  = 0
    combined_slip_distribution['n_slips_over50'] = 0
    combined_slip_distribution['n_slips_NaN']    = 0
    combined_slip_distribution['n_slips_Tot']    = 0
    

        
    # for i = 1:nSat
    #     # create struct for current sat
    #    range1_slip_distribution_per_sat{i} = struct;
    #    slip_epochs = [];
    #    [nSlipPeriods, ~] = size(range1_slip_periods{i});
       
    for i in range(0,nSat):
        ## Create struct for current sat
        range1_slip_distribution_per_sat[i] = {}
        slip_epochs = []
        nSlipPeriods, _= range1_slip_periods[i].shape
        
       ## Get all slip epochs of current sat.
        # for j = 1:nSlipPeriods
        #     # if slip period is shorter than 60 seconds
        #     slip_epochs = [slip_epochs; range1_slip_periods{i}(j,1)];
        for j in range(0,nSlipPeriods):
            ## if slip period is shorter than 60 seconds
            slip_epochs = concatenate((slip_epochs, range1_slip_periods[i][j,0]))
        
        ## Get elevation angles for every slip of current sat
        # slip_epoch_elevation_angles = current_sat_elevation_angles(slip_epochs, i);
        slip_epoch_elevation_angles = current_sat_elevation_angles[slip_epochs, i]
        
        ## -- Store number of slips into groups of their elevation angles. Groups are: 0-10, 10-20, 20-30, 30-40, 40-50, >50, and NaN
        # range1_slip_distribution_per_sat{i}.n_slips_0_10    = length(slip_epoch_elevation_angles(slip_epoch_elevation_angles>=0 & slip_epoch_elevation_angles<10));
        # range1_slip_distribution_per_sat{i}.n_slips_10_20   = length(slip_epoch_elevation_angles(slip_epoch_elevation_angles>=10 & slip_epoch_elevation_angles<20));
        # range1_slip_distribution_per_sat{i}.n_slips_20_30   = length(slip_epoch_elevation_angles(slip_epoch_elevation_angles>=20 & slip_epoch_elevation_angles<30));
        # range1_slip_distribution_per_sat{i}.n_slips_30_40   = length(slip_epoch_elevation_angles(slip_epoch_elevation_angles>=30 & slip_epoch_elevation_angles<40));
        # range1_slip_distribution_per_sat{i}.n_slips_40_50   = length(slip_epoch_elevation_angles(slip_epoch_elevation_angles>=40 & slip_epoch_elevation_angles<50));
        # range1_slip_distribution_per_sat{i}.n_slips_over50  = length(slip_epoch_elevation_angles(slip_epoch_elevation_angles>=50));
        # range1_slip_distribution_per_sat{i}.n_slips_NaN     = length(slip_epoch_elevation_angles(isnan(slip_epoch_elevation_angles)));
        # range1_slip_distribution_per_sat{i}.n_slips_Tot     = length(slip_epoch_elevation_angles);
        
        # range1_slip_distribution.n_slips_0_10   = range1_slip_distribution.n_slips_0_10     + range1_slip_distribution_per_sat{i}.n_slips_0_10;
        # range1_slip_distribution.n_slips_10_20  = range1_slip_distribution.n_slips_10_20    + range1_slip_distribution_per_sat{i}.n_slips_10_20;
        # range1_slip_distribution.n_slips_20_30  = range1_slip_distribution.n_slips_20_30    + range1_slip_distribution_per_sat{i}.n_slips_20_30;
        # range1_slip_distribution.n_slips_30_40  = range1_slip_distribution.n_slips_30_40    + range1_slip_distribution_per_sat{i}.n_slips_30_40;
        # range1_slip_distribution.n_slips_40_50  = range1_slip_distribution.n_slips_40_50    + range1_slip_distribution_per_sat{i}.n_slips_40_50;
        # range1_slip_distribution.n_slips_over50 = range1_slip_distribution.n_slips_over50   + range1_slip_distribution_per_sat{i}.n_slips_over50;
        # range1_slip_distribution.n_slips_NaN    = range1_slip_distribution.n_slips_NaN      + range1_slip_distribution_per_sat{i}.n_slips_NaN;
        # range1_slip_distribution.n_slips_Tot    = range1_slip_distribution.n_slips_Tot      + range1_slip_distribution_per_sat{i}.n_slips_Tot;
        range1_slip_distribution_per_sat[i]['n_slips_0_10']    = len(slip_epoch_elevation_angles[slip_epoch_elevation_angles >=0  and slip_epoch_elevation_angles <10])
        range1_slip_distribution_per_sat[i]['n_slips_10_20']   = len(slip_epoch_elevation_angles[slip_epoch_elevation_angles >=10 and slip_epoch_elevation_angles <20])
        range1_slip_distribution_per_sat[i]['n_slips_20_30']   = len(slip_epoch_elevation_angles[slip_epoch_elevation_angles >=20 and slip_epoch_elevation_angles <30])
        range1_slip_distribution_per_sat[i]['n_slips_30_40']   = len(slip_epoch_elevation_angles[slip_epoch_elevation_angles >=30 and slip_epoch_elevation_angles <40])
        range1_slip_distribution_per_sat[i]['n_slips_40_50']   = len(slip_epoch_elevation_angles[slip_epoch_elevation_angles >=40 and slip_epoch_elevation_angles <50])
        range1_slip_distribution_per_sat[i]['n_slips_over50']  = len(slip_epoch_elevation_angles[slip_epoch_elevation_angles >=50])
        range1_slip_distribution_per_sat[i]['n_slips_NaN']     = len(slip_epoch_elevation_angles[isnan(slip_epoch_elevation_angles)])
        range1_slip_distribution_per_sat[i]['n_slips_Tot']     = len(slip_epoch_elevation_angles)
        
        range1_slip_distribution['n_slips_0_10']   = range1_slip_distribution['n_slips_0_10']     + range1_slip_distribution_per_sat[i]['n_slips_0_10']
        range1_slip_distribution['n_slips_10_20']  = range1_slip_distribution['n_slips_10_20']    + range1_slip_distribution_per_sat[i]['n_slips_10_20']
        range1_slip_distribution['n_slips_20_30']  = range1_slip_distribution['n_slips_20_30']    + range1_slip_distribution_per_sat[i]['n_slips_20_30']
        range1_slip_distribution['n_slips_30_40']  = range1_slip_distribution['n_slips_30_40']    + range1_slip_distribution_per_sat[i]['n_slips_30_40']
        range1_slip_distribution['n_slips_40_50']  = range1_slip_distribution['n_slips_40_50']    + range1_slip_distribution_per_sat[i]['n_slips_40_50']
        range1_slip_distribution['n_slips_over50'] = range1_slip_distribution['n_slips_over50']   + range1_slip_distribution_per_sat[i]['n_slips_over50']
        range1_slip_distribution['n_slips_NaN']    = range1_slip_distribution['n_slips_NaN']      + range1_slip_distribution_per_sat[i]['n_slips_NaN']
        range1_slip_distribution['n_slips_Tot']    = range1_slip_distribution['n_slips_Tot']      + range1_slip_distribution_per_sat[i]['n_slips_Tot']        
        
        ## Create a dictionary for current sat
        # LLI_slip_distribution_per_sat{i} = struct;
        # LLI_slip_epochs = [];
        # [nSlipPeriods, ~] = size(LLI_slip_periods{i});
        
        LLI_slip_distribution_per_sat[i] = {}
        LLI_slip_epochs = []
        nSlipPeriods, _ = LLI_slip_periods[i].shape
        
        ## -- Get all LLI slip epochs of current sat.
        # for j = 1:nSlipPeriods
        #     LLI_slip_epochs = [LLI_slip_epochs; LLI_slip_periods{i}(j,1)];
        for j in range(0,nSlipPeriods):
            LLI_slip_epochs = concatenate((LLI_slip_epochs, LLI_slip_periods[i][j,0]))
            
        ## Get elevation angles for every LLI slip of current sat
        # LLI_slip_epoch_elevation_angles = current_sat_elevation_angles(LLI_slip_epochs, i);
        LLI_slip_epoch_elevation_angles = current_sat_elevation_angles[LLI_slip_epochs, i]
        
        ## -- Store number of LLI slips into groups of their elevation angles. Groups are: 0-10, 10-20, 20-30, 30-40, 40-50, >50, and NaN
        # LLI_slip_distribution_per_sat{i}.n_slips_0_10    = length(LLI_slip_epoch_elevation_angles(LLI_slip_epoch_elevation_angles>=0 & LLI_slip_epoch_elevation_angles<10));
        # LLI_slip_distribution_per_sat{i}.n_slips_10_20   = length(LLI_slip_epoch_elevation_angles(LLI_slip_epoch_elevation_angles>=10 & LLI_slip_epoch_elevation_angles<20));
        # LLI_slip_distribution_per_sat{i}.n_slips_20_30   = length(LLI_slip_epoch_elevation_angles(LLI_slip_epoch_elevation_angles>=20 & LLI_slip_epoch_elevation_angles<30));
        # LLI_slip_distribution_per_sat{i}.n_slips_30_40   = length(LLI_slip_epoch_elevation_angles(LLI_slip_epoch_elevation_angles>=30 & LLI_slip_epoch_elevation_angles<40));
        # LLI_slip_distribution_per_sat{i}.n_slips_40_50   = length(LLI_slip_epoch_elevation_angles(LLI_slip_epoch_elevation_angles>=40 & LLI_slip_epoch_elevation_angles<50));
        # LLI_slip_distribution_per_sat{i}.n_slips_over50  = length(LLI_slip_epoch_elevation_angles(LLI_slip_epoch_elevation_angles>=50));
        # LLI_slip_distribution_per_sat{i}.n_slips_NaN     = length(LLI_slip_epoch_elevation_angles(isnan(LLI_slip_epoch_elevation_angles)));
        # LLI_slip_distribution_per_sat{i}.n_slips_Tot     = length(LLI_slip_epoch_elevation_angles);
        
        # LLI_slip_distribution.n_slips_0_10   = LLI_slip_distribution.n_slips_0_10     + LLI_slip_distribution_per_sat{i}.n_slips_0_10;
        # LLI_slip_distribution.n_slips_10_20  = LLI_slip_distribution.n_slips_10_20    + LLI_slip_distribution_per_sat{i}.n_slips_10_20;
        # LLI_slip_distribution.n_slips_20_30  = LLI_slip_distribution.n_slips_20_30    + LLI_slip_distribution_per_sat{i}.n_slips_20_30;
        # LLI_slip_distribution.n_slips_30_40  = LLI_slip_distribution.n_slips_30_40    + LLI_slip_distribution_per_sat{i}.n_slips_30_40;
        # LLI_slip_distribution.n_slips_40_50  = LLI_slip_distribution.n_slips_40_50    + LLI_slip_distribution_per_sat{i}.n_slips_40_50;
        # LLI_slip_distribution.n_slips_over50 = LLI_slip_distribution.n_slips_over50   + LLI_slip_distribution_per_sat{i}.n_slips_over50;
        # LLI_slip_distribution.n_slips_NaN    = LLI_slip_distribution.n_slips_NaN      + LLI_slip_distribution_per_sat{i}.n_slips_NaN;
        # LLI_slip_distribution.n_slips_Tot    = LLI_slip_distribution.n_slips_Tot      + LLI_slip_distribution_per_sat{i}.n_slips_Tot;
        LLI_slip_distribution_per_sat[i]['n_slips_0_10']    = len(LLI_slip_epoch_elevation_angles(LLI_slip_epoch_elevation_angles >=0  and LLI_slip_epoch_elevation_angles <10))
        LLI_slip_distribution_per_sat[i]['n_slips_10_20']   = len(LLI_slip_epoch_elevation_angles(LLI_slip_epoch_elevation_angles >=10 and LLI_slip_epoch_elevation_angles <20))
        LLI_slip_distribution_per_sat[i]['n_slips_20_30']   = len(LLI_slip_epoch_elevation_angles(LLI_slip_epoch_elevation_angles >=20 and LLI_slip_epoch_elevation_angles <30))
        LLI_slip_distribution_per_sat[i]['n_slips_30_40']   = len(LLI_slip_epoch_elevation_angles(LLI_slip_epoch_elevation_angles >=30 and LLI_slip_epoch_elevation_angles <40))
        LLI_slip_distribution_per_sat[i]['n_slips_40_50']   = len(LLI_slip_epoch_elevation_angles(LLI_slip_epoch_elevation_angles >=40 and LLI_slip_epoch_elevation_angles <50))
        LLI_slip_distribution_per_sat[i]['n_slips_over50']  = len(LLI_slip_epoch_elevation_angles(LLI_slip_epoch_elevation_angles >=50))
        LLI_slip_distribution_per_sat[i]['n_slips_NaN']     = len(LLI_slip_epoch_elevation_angles(isnan(LLI_slip_epoch_elevation_angles)))
        LLI_slip_distribution_per_sat[i]['n_slips_Tot']     = len(LLI_slip_epoch_elevation_angles)
        
        LLI_slip_distribution['n_slips_0_10']   = LLI_slip_distribution['n_slips_0_10']     + LLI_slip_distribution_per_sat[i]['n_slips_0_10']
        LLI_slip_distribution['n_slips_10_20']  = LLI_slip_distribution['n_slips_10_20']    + LLI_slip_distribution_per_sat[i]['n_slips_10_20']
        LLI_slip_distribution['n_slips_20_30']  = LLI_slip_distribution['n_slips_20_30']    + LLI_slip_distribution_per_sat[i]['n_slips_20_30']
        LLI_slip_distribution['n_slips_30_40']  = LLI_slip_distribution['n_slips_30_40']    + LLI_slip_distribution_per_sat[i]['n_slips_30_40']
        LLI_slip_distribution['n_slips_40_50']  = LLI_slip_distribution['n_slips_40_50']    + LLI_slip_distribution_per_sat[i]['n_slips_40_50']
        LLI_slip_distribution['n_slips_over50'] = LLI_slip_distribution['n_slips_over50']   + LLI_slip_distribution_per_sat[i]['n_slips_over50']
        LLI_slip_distribution['n_slips_NaN']    = LLI_slip_distribution['n_slips_NaN']      + LLI_slip_distribution_per_sat[i]['n_slips_NaN']
        LLI_slip_distribution['n_slips_Tot']    = LLI_slip_distribution['n_slips_Tot']      + LLI_slip_distribution_per_sat[i]['n_slips_Tot']
        
         ## create dict for current sat
        combined_slip_distribution_per_sat[i] = {}
        
        ## -- Get all combined slips for current satellite
        # combined_slip_epochs = intersect(slip_epochs, LLI_slip_epochs);
        combined_slip_epochs = intersect1d(slip_epochs, LLI_slip_epochs)
        
        # get elevation angles for every combined slip of current sat
        combined_slip_epoch_elevation_angles = current_sat_elevation_angles[combined_slip_epochs, i]
        
        # Store number of LLI slips into groups of their elevation angles. Groups are: 0-10, 10-20, 20-30, 30-40, 40-50, >50, and NaN
        combined_slip_distribution_per_sat[i]['n_slips_0_10']    = len(combined_slip_epoch_elevation_angles(combined_slip_epoch_elevation_angles >=0 and combined_slip_epoch_elevation_angles <10))
        combined_slip_distribution_per_sat[i]['n_slips_10_20']   = len(combined_slip_epoch_elevation_angles(combined_slip_epoch_elevation_angles >=10 and combined_slip_epoch_elevation_angles <20))
        combined_slip_distribution_per_sat[i]['n_slips_20_30']   = len(combined_slip_epoch_elevation_angles(combined_slip_epoch_elevation_angles >=20 and combined_slip_epoch_elevation_angles <30))
        combined_slip_distribution_per_sat[i]['n_slips_30_40']   = len(combined_slip_epoch_elevation_angles(combined_slip_epoch_elevation_angles >=30 and combined_slip_epoch_elevation_angles <40))
        combined_slip_distribution_per_sat[i]['n_slips_40_50']   = len(combined_slip_epoch_elevation_angles(combined_slip_epoch_elevation_angles >=40 and combined_slip_epoch_elevation_angles <50))
        combined_slip_distribution_per_sat[i]['n_slips_over50']  = len(combined_slip_epoch_elevation_angles(combined_slip_epoch_elevation_angles >=50))
        combined_slip_distribution_per_sat[i]['n_slips_NaN']     = len(combined_slip_epoch_elevation_angles(isnan(combined_slip_epoch_elevation_angles)))
        combined_slip_distribution_per_sat[i]['n_slips_Tot']     = len(combined_slip_epoch_elevation_angles)
        
        combined_slip_distribution['n_slips_0_10']   = combined_slip_distribution['n_slips_0_10']     + combined_slip_distribution_per_sat[i]['n_slips_0_10']
        combined_slip_distribution['n_slips_10_20']  = combined_slip_distribution['n_slips_10_20']    + combined_slip_distribution_per_sat[i]['n_slips_10_20']
        combined_slip_distribution['n_slips_20_30']  = combined_slip_distribution['n_slips_20_30']    + combined_slip_distribution_per_sat[i]['n_slips_20_30']
        combined_slip_distribution['n_slips_30_40']  = combined_slip_distribution['n_slips_30_40']    + combined_slip_distribution_per_sat[i]['n_slips_30_40']
        combined_slip_distribution['n_slips_40_50']  = combined_slip_distribution['n_slips_40_50']    + combined_slip_distribution_per_sat[i]['n_slips_40_50']
        combined_slip_distribution['n_slips_over50'] = combined_slip_distribution['n_slips_over50']   + combined_slip_distribution_per_sat[i]['n_slips_over50']
        combined_slip_distribution['n_slips_NaN']    = combined_slip_distribution['n_slips_NaN']      + combined_slip_distribution_per_sat[i]['n_slips_NaN']
        combined_slip_distribution['n_slips_Tot']    = combined_slip_distribution['n_slips_Tot']      + combined_slip_distribution_per_sat[i]['n_slips_Tot']

    
    ## Amount of range1 observations
    # nRange1Obs_Per_Sat = zeros(1, nSat);
    # for i =1:nSat
    #     first_obs_epoch = find(range1_observations(:, i) ~= 0, 1, 'first');
    #     last_obs_epoch = find(range1_observations(:, i) ~= 0, 1, 'last');
    #    nRange1Obs_Per_Sat(i) = length(range1_observations(first_obs_epoch:last_obs_epoch, i));
    # end
    # nRange1Obs = sum(nRange1Obs_Per_Sat);
    nRange1Obs_Per_Sat = zeros([1, nSat])
    for i in range(0,nSat):
        try:
            first_obs_epoch = where(range1_observations[:, i] != 0)[0][0]
        except:
            first_obs_epoch = where(range1_observations[:, i] != 0)[0] # if empty
            
        try:
            last_obs_epoch = where(range1_observations[:, i] != 0)[0][-1]
        except:
            last_obs_epoch = where(range1_observations[:, i] != 0)[0] # if empty
        
        nRange1Obs_Per_Sat[i] = len(range1_observations[first_obs_epoch:last_obs_epoch, i])
    
    nRange1Obs = sum(nRange1Obs_Per_Sat)

    return mean_multipath_range1, overall_mean_multipath_range1, rms_multipath_range1, average_rms_multipath_range1,\
        mean_ion_delay_phase1, overall_mean_ion_delay_phase1, mean_sat_elevation_angles, nEstimates, nEstimates_per_sat,\
        nRange1Obs_Per_Sat, nRange1Obs, range1_slip_distribution_per_sat, range1_slip_distribution, LLI_slip_distribution_per_sat, LLI_slip_distribution,\
        combined_slip_distribution_per_sat, combined_slip_distribution, elevation_weighted_rms_multipath_range1, elevation_weighted_average_rms_multipath_range1