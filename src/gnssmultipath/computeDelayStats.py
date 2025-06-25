"""
This module is computing the stats.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""

import numpy as np
import warnings
warnings.filterwarnings("ignore")

def computeDelayStats(ion_delay_phase1, multipath_range1, current_sat_elevation_angles, range1_slip_periods, ambiguity_slip_periods ,LLI_slip_periods, range1_observations, tInterval):
    """
    Function that computes statistical values on estimates of multipath delay, ionospheric delay and satellite elevation angles


     INPUTS:
    -------

    ion_delay_phase1:      matrix containing estimates of ionospheric delays
                           on the first phase signal for each PRN, at each epoch.

                           ion_delay_phase1(epoch, PRN)

     multipath_range1:     matrix containing estimates of multipath delays
                           on the first range signal for each PRN, at each epoch.

                           multipath_range1(epoch, PRN)

     current_sat_elevation_angles: Array contaning satellite elevation angles at each
                                   epoch, for current GNSS system.

                                   sat_elevation_angles(epoch, PRN)

    range1_slip_periods:           dict, each cell element contains one matrix for
                                   every PRN. Each matrix contains epochs of range1
                                   slip starts in first column and slip ends in second
                                   column.

    ambiguity_slip_periods:       dict, each element contains one array for
                                  every PRN. Each array contains epochs of detected cycle slips (ion & code phase)
                                  starts in first column and slip ends in second
                                  column.



     LLI_slip_periods_per_sat:     dicct, each cell element contains one matrix for
                                   every PRN. Each matrix contains epochs of LLI
                                   slip starts in first column and LLI slip ends in second
                                   column. These may be often be empty depending on if
                                   RINEX observation files include LLI indicators or
                                   not.


     range1_observations:  matrix. Contains all range1 observations for all
                           epochs and for all SVs.

                           range1_observations(epoch, PRN)

     tInterval:            observation interval in seconds




     OUTPUTS:
     -------


     mean_multipath_range1:        array. contains mean values of estimates of
                                   multipath delay of the first range signal,
                                   for each satellite

                                   mean_multipath_range1(PRN)


     overall_mean_multipath_range1: overall mean of all estimates of
                                    multipath delay of the first range signal

     rms_multipath_range1:         rms values of estimates of multipath delay
                                   of the first range signal, for each
                                   satellite.

                                   rms_multipath_range1(PRN)

     average_rms_multipath_range1: average rms of estimates of multipath delay
                                   of the first range signal

     elevation_weighted_rms_multipath_range1:  rms_multipath_range1, but
                                               weighted based on elevation
                                               angles

     elevation_weighted_average_rms_multipath_range1:  average_rms_multipath_range1, but
                                                       weighted based on elevation
                                                       angles

     mean_ion_delay_phase1:        array. contains mean values of estimates of
                                   ionospheric delay on the first phase signal,
                                   for each satellite

                                   mean_ion_delay_phase1(PRN)

     overall_mean_ion_delay_phase1: overall mean of all estimates of
                                    ionospheric delay on the first phase signal

     mean_sat_elevation_angles:    array. contains mean elevation angle, for
                                   each satellite.

                                   mean_sat_elevation_angles(PRN)

     nEstimates:                   total amount of epochs with multipath estimates, double.

     nEstimates_per_sat:           amount of epochs with multipath estimates
                                   per sat.

     nRange1Obs_Per_Sat:           array. each elements says how many range1
                                   observations for that SV

                                   nRange1Obs_Per_Sat(PRN)

     nRange1Obs:                   total number of range1 observations, all SVs

     range1_slip_distribution_per_sat:     dict. each element contains a structure for
                                   each satellite. each structure has
                                   information on number of range1 slips
                                   distributed over groups depending on the
                                   elevation angle of the satellite at the
                                   time of the slip.

     range1_slip_distribution:     dict. Stores information on number of
                                   range1 slips distributed over groups depending
                                   on the elevation angle of the satellite at the
                                   time of the slip.

     LLI_slip_distribution_per_sat:     dict. each element contains a structure for
                                   each satellite. each structure has
                                   information on number of LLI detected slips
                                   distributed over groups depending on the
                                   elevation angle of the satellite at the
                                   time of the slip.

     LLI_slip_distribution:        dict. Stores information on number of
                                   LLI detected slips distributed over groups depending
                                   on the elevation angle of the satellite at the
                                   time of the slip.

     combined_slip_distribution_per_sat:     dict. each element contains a structure for
                                   each satellite. each structure has
                                   information on number of slips detected
                                   both by LLI and this siftware analysis,
                                   distributed over groups depending on the
                                   elevation angle of the satellite at the
                                   time of the slip.

     combined_slip_distribution:   dict. Stores information on number of
                                   slips detected both ny LLI and this software
                                   analysis, distributed over groups depending
                                   on the elevation angle of the satellite at the
                                   time of the slip.
    --------------------------------------------------------------------------------------------------------------------------
    """
    nSat = len(range1_slip_periods)

    # set all 0 values to NaN so they are excluded from stats calculation
    ion_delay_phase1[ion_delay_phase1==0] = np.nan
    multipath_range1[multipath_range1==0] = np.nan
    current_sat_elevation_angles[current_sat_elevation_angles==0] = np.nan

    mean_multipath_range1 = np.nanmean(multipath_range1,axis=0)              # Compute mean
    overall_mean_multipath_range1 = np.nanmean(mean_multipath_range1,axis=0) # overall mean multipath, excluding NaN values
    rms_multipath_range1 = np.nanstd(multipath_range1,axis=0, ddof=0)        # RMS multipath of each satellite, excluding NaN
    average_rms_multipath_range1 = np.nanstd(multipath_range1, ddof=0)       # Average RMS multipath, excluding NaN

    #  Weighted RMS multipath
    weights     = current_sat_elevation_angles.copy()
    crit_weight = 4*np.sin(30*np.pi/180)**2
    weights     = 4*np.sin(weights*np.pi/180)**2
    weights[weights > crit_weight] = 1
    elevation_weighted_multipath_range1 = multipath_range1*weights

    # RMS multipath of each satellite, excluding NaN
    elevation_weighted_rms_multipath_range1 = np.sqrt(np.nanmean(elevation_weighted_multipath_range1*elevation_weighted_multipath_range1,axis=0))
    # Average RMS multipath, excluding NaN
    elevation_weighted_average_rms_multipath_range1 = np.sqrt(np.nanmean(elevation_weighted_multipath_range1*elevation_weighted_multipath_range1))

    # Ionosphere
    mean_ion_delay_phase1 = np.nanmean(ion_delay_phase1,axis=0) # mean ionospheric delay for each satellite, excluding NaN
    overall_mean_ion_delay_phase1 = np.nanmean(mean_ion_delay_phase1) # Overall mean ionospheric delay, excluding NaN

    ## Average elevation angle for each satellite, excluding NaN
    dumm1 = (range1_observations!= 0)*1 # multiplying with 1 to get from True/False -> 1/0 # commented out 25.11.2023
    dumm2 = (~np.isnan(range1_observations))*1
    dummy = (dumm1 & dumm2)
    obs_elevations = current_sat_elevation_angles*dummy
    obs_elevations[obs_elevations==0]=np.nan

    # Compute mean_sat_elevation_angles safely for current_sat_elevation_angles
    mean_sat_elevation_angles = np.nanmean(current_sat_elevation_angles, axis=0)

    ## -- Amount of epochs with estimates
    nEstimates = np.sum(~np.isnan(multipath_range1))
    nEstimates_per_sat = np.sum(~np.isnan(multipath_range1),axis=0)

    ## -- Cycle slip distribution based on code-phase difference ---
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

    ## -- Cycle slip distribution based on Loss of lock indicators---
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


    ## -- Combining cycle slip based on code-phase difference and LLI ---
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


    ## -- Detected Cycle slip distribution based on both code-phase difference and ionosphere residuals ---
    ambiguity_slip_distribution_per_sat = {}
    ambiguity_slip_distribution         = {}

    ambiguity_slip_distribution['n_slips_0_10']   = 0
    ambiguity_slip_distribution['n_slips_10_20']  = 0
    ambiguity_slip_distribution['n_slips_20_30']  = 0
    ambiguity_slip_distribution['n_slips_30_40']  = 0
    ambiguity_slip_distribution['n_slips_40_50']  = 0
    ambiguity_slip_distribution['n_slips_over50'] = 0
    ambiguity_slip_distribution['n_slips_NaN']    = 0
    ambiguity_slip_distribution['n_slips_Tot']    = 0

    for i in np.arange(0,nSat):
        ## Create struct for current sat
        range1_slip_distribution_per_sat[i] = {}
        slip_epochs = []
        nSlipPeriods = len(range1_slip_periods[i+1])

       ## -- Get all slip epochs of current sat.
        for j in np.arange(0,nSlipPeriods):
            ## if slip period is shorter than 60 seconds
            slip_epochs.append(int(range1_slip_periods[i+1][j,0]))

        ## Get elevation angles for every slip of current sat
        slip_epoch_elevation_angles = current_sat_elevation_angles[slip_epochs, i+1]

        ## -- Store number of slips into groups of their elevation angles. Groups are: 0-10, 10-20, 20-30, 30-40, 40-50, >50, and NaN
        range1_slip_distribution_per_sat[i]['n_slips_0_10']    = len(slip_epoch_elevation_angles[(slip_epoch_elevation_angles >=0)  & (slip_epoch_elevation_angles <10)])
        range1_slip_distribution_per_sat[i]['n_slips_10_20']   = len(slip_epoch_elevation_angles[(slip_epoch_elevation_angles >=10) & (slip_epoch_elevation_angles <20)])
        range1_slip_distribution_per_sat[i]['n_slips_20_30']   = len(slip_epoch_elevation_angles[(slip_epoch_elevation_angles >=20) & (slip_epoch_elevation_angles <30)])
        range1_slip_distribution_per_sat[i]['n_slips_30_40']   = len(slip_epoch_elevation_angles[(slip_epoch_elevation_angles >=30) & (slip_epoch_elevation_angles <40)])
        range1_slip_distribution_per_sat[i]['n_slips_40_50']   = len(slip_epoch_elevation_angles[(slip_epoch_elevation_angles >=40) & (slip_epoch_elevation_angles <50)])
        range1_slip_distribution_per_sat[i]['n_slips_over50']  = len(slip_epoch_elevation_angles[(slip_epoch_elevation_angles >=50)])
        range1_slip_distribution_per_sat[i]['n_slips_NaN']     = len(slip_epoch_elevation_angles[np.isnan(slip_epoch_elevation_angles)])
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
        LLI_slip_distribution_per_sat[i] = {}
        LLI_slip_epochs = []
        nSlipPeriods = len(LLI_slip_periods[i])

        ## -- Get all LLI slip epochs of current sat.
        for j in np.arange(0,nSlipPeriods):
            LLI_slip_epochs.append(int(LLI_slip_periods[i][j,0]))

        ## Get elevation angles for every LLI slip of current sat
        if not np.isnan(current_sat_elevation_angles[LLI_slip_epochs, i+1]).any():
            LLI_slip_epoch_elevation_angles = current_sat_elevation_angles[LLI_slip_epochs, i+1].astype(int)
        else:
            LLI_slip_epoch_elevation_angles = current_sat_elevation_angles[LLI_slip_epochs, i+1]

        ## -- Store number of LLI slips into groups of their elevation angles. Groups are: 0-10, 10-20, 20-30, 30-40, 40-50, >50, and NaN
        LLI_slip_distribution_per_sat[i]['n_slips_0_10']    = len(LLI_slip_epoch_elevation_angles[(LLI_slip_epoch_elevation_angles >=0)  & (LLI_slip_epoch_elevation_angles <10)])
        LLI_slip_distribution_per_sat[i]['n_slips_10_20']   = len(LLI_slip_epoch_elevation_angles[(LLI_slip_epoch_elevation_angles >=10) & (LLI_slip_epoch_elevation_angles <20)])
        LLI_slip_distribution_per_sat[i]['n_slips_20_30']   = len(LLI_slip_epoch_elevation_angles[(LLI_slip_epoch_elevation_angles >=20) & (LLI_slip_epoch_elevation_angles <30)])
        LLI_slip_distribution_per_sat[i]['n_slips_30_40']   = len(LLI_slip_epoch_elevation_angles[(LLI_slip_epoch_elevation_angles >=30) & (LLI_slip_epoch_elevation_angles <40)])
        LLI_slip_distribution_per_sat[i]['n_slips_40_50']   = len(LLI_slip_epoch_elevation_angles[(LLI_slip_epoch_elevation_angles >=40) & (LLI_slip_epoch_elevation_angles <50)])
        LLI_slip_distribution_per_sat[i]['n_slips_over50']  = len(LLI_slip_epoch_elevation_angles[(LLI_slip_epoch_elevation_angles >=50)])
        LLI_slip_distribution_per_sat[i]['n_slips_NaN']     = len(LLI_slip_epoch_elevation_angles[np.isnan(LLI_slip_epoch_elevation_angles)])
        LLI_slip_distribution_per_sat[i]['n_slips_Tot']     = len(LLI_slip_epoch_elevation_angles)

        LLI_slip_distribution['n_slips_0_10']   = LLI_slip_distribution['n_slips_0_10']     + LLI_slip_distribution_per_sat[i]['n_slips_0_10']
        LLI_slip_distribution['n_slips_10_20']  = LLI_slip_distribution['n_slips_10_20']    + LLI_slip_distribution_per_sat[i]['n_slips_10_20']
        LLI_slip_distribution['n_slips_20_30']  = LLI_slip_distribution['n_slips_20_30']    + LLI_slip_distribution_per_sat[i]['n_slips_20_30']
        LLI_slip_distribution['n_slips_30_40']  = LLI_slip_distribution['n_slips_30_40']    + LLI_slip_distribution_per_sat[i]['n_slips_30_40']
        LLI_slip_distribution['n_slips_40_50']  = LLI_slip_distribution['n_slips_40_50']    + LLI_slip_distribution_per_sat[i]['n_slips_40_50']
        LLI_slip_distribution['n_slips_over50'] = LLI_slip_distribution['n_slips_over50']   + LLI_slip_distribution_per_sat[i]['n_slips_over50']
        LLI_slip_distribution['n_slips_NaN']    = LLI_slip_distribution['n_slips_NaN']      + LLI_slip_distribution_per_sat[i]['n_slips_NaN']
        LLI_slip_distribution['n_slips_Tot']    = LLI_slip_distribution['n_slips_Tot']      + LLI_slip_distribution_per_sat[i]['n_slips_Tot']


        ## -- For ambiguity
        ambiguity_slip_distribution_per_sat[i] = {}
        ambiguity_slip_epochs = []
        ambiguity_nSlipPeriods = len(ambiguity_slip_periods[i+1])

       ## -- Get all slip epochs of current sat.
        for j in np.arange(0,ambiguity_nSlipPeriods):
            ## if slip period is shorter than 60 seconds
            ambiguity_slip_epochs.append(int(ambiguity_slip_periods[i+1][j,0]))

        ## Get elevation angles for every slip of current sat
        slip_epoch_elevation_angles = current_sat_elevation_angles[ambiguity_slip_epochs, i+1]

        ## -- Store number of slips into groups of their elevation angles. Groups are: 0-10, 10-20, 20-30, 30-40, 40-50, >50, and NaN
        ambiguity_slip_distribution_per_sat[i]['n_slips_0_10']    = len(slip_epoch_elevation_angles[(slip_epoch_elevation_angles >=0)  & (slip_epoch_elevation_angles <10)])
        ambiguity_slip_distribution_per_sat[i]['n_slips_10_20']   = len(slip_epoch_elevation_angles[(slip_epoch_elevation_angles >=10) & (slip_epoch_elevation_angles <20)])
        ambiguity_slip_distribution_per_sat[i]['n_slips_20_30']   = len(slip_epoch_elevation_angles[(slip_epoch_elevation_angles >=20) & (slip_epoch_elevation_angles <30)])
        ambiguity_slip_distribution_per_sat[i]['n_slips_30_40']   = len(slip_epoch_elevation_angles[(slip_epoch_elevation_angles >=30) & (slip_epoch_elevation_angles <40)])
        ambiguity_slip_distribution_per_sat[i]['n_slips_40_50']   = len(slip_epoch_elevation_angles[(slip_epoch_elevation_angles >=40) & (slip_epoch_elevation_angles <50)])
        ambiguity_slip_distribution_per_sat[i]['n_slips_over50']  = len(slip_epoch_elevation_angles[(slip_epoch_elevation_angles >=50)])
        ambiguity_slip_distribution_per_sat[i]['n_slips_NaN']     = len(slip_epoch_elevation_angles[np.isnan(slip_epoch_elevation_angles)])
        ambiguity_slip_distribution_per_sat[i]['n_slips_Tot']     = len(slip_epoch_elevation_angles)

        ambiguity_slip_distribution['n_slips_0_10']   = ambiguity_slip_distribution['n_slips_0_10']     + ambiguity_slip_distribution_per_sat[i]['n_slips_0_10']
        ambiguity_slip_distribution['n_slips_10_20']  = ambiguity_slip_distribution['n_slips_10_20']    + ambiguity_slip_distribution_per_sat[i]['n_slips_10_20']
        ambiguity_slip_distribution['n_slips_20_30']  = ambiguity_slip_distribution['n_slips_20_30']    + ambiguity_slip_distribution_per_sat[i]['n_slips_20_30']
        ambiguity_slip_distribution['n_slips_30_40']  = ambiguity_slip_distribution['n_slips_30_40']    + ambiguity_slip_distribution_per_sat[i]['n_slips_30_40']
        ambiguity_slip_distribution['n_slips_40_50']  = ambiguity_slip_distribution['n_slips_40_50']    + ambiguity_slip_distribution_per_sat[i]['n_slips_40_50']
        ambiguity_slip_distribution['n_slips_over50'] = ambiguity_slip_distribution['n_slips_over50']   + ambiguity_slip_distribution_per_sat[i]['n_slips_over50']
        ambiguity_slip_distribution['n_slips_NaN']    = ambiguity_slip_distribution['n_slips_NaN']      + ambiguity_slip_distribution_per_sat[i]['n_slips_NaN']
        ambiguity_slip_distribution['n_slips_Tot']    = ambiguity_slip_distribution['n_slips_Tot']      + ambiguity_slip_distribution_per_sat[i]['n_slips_Tot']




        ## create dict for current sat
        combined_slip_distribution_per_sat[i] = {}

        ## -- Get all combined slips for current satellite
        combined_slip_epochs = np.union1d(slip_epochs, LLI_slip_epochs) ## må være union her??

        # get elevation angles for every combined slip of current sat
        if len(combined_slip_epochs) != 0:
            combined_slip_epochs = combined_slip_epochs.astype(int)
            combined_slip_epoch_elevation_angles = current_sat_elevation_angles[combined_slip_epochs, i+1].astype(int)

            # Store number of LLI slips into groups of their elevation angles. Groups are: 0-10, 10-20, 20-30, 30-40, 40-50, >50, and NaN
            combined_slip_distribution_per_sat[i]['n_slips_0_10']    = len(combined_slip_epoch_elevation_angles[(combined_slip_epoch_elevation_angles >=0)  & (combined_slip_epoch_elevation_angles <10)])
            combined_slip_distribution_per_sat[i]['n_slips_10_20']   = len(combined_slip_epoch_elevation_angles[(combined_slip_epoch_elevation_angles >=10) & (combined_slip_epoch_elevation_angles <20)])
            combined_slip_distribution_per_sat[i]['n_slips_20_30']   = len(combined_slip_epoch_elevation_angles[(combined_slip_epoch_elevation_angles >=20) & (combined_slip_epoch_elevation_angles <30)])
            combined_slip_distribution_per_sat[i]['n_slips_30_40']   = len(combined_slip_epoch_elevation_angles[(combined_slip_epoch_elevation_angles >=30) & (combined_slip_epoch_elevation_angles <40)])
            combined_slip_distribution_per_sat[i]['n_slips_40_50']   = len(combined_slip_epoch_elevation_angles[(combined_slip_epoch_elevation_angles >=40) & (combined_slip_epoch_elevation_angles <50)])
            combined_slip_distribution_per_sat[i]['n_slips_over50']  = len(combined_slip_epoch_elevation_angles[(combined_slip_epoch_elevation_angles >=50)])
            combined_slip_distribution_per_sat[i]['n_slips_NaN']     = len(combined_slip_epoch_elevation_angles[np.isnan(combined_slip_epoch_elevation_angles)])
            combined_slip_distribution_per_sat[i]['n_slips_Tot']     = len(combined_slip_epoch_elevation_angles)

            combined_slip_distribution['n_slips_0_10']   = combined_slip_distribution['n_slips_0_10']     + combined_slip_distribution_per_sat[i]['n_slips_0_10']
            combined_slip_distribution['n_slips_10_20']  = combined_slip_distribution['n_slips_10_20']    + combined_slip_distribution_per_sat[i]['n_slips_10_20']
            combined_slip_distribution['n_slips_20_30']  = combined_slip_distribution['n_slips_20_30']    + combined_slip_distribution_per_sat[i]['n_slips_20_30']
            combined_slip_distribution['n_slips_30_40']  = combined_slip_distribution['n_slips_30_40']    + combined_slip_distribution_per_sat[i]['n_slips_30_40']
            combined_slip_distribution['n_slips_40_50']  = combined_slip_distribution['n_slips_40_50']    + combined_slip_distribution_per_sat[i]['n_slips_40_50']
            combined_slip_distribution['n_slips_over50'] = combined_slip_distribution['n_slips_over50']   + combined_slip_distribution_per_sat[i]['n_slips_over50']
            combined_slip_distribution['n_slips_NaN']    = combined_slip_distribution['n_slips_NaN']      + combined_slip_distribution_per_sat[i]['n_slips_NaN']
            combined_slip_distribution['n_slips_Tot']    = combined_slip_distribution['n_slips_Tot']      + combined_slip_distribution_per_sat[i]['n_slips_Tot']

        else:

            ## Set to zero
            combined_slip_distribution_per_sat[i]['n_slips_0_10']    = 0
            combined_slip_distribution_per_sat[i]['n_slips_10_20']   = 0
            combined_slip_distribution_per_sat[i]['n_slips_20_30']   = 0
            combined_slip_distribution_per_sat[i]['n_slips_30_40']   = 0
            combined_slip_distribution_per_sat[i]['n_slips_40_50']   = 0
            combined_slip_distribution_per_sat[i]['n_slips_over50']  = 0
            combined_slip_distribution_per_sat[i]['n_slips_NaN']     = 0
            combined_slip_distribution_per_sat[i]['n_slips_Tot']     = 0

            combined_slip_distribution['n_slips_0_10']   = 0
            combined_slip_distribution['n_slips_10_20']  = 0
            combined_slip_distribution['n_slips_20_30']  = 0
            combined_slip_distribution['n_slips_30_40']  = 0
            combined_slip_distribution['n_slips_40_50']  = 0
            combined_slip_distribution['n_slips_over50'] = 0
            combined_slip_distribution['n_slips_NaN']    = 0
            combined_slip_distribution['n_slips_Tot']    = 0


    ## --Amount of range1 observations
    nRange1Obs_Per_Sat = np.sum(~np.isnan(range1_observations),axis=0).reshape(1, -1) # reshape to get 2D array
    nRange1Obs =  np.sum(~np.isnan(range1_observations))


    return mean_multipath_range1, overall_mean_multipath_range1, rms_multipath_range1, average_rms_multipath_range1,\
        mean_ion_delay_phase1, overall_mean_ion_delay_phase1, mean_sat_elevation_angles, nEstimates, nEstimates_per_sat,\
        nRange1Obs_Per_Sat, nRange1Obs, range1_slip_distribution_per_sat, range1_slip_distribution, ambiguity_slip_distribution_per_sat, ambiguity_slip_distribution,LLI_slip_distribution_per_sat, LLI_slip_distribution,\
        combined_slip_distribution_per_sat, combined_slip_distribution, elevation_weighted_rms_multipath_range1, elevation_weighted_average_rms_multipath_range1