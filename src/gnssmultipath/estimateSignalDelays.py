"""
This module is computing the multipath effect, ionospheric delay and other
relevant parameters.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""

import warnings
import numpy as np
from gnssmultipath.detectCycleSlips import detectCycleSlips, orgSlipEpochs
warnings.filterwarnings("ignore")

def estimateSignalDelays(range1_Code, range2_Code,phase1_Code, phase2_Code, carrier_freq1, \
                         carrier_freq2, nepochs, max_sat, GNSS_SVs, obsCodes, GNSS_obs, \
                             currentGNSSsystem, tInterval, phaseCodeLimit, ionLimit):
    """
     Function that takes observations of from the observation period and
     estimates the following delays on the signal:
                                                   ionospheric delay on the
                                                   first phase signal

                                                   multipath delay on the
                                                   first range signal

                                                   multipath delay on the
                                                   second range signal

     Also estimates epochs with ambiguity slip for the first signal.

     The estimates are all relative estimates. The multipath estimates are
     relative to a mean values. The ionosphere delay estimates are relative to
     the first estimate. These estimates must be reduced again by the
     releative value at ambiguity slips. The function calls on another
     function to estimate periods of ambiguity slips. It then corrects the
     delay estimates for these slips.
    --------------------------------------------------------------------------------------------------------------------------
     INPUTS

     range1_Code:          string, defines the observation type that will be
                           used as the first range observation. ex. "C1X", "C5X"

     range2_Code:          string, defines the observation type that will be
                           used as the second range observation. ex. "C1X", "C5X"

     phase1_Code:          string, defines the observation type that will be
                           used as the first phase observation. ex. "L1X", "L5X"

     phase2_Code:          string, defines the observation type that will be
                           used as the second phase observation. ex. "L1X", "L5X"

     carrier_freq1:        carrier frequency of first phase signal. unit Hz

     carrier_freq2:        carrier frequency of second phase signal. unit Hz

     nepochs:              number of epochs with observations in rinex observation file.

     max_sat:              max PRN number of current GNSS system.

     GNSS_SVs:             matrix containing number of satellites with
                           obsevations for each epoch, and PRN for those satellites

                           GNSS_SVs(epoch, j)  j=1: number of observed satellites
                                               j>1: PRN of observed satellites

     obsCodes:             cell contaning strings. Cell defines the observation
                           codes available from the current GNSS system. Each
                           string is a three-character string, the first
                           character (a capital letter) is an observation code
                           ex. "L" or "C". The second character (a digit)
                           is a frequency code. The third character(a Capital letter)
                           is the attribute, ex. "P" or "X"

     GNSS_obs:             3D matrix containing all observation of current
                           GNSS system for all epochs. Order of obsType index
                           is same order as in obsCodes cell

                           GNSS_obs(PRN, obsType, epoch)
                                               PRN: double
                                               ObsType: double: 1,2,...,numObsTypes
                                               epoch: double

     currentGNSSsystem:    string, code to indicate which GNSS system is the
                           current observations are coming from. Only used to
                           decide if system uses CDMA or FDMA.

     tInterval:            observations interval seconds.

     phaseCodeLimit:       critical limit that indicates cycle slip for
                           phase-code combination. Unit: m/s. If set to 0,
                           default value will be used

     ionLimit:             critical limit that indicates cycle slip for
                           the rate of change of the ionopheric delay.
                           Unit: m/s. If set to 0, default value will be used
    --------------------------------------------------------------------------------------------------------------------------
     OUTPUTS

     ion_delay_phase1:     2D array containing estimates of ionospheric delays
                           on the first phase signal for each PRN, at each epoch.

                           ion_delay_phase1[epoch, PRN]

     multipath_range1:     2D array containing estimates of multipath delays
                           on the first range signal for each PRN, at each epoch.

                           multipath_range2(epoch, PRN)

     multipath_range2:     2D array containing estimates of multipath delays
                           on the second range signal for each PRN, at each epoch.

                           multipath_range2[epoch, PRN]

     range1_slip_periods:  array, each cell element contains one matrix for
                           every PRN. Each matrix contains epochs of range1
                           slip starts in first column and slip ends in second
                           column.

     range1_observations:  2D array. Contains all range1 observations for all
                           epochs and for all SVs.

                           range1_observations[epoch, PRN]

     phase1_observations:  2D array. Contains all phase1 observations for all
                           epochs and for all SVs.

                           phase1_observations[epoch, PRN]

     success:              boolean, 1 if no error occurs, 0 otherwise

    --------------------------------------------------------------------------------------------------------------------------
    """

    FDMA_used = 0
    success = 1

    if 'R' in currentGNSSsystem:
        FDMA_used = 1
        carrier_freq1_list = carrier_freq1
        carrier_freq2_list = carrier_freq2
    else:
        alpha = carrier_freq1**2/carrier_freq2**2 # amplfication factor

    # Define parameters
    c = 299792458 # speed of light

    if ionLimit ==0:
        ionLimit = 4/60   # critical rate of change of ionosphere delay  to indicate ambiguity slip on either
                          # range1/phase1 signal, range2/phase2 signal, or both


    if phaseCodeLimit == 0:
        phaseCodeLimit = 4/60*100   # critical rate of change of
                                    # N1_pseudo_estimate to indicate ambiguity slip on
                                    # the range1/phase1 signal

    ## -- Initialize data matrices
    ion_delay_phase1        = np.zeros([nepochs, max_sat+1])
    multipath_range1        = np.zeros([nepochs, max_sat+1])
    # multipath_range2        = np.zeros([nepochs, max_sat+1]) # kommenterte bort 23.01.2023
    N1_pseudo_estimate      = np.zeros([nepochs, max_sat+1])
    missing_obs_overview    = np.zeros([nepochs, max_sat+1])
    missing_range1_overview = np.zeros([nepochs, max_sat+1])

    ambiguity_slip_periods = {int(key): [] for key in range(1,  max_sat+1)} # Initialize cell for storing phase slip periods
    range1_slip_periods = {int(key): [] for key in range(1,  max_sat+1)} # Initialize cell for storing slip periods for only range1/phase1

    if FDMA_used: # denne er inødnendig etter vektoriseringen
        carrier_freq1 = carrier_freq1_list
        carrier_freq2 = carrier_freq2_list
        alpha = carrier_freq1**2/carrier_freq2**2 # amplfication factor


    #  Get observations
    range1 = create_array_for_current_obscode(GNSS_obs, ismember(obsCodes[currentGNSSsystem],range1_Code))
    range2 = create_array_for_current_obscode(GNSS_obs, ismember(obsCodes[currentGNSSsystem],range2_Code))

    phase1 = create_array_for_current_obscode(GNSS_obs, ismember(obsCodes[currentGNSSsystem],phase1_Code))*c/carrier_freq1
    phase2 = create_array_for_current_obscode(GNSS_obs, ismember(obsCodes[currentGNSSsystem],phase2_Code))*c/carrier_freq2

    if all(v is None for v in [range1, range2, phase1, phase2]):
        print('ERROR(estimateSignalDelays): Some observations are missing. Check for missing data in RINEX observation file!')
        success = 0
        return success


    ## -- Calculate estimate for Ionospheric delay on phase 1 signal
    ion_delay_phase1 = 1/(alpha-1)*(phase1-phase2)
    ## -- Calculate multipath on both code signals
    multipath_range1 = range1 - (1 + 2/(alpha-1))*phase1 + (2/(alpha-1))*phase2
    # multipath_range2[epoch,PRN] = range2 - (2*alpha/(alpha-1))*phase1 + ((2*alpha)/(alpha-1) - 1)*phase2 ## kommenterte bort 23.01.2023
    N1_pseudo_estimate = phase1 - range1

    # If any of the four observations are missing, ie np.nan, estimate remains 0 for that epoch and satellite
    missing_obs_overview = find_missing_observation(range1, range2, array3=phase1, array4=phase2)
    missing_range1_overview = find_missing_observation(range1,phase1)

    ## -- Get first and last epoch with observations for current PRN
    epoch_first_obs = find_epoch_of_first_obs(ion_delay_phase1)
    epoch_last_obs = find_epoch_of_last_obs(ion_delay_phase1)

    ## -- Get first and last epoch with range 1 observations for current PRN
    epoch_first_range1obs = find_epoch_of_first_obs(N1_pseudo_estimate)
    epoch_last_range1obs =  find_epoch_of_last_obs(N1_pseudo_estimate)

    range1_slip_epochs = detectCycleSlips(N1_pseudo_estimate, missing_range1_overview, epoch_first_range1obs, epoch_last_range1obs, tInterval, phaseCodeLimit) # detect cycle slips for current epoch for range1/phase1 only
    ionosphere_slip_epochs = detectCycleSlips(ion_delay_phase1, missing_obs_overview, epoch_first_obs, epoch_last_obs, tInterval, ionLimit) # detect cycle slips for current epoch for either range1/phase1 signal, range2/phase2 signal, or both


    ## -- Detect and correct for ambiguity slips
    for PRN in np.arange(0,max_sat):
        PRN = PRN + 1
        if np.isnan(epoch_first_obs[PRN]):
            continue  # Skip the current iteration if it's np.nan
        epoch_first_obs_prn = int(epoch_first_obs[PRN])
        epoch_last_obs_prn = int(epoch_last_obs[PRN])

        ambiguity_slip_epochs = np.union1d(np.array(range1_slip_epochs[str(PRN)]),np.array(ionosphere_slip_epochs[str(PRN)])) # Make combined array of slip epochs from both lin. combinations used to detects slips
        range1_slip_periods[int(PRN)],_ = orgSlipEpochs(np.array(range1_slip_epochs[str(PRN)])) # Organize slips detected on range1/phase1 signal only
        ambiguity_slip_periods[int(PRN)], n_slip_periods = orgSlipEpochs(ambiguity_slip_epochs) # Orginize combined slips detected on range1/phase1 signal only

        ## If there are no slips then there is only one "ambiguity period". All estimates are therefore reduced by the same relative value
        if len(ambiguity_slip_periods[int(PRN)]) == 0:
            ion_delay_phase1[epoch_first_obs_prn::, PRN] = ion_delay_phase1[epoch_first_obs_prn::, PRN] - ion_delay_phase1[epoch_first_obs_prn, PRN]
            multipath_range1[epoch_first_obs_prn::, PRN] = multipath_range1[epoch_first_obs_prn::, PRN] - np.nanmean(multipath_range1[epoch_first_obs_prn::, PRN])
            # multipath_range2[epoch_first_obs::, PRN] = multipath_range2[epoch_first_obs::, PRN] - np.mean(np.nonzero(multipath_range2[epoch_first_obs::, PRN]))
        else:
            ## -- Set all estimates of epochs with cycle slips to nan
            for slip_period in np.arange(0,n_slip_periods):
                slip_start   = int(ambiguity_slip_periods[PRN][slip_period,0])
                slip_end     = int(ambiguity_slip_periods[PRN][slip_period,1])
                if slip_start == slip_end:  # need a if test bacause if there equal, python dont set to nan
                    ion_delay_phase1[slip_start, PRN] = np.nan
                    multipath_range1[slip_start, PRN] = np.nan
                    # multipath_range2[slip_start, PRN] = np.nan
                else:
                    ion_delay_phase1[slip_start:slip_end+1, PRN] = np.nan
                    multipath_range1[slip_start:slip_end+1, PRN] = np.nan
                    # multipath_range2[slip_start:slip_end+1, PRN] = np.nan

            ## Extract start and end of each segment and correct multipath and ionosphere estimates for each segment
            for ambiguity_period in np.arange(0,n_slip_periods+1):
                if ambiguity_period == 0:
                    ambiguity_period_start  = epoch_first_obs_prn
                    ambiguity_period_end    = int(ambiguity_slip_periods[PRN][0,0])
                ## -- If last ambiguity period
                elif ambiguity_period == n_slip_periods: # removed + 1
                    ambiguity_period_start       = int(ambiguity_slip_periods[PRN][-1,1] + 1)
                    ambiguity_period_end         = epoch_last_obs_prn +1

                    ## If last epoch with observation is a slip, then there is no  last ambiguity period
                    if ambiguity_period_start > epoch_last_obs_prn:
                        ambiguity_period_start = []
                        ambiguity_period_end = []

                else:
                    ambiguity_period_start = int(ambiguity_slip_periods[PRN][ambiguity_period-1, 1] + 1)
                    ambiguity_period_end   = int(ambiguity_slip_periods[PRN][ambiguity_period, 0])


                #  Ionosphere delay estimates of current ambiguity period is reduced by first estimate of ambiguity period
                if ambiguity_period_start != ambiguity_period_end:
                    ion_delay_phase1[ambiguity_period_start:ambiguity_period_end, PRN] = ion_delay_phase1[ambiguity_period_start:ambiguity_period_end, PRN] - \
                        ion_delay_phase1[ambiguity_period_start, PRN]

                    # Multipath delays of current ambiguity period are reduced by mean of estimates in ambiguity period, excluding NaN and
                    slice_data = multipath_range1[ambiguity_period_start:ambiguity_period_end, PRN]
                    if slice_data.size > 0 and not np.all(np.isnan(slice_data)):
                        multipath_range1[ambiguity_period_start:ambiguity_period_end, PRN] = slice_data - np.nanmean(slice_data)
                    else:
                        multipath_range1[ambiguity_period_start:ambiguity_period_end, PRN] = np.nan

                    # multipath_range2[ambiguity_period_start:ambiguity_period_end, PRN] = multipath_range2[ambiguity_period_start:ambiguity_period_end, PRN] -\
                        # np.nanmean(multipath_range2[ambiguity_period_start:ambiguity_period_end, PRN])
                else:

                    ion_delay_phase1[ambiguity_period_start, PRN] = ion_delay_phase1[ambiguity_period_start, PRN] - \
                        ion_delay_phase1[ambiguity_period_start, PRN]
                    #  Multipath delays of current ambiguity period are reduced by mean of estimates in ambiguity period, excluding NaN and
                    multipath_range1[ambiguity_period_start, PRN] = multipath_range1[ambiguity_period_start, PRN] - \
                        np.nanmean(multipath_range1[ambiguity_period_start, PRN])

                    # multipath_range2[ambiguity_period_start, PRN] = multipath_range2[ambiguity_period_start, PRN] -\
                        # np.nanmean(multipath_range2[ambiguity_period_start, PRN])

    # Create array of all code and phase observation for the current obscodes
    range1_observations = create_array_for_current_obscode(GNSS_obs, ismember(obsCodes[currentGNSSsystem],range1_Code))
    phase1_observations = create_array_for_current_obscode(GNSS_obs, ismember(obsCodes[currentGNSSsystem],phase1_Code))

    # return ion_delay_phase1, multipath_range1, multipath_range2, ambiguity_slip_periods, range1_observations, phase1_observations, success # changeing from range1slip to amgiguity
    return ion_delay_phase1, multipath_range1, range1_slip_periods, ambiguity_slip_periods, range1_observations, phase1_observations, success #removed multipath_range2



def ismember(list_,code):
    """
    The function takes in a string and a list, and finds the index of
    """
    indx = [idx for idx, val in enumerate(list_) if val == code]
    if indx != []:
        indx = indx[0]
    return indx


def create_array_for_current_obscode(GNSS_obs, obscode_idx):
    """
    Takes inn a dict of numpy arrays and returns a
    numpy array where the keys are rows in the array.
    Replaces all zeros with np.nan

    """
    try:
        code_array = np.stack(list(GNSS_obs.values()))[:, :, obscode_idx]
        code_array = np.squeeze(code_array)
        code_array[code_array == 0] = np.nan
    except:
        code_array = None
    return code_array



def find_epoch_of_first_obs(ion_delay_phase1):
    """
    Finds the first epoch with observation and
    stores the epoch/row index for each satellite as a single-column array.
    A column with no observation at all, is stored as  np.nan
    """
    not_nan_indices = np.isfinite(ion_delay_phase1)
    first_epochs_array = np.argmax(not_nan_indices, axis=0).astype(float)
    no_true_columns = np.all(~not_nan_indices, axis=0) # Find columns with no True values and set their indices to np.nan
    first_epochs_array[no_true_columns] = np.nan
    return first_epochs_array


def find_epoch_of_last_obs(ion_delay_phase1):
    """
    Finds the last epoch with observations and
    stores the epoch/row index for each satellite as a single-column array.
    A column with no observation at all, is stored as  np.nan
    """
    not_nan_indices = np.isfinite(ion_delay_phase1)
    reversed_array = np.flip(not_nan_indices, axis=0) # Reverse the boolean array along the rows
    last_epochs_array = not_nan_indices.shape[0] - 1 - np.argmax(reversed_array, axis=0).astype(float) # Find the indices of the first True element in each column
    last_epochs_array[~np.any(not_nan_indices, axis=0)] =  np.nan # Mask columns with no True values with np.nan
    return last_epochs_array



def find_missing_observation(array1,array2,array3=None,array4=None):
    """
    Creates a array with overview of epoch and sats
    that are missing observations.
    """
    if array3 is not None:
        mask1 = np.isnan(array1).astype(int)
        mask2 = np.isnan(array2).astype(int)
        mask3 = np.isnan(array3).astype(int)
        mask4 = np.isnan(array4).astype(int)
        # Combine the masks to create the final result
        missing_obs_overview = np.maximum.reduce([mask1, mask2, mask3, mask4])
    else:
        mask1 = np.isnan(array1).astype(int)
        mask2 = np.isnan(array2).astype(int)
        # Combine the masks to create the final result
        missing_obs_overview = np.maximum.reduce([mask1, mask2])
    return missing_obs_overview




