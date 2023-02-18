import numpy as np
from detectCycleSlips import detectCycleSlips, orgSlipEpochs
import warnings
warnings.filterwarnings(action='ignore', message='Mean of empty slice')

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
    
     tInterval:            observations interval; seconds. 
    
     phaseCodeLimit:       critical limit that indicates cycle slip for
                           phase-code combination. Unit: m/s. If set to 0,
                           default value will be used
    
     ionLimit:             critical limit that indicates cycle slip for
                           the rate of change of the ionopheric delay. 
                           Unit: m/s. If set to 0, default value will be used
    --------------------------------------------------------------------------------------------------------------------------
     OUTPUTS
    
     ion_delay_phase1:     matrix containing estimates of ionospheric delays 
                           on the first phase signal for each PRN, at each epoch.
    
                           ion_delay_phase1(epoch, PRN)
    
     multipath_range1:     matrix containing estimates of multipath delays 
                           on the first range signal for each PRN, at each epoch.
    
                           multipath_range2(epoch, PRN)
    
     multipath_range2:     matrix containing estimates of multipath delays 
                           on the second range signal for each PRN, at each epoch.
    
                           multipath_range2(epoch, PRN)
    
     range1_slip_periods:  cell, each cell element contains one matrix for
                           every PRN. Each matrix contains epochs of range1
                           slip starts in first column and slip ends in second
                           column.
    
     range1_observations:  matrix. Contains all range1 observations for all
                           epochs and for all SVs. 
    
                           range1_observations(epoch, PRN)
    
     phase1_observations:  matrix. Contains all phase1 observations for all
                           epochs and for all SVs. 
    
                           range1_observations(epoch, PRN)
    
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
    ion_delay_phase1        = np.zeros([nepochs, max_sat+1]) # Addin 1 to since python is nullindexed
    multipath_range1        = np.zeros([nepochs, max_sat+1])
    # multipath_range2        = np.zeros([nepochs, max_sat+1]) # kommenterte bort 23.01.2023
    N1_pseudo_estimate      = np.zeros([nepochs, max_sat+1])
    missing_obs_overview    = np.zeros([nepochs, max_sat+1])
    missing_range1_overview = np.zeros([nepochs, max_sat+1])
    
    ## -- Initialize cell for storing phase slip periods
    ambiguity_slip_periods = {}
    ## -- Initialize cell for storing slip periods for only range1/phase1

    range1_slip_periods = {}
    
    ## -- Make estimates of delays
    for epoch in np.arange(0,nepochs):
        ## -- Calculate how many satelites with observations in current epoch
        n_sat = int(GNSS_SVs[epoch,0])
        ## - Iterate over all satellites of current epoch
        for i in np.arange(0,n_sat):  
            PRN = int(GNSS_SVs[epoch, 1+i]) # Get PRN ## KANSJE FJERNE +1??? for å få med PRN 11
            
            if FDMA_used:
                carrier_freq1 = carrier_freq1_list[PRN]
                carrier_freq2 = carrier_freq2_list[PRN]
                alpha = carrier_freq1**2/carrier_freq2**2 # amplfication factor
            
            
            ## -- Get observations            
            range1 = GNSS_obs[epoch+1][PRN, ismember(obsCodes[currentGNSSsystem],range1_Code)] 
            range2 = GNSS_obs[epoch+1][PRN, ismember(obsCodes[currentGNSSsystem],range2_Code)]
            
            phase1 = GNSS_obs[epoch+1][PRN, ismember(obsCodes[currentGNSSsystem],phase1_Code)]*c/carrier_freq1;
            phase2 = GNSS_obs[epoch+1][PRN, ismember(obsCodes[currentGNSSsystem],phase2_Code)]*c/carrier_freq2;
            
            if any([str(range1), str(range2), str(phase1), str(phase2)]) == '[]':
                print('ERROR(estimateSignalDelays): There is no observation type #s. Check for missing data in RINEX observation file!' % (phase1_Code))
                success = 0
                return success 
            
            
            ## -- If any of the four observations are missing, ie 0, estimate remains 0 for that epoch and satellite
            if [range1, range2, phase1, phase2].count(0) == 0:
                ## -- Calculate estimate for Ionospheric delay on phase 1 signal
                ion_delay_phase1[epoch, PRN] = 1/(alpha-1)*(phase1-phase2)
                
                ## -- Calculate multipath on both code signals
                multipath_range1[epoch,PRN] = range1 - (1 + 2/(alpha-1))*phase1 + (2/(alpha-1))*phase2
                # multipath_range2[epoch,PRN] = range2 - (2*alpha/(alpha-1))*phase1 + ((2*alpha)/(alpha-1) - 1)*phase2 ## kommenterte bort 23.01.2023
                N1_pseudo_estimate[epoch, PRN] = phase1 - range1
            else:
                ## Flag epoch and PRN as missing obs
                missing_obs_overview[epoch, PRN] = 1
                if range1 == 0 or phase1 == 0:
                    missing_range1_overview[epoch, PRN] = 1 
                else:
                    N1_pseudo_estimate[epoch, PRN] = phase1 - range1

    
    ## -- Detect and correct for ambiguity slips
    for PRN in np.arange(0,max_sat):
        PRN = PRN + 1
        ## -- Get first and last epoch with observations for current PRN
        if len(np.nonzero(ion_delay_phase1[:,PRN])[0]) != 0:
            epoch_first_obs = np.nonzero(ion_delay_phase1[:,PRN])[0][0]
            epoch_last_obs = np.nonzero(ion_delay_phase1[:,PRN])[0][-1]
        else:
            epoch_first_obs = np.nonzero(ion_delay_phase1[:,PRN])[0]
            epoch_last_obs = np.nonzero(ion_delay_phase1[:,PRN])[0]
                
               
        ## -- Get first and last epoch with range 1 observations for current PRN
        if len(np.nonzero(N1_pseudo_estimate[:,PRN])[0]) != 0:
            epoch_first_range1obs =  np.nonzero(N1_pseudo_estimate[:,PRN])[0][0]
            epoch_last_range1obs =  np.nonzero(N1_pseudo_estimate[:,PRN])[0][-1]
        else:
            epoch_first_range1obs =  np.nonzero(N1_pseudo_estimate[:,PRN])[0]
            epoch_last_range1obs =  np.nonzero(N1_pseudo_estimate[:,PRN])[0]
        
        
        ## -- Run function to detect cycle slips for current epoch for range1/phase1 only 
        range1_slip_epochs = detectCycleSlips(N1_pseudo_estimate[:, PRN], missing_range1_overview[:, PRN],epoch_first_range1obs, epoch_last_range1obs, tInterval, phaseCodeLimit)
       
        # Run function to detect cycle slips for current epoch for either range1/phase1 signal, range2/phase2 signal, or both
        ionosphere_slip_epochs = detectCycleSlips(ion_delay_phase1[:, PRN], missing_obs_overview[:, PRN],epoch_first_obs, epoch_last_obs, tInterval, ionLimit)
           
        ## -- Make combined array of slip epochs from both lin. combinations used to detects slips
        ambiguity_slip_epochs = np.union1d(range1_slip_epochs, ionosphere_slip_epochs) 
        # range1_slip_epochs = np.intersect1d(range1_slip_epochs, ionosphere_slip_epochs) # hvorfor intersect her?
        range1_slip_epochs = range1_slip_epochs #tester om det gir stor forskjell uten intersect

        ## -- Organize slips detected on range1/phase1 signal only
        range1_slip_periods[PRN],_ = orgSlipEpochs(range1_slip_epochs)
        
        ## -- Orginize combined slips detected on range1/phase1 signal only
        ambiguity_slip_periods[PRN], n_slip_periods = orgSlipEpochs(ambiguity_slip_epochs)
        
        # Set all zero estimates to NaN so that epochs with missing observations are not corrected        
        ion_delay_phase1[ion_delay_phase1[:, PRN]==0, PRN] = np.nan
        multipath_range1[multipath_range1[:, PRN]==0, PRN] = np.nan
        # multipath_range2[multipath_range2[:, PRN]==0, PRN] = np.nan ## 23.01.2023 kommenterer bort
        
        ## If there are no slips then there is only one "ambiguity period". All estimates are therefore reduced by the same relative value
        if len(ambiguity_slip_periods[PRN]) == 0:
            # if len(list(epoch_first_obs)) == 0:
            if epoch_first_obs.size == 0: # added 18.02.2023 because of error when running on RINEX v2 (should be like this either way!)
                pass
            else:    
                ion_delay_phase1[epoch_first_obs::, PRN] = ion_delay_phase1[epoch_first_obs::, PRN] - ion_delay_phase1[epoch_first_obs, PRN]  
                # multipath_range1[epoch_first_obs::, PRN] = multipath_range1[epoch_first_obs::, PRN] - np.mean(np.nonzero(multipath_range1[epoch_first_obs::, PRN]))
                ### removed np.nonzero and use np.nanmean instead since every zero is replaced by nan in line 247 and 248
                multipath_range1[epoch_first_obs::, PRN] = multipath_range1[epoch_first_obs::, PRN] - np.nanmean(multipath_range1[epoch_first_obs::, PRN])  

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
                    ion_delay_phase1[slip_start:slip_end+1, PRN] = np.nan # + 1 because a[2:3] gives one element. Matlab a(2:3) gives 2 element.
                    multipath_range1[slip_start:slip_end+1, PRN] = np.nan # if error msg here, try add "if slip_end != epoch_last_obs else epoch_last_obs" in the slicing (oneliner)
                    # multipath_range2[slip_start:slip_end+1, PRN] = np.nan
            
              
            ## Extract start and end of each segment and correct multipath and ionosphere estimates for each segment
            for ambiguity_period in np.arange(0,n_slip_periods+1): # removed + 1 cause of indexproblem 29.11
                if ambiguity_period == 0: 
                    ambiguity_period_start  = epoch_first_obs
                    ambiguity_period_end    = int(ambiguity_slip_periods[PRN][0,0])   #INK -1 igjen???
                ## -- If last ambiguity period
                elif ambiguity_period == n_slip_periods: # removed + 1 
                    ambiguity_period_start       = int(ambiguity_slip_periods[PRN][-1,1] + 1) 
                    ambiguity_period_end         = epoch_last_obs +1
                   
                    ## If last epoch with observation is a slip, then there is no  last ambiguity period
                    if ambiguity_period_start > epoch_last_obs:
                        ambiguity_period_start = []
                        ambiguity_period_end = []
                    
                else:
                    ambiguity_period_start = int(ambiguity_slip_periods[PRN][ambiguity_period-1, 1] + 1)
                    ambiguity_period_end   = int(ambiguity_slip_periods[PRN][ambiguity_period, 0])
                
               
                ## -- Ionosphere delay estimates of current ambiguity period is reduced by first estimate of ambiguity period
                if ambiguity_period_start != ambiguity_period_end:
                    ion_delay_phase1[ambiguity_period_start:ambiguity_period_end, PRN] = ion_delay_phase1[ambiguity_period_start:ambiguity_period_end, PRN] - \
                        ion_delay_phase1[ambiguity_period_start, PRN]
                   
                    ## -- Multipath delays of current ambiguity period are reduced by mean of estimates in ambiguity period, excluding NaN and
                    multipath_range1[ambiguity_period_start:ambiguity_period_end, PRN] = multipath_range1[ambiguity_period_start:ambiguity_period_end, PRN] - \
                        np.nanmean(multipath_range1[ambiguity_period_start:ambiguity_period_end, PRN]) # added nanmean 30.11
                        
                    # multipath_range2[ambiguity_period_start:ambiguity_period_end, PRN] = multipath_range2[ambiguity_period_start:ambiguity_period_end, PRN] -\
                        # np.nanmean(multipath_range2[ambiguity_period_start:ambiguity_period_end, PRN])
                else:
                    
                    ion_delay_phase1[ambiguity_period_start, PRN] = ion_delay_phase1[ambiguity_period_start, PRN] - \
                        ion_delay_phase1[ambiguity_period_start, PRN]
                    ## -- Multipath delays of current ambiguity period are reduced by mean of estimates in ambiguity period, excluding NaN and
                    multipath_range1[ambiguity_period_start, PRN] = multipath_range1[ambiguity_period_start, PRN] - \
                        np.nanmean(multipath_range1[ambiguity_period_start, PRN])
                       
                    # multipath_range2[ambiguity_period_start, PRN] = multipath_range2[ambiguity_period_start, PRN] -\
                        # np.nanmean(multipath_range2[ambiguity_period_start, PRN])
                        
        ## -- Get range1 and phase 1 observations for all epochs and PRN
        range1_observations =  np.zeros([nepochs, max_sat+1]) 
        phase1_observations =  np.zeros([nepochs, max_sat+1]) 
        for ep in np.arange(0, len(GNSS_obs)):
            for PRN in np.arange(0,max_sat):
                range1_observations[ep,PRN] = GNSS_obs[ep+1][PRN,ismember(obsCodes[currentGNSSsystem],range1_Code)] 
                phase1_observations[ep,PRN] = GNSS_obs[ep+1][PRN,ismember(obsCodes[currentGNSSsystem],phase1_Code)]

    # return ion_delay_phase1, multipath_range1, multipath_range2, range1_slip_periods, range1_observations, phase1_observations, success
    # return ion_delay_phase1, multipath_range1, multipath_range2, ambiguity_slip_periods, range1_observations, phase1_observations, success # changeing from range1slip to amgiguity
    return ion_delay_phase1, multipath_range1, range1_slip_periods,ambiguity_slip_periods, range1_observations, phase1_observations, success #removed multipath_range2


def ismember(list_,code):
    """
    The function takes in a string and a list, and finds the index of 
    """
    indx = [idx for idx, val in enumerate(list_) if val == code]
    if indx != []:
        indx = indx[0]
    return indx


