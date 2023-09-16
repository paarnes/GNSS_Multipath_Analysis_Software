import numpy as np
def detectCycleSlips(estimates, missing_obs_overview,epoch_first_obs, epoch_last_obs, tInterval, crit_slip_rate):
    """
     Function that detects epochs with cycle slips, given test estimates
     and a critical rate of change
    --------------------------------------------------------------------------------------------------------------------------
    
     INPUTS:            
    --------------------------------------------------------------------------------------------------------------------------
     estimates:            matrix containing estimates from a linear combination
                           for a specific PRN estimates(epoch, PRN)
    
     missing_obs_overview: matrix of size nepochs x nPRN containing 1 or 0.
                           1 indicates that the satellite with this PRN has no
                           estimate at this epoch. There are no 1s before
                           first observation epoch or after last observation epoch            
                           as lack of estimates are implied here.
    
                           missing_obs_overview(epoch, PRN)
    
    
     epoch_first_obs:      epoch of first estimate for the current satellite
    
     epoch_last_obs:       epoch of last estimate for the current satellite
    
     tInterval:            observations interval in seconds. 
    
     crit_slip_rate:       critical rate of change of estimate to
                           indicate an cycle slip. [m/seconds]. 
    
    --------------------------------------------------------------------------------------------------------------------------
    
     OUTPUTS:
    
     slip_epochs:          array, contains epochs with detected cycle slip
    --------------------------------------------------------------------------------------------------------------------------
    """
    ## -- Preallocating
    slips_from_missing_obs = []
    ## -- Calculate rate of change of estimates of ionospheric delay (time derivative)
    estimates_rate_of_change = np.diff(estimates)/tInterval

    ## -- Detect slips for current sat. as epochs with either estimates_rate_of_change
    ## higher than critical value, or epochs with missing estimates
    slips_from_crit_rate =  [idx for idx,val in enumerate(estimates_rate_of_change) if abs(val) > crit_slip_rate] # use listcomp instead

    if type(epoch_first_obs) !=np.ndarray and type(epoch_last_obs) != np.ndarray:
        slips_from_missing_obs = (np.where(missing_obs_overview[epoch_first_obs:epoch_last_obs] == 1) + epoch_first_obs)[0].tolist() 

    if len(slips_from_missing_obs) != 0:
        slip_epochs = np.array(sorted(set(slips_from_crit_rate + slips_from_missing_obs)))
    else:
        slip_epochs = np.array(sorted(set(slips_from_crit_rate)))

    return slip_epochs 




def orgSlipEpochs(slip_epochs):
    """
     Function that takes array of epochs with detected cycle slips and for
     one satellite and organizing them into slip periods
    --------------------------------------------------------------------------------------------------------------------------
    
     INPUTS
    
     slip_epochs:          array, contains epochs with detected cycle slip
    --------------------------------------------------------------------------------------------------------------------------
    
     OUTPUTS
    
     slip_periods:         Matrix that contains the start of periods with cycle slips in
                           the first column and the ends of the same periods in
                           the second column.
    
                           cycle_slip_periods(ambiguity_slip_priod, j),
                           j = 1, ambiguity period start
                           j = 2, ambiguity period ends
    
     n_slip_periods:       amount of slip periods for current satellite
    
    --------------------------------------------------------------------------------------------------------------------------
    """

    ## -- If no slips occurs there are no slip periods for this sat.
    if len(slip_epochs) != 0:    
        # dummy is logical. It will be 1 at indices where the following
        # slip epoch is NOT the epoch following the current slip epoch.
        # These will therefor be the indices where slip periods end.
        # The last slip end is not detected this way and is inserted manually
        dummy = (np.diff(slip_epochs) != 1) * 1 # multiply with one to get from "False" and "True" to "0" and "1"
        dummy2 = np.where(dummy==1)
        slip_period_ends = np.append(slip_epochs[dummy2], np.array([slip_epochs[-1]]))
        n_slip_periods = np.sum(dummy) + 1 
        
        ## -- Slip_periods = zeros(n_slip_periods,2); 
        slip_periods = np.zeros([n_slip_periods,2])
        ## -- Store slip ends
        slip_periods[:,1] = slip_period_ends
        ## -- Store first slip start manually
        slip_periods[0,0] = slip_epochs[0]
       
        ## Insert remaining slip period starts
        for k in range(1,n_slip_periods):
            indx = [x+1 for x, val in enumerate(slip_epochs) if val == slip_periods[k-1, 1]]
            if len(indx) != 0:
                indx = indx[0]
                slip_periods[k, 0] = slip_epochs[indx]
             
    else:
        slip_periods = []
        n_slip_periods = 0

    return slip_periods, n_slip_periods


