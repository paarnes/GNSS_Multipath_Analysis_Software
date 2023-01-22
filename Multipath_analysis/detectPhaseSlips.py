def detectPhaseSlips(estimates, missing_obs_overview,epoch_first_obs, epoch_last_obs, tInterval, crit_slip_rate):
    """
     Function that detects epochs with cycle slips, given test estimates
     and a critical rate of change
    --------------------------------------------------------------------------------------------------------------------------
    
     INPUTS:            
    --------------------------------------------------------------------------------------------------------------------------
     estimates:     matrix containing estimates from a lin. combination
                           for a specific PRN
    
                           estimates(epoch, PRN)
    
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
    
     OUTPUTS
    
     slip_epochs:          array, contains epochs with detected cycle slip
    --------------------------------------------------------------------------------------------------------------------------
    """
    import numpy as np
    
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