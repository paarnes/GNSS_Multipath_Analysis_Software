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

    ## -- Calculate rate of change of estimates of ionospheric delay (time derivative)
    estimates_rate_of_change = np.diff(estimates,axis=0)/tInterval

    ## -- Detect slips for epochs with either estimates_rate_of_change
    ## higher than critical value, or epochs with missing estimates
    condition_met = np.abs(estimates_rate_of_change) > crit_slip_rate  # Create a boolean array where True indicates the condition is met
    slips_from_crit_rate = np.where(condition_met)


    slip_epochs = {str(key): [] for key in range(1, estimates.shape[1])}
    if slips_from_crit_rate[0].size > 0:
        for key, value in zip(slips_from_crit_rate[1].astype(str).tolist(), slips_from_crit_rate[0]):
            slip_epochs[key].append(value)


    slips_from_missing_obs_curr_sat = []
    if epoch_first_obs.size > 0:
        for PRN in np.arange(len(epoch_first_obs)):
            epoch_first_obs_temp = epoch_first_obs[PRN].astype(int) if not np.isnan(epoch_first_obs[PRN]) else None
            epoch_last_obs_temp = epoch_last_obs[PRN].astype(int) if not np.isnan(epoch_last_obs[PRN]) else None

            if any(item is None for item in [epoch_first_obs_temp, epoch_last_obs_temp]):
                continue
            else:
                slips_from_missing_obs_curr_sat = (np.where(missing_obs_overview[epoch_first_obs_temp:epoch_last_obs_temp,PRN] == 1) + epoch_first_obs_temp)[0].tolist()

            slips_from_crit_rate_curr = slip_epochs[str(PRN)]
            if len(slips_from_missing_obs_curr_sat) != 0:
                slip_epochs[str(PRN)].extend(list(sorted(set(slips_from_crit_rate_curr + slips_from_missing_obs_curr_sat))))
            else:
                slip_epochs[str(PRN)].extend(list(sorted(set(slips_from_missing_obs_curr_sat))))

    # Ensure no duplicates
    for key, value in slip_epochs.items():
        slip_epochs[key] = list(sorted(set(value)))

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


