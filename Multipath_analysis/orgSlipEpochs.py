def orgSlipEpochs(slip_epochs):
    """
    # Function that takes array of epochs with detected ambiguity slips and for
    # one satellite and orginses them into slip periods
    #--------------------------------------------------------------------------------------------------------------------------
    
    # INPUTS
    
    # slip_epochs:          array, contains epochs with detected ambiguity slip
    #--------------------------------------------------------------------------------------------------------------------------
    
    # OUTPUTS
    
    # slip_periods:         Matrix that contains the start of periods with ambiguity slips in
    #                       the first column and the ends of the same periods in
    #                       the second column.
    
    #                       phase_slip_periods(ambiguity_slip_priod, j),
    #                       j = 1, ambiguity period start
    #                       j = 2, ambiguity period ends
    
    # n_slip_periods:       amount of slip periods for current satellite
    
    #--------------------------------------------------------------------------------------------------------------------------
    """
    import numpy as np
    ## -- If no slips occurs there are no slip periods for this sat.
    # if ~isempty(slip_epochs)
    if len(slip_epochs) != 0:
       
        # dummy is logical. It will be 1 at indices where the following
        # slip epoch is NOT the epoch following the current slip epoch.
        # These will therefor be the indices where slip periods end.
        # The last slip end is not detected this way and is inserted
        # manually
        # dummy = diff(slip_epochs)~=1;
        dummy = (np.diff(slip_epochs) != 1) * 1 # multiply with one to get from "False" and "True" to "0" and "1"
        dummy2 = np.where(dummy==1)
        # slip_period_ends = [slip_epochs(dummy); slip_epochs(end)];
        # slip_period_ends = np.concatenate((slip_epochs[dummy], slip_epochs[-1]))
        # slip_period_ends = slip_epochs[dummy]
        # slip_period_ends = slip_epochs[dummy2]
        # slip_period_ends = np.append(slip_period_ends, np.array([slip_epochs[-1]]))
        # slip_period_ends = slip_epochs[dummy2]
        slip_period_ends = np.append(slip_epochs[dummy2], np.array([slip_epochs[-1]]))
        n_slip_periods = np.sum(dummy) + 1 
        
        ## --Slip_periods = zeros(n_slip_periods,2); 
        slip_periods = np.zeros([n_slip_periods,2])
        ## -- Store slip ends
        # slip_periods(:,2) = slip_period_ends
        slip_periods[:,1] = slip_period_ends
        # store first slip start manually
        # slip_periods(1,1) = slip_epochs(1);
        slip_periods[0,0] = slip_epochs[0]
       
        ## Insert remaining slip period starts
        #     for k = 2:n_slip_periods
        #        slip_periods(k, 1) = slip_epochs(find(slip_epochs == slip_periods(k-1, 2)) + 1)
        #     end
        # else
        #     slip_periods = []
        #     n_slip_periods = 0
        for k in range(1,n_slip_periods):
            # indx = next(x for x, val in enumerate(slip_epochs) if val == slip_periods[k-1, 1] + 1)
            indx = [x+1 for x, val in enumerate(slip_epochs) if val == slip_periods[k-1, 1]]
            if len(indx) != 0:
                indx = indx[0]
                slip_periods[k, 0] = slip_epochs[indx]
          
       
    else:
        slip_periods = []
        n_slip_periods = 0

    return slip_periods, n_slip_periods


