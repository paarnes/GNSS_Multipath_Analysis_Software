import warnings
warnings.filterwarnings("ignore")

def getLLISlipPeriods(LLI_current_phase):
    """
    Function that sorts all ambiguity slips indicated by LLI in RINEX
    observation file.

     INPUTS:
     -------
     LLI_current_phase:        matrix. contains LLI indicators for all epochs
                               for current GNSS system and phase pbservation,
                               for all satellites.

                               LLI_current_phase(epoch, satID)

     OUTPUTS:
     --------
     LLI_slip_periods:         dict. One element for each satellite. Each
                               cell is a matrix. The matrix contains indicated
                               slip period start epochs in first column, and
                               end epochs in second column.
    """
    import numpy as np
    _, nSat = LLI_current_phase.shape
    LLI_slip_periods = {}

    # Itterrate through all satellites in current GNSS system
    for sat in range(0,nSat-1):
        LLI_current_sat = LLI_current_phase[:, sat+1]
        ## Get epochs where LLI indicate slip, for current satellite
        LLI_slips = np.array(ismember2([1, 2, 3, 5, 6, 7], LLI_current_sat)).reshape(len(ismember2([1, 2, 3, 5, 6, 7], LLI_current_sat)),1)
        # if there are slips
        if not len(LLI_slips) == 0:
            ## dummy is logical. It will be 1 at indices where the following
            ## slip epoch is NOT the epoch following the current slip epoch.
            ## These will therefor be the indices where slip periods end.
            ## The last slip end is not detected this way and is inserted manually
            dummy = np.diff(LLI_slips) !=1 * 1
            slip_period_ends = np.concatenate((LLI_slips[dummy], LLI_slips[-1]))
            n_slip_periods = np.sum(dummy) + 1
            current_slip_periods = np.zeros([n_slip_periods,2])
            # store slip ends
            current_slip_periods[:,1] = slip_period_ends
            # store first slip start manually
            # store first slip start manually
            current_slip_periods[0,0] = LLI_slips[0].item() if hasattr(LLI_slips[0], "item") else LLI_slips[0]
            # Insert remaining slip period starts
            for k in range(1,n_slip_periods):
                indx = next(x for x, val in enumerate(LLI_slips) if abs(val) ==  current_slip_periods[k-1, 1])
                current_slip_periods[k, 0] = LLI_slips[indx + 1]
        else:
            current_slip_periods = []

        LLI_slip_periods[sat] = current_slip_periods

    return LLI_slip_periods

def ismember2(LLI_codes,LLI_current_sat):
    """The function takes in arrays, and finds the indencies of where the LLI_current_sat matches LLI_codes"""
    indx = [i for i, e in enumerate(LLI_current_sat) if e in LLI_codes]
    return indx