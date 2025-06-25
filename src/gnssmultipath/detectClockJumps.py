"""
This moduling is for detecting reciever clock jumps.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""

import warnings
import numpy as np
warnings.filterwarnings(action='ignore')


def detectClockJumps(GNSS_obs, nGNSSsystems, obsCodes, time_epochs, tInterval, GNSSsystems):
    """
    Function that detects receiver clock jumps from GNSS observations.

    INPUTS:
    -------

    GNSS_obs:           Dict containing a arrays for each GNSS system.
                        Each matrix is a 3D matrix containing all
                        observation of current GNSS system for all epochs.
                        Order of obsType index is same order as in
                        obsCodes dict

    nGNSSsystems:       number of GNSS systems present

    obsCodes:           Dict that defines the observation
                        codes available for all GNSS system. Each cell
                        element is another cell containing the codes for
                        that GNSS system. Each element in this cell is a
                        string with three-characters. The first
                        character (a capital letter) is an observation code
                        ex. "L" or "C". The second character (a digit)
                        is a frequency code. The third character(a Capital letter)
                        is the attribute, ex. "P" or "X"

    time_epochs:        matrix containing gps-week and "time of week"
                        for each epoch
                        time_epochs(epoch,i),   i=1: week
                                                i=2: time-of-week in seconds (tow)

    tInterval:          observations interval; seconds.



    OUTPUTS:
    --------

    nClockJumps:            number of receiver clock jumps detected

    meanClockJumpInterval:  average time between receiver clock jumps,
                            seconds

    stdClockJumpInterval:   standard deviation of receiver clock
                            jump intervals

    """


    for i in range(0,nGNSSsystems):
        curr_sys = GNSSsystems[i+1]
        obsTypes = obsCodes[i+1][curr_sys]
        # curr_sys = GNSSsystems[1]
        # obsTypes = obsCodes[1][curr_sys]
        codeIndices = [idx for idx ,obstype in enumerate(obsTypes) if 'C' in obstype[0]]
        nCodeObs = len(codeIndices)
        max_sat_sys,max_ncodes = GNSS_obs[curr_sys][1].shape
        nepochs = len(GNSS_obs[curr_sys])
        current_obs = np.zeros([nepochs,max_ncodes,max_sat_sys-1])
        # current_obs = permute(GNSS_obs{i},[3 2 1]);
        for ep in range(0,nepochs):
            current_obs[ep,:,:] = np.transpose(GNSS_obs[curr_sys][ep+1][1::, None], (1,2,0))

        nepochs, _, nSat = current_obs.shape

        if i == 0:
            reshaped_obs = np.reshape(current_obs[:,codeIndices,:], (nepochs, nSat*nCodeObs), order='F') # 'F' is adding obscodes together (Fortrans-like indexing)
        else:
            reshaped_obs = np.append(
                reshaped_obs,
                np.reshape(current_obs[:,codeIndices,:], (nepochs, nSat*nCodeObs), order='F'),
                axis=1
            )

    obsChange =np.diff(reshaped_obs,axis=0)

    # code_jump_epochs = find(all((abs(obsChange)> 2e5 )| obsChange == 0,2));
    code_jump_epochs = []
    for ep in range(0,len(obsChange)):
        non_zero_columns = np.nonzero(obsChange[ep,:])[0]
        if np.all(abs(obsChange[ep,non_zero_columns]) > 2e5):
            code_jump_epochs.append(ep)

    time_diff = np.diff(time_epochs,axis=0)[:,1].reshape(len(np.diff(time_epochs,axis=0)[:,1]),1)


    if np.all(time_diff[:,0] == 0):
        time_diff = time_diff[:,1]

    # time_jump_epochs = find(abs(time_diff-tInterval)> 1e-4);
    time_jump_epochs = np.argwhere(abs(time_diff - tInterval) > 1e-4)

    jump_epochs = np.union1d(time_jump_epochs, code_jump_epochs)
    jump_epochs = np.unique(jump_epochs)

    nClockJumps = len(jump_epochs)
    meanClockJumpInterval = np.mean(np.diff(jump_epochs)*tInterval)
    stdClockJumpInterval = np.std(np.diff(jump_epochs)*tInterval)
    if np.isnan(meanClockJumpInterval):
        meanClockJumpInterval = 0
        stdClockJumpInterval = 0


    return nClockJumps, meanClockJumpInterval, stdClockJumpInterval
