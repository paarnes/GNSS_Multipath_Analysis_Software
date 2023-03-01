from datetime import date,timedelta,datetime
from barylag import barylag
import numpy as np
from scipy.interpolate import BarycentricInterpolator
from Geodetic_functions import date2gpstime


def preciseOrbits2ECEF(sys, PRN, date_, dates, epochInterval, nEpochs, sat_positions, navGNSSsystems):
    """
      Function that finds positions of speficied satelite at nearest epochs.
      Then interpolates position at specified time from these epoch positions
    --------------------------------------------------------------------------------------------------------------------------
      INPUTS
    
      sys:              Satellite system, string, 
                          ex. "E" or "G"
    
      PRN:              Satellite identification number, integer
    
      date_:            array, date of specified epoch in form of [year, month, day, hour, min, sec]
    
      dates:            matrix. Each row contains date of one of the epochs in 
                        the SP3 orbit file. 
                        [nEpochs x 6]
    
    
      epochInterval:    interval of position epochs in SP3 file, seconds
    
      nEpochs:          number of position epochs in SP3 file, integer
    
      sat_positions:    Dict. Each cell elements contains position data for a
                        specific GNSS system. Order is defined by order of 
                        navGNSSsystems. Each cell element is another cell that 
                        stores position data of specific satellites of that 
                        GNSS system. Each of these cell elements is a matrix 
                        with [X, Y, Z] position of a epoch in each row.
    
                        sat_positions[GNSSsystemIndex][PRN][epoch, :] = [X, Y, Z]
    
      navGNSSsystems:   Dict. Contains char. Each string is a code for a
                        GNSS system with position data stored in sat_positions.
                        Must be one of: 'G', 'R', 'E', 'C'
    --------------------------------------------------------------------------------------------------------------------------
      OUTPUTS
    
      X, Y, Z:          ECEF coordinates of satellite at desired time computed
                        by interpolation 
    --------------------------------------------------------------------------------------------------------------------------
    """

    
    ## Formating the date 
    date_ = str(date_[0]) + "/" + str(date_[1]) + "/" + str(date_[2]) + " " + str(int(date_[3])) \
        + ":" + str(int(date_[4])) + ":" + str(date_[5])[0:9]
    
    ## Need at tray/except for handelig date that have more decimals in the seconds. 
    ## Ex date_ should look like 2020/10/30 13:22:59.000000, but sometimes i could look like 2020/10/30 13:23:3.9262886 (RINEX from RINEX ON Andriod app)
    try:
        date_ = datetime.strptime(date_, "%Y/%m/%d %H:%M:%S.%f")
    except:
        sec_dum = date_.split(':')[-1]
        date_ = date_.replace(sec_dum, str(format(float(sec_dum), '.6f'))) # removing decimals if more than 3
        date_ = datetime.strptime(date_, "%Y/%m/%d %H:%M:%S.%f")
    
    date_ = datetime.timestamp(date_) # convert to timestamp
    
    ## --Degree of lagrange polynomial to be used
    # lagrangeDegree = 3
    lagrangeDegree = 7 # changed 27.01.2023

    ## -- Amount of nodes
    nNodes = lagrangeDegree + 1
    
    ## -- GNSSsystemIndex = find(navGNSSsystems == sys, 1);
    GNSSsystemIndex = [idx for idx,val in enumerate(navGNSSsystems) if navGNSSsystems[idx]==sys][0] 
    curr_sys = navGNSSsystems[GNSSsystemIndex]

    ## --Date of first epoch
    tFirstEpoch = np.array(dates[0, :]).astype(float)
    tFirstEpoch = str(int(tFirstEpoch[0])) + "/" + str(int(tFirstEpoch[1])) + "/" + str(int(tFirstEpoch[2])) + " " + str(int(tFirstEpoch[3])) \
        + ":" + str(int(tFirstEpoch[4])) + ":" + str(tFirstEpoch[5])[0:9] 
        
    tFirstEpoch = datetime.strptime(tFirstEpoch, "%Y/%m/%d %H:%M:%S.%f")
    tFirstEpoch = datetime.timestamp(tFirstEpoch) # convert to timestamp
    
    tk = date_ - tFirstEpoch 
    
    ## -- Closest node before desired time
    closestEpochBeforeIndex = np.floor(tk/epochInterval) + 1
    
    # Get the index of the first and last node. If there is not enough epochs
    # before or after desired time the degree of lagrange polynomial i reduces,
    # as well as number of nodes
    node1EpochIndex = closestEpochBeforeIndex - min(nNodes/2 - 1, closestEpochBeforeIndex-1);
    diff1 = node1EpochIndex - (closestEpochBeforeIndex - (nNodes/2 - 1));
    
    node8EpochIndex = closestEpochBeforeIndex + min(nNodes/2, nEpochs - closestEpochBeforeIndex);
    diff2 = node8EpochIndex - (closestEpochBeforeIndex + nNodes/2);
    
    node1EpochIndex = int(node1EpochIndex - diff2) -1 # - 1 because null-indexed
    node8EpochIndex = int(node8EpochIndex - diff1) -1 
    
    ## Indices of node epochs
    nodeEpochs = np.array(range(node1EpochIndex,node8EpochIndex+1))
    
    ## Reduce number of nodes if necessary
    nNodes = int(nNodes - diff1*2 + diff2*2)
    
    # Get positions at each node and relative time
    nodePositions = np.zeros([nNodes, 3])
    nodeTimes = np.zeros([nNodes, 1])
    for i in np.arange(0,nNodes):
        try: 
            nodePositions[i, :] = sat_positions[curr_sys][nodeEpochs[i]][PRN]
        except:
            nodePositions[i, :] = np.nan
        date_dum = datetime(year=int(dates[nodeEpochs[i], :][0]), month= int(dates[nodeEpochs[i], :][1]), 
                            day = int(dates[nodeEpochs[i], :][2]), hour = int(dates[nodeEpochs[i], :][3]), 
                            minute = int(dates[nodeEpochs[i], :][4]),second = int(dates[nodeEpochs[i], :][5][0:1]))

        
        date_dum = datetime.timestamp(date_dum) # convert to timestamp
        nodeTimes[i] = date_dum - tFirstEpoch
    
    ## -- Interpolate new posistion of satellite using a lagrange polynomial
    try: 
        X=barylag(np.hstack([nodeTimes, nodePositions[:, 0].reshape(len(nodePositions[:, 0]),1)]), tk)
        Y=barylag(np.hstack([nodeTimes, nodePositions[:, 1].reshape(len(nodePositions[:, 1]),1)]), tk)
        Z=barylag(np.hstack([nodeTimes, nodePositions[:, 2].reshape(len(nodePositions[:, 2]),1)]), tk)
    except:
        X = np.nan
        Y = np.nan
        Z = np.nan
    return X, Y, Z