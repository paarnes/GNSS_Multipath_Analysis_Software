"""
This module is performing an analysis of the different GNSS signals.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""

import numpy as np
from gnssmultipath.estimateSignalDelays import estimateSignalDelays
from gnssmultipath.getLLISlipPeriods import getLLISlipPeriods
from gnssmultipath.computeDelayStats import computeDelayStats


def signalAnalysis(currentGNSSsystem, range1_Code, range2_Code, GNSSsystems, frequencyOverview, nepochs, \
    tInterval, current_max_sat, current_GNSS_SVs, current_obsCodes, current_GNSS_obs, current_GNSS_LLI, current_sat_elevation_angles,\
    phaseCodeLimit, ionLimit, cutoff_elevation_angle):
    """
    Function that executes a signal analysis on a specific GNSS code range
    signal for a specific GNSS system.


    INPUTS:
    -------

    currentGNSSsystem:          string. Code that gives current GNSS
                                system. Ex: "G" or "E"

    range1_Code:                string. obs code for first code pseudorange
                                observation

    range2_Code:                string. obs code for second code pseudorange
                                observation

    GNSSsystems:                List containing codes of GNSS systems.
                                Elements are strings.
                                ex. "G" or "E"

    frequencyOverview:          dict. each elements contains carrier band
                                frequenies for a specific GNSS system.
                                Order is set by GNSSsystems. If system is
                                GLONASS then element is matrix where each
                                row i gives carrier band i frequencies for
                                all GLONASS SVs. If not GLONASS, element is
                                array with one frequency for each carrier
                                band.

    nepochs:                    nepochs of observations in RINEX observation file.

    tInterval:                  observation interval in seconds

    current_max_sat:            max PRN number of current GNSS system

    current_GNSS_SVs:           Dict containing a matrix for current GNSS system.
                                Each matrix contains number of satellites with
                                obsevations for each epoch, and PRN for those
                                satellites

                                GNSS_SVs(epoch, j)
                                        j=1: number of observed satellites
                                        j>1: PRN of observed satellites

    current_obsCodes:           Dict that defines the observation
                                codes available in current GNSS system.
                                Each element in this cell is a
                                string with three-characters. The first
                                character (a capital letter) is an observation code
                                ex. "L" or "C". The second character (a digit)
                                is a frequency code. The third character(a Capital letter)
                                is the attribute, ex. "P" or "X"

    current_GNSS_obs:           3D matrix  containing all observation of
                                current GNSS system for all epochs.
                                Order of obsType index is same order as in
                                current_obsCodes

                                current_GNSS_obs(PRN, obsType, epoch)
                                            PRN: int
                                            ObsType: int: 1,2,...,numObsTypes
                                            epoch: int

    current_GNSS_LLI:            matrix  containing all Loss of Lock
                                 indicators of current GNSS system for all epochs.
                                 Order of obsType index is same order as in
                                 current_obsCodes

                                 current_GNSS_LLI(PRN, obsType, epoch)
                                            PRN: int
                                            ObsType: int: 1,2,...,numObsTypes
                                            epoch: int

    current_sat_elevation_angles: matrix contaning satellite elevation angles
                                  at each epoch, for current GNSS system.

                                  sat_elevation_angles(epoch, PRN)

    phaseCodeLimit:               critical limit that indicates cycle slip for
                                  phase-code combination. Unit: m/s. If set to 0,
                                  default value will be used

    ionLimit:                     critical limit that indicates cycle slip for
                                  the rate of change of the ionopheric delay.
                                  Unit: m/s. If set to 0, default value will be used

    cutoff_elevation_angle        Critical cutoff angle for satellite elevation angles, degrees
                                  Estimates where satellite elevation angle
                                  is lower than cutoff are removed, so are
                                  estimated slip periods

    OUTPUTS:
    --------
    currentStats:                 dict. Contains statitics from analysis
                                  executed. More detail of each stattistic
                                  is given in function computeDelayStats.m

    success:                      boolean. 1 if no errors thrown, 0 otherwise

    """

    ## --Get corrosponding phase codes to the range codes
    phase1_Code = "L" + range1_Code[1::]
    phase2_Code = "L" + range2_Code[1::]
    # Get current GNSS system index
    GNSSsystemIndex = [k for k in GNSSsystems if GNSSsystems[k]==currentGNSSsystem][0]
    ## -- Get frequencies of carrier bands. These are arrays if current GNSS system is GLONASS. Each element for a specific GLONASS SV.
    if currentGNSSsystem == 'R':
        carrier_freq1 = frequencyOverview[GNSSsystemIndex][int(range1_Code[1])-1, :]
        carrier_freq2 = frequencyOverview[GNSSsystemIndex][int(range2_Code[1])-1, :]
    else:
        carrier_freq1 = frequencyOverview[GNSSsystemIndex][int(range1_Code[1])-1, :][0]
        carrier_freq2 = frequencyOverview[GNSSsystemIndex][int(range2_Code[1])-1, :][0]
    # Run function to compute estimates of ionospheric delay, multipath delays slip periods of range1 signal.
    ion_delay_phase1, multipath_range1, range1_slip_periods,ambiguity_slip_periods ,range1_observations, phase1_observations, success = estimateSignalDelays(range1_Code, range2_Code, \
        phase1_Code, phase2_Code, carrier_freq1, carrier_freq2,nepochs, current_max_sat,\
          current_GNSS_SVs, current_obsCodes, current_GNSS_obs, currentGNSSsystem, tInterval, phaseCodeLimit, ionLimit)
    # Create logical mask for epochs where sat elevation is lower than cutoff or missing
    cutoff_elevation_mask = current_sat_elevation_angles.copy()
    cutoff_elevation_mask[cutoff_elevation_mask < cutoff_elevation_angle] = 0
    cutoff_elevation_mask[cutoff_elevation_mask >= cutoff_elevation_angle] = 1
    #  Apply satellite elevation cutoff mask to estimates
    ion_delay_phase1 = ion_delay_phase1 * cutoff_elevation_mask
    multipath_range1 = multipath_range1 * cutoff_elevation_mask
    range1_observations = range1_observations * cutoff_elevation_mask
    phase1_observations = phase1_observations * cutoff_elevation_mask

    # Remove estimated slip periods (range_1 slips) if satellite elevation angle was lower than cutoff or missing.
    for sat in np.arange(0,len(range1_slip_periods)):
        current_sat_slip_periods = np.array(range1_slip_periods[sat+1]).astype(int)
        if len(current_sat_slip_periods) > 0:
            n_slip_periods,_ = current_sat_slip_periods.shape
            n_slips_removed = 0
            for slip_period in np.arange(0,n_slip_periods):
                if cutoff_elevation_mask[current_sat_slip_periods[slip_period - n_slips_removed, 0], sat] == 0 \
                    or cutoff_elevation_mask[current_sat_slip_periods[slip_period - n_slips_removed, 1], sat] == 0:

                    current_sat_slip_periods = np.delete(current_sat_slip_periods, slip_period - n_slips_removed, axis=0)
                    n_slips_removed = n_slips_removed + 1

            range1_slip_periods[sat+1] = current_sat_slip_periods

    ## -- Remove estimated slip periods (both ionspher residuals and code phase) if satellite elevation angle was lower than cutoff or missing.
    for sat in np.arange(0,len(ambiguity_slip_periods)):
        current_sat_slip_periods = np.array(ambiguity_slip_periods[sat+1]).astype(int)
        if len(current_sat_slip_periods) > 0:
            n_slip_periods,_ = current_sat_slip_periods.shape
            n_slips_removed = 0
            for slip_period in np.arange(0,n_slip_periods):
                if cutoff_elevation_mask[current_sat_slip_periods[slip_period - n_slips_removed, 0], sat] == 0 \
                    or cutoff_elevation_mask[current_sat_slip_periods[slip_period - n_slips_removed, 1], sat] == 0:

                    current_sat_slip_periods = np.delete(current_sat_slip_periods, slip_period - n_slips_removed, axis=0)
                    n_slips_removed = n_slips_removed + 1

            ambiguity_slip_periods[sat+1] = current_sat_slip_periods

    if not success:
        currentStats = np.nan
        return currentStats

    # Compute slips from LLI in rinex file
    max_sat = len(current_GNSS_obs[1])
    LLI_current_phase =  np.zeros([nepochs,max_sat])
    for ep in np.arange(0, nepochs):
        LLI_current_dum = np.array(current_GNSS_LLI[ep+1][:,ismember(current_obsCodes[currentGNSSsystem],phase1_Code)]).reshape(1, len(current_GNSS_LLI[ep+1][:,ismember(current_obsCodes[currentGNSSsystem],phase1_Code)]))
        LLI_current_phase[ep,:] = np.squeeze(LLI_current_dum)
    LLI_slip_periods = getLLISlipPeriods(LLI_current_phase)

    # Compute statistics of estimates
    mean_multipath_range1, overall_mean_multipath_range1,\
        rms_multipath_range1, average_rms_multipath_range1,\
        mean_ion_delay_phase1, overall_mean_ion_delay_phase1, mean_sat_elevation_angles, nEstimates, nEstimates_per_sat,\
        nRange1Obs_Per_Sat, nRange1Obs, range1_slip_distribution_per_sat, range1_slip_distribution,ambiguity_slip_distribution_per_sat, ambiguity_slip_distribution,LLI_slip_distribution_per_sat, LLI_slip_distribution, \
        combined_slip_distribution_per_sat, combined_slip_distribution, elevation_weighted_rms_multipath_range1, \
        elevation_weighted_average_rms_multipath_range1 = \
        computeDelayStats(ion_delay_phase1, multipath_range1, current_sat_elevation_angles,range1_slip_periods,ambiguity_slip_periods,LLI_slip_periods, range1_observations, tInterval)

    currentStats = {'mean_multipath_range1_satellitewise' : mean_multipath_range1,
                    'mean_multipath_range1_overall' : overall_mean_multipath_range1,
                    'rms_multipath_range1_satellitewise' : rms_multipath_range1,
                    'rms_multipath_range1_averaged' : average_rms_multipath_range1,
                    'mean_ion_delay_phase1_satellitewise' : mean_ion_delay_phase1,
                    'mean_ion_delay_phase1_overall' : overall_mean_ion_delay_phase1,
                    'mean_sat_elevation_angles' : mean_sat_elevation_angles,
                    'nEstimates' : nEstimates,
                    'nEstimates_per_sat' : nEstimates_per_sat,
                    'n_range1_obs_per_sat': nRange1Obs_Per_Sat,
                    'nRange1Obs' : nRange1Obs,
                    'range1_slip_distribution_per_sat' : range1_slip_distribution_per_sat,
                    'range1_slip_distribution' : range1_slip_distribution,
                    'cycle_slip_distribution_per_sat' : ambiguity_slip_distribution_per_sat,
                    'cycle_slip_distribution' : ambiguity_slip_distribution,
                    'LLI_slip_distribution_per_sat' : LLI_slip_distribution_per_sat,
                    'LLI_slip_distribution' : LLI_slip_distribution,
                    'slip_distribution_per_sat_LLI_fusion' : combined_slip_distribution_per_sat,
                    'slip_distribution_LLI_fusion' : combined_slip_distribution,
                    'elevation_weighted_rms_multipath_range1_satellitewise' : elevation_weighted_rms_multipath_range1,
                    'elevation_weighted_average_rms_multipath_range1' :  elevation_weighted_average_rms_multipath_range1,
                    'range1_observations' : range1_observations,
                    'phase1_observations' : phase1_observations
                    }

    ## -- Store estimates needed for plotting
    currentStats['ion_delay_phase1'] = ion_delay_phase1
    currentStats['multipath_range1'] = multipath_range1
    currentStats['sat_elevation_angles'] = current_sat_elevation_angles
    ## -- Store codes
    currentStats['range1_Code'] = range1_Code
    currentStats['range2_Code'] = range2_Code
    currentStats['phase1_Code'] = phase1_Code
    currentStats['phase2_Code'] = phase2_Code
    ## -- Store slips
    currentStats['range1_slip_periods'] = range1_slip_periods
    currentStats['cycle_slip_periods']  = ambiguity_slip_periods

    return currentStats, success


def ismember(list_,code):
    """
    The function takes in a string and a list, and finds the index of
    """
    indx = [idx for idx, val in enumerate(list_) if val == code]
    if indx != []:
        indx = indx[0]
    return indx
