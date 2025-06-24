"""
Module for reading RINEX observation files in v2 and v3.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""
import time
import os
import re
from datetime import datetime
import numpy as np
from tqdm import tqdm
from gnssmultipath.Geodetic_functions import date2gpstime


global tFirstObs

def readRinexObs(filename, readSS=None, readLLI=None, includeAllGNSSsystems=None,includeAllObsCodes=None, \
                                      desiredGNSSsystems=None, desiredObsCodes=None, desiredObsBands=None):
    """
    Function that chooses which function to use based on header info.
    """

    fid = open(filename,'r')
    if os.stat(filename).st_size == 0:
        raise ValueError('ERROR: This file seems to be empty')
    line = fid.readline().rstrip()
    rinexVersion = line[0:9].strip()
    if '2' in rinexVersion.split('.')[0]:
        GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
            obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
            rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success=  readRinexObs211(filename, readSS=None, readLLI=None, includeAllGNSSsystems=None,includeAllObsCodes=None, \
                            desiredGNSSsystems=desiredGNSSsystems, desiredObsCodes=None, desiredObsBands=None) ## WHEN A SOUTION IS FOUND ON desiredGNSSsystems, =None must be removed.
    else:
        GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
            obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
            rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success =  readRinexObs304(filename, readSS, readLLI, includeAllGNSSsystems,includeAllObsCodes, \
                                desiredGNSSsystems, desiredObsCodes, desiredObsBands)

    return GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
        obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
        rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success



def readRinexObs304(filename, readSS=None, readLLI=None, includeAllGNSSsystems=None,includeAllObsCodes=None, \
                    desiredGNSSsystems=None, desiredObsCodes=None, desiredObsBands=None):
    """
    Program/function to read GNSS observations in RINEX 3.04 observation files
    The main core of the program is 4 functions:
                                  rinexReadObsFileHeader304
                                  rinexReadObsBlockHead304
                                  rinexReadObsBlock304
                                  rinexFindNEpochs304


    To export every parameter use this code:

    GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
         obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
         rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success = \
         readRinexObs304(rinObsFilename)

    Tip: Unwanted outputs can be ignored using _

    --------------------------------------------------------------------------------------------------------------------------
    INPUTS

    filename:                 path and name of RINEX 3.04 observation file,
                              string

    readSS:                   Boolean, 0 or 1.
                              1 = read "Signal Strength" Indicators
                              0 = do not read "Signal Strength" Indicators

    readLLI:                  Boolean, 0 or 1.
                              1 = read "Loss-Of-Lock Indicators"
                              0 = do not read "Loss-Of-Lock Indicators"

    includeAllGNSSsystems:    Boolean, 0 or 1.
                              1 = include alle GNSS systems(GPS, GLONASS, Galieo, BeiDou)
                              0 = include only GNSSsystems specified in desiredGNSSsystems

    includeAllObsTypes:       Boolean, 0 or 1.
                              1 = include all valid ObsTypes
                              0 = include only ObsTypes specified in desiredObsTypes

    desiredGNSSsystems:       array og strings containing  codes of desired
                              GNSSsystems to be included,
                              ex. ["G", "E"]
                              OBS: Must be string array, NOT char vector

    desiredObsTypes:          array of strings containing desired ObsTypes to be
                              included, ex. ["C", "L", "S", "D"]
                              OBS: Must be string array, NOT char vector

    desiredObsBands:          array of desired obs Bands to be included,
                              ex [1, 5]
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS

    GNSS_obs:                 cell containing a matrix for each GNSS system.
                              Each matrix is a 3D matrix containing all
                              observation of current GNSS system for all epochs.
                              Order of obsType index is same order as in
                              obsCodes cell

                              GNSS_obs{GNSSsystemIndex}(PRN, obsType, epoch)
                                              GNSSsystemIndex: double,
                                              1,2,...,numGNSSsystems
                                              PRN: double
                                              ObsType: double: 1,2,...,numObsTypes
                                              epoch: double

    GNSS_LLI:                 cell containing a matrix for each GNSS system
                              Each matrix stores loss of lock indicators for
                              each epoch for that GNSS system

    GNSS_SS:                  cell containing a matrix for each GNSS system
                              Each matrix stores signal strength indicators
                              for each epoch for that GNSS system

    GNSS_SVs:                 cell containing a matrix for each GNSS system.
                              Each matrix contains number of satellites with
                              obsevations for each epoch, and PRN for those
                              satellites

                              GNSS_SVs{GNSSsystemIndex}(epoch, j)
                                              j=1: number of observed satellites
                                              j>1: PRN of observed satellites

    time_epochs:              matrix conatining gps-week and "time of week"
                              for each epoch
                              time_epochs(epoch,i),   i=1: week
                                                      i=2: time-of-week in seconds (tow)

    nepochs:                  number of epochs with observations in rinex observation file.

    GNSSsystems:              cell array containing codes of GNSS systems included
                              in RINEX observationfile. Elements are strings.
                              ex. "G" or "E"

    obsCodes:                 Cell that defines the observation
                              codes available for all GNSS system. Each cell
                              element is another cell containing the codes for
                              that GNSS system. Each element in this cell is a
                              string with three-characters. The first
                              character (a capital letter) is an observation code
                              ex. "L" or "C". The second character (a digit)
                              is a frequency code. The third character(a Capital letter)
                              is the attribute, ex. "P" or "X"

    approxPosition:           array containing approximate position from rinex
                              observation file header. [X, Y, Z]

    max_sat:                  array conataining max PRN number for each GNSS
                              system. Follows same order as GNSSsystems

    tInterval:                observations interval; seconds.

    markerName:               name of the antenna marker; '' if not specified

    rinexVersion:             string. rinex observation file version

    recType:                  receiver type, char vector

    timeSystem:               three-character code string of the time system
                              used for expressing tfirstObs;
                              can be GPS, GLO or GAL;

    leapSec:                  number of leap seconds since 6-Jan-1980.
                              UTC=GPST-leapSec. NaN by default.
                              THIS IS RINEX 3.04 OPTIONAL DATA

                              rinexHeader: cell column-vector containing the
                              following data:
                              rinexVersion:   RINEX version number; string.
                                              '' if not specified
                              rinexType:      RINEX file type; char

    gnssType:                 GNSS system of the satellites observed; can be
                              'G', 'R', 'E', 'C' or 'M' that stand for
                              GPS, GLONASS, GALILEO, BeiDou or Mixed ; char

    rinexProgr:               name of the software used to produce de RINEX
                              GPS obs file; '' if not specified


    rinexDate:                date/time of the RINEX file creation; '' if not
                              specified; char


    antDelta:                 column vector ot the three components of the
                              distance from the marker to the antenna,
                              in the following order - up, east and north;
                              null vector by default

    tFirstObs:                time stamp of the first observation record in the RINEX
                              observations file; column vector
                              [YYYY; MM; DD; hh; mm; ss.sssssss];
                              THIS IS CRITICAL DATA

    tLastObs:                 time stamp of the last observation record in the RINEX
                              observations file; column vector
                              [YYYY, MM, DD, hh, mm,ss.sssssss]. NaN by default.
                              THIS IS RINEX 3.04 OPTIONAL DATA

    clockOffsetsON:           receiver clock offsets flag. O if no realtime-derived
                              receiver clock offset was applied to epoch,
                              code and phase data (in other words, if the
                              file only has raw data), 1 otherwise. 0 by default.
                              THIS IS RINEX 3.04 OPTIONAL DATA

    GLO_Slot2ChannelMap:      map container that maps GLONASS slot numbers to
                              their respective channel number.
                              GLO_Slot2ChannelMap(slotnumber)

    success:                  Boolean. 1 if the reading of the RINEX
                              observations file seems to be successful,
                              0 otherwise
    --------------------------------------------------------------------------------------------------------------------------

    ADVICE: The function rinexFindNEpochs() calculates the amount of observation epochs in
    advance. This calculation will be incredibly more effective if TIME OF
    LAST OBS is included in the header of the observation file. It is
    strongly advized to manually add this information to the header if it is
    not included by default.
    --------------------------------------------------------------------------------------------------------------------------

    According to RINEX 3.04 the observation type codes are:
    Observation code
            C: Pseudorange
    GPS: C/A, L2C
    Glonass: C/A
    Galileo: All
    L: Carrier phase
    D: Doppler frequency
    S: Raw signal strengths or SNR values as given by the receiver for the
    respective phase observations
    I: Ionosphere phase delay
    X: Receiver channel numbers

    Frequency code
    GPS Glonass Galileo SBAS
    1: L1 G1 E1 B1    (GPS,QZSS,SBAS,BDS)
    2: L2 G2 B1-2     (GLONASS)
    4: G1a            (Galileo)
    5: L5 E5a B2/B2a  (GPS, QZSS, SBAS, IRNSS)
    6: L6 E6 B3 G2a   (Galileo, QZSS, BDS, GLONASS)
    7: E5b B2/B2b     (Galileo)
    8: E5a+b E5a+b    (Galileo, BDS)
    9: S              (IRNSS)
    0: for type X     (all)

    Attribute:
    A = A channel     (Galileo,IRNSS,GLONASS)
    B = B channel     (Galileo,IRNSS,GLONASS)
    C = C channel     (Galiloe, IRNSS)
    C code-based  (SBAS, GPS, GLONASS, QZSS)
    D = Semi-codelss  (GPS)

    I = I channel     (GPS, Galileo, QZSS, BDS)
    L = L channel     (L2C GPS, QZSS)
    P channel     (GPS. QZSS)
    M = M code-based  (GPS)
    N = Codeless      (GPS)
    P = P code-based  (GPS, GLONASS)
    Pilot channel (BDS)

    Q = Q channel     (GPS, Galileo, QZSS, BDS)
    S = D channel     (GPS, Galileo, QZSS, BDS)
    M channel     (L2C GPS, QZSS)

    W = Based on Z-tracking (GPS)
    X = B+C channels  (Galileo, IRNSS)
    I+Q channels  (GPS, IRNSS)
    M+L channels  (GPS, QZSS)
    D+P channels  (GPS, QZSS, BDS)

    Y = Y code based  (GPS)
    Z = A+B+C channels(Galileo)
    D+P channels  (BDS)
    --------------------------------------------------------------------------------------------------------------------------
    """
    ## -- Setting None arguments
    if readSS is None:
        readSS = 1
    if readLLI is None:
        readLLI = 1
    if includeAllGNSSsystems is None and desiredGNSSsystems is None:
        includeAllGNSSsystems = 1
    if includeAllObsCodes is None and desiredObsCodes is None:
        includeAllObsCodes = 1
    if desiredGNSSsystems is None:
        desiredGNSSsystems = ['G','R','E','C']
    if desiredObsCodes is None:
        desiredObsCodes = ['C','L','S','D']
    if desiredObsBands is None:
        desiredObsBands = list(np.arange(1,10))

    ## Get the start time
    t = time.process_time()

    ## - Initialize variables in case of input error
    GNSS_obs       = np.nan
    GNSS_LLI       = np.nan
    GNSS_SS        = np.nan
    GNSS_SVs       = np.nan
    time_epochs    = np.nan
    nepochs        = np.nan
    GNSSsystems    = np.nan
    obsCodes       = np.nan
    approxPosition = np.nan
    max_sat        = np.nan
    tInterval      = np.nan
    markerName     = np.nan
    rinexVersion   = np.nan
    recType        = np.nan
    timeSystem     = np.nan
    leapSec        = np.nan
    gnssType       = np.nan
    rinexProgr     = np.nan
    rinexDate      = np.nan
    antDelta       = np.nan
    tFirstObs      = np.nan
    tLastObs       = np.nan
    clockOffsetsON = np.nan
    GLO_Slot2ChannelMap = np.nan

    ### --- Dict for storing data
    GNSS_obs = {}
    GPS = {}
    GLONASS = {}
    Galileo = {}
    BeiDou = {}
    GPS_LLI = {}
    GLONASS_LLI = {}
    Galileo_LLI = {}
    BeiDou_LLI = {}
    GPS_SS = {}
    GLONASS_SS = {}
    Galileo_SS = {}
    BeiDou_SS = {}

    ## -- Test if readSS is boolean
    if readSS!=1 and readSS!=0:
        print('INPUT ERROR(readRinexObs304): The input argument readSS must be either 1 or 0')
        success = 0
        return


    ## -- Test if readLLI is boolean
    if readLLI!=1 and readLLI!=0:
        print('INPUT ERROR(readRinexObs304): The input argument readLLI must be either 1 or 0')
        success = 0
        return

    max_GPS_PRN     = 36 # Max number of GPS PRN in constellation
    max_GLONASS_PRN = 36 # Max number of GLONASS PRN in constellation
    max_Galileo_PRN = 36 # Max number of Galileo PRN in constellation
    max_Beidou_PRN  = 100 # Max number of BeiDou PRN in constellation

    ## -- Read header of observation file
    [success, rinexVersion, gnssType, markerName, recType, antDelta,\
    GNSSsystems,numOfObsCodes, obsCodes, obsCodeIndex,tFirstObs, tLastObs, tInterval, \
    timeSystem, _, clockOffsetsON, rinexProgr, rinexDate,leapSec, approxPosition, GLO_Slot2ChannelMap, _, fid] = \
    rinexReadObsFileHeader304(filename, includeAllGNSSsystems, includeAllObsCodes,desiredGNSSsystems, desiredObsCodes, desiredObsBands)

    if success==0:
        return

    ## -- Compute number of epochs with observations
    nepochs, tLastObs, tInterval, success =\
        rinexFindNEpochs304(filename, tFirstObs, tLastObs, tInterval) #computes number of epochs in observation file

    if success==0:
        return

    nGNSSsystems = len(GNSSsystems)
    GNSS_SVs = {}
    max_sat  =  np.zeros([nGNSSsystems,1])
    t_week = []
    t_tow = []
    GNSSsystems_full_names =  [""]*nGNSSsystems
    GNSS_LLI = {}
    GNSS_SS = {}

    # Create array for max_sat. Initialize cell elements in dicts
    for k in np.arange(0,nGNSSsystems):
        if GNSSsystems[k+1] == 'G':
            max_sat[k] = max_GPS_PRN
            GNSS_SVs['G'] = np.zeros([nepochs, int(max_sat[k].item() + 1)])
            GNSSsystems_full_names[k] = "GPS"
        elif GNSSsystems[k+1] == 'R':
            max_sat[k] = max_GLONASS_PRN
            GNSS_SVs['R'] = np.zeros([nepochs, int(max_sat[k].item() + 1)])
            GNSSsystems_full_names[k] = "GLONASS"

        elif GNSSsystems[k+1] == 'E':
            max_sat[k] = max_Galileo_PRN
            GNSS_SVs['E'] = np.zeros([nepochs, int(max_sat[k].item() + 1)])
            GNSSsystems_full_names[k] = "Galileo"

        elif GNSSsystems[k+1] == 'C':
            max_sat[k] = max_Beidou_PRN
            GNSS_SVs['C'] = np.zeros([nepochs, int(max_sat[k].item() + 1)])
            GNSSsystems_full_names[k] = "BeiDou"
        else:
            print(f'ERROR(readRinexObs304): Only following GNSS systems are compatible with this program: GPS, GLONASS, Galileo, Beidou. {GNSSsystems[k]} is not valid')
            return


        curr_sys = GNSSsystems[k+1]
        max_sat_k = int(max_sat[k].item()) if hasattr(max_sat[k], "item") else int(max_sat[k])
        GNSS_obs[curr_sys] = np.zeros([max_sat_k, numOfObsCodes[k], nepochs])


        # Preallocation LLI and SS
        if readLLI:
            GNSS_LLI[curr_sys] = np.zeros([nGNSSsystems,1])
        else:
            GNSS_LLI[curr_sys] = np.nan
        if readSS:
            GNSS_SS[curr_sys] = np.zeros([nGNSSsystems,1])
        else:
            GNSS_SS[curr_sys] = np.nan


    GNSS_names = dict(zip(['G', 'R', 'E', 'C'],['GPS','GLONASS','Galileo','Beidou']))
    current_epoch      = 0
    ## -- Initialize progress bar
    n_update_break = int(np.floor(nepochs/10)) #number of epoch before updating progressbar
    bar_format = '{desc}: {percentage:3.0f}%|{bar}| ({n_fmt}/{total_fmt})'
    # with tqdm(total=100,desc ="Rinex observations are being read" , position=0, leave=True) as pbar:
    with tqdm(total=100,desc ="Rinex observations are being read" , position=0, leave=True, bar_format=bar_format) as pbar:
        while 1:
            success, _, _, date, numSV, eof = rinexReadObsBlockHead304(fid) # Read Obs Block Header
            if success==0 or eof==1:
                break
            success, Obs,SVlist, numSV, LLI, SS, eof = rinexReadObsBlock304(fid, numSV, numOfObsCodes, GNSSsystems, obsCodeIndex, readSS, readLLI) # Read current block of observations
            if success ==0 or eof==1:
                break
            current_epoch += 1
            if np.mod(current_epoch, n_update_break) == 0:  # Update progress bar every n_update_break epochs
                pbar.update(10)

            ## Convert date to GPS-week and "time-of-week"
            week, tow = date2gpstime(int(date[0]), int(date[1]), int(date[2]), int(date[3]), int(date[4]), int(date[5]))
            ## Store GPS-week and "time-of-week" of current epoch
            t_week.append(week)
            t_tow.append(tow)
            time_epochs = np.column_stack((t_week,t_tow))
            nGNSS_sat_current_epoch = np.zeros([nGNSSsystems,1])
            ## Initialize dummy variables
            GNSS_obs_dum = {}
            GNSS_LLI_dum = {}
            GNSS_SS_dum  = {}
            for k in np.arange(0, nGNSSsystems):
                max_sat_k = int(max_sat[k].item()) if hasattr(max_sat[k], "item") else int(max_sat[k])
                GNSS_obs_dum[k+1] = np.zeros([max_sat_k + 1, numOfObsCodes[k]])
                GNSS_LLI_dum[k+1] = np.zeros([max_sat_k + 1, numOfObsCodes[k]])
                GNSS_SS_dum[k+1]  = np.zeros([max_sat_k + 1, numOfObsCodes[k]])

            ## -- Iterate through satellites of epoch and store obs, LLI and SS
            for sat in np.arange(0,numSV):
                curr_sys = SVlist[sat][0]
                GNSSsystemIndex = [i for i in GNSSsystems if GNSSsystems[i]==curr_sys][0]
                nGNSS_sat_current_epoch[GNSSsystemIndex-1] +=  1 # Increment amount of satellites this epoch for this GNSS system
                SV = int(SVlist[sat][1:3]) # Get just PRN number
                nObsTypes_current_sat = numOfObsCodes[GNSSsystemIndex-1] # Number of obs types for current satellite
                GNSS_obs_dum[GNSSsystemIndex][SV][0:nObsTypes_current_sat] = Obs[sat,0:nObsTypes_current_sat] # Store observations, LLI, and SS of current satellite this epoch

                if readLLI:
                    GNSS_LLI_dum[GNSSsystemIndex][SV][0:nObsTypes_current_sat] = LLI[sat, 0:nObsTypes_current_sat]
                if readSS:
                    GNSS_SS_dum[GNSSsystemIndex][SV][0:nObsTypes_current_sat] = SS[sat, 0:nObsTypes_current_sat]

                # GNSS_SVs[curr_sys][current_epoch-1, int(nGNSS_sat_current_epoch[GNSSsystemIndex-1])] = SV  # Store PRN number of current sat to PRNs of this epoch
                GNSS_SVs[curr_sys][current_epoch-1, int(nGNSS_sat_current_epoch[GNSSsystemIndex-1].item()) if hasattr(nGNSS_sat_current_epoch[GNSSsystemIndex-1], "item") else int(nGNSS_sat_current_epoch[GNSSsystemIndex-1])] = SV


            for k in np.arange(0,nGNSSsystems):
                curr_sys = GNSSsystems[k+1]
                # GNSS_SVs[curr_sys][current_epoch-1, 0]  = nGNSS_sat_current_epoch[k] # Set number of satellites with obs for each GNSS system this epoch
                GNSS_SVs[curr_sys][current_epoch-1, 0] = float(nGNSS_sat_current_epoch[k].item()) if hasattr(nGNSS_sat_current_epoch[k], "item") else float(nGNSS_sat_current_epoch[k])

                if curr_sys == 'G':
                    GPS[current_epoch] = GNSS_obs_dum[k+1]
                elif curr_sys == 'R':
                    GLONASS[current_epoch] = GNSS_obs_dum[k+1]
                elif curr_sys == 'E':
                    Galileo[current_epoch] = GNSS_obs_dum[k+1]
                    Galileo_LLI[current_epoch]  =GNSS_LLI_dum[k+1]
                elif curr_sys == 'C':
                    BeiDou[current_epoch] = GNSS_obs_dum[k+1]
                    BeiDou_LLI[current_epoch]  =GNSS_LLI_dum[k+1]


                if readLLI and curr_sys == 'G':
                    GPS_LLI[current_epoch]  =GNSS_LLI_dum[k+1]
                elif readLLI and curr_sys == 'R':
                    GLONASS_LLI[current_epoch] = GNSS_LLI_dum[k+1]
                elif readLLI and curr_sys == 'E':
                    Galileo_LLI[current_epoch]  = GNSS_LLI_dum[k+1]
                elif readLLI and curr_sys == 'C':
                    BeiDou_LLI[current_epoch]  = GNSS_LLI_dum[k+1]

                if readSS and curr_sys =='G':
                    GPS_SS[current_epoch] = GNSS_SS_dum[k+1]
                elif readSS and curr_sys =='R':
                    GLONASS_SS[current_epoch] = GNSS_SS_dum[k+1]
                elif readSS and curr_sys =='E':
                    Galileo_SS[current_epoch] = GNSS_SS_dum[k+1]
                elif readSS and curr_sys =='C':
                    BeiDou_SS[current_epoch]  = GNSS_SS_dum[k+1]



    ## -- Storing observation in dictionary
    GNSS_obs['G'] = GPS
    GNSS_obs['R'] = GLONASS
    GNSS_obs['E'] = Galileo
    GNSS_obs['C'] = BeiDou
    ## -- Storing "loss of lock indicaors"  in dict
    GNSS_LLI['G'] = GPS_LLI
    GNSS_LLI['R'] = GLONASS_LLI
    GNSS_LLI['E'] = Galileo_LLI
    GNSS_LLI['C'] = BeiDou_LLI
    ## -- Storing  SS in dict
    GNSS_SS['G'] = GPS_SS
    GNSS_SS['R'] = GLONASS_SS
    GNSS_SS['E'] = Galileo_SS
    GNSS_SS['C'] = BeiDou_SS
    # Deleting systems with no observations
    del_sys = list(GNSS_obs.keys())
    for sys in del_sys:
        if not GNSS_obs[sys]:
            del GNSS_obs[sys]

    if current_epoch!= nepochs and success == 1:
        print('ERROR(readRinexObs304): The amount of epochs calculated in advance(nepochs = %d) does not equal number og epochs prossesed(current_epoch = %d).\nCheck that header information concerning TIME OF FIRST OBS and TIME OF LAST OBS is correct.\n' %(nepochs, current_epoch))

    messages = {}
    if success == 1:
        messages[0]= 'INFO(readRinexObs304): The following GNSS systems have been read into the data:'
        for k in np.arange(0,nGNSSsystems):
            messages[k+1]= 'INFO(readRinexObs304): The following %s observation types have been registered:' % (GNSS_names[GNSSsystems[k+1]])
            curr_sys = GNSSsystems[k+1]
            for obs in np.arange(0, len(obsCodes[k+1][curr_sys])):
                if obs == 0:
                    messages[k+1]= messages[k+1] + ' %s' % (obsCodes[k+1][curr_sys][obs])
                else:
                    messages[k+1]= messages[k+1] + ', %s' % (obsCodes[k+1][curr_sys][obs])
            if k == 0:
                messages[0]= messages[0] + ' %s' % GNSS_names[GNSSsystems[1]]
            else:
                messages[0]= messages[0] + ', %s' % GNSS_names[GNSSsystems[k+1]]

        for msg in np.arange(0,len(messages)):
            print(messages[msg])

    if readLLI:
        print('INFO(readRinexObs304): LLI have been read (if present in observation file)')
    else:
        print('INFO(readRinexObs304): LLI have not been read')


    if readSS:
        print('INFO(readRinexObs304): SS have been read (if present in observation file)')
    else:
        print('INFO(readRinexObs304): SS have not been read')


    ## --  Finding processing time
    et = time.process_time()  # get the end time
    e = et - t                # get execution time

    if e >= 3600:
        hours = np.floor(e/3600)
        minutes = np.floor((e-hours*3600)/60)
        seconds = e-hours*3600-minutes*60
        print('INFO(readRinexObs304): Total processing time: %d hours, %d minutes, %f seconds\n' % (hours, minutes, seconds))
    elif e>60:
        minutes = np.floor(e/60)
        seconds = e-minutes*60
        print('INFO(readRinexObs304): Total processing time: %d minutes, %f seconds\n' % (minutes, seconds))
    else:
        print('INFO(readRinexObs304): Total processing time: %f seconds\n\n' % (e))


    return GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
        obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
        rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success



def rinexFindNEpochs304(filename, tFirstObs, tLastObs, tInterval):
    """
    Function that computes number of epochs in Rinex 3.xx observation file.
    --------------------------------------------------------------------------------------------------------------------------
    INPUTS

    filename:         RINEX observation filename

    tFirstObs:        time stamp of the first observation record in the RINEX
                      observations file; column vector
                      [YYYY; MM; DD; hh; mm; ss.sssssss];

    tLastObs:         time stamp of the last observation record in the RINEX
                      observations file; column vector
                      [YYYY; MM; DD; hh; mm; ss.sssssss]. If this information
                      was not available in rinex observation header the
                      default value is Nan. In this case the variable is
                      determined in this function

    tInterval:        observations interval; seconds. If this information
                      was not available in rinex observation header the
                      default value is Nan. In this case the variable is
                      determined in this function.
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS:
    -------

    nepochs:          number of epochs in Rinex observation file with
                      observations

    tLastObs:         time stamp of the last observation record in the RINEX
                      observations file; column vector
                      [YYYY, MM, DD, hh, mm, ss.sssssss]. If this information
                      was not available in rinex observation header the
                      default value is Nan. In this case the variable is
                      determined in this function

    tInterval:        observations interval; seconds. If this information
                      was not available in rinex observation header the
                      default value is Nan.

    success:                  Boolean. 1 if the function seems to be successful,
                              0 otherwise
    --------------------------------------------------------------------------------------------------------------------------

    ADVICE: The function rinexFindNEpochs() calculates the amount of observation epochs in
    advance. This calculation will be incredibly more effective if TIME OF
    LAST OBS is included in the header of the observation file. It is
    strongly advized to manually add this information to the header if it is
    not included by default.
    --------------------------------------------------------------------------------------------------------------------------
    """
    # #  Testing input arguments
    success = 1
    nepochs = 0

    ## --Test if filename is valid format
    if type(filename) is not str:
        raise TypeError('INPUT ERROR(rinexFindNEpoch): The input argument filename'\
            'is of type %s. Must be of type string' % type(filename))


    ## --Open observation file
    fid = open(filename, 'rt')
    #  tLastObs is in header
    if ~np.all(np.isnan(tLastObs)):
        #  tInterval is not in header
        if np.isnan(tInterval):
            tInterval_found = 0
            first_epoch_found = 0
            #  calculate tInterval
            while not tInterval_found: #  calculate tInterval
                line = fid.readline().rstrip()
                #  start of new epoch
                if '>' in line:
                    if not first_epoch_found: #  first epoch
                        first_epoch_time = line[1::]
                        first_epoch_time = [float(el) for el in line[1::].split(" ") if el != ""]
                        first_epoch_time = first_epoch_time[:6]
                        first_epoch_found = 1
                    else: #  seconds epoch
                        second_epoch_time = line[1::]
                        second_epoch_time = [float(el) for el in line[1::].split(" ") if el != ""]
                        second_epoch_time = second_epoch_time[:6]
                        tInterval = second_epoch_time[5]-first_epoch_time[5]
                        tInterval_found = 1

        fid.close(); fid = open(filename, 'rt')
        tFirstObs = tFirstObs.astype(int)
        tLastObs = tLastObs.astype(int)
        rinex_lines = fid.readlines()
        epoch_line = [line for line in rinex_lines if line.startswith('>')] # list with all the line thats defines a epoch
        nepochs = len(epoch_line)
    #  if tLastObs is not in header. Function counts number of epochs manually
    else:
    ## New code for finding last
        print('INFO(rinexFindEpochs304): The header of the rinex observation file does not contain TIME OF LAST OBS.\n' \
            'This will be calculated, but consider editing rinex header to include TIME OF LAST HEADER')

        fid.close(); fid = open(filename, 'rt')
        rinex_lines = fid.readlines()
        epoch_lines = [line for line in rinex_lines if '>' in line] # list with all the line thats defines a epoch
        nepochs = len(epoch_lines)
        ## Computing the tInterval if not present in the header
        if np.isnan(tInterval):
            first_epoch_line  = epoch_lines[0][1::]
            second_epoch_line = epoch_lines[1][1::]
            first_epoch_time  = [float(el) for el in first_epoch_line[1::].split(" ") if el != ""]
            second_epoch_time = [float(el) for el in second_epoch_line[1::].split(" ") if el != ""]
            tInterval = second_epoch_time[5]-first_epoch_time[5]

        line = epoch_lines[-1]
        line = line[1:60]     #  deletes 'TIME OF LAST OBS'
        line_ = [el for el in line.split(" ") if el != ""]
        for k in np.arange(0,6):
            tok = line_.pop(0)
            if k ==0:
                yyyy = int(tok)
            elif k ==1:
                mm = int(tok)
            elif k ==2:
                dd = int(tok)
            elif k ==3:
                hh = int(tok)
            elif k ==4:
                mnt = int(tok)
            elif k ==5:
                ss = float(tok)

        tLastObs = np.array([[yyyy],[mm],[dd],[hh],[mnt],[ss]]).astype(int)
        print('INFO(rinexFindNEpochs304): TIME OF LAST OBS has been found and amount of epochs have been computed')
        fid.close()
    return int(nepochs), tLastObs, tInterval, success



def rinexReadObsFileHeader304(filename, includeAllGNSSsystems, includeAllObsCodes,desiredGNSSsystems,desiredObsCodes, desiredObsBands):
    """
    Extracts relevant data from the header of a RINEX 3.xx GNSS observations
    file. Excludes undesired GNSS systems, obsevation codes and/or frequency
    bands.

    --------------------------------------------------------------------------------------------------------------------------
    INPUTS:
    ------

    filename:                     RINEX observation filename and path

    includeAllGNSSsystems:        Boolean, 0 or 1.
                                      1 = include alle GNSS systems
                                          (GPS, GLONASS, Galieo, BeiDou)
                                      0 = include only GNSSsystems
                                          specified in desiredGNSSsystems

    includeAllobsCodes:           Boolean, 0 or 1.
                                      1 = include all valid obsCodes
                                      0 = include only obsCodes
                                          specified in desiredobsCodes

    desiredGNSSsystems:           string array containing desired GNSSsystems
                                  to be included, ex. ["G", "E", "C"]

    desiredobsCodes:              string array containing desired obsCodes to
                                  be included, ex. ["C", "L", "S", "D"]

    desiredObsBands:              array of desired obs Bands to be included,
                                  ex [1, 5]

    NOTE: If both includeAllGNSSsystems and includeAllobsCodes Boolean are 1
          then the last three input arguments are optional to include and may
          be left blank without en error.
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS:
    -------

    success:                      1 if the reading of the RINEX observations
                                  file seems to be successful, 0 otherwise

    rinexVersion:                 string. rinex observation file version

    gnssType:                     GNSS system of the satellites observed; can
                                  be 'G', 'R', 'E', 'C' or 'M' that stand for
                                  GPS, GLONASS, GALILEO, BeiDou or Mixed; char

    markerName:                   name of the antenna marker; '' if not
                                  specified

    recType:                      Receiver type, char vector

    antDelta:                     column vector ot the three components of
                                  the distance from the marker to the antenna,
                                  in the following order - up, east and north;
                                  null vector by default

    GNSSsystems:                  cell array containing codes of GNSS systems
                                  included in RINEX observationfile. Elements
                                  are strings. ex. "G" or "E"

    numOfObsCodes:                column vector containing number of observation
                                  types for each GNSS system. Order is the same
                                  as GNSSsystems

    obsCodes:                     Cell that defines the observation
                                  codes available for all GNSS system. Each
                                  cell element is another cell containing the
                                  codes for that GNSS system. Each element in
                                  this cell is a string with three-characters.
                                  The first character (a capital letter) is
                                  an observation code ex. "L" or "C". The
                                  second character (a digit) is a frequency
                                  code. The third character(a Capital letter)
                                  is the attribute, ex. "P" or "X"

    obsCodeIndex:                 cell with one cell element for each GNSS
                                  system. Order is the same as GNSSsystems.
                                  Each cell element contains an array of
                                  indices. These indices indicate the
                                  observation types that should be read
                                  for each GNSS system. ex. If one index for
                                  GPS is 1 then the first observation type
                                  for GPS should  be read.

    tFirstObs:                    time stamp of the first observation record
                                  in the RINEX observations file; column vector
                                  [YYYY; MM; DD; hh; mm; ss.sssssss];
                                  THIS IS CRITICAL DATA

    tLastObs:                     time stamp of the last observation record
                                  in the RINEX observations file; column vector
                                  [YYYY; MM; DD; hh; mm;ss.sssssss].
                                  NaN by default.
                                  THIS IS RINEX 3.04 OPTIONAL DATA

    tInterval:                    observations interval; seconds.

    timeSystem:                   three-character code string of the time
                                  system used for expressing tfirstObs;
                                  can be GPS, GLO or GAL;

    numHeaderLines:               number of lines in header

    rinexProgr:                   name of the software used to produce de
                                  RINEX GPS obs file; '' if not specified

    rinexDate:                    date/time of the RINEX file creation; ''
                                  if not specified; char

    leapSec:                      number of leap seconds since 6-Jan-1980.
                                  UTC=GPST-leapSec. NaN by default.
                                  THIS IS RINEX 3.04 OPTIONAL DATA

    approxPosition:               array containing approximate position from
                                  rinex observation file header. [X, Y, Z]

    GLO_Slot2ChannelMap:          map container that maps GLONASS slot
                                  numbers to their respective channel number.
                                  GLO_Slot2ChannelMap(slotnumber)

    eof:                          end-of-file flag; 1 if end-of-file was reached,
                                  0 otherwise

    fid:                          Matlab file identifier of a Rinex
                                  observations text file
    --------------------------------------------------------------------------------------------------------------------------

    According to RINEX 3.04 these codes are:

       Observation codes:
       ------------------
       C: Pseudorange
          GPS: C/A, L2C
          Glonass: C/A
          Galileo: All
       L: Carrier phase
       D: Doppler frequency
       S: Raw signal strengths or SNR values as given by the receiver for the
          respective phase observations
       I: Ionosphere phase delay
       X: Receiver channel numbers

       Frequency code:
       ---------------
       GPS Glonass Galileo SBAS
       1: L1 G1 E1 B1    (GPS,QZSS,SBAS,BDS)
       2: L2 G2 B1-2     (GLONASS)
       4: G1a            (Galileo)
       5: L5 E5a B2/B2a  (GPS, QZSS, SBAS, IRNSS)
       6: L6 E6 B3 G2a   (Galileo, QZSS, BDS, GLONASS)
       7: E5b B2/B2b     (Galileo)
       8: E5a+b E5a+b    (Galileo, BDS)
       9: S              (IRNSS)
       0: for type X     (all)

       Attribute:
       ----------
       A = A channel     (Galileo,IRNSS,GLONASS)
       B = B channel     (Galileo,IRNSS,GLONASS)
       C = C channel     (Galiloe, IRNSS)
           C code-based  (SBAS, GPS, GLONASS, QZSS)
       D = Semi-codelss  (GPS)

       I = I channel     (GPS, Galileo, QZSS, BDS)
       L = L channel     (L2C GPS, QZSS)
           P channel     (GPS. QZSS)
       M = M code-based  (GPS)
       N = Codeless      (GPS)
       P = P code-based  (GPS, GLONASS)
           Pilot channel (BDS)

       Q = Q channel     (GPS, Galileo, QZSS, BDS)
       S = D channel     (GPS, Galileo, QZSS, BDS)
           M channel     (L2C GPS, QZSS)

       W = Based on Z-tracking (GPS)
       X = B+C channels  (Galileo, IRNSS)
           I+Q channels  (GPS, IRNSS)
           M+L channels  (GPS, QZSS)
           D+P channels  (GPS, QZSS, BDS)

       Y = Y code based  (GPS)
       Z = A+B+C channels(Galileo)
           D+P channels  (BDS)
    -------------------------------------------------------------------------------------------------------------------------
    """
    eof         = 0
    success     = 1
    warnings    = 0
    antDelta    = []
    timeSystem  = ''
    tFirstObs   = []
    tLastObs    = np.nan
    tInterval   = np.nan
    rinexProgr  = np.nan
    rinexDate   = np.nan
    obsCodes    = {}
    GNSSsystems = {}
    gnssType    = ""
    markerName  = ""
    numHeaderLines  = 0
    clockOffsetsON  = 0
    numGNSSsystems  = 0
    leapSec         = np.nan
    numOfObsCodes   = []
    rinexHeader     = {}
    approxPosition  = [0, 0, 0]
    obsCodeIndex = {}
    rinexVersion = np.nan
    recType = np.nan
    GLO_Slot2ChannelMap = np.nan

    ## -------Testing input arguments
    # Test if filename is valid format
    if type(filename) != str:
        print('INPUT ERROR(rinexReadsObsHeader304): The input argument filename is of type %s.\n Must be of type string or char' %(type(filename)))
        success = 0
        fid     = 0
        return success
    ## -- Open rinex observation file
    fid = open(filename,'r')
    if os.stat(filename).st_size == 0:
        raise ValueError('ERROR: This file seems to be empty')

    while 1: # Gobbling the header
        numHeaderLines = numHeaderLines + 1
        line = fid.readline().rstrip()
        if 'END OF HEADER' in line:
            break
        if numHeaderLines == 1: # if first line of header
            rinexVersion = line[0:9]
            # store rinex type, ex. "N" or "O"
            rinexType = line[20]
            # if rinex file is not an observation file
            if rinexType != 'O':  # Rinex file is oservation file
                print('ERROR(rinexReadObsFileHeader304): the file is not a RINEX observations data file!')
                success = 0
                fid.close()
                return

            ## -- Check gnss type  ## Changend indent here 09.12.2022 (was apart of the if test above earlier, and thats wrong)
            gnssType = line[40] # reads the GNSS system type
            if gnssType not in [' ', 'G', 'R', 'C', 'E', 'M' ]:
                if gnssType in ['J', 'I', 'S']:
                    print('ERROR(rinexReadObsFileHeader304): This software is meant for reading GNSS data only.\
                           %s is an invalid satellite system type.' %(gnssType))
                else:
                    print('ERROR(rinexReadObsFileHeader304): %s is an unrecognized satellite system type.' %(gnssType))

                success = 0
                fid.close()
            ## -- If no system type, set G
            if gnssType == ' ':
                gnssType = 'G'

        if 'PGM / RUN BY / DATE' in line:
            rinexProgr = line[0:20] # rinex program
            rinexDate = line[40:60] # rinex date

        if 'MARKER NAME' in line:
            markerName = line.strip() # markername

        ## if no marker name, "MARKER" is read, so set to blank
        if 'Marker' in markerName:
            markerName = ''

        if 'ANTENNA: DELTA H/E/N' in line:
            for k in np.arange(0,3):
                line_ = [el for el in line.split(" ") if el != ""]
                antDelta = [line_[0],line_[1],line_[2]]

        ## Section describing what GNSS systems are present, and their obs types
        if 'SYS / # / OBS TYPES' in line:
            line = line[0:60]     # deletes 'SYS / # / OBS TYPES'
            line_ = [el for el in line.split(" ") if el != ""]
            Sys = line_.pop(0) # assingning system to variable and removing it from the list
            if Sys not in ["G","R","E","C"]: # added this line 29.01.2023 to fix bug where Only one system and several lines with Obscodes in rinex file
                continue
            nObs = int(line_.pop(0))
            ## array for storing indeces of undesired ObsCodes for this GNSS system
            undesiredobsCodeIndex = []
            desiredObsCodeIndex = []
            ## is Sys amoung desired GNSS systems
            if (includeAllGNSSsystems and Sys in ["G", "R", "E", "C"] or Sys in desiredGNSSsystems):
                numGNSSsystems  = numGNSSsystems + 1 # increment number of GNSS systems
                GNSSsystems[numGNSSsystems] = str(Sys) # Store current GNSS system
                GNSSSystemObsCodes = {}  # Reset cell of obsCodes for this GNSS system
                obsCode_list = []
                for k in np.arange(0,nObs):
                    obsCode = line_.pop(0)
                    # Checking if obsCode is valid
                    if len(obsCode) != 3 or obsCode[0] not in ['C', 'L', 'D','S', 'I', 'X'] or  \
                              obsCode[1] not in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'] or \
                               obsCode[2] not in ['A', 'B', 'C', 'D', 'I', 'L', 'M', 'N', 'P', 'Q', 'S', 'W', 'X', 'Y', 'Z']:
                        print('ERROR (rinexReadsObsHeader304):  obsCode %s is a not a standard RINEX 3.04 observation type!' %(obsCode))

                    ## is obsCode amoung desired obscodes and frequency bands
                    if includeAllObsCodes or obsCode[0] in desiredObsCodes and int(obsCode[1]) in desiredObsBands:
                         ## store obsCode if amoung desire obsCodes
                        obsCode_list.append(obsCode)
                        GNSSSystemObsCodes[Sys] =  obsCode_list
                        desiredObsCodeIndex.append(k)
                    else:
                        # store index of discareded obsCode
                        undesiredobsCodeIndex.append(k)

                    # Every 13 obsCodes is at end of line. In this case read next line and continue
                    if np.mod(k+1, 13) == 0 and nObs != 13:
                        numHeaderLines = numHeaderLines + 1
                        line = fid.readline().rstrip()
                        line = line[0:60]     # deletes 'SYS / # / OBS TYPES'
                        line_ = [el for el in line.split(" ") if el != ""]

                numOfObsCodes.append(len(GNSSSystemObsCodes[Sys]))
                obsCodes[numGNSSsystems] = GNSSSystemObsCodes
                obsCodeIndex[numGNSSsystems] = desiredObsCodeIndex # Store indices of desired obsCodes


        if 'TIME OF FIRST OBS' in line:
            line = line[0:60]     #  deletes 'TIME OF FIRST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            for k in np.arange(0,6):
                tok = line_.pop(0)
                if k ==0:
                    yyyy = int(tok)
                elif k ==1:
                    mm = int(tok)
                elif k ==2:
                    dd = int(tok)
                elif k ==3:
                    hh = int(tok)
                elif k ==4:
                    mnt = int(tok)
                elif k ==5:
                    ss = float(tok)


            tFirstObs = np.array([[yyyy],[mm],[dd],[hh],[mnt],[ss]])

            # Get Time system
            if len(line_) != 0:
                aux = line_.pop(0)
                if aux == 'GPS':
                    timeSystem = 'GPS'
                elif aux == 'GLO':
                    timeSystem = 'GLO'
                elif aux == 'GAL':
                    timeSystem = 'GAL'
                elif aux == 'BDT':
                    timeSystem = 'BDT'

            else:
                try:
                    if gnssType == 'G':
                        timeSystem = 'GPST'
                    elif gnssType == 'R':
                        timeSystem = 'GLOT'
                    elif gnssType == 'E':
                        timeSystem = 'GALT'
                    elif gnssType == 'C':
                        timeSystem = 'BDT'
                except:
                    timeSystem = "GPS"



        if 'TIME OF LAST OBS' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            for k in np.arange(0,6):
                tok = line_.pop(0)
                if k ==0:
                    yyyy = int(tok)
                elif k ==1:
                    mm = int(tok)
                elif k ==2:
                    dd = int(tok)
                elif k ==3:
                    hh = int(tok)
                elif k ==4:
                    mnt = int(tok)
                elif k ==5:
                    ss = float(tok)

            tLastObs = np.array([[yyyy],[mm],[dd],[hh],[mnt],[ss]])

        if 'INTERVAL' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            tInterval = float(line_.pop(0))



          ## -- This is an optional record!
          # if 'RCV CLOCK OFFS APPL' in line:
          #     if (strtok(line)=='0'):
          #         clockOffsetsON = 0;
          #     elif (strtok(line)=='1'):
          #         clockOffsetsON = 1;
          #     else:
          #         success = 0;
          #         print('ERROR (rinexReadsObsHeader304): unrecognized receiver clock offsets flag!')
          #         fid.close()


           ## This is an optional record
        if 'LEAP SECONDS' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            leapSec = int(line_.pop(0))



           ## -- store approximate receiver position
        if 'APPROX POSITION XYZ' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            approxPosition = np.array([[float(line_[0])],[float(line_[1])],[float(line_[2])]])


         ## GLOANSS SLOTS
        if 'GLONASS SLOT / FRQ #' in line:
            line = line[0:60]     #  deletes 'GLONASS SLOT / FRQ #'
            line_ = [el for el in line.split(" ") if el != ""]
            nGLOSat = int(line_.pop(0))
            slotNumbers = np.array([])
            channels = np.array([])
            for k in np.arange(0,nGLOSat):
                slotNumber = line_.pop(0)[1::]
                channel = int(line_.pop(0))
                slotNumbers = np.append(slotNumbers,slotNumber)
                channels = np.append(channels,channel)

                # if np.mod(k+1, 8) == 0:
                if np.mod(k+1, 8) == 0 and k+1 != 24: # add test if  k ==24 to prevnet skip extra line
                    # line = fgetl(fid); # end of line is reached so read next line
                    line = fid.readline().rstrip()
                    numHeaderLines = numHeaderLines + 1
                    line = line[0:60]     #  deletes 'TIME OF LAST OBS'
                    line_ = [el for el in line.split(" ") if el != ""]
                elif np.mod(k+1, 8) == 0 and k+1 == 24:
                    break

            GLO_Slot2ChannelMap = dict(zip(slotNumbers.astype(int),channels.astype(int)))

        if 'REC # / TYPE / VERS' in line:
            recType = line[20:40]



     # End of Gobbling Header Loop
    for k in np.arange(0,numGNSSsystems):
        # Give info if any of GNSS systems had zero of desired obscodes.
        if numOfObsCodes[k] == 0 or sum(tFirstObs) == 0:
            if GNSSsystems[k] == 'G':
                print('INFO: (rinexReadsObsHeader304)\nNone of the GPS satellites had any of the desired obsCodes\n\n')
            elif GNSSsystems[k] == 'R':
                print('INFO: (rinexReadsObsHeader304)\nNone of the GLONASS satellites had any of the desired obsCodes\n\n')
            elif GNSSsystems[k] == 'E':
                print('INFO: (rinexReadsObsHeader304)\nNone of the Galileo satellites had any of the desired obsCodes\n\n')
            elif GNSSsystems[k] == 'C':
                print('INFO: (rinexReadsObsHeader304)\nNone of the BeiDou satellites had any of the desired obsCodes\n\n')

    ## store rinex header info
    rinexHeader['rinexVersion'] =rinexVersion
    rinexHeader['rinexType'] = rinexType
    rinexHeader['gnssType'] =gnssType
    rinexHeader['rinexProgr'] =rinexProgr
    rinexHeader['rinexDate'] =rinexDate


    print('INFO(rinexReadObsFileHeader304): Rinex header has been read')

    return success, rinexVersion, gnssType, markerName, recType, antDelta,GNSSsystems,numOfObsCodes, \
    obsCodes, obsCodeIndex,tFirstObs, tLastObs, tInterval,timeSystem, numHeaderLines, clockOffsetsON, \
    rinexProgr, rinexDate,leapSec, approxPosition, GLO_Slot2ChannelMap, eof, fid



def rinexReadObsBlock304(fid, numSV, nObsCodes, GNSSsystems, obsCodeIndex, readSS, readLLI):
    """
    Reads all the observations from a RINEX observation block.

    Positioned at the beginning of the line immediately after the header of the
    observations block, reads all the observations in this block of a RINEX
    observations file. This function is meant to be used after using function
    rinexReadObsFileHeader304

    Based in the work of Antonio Pestana, rinexReadObsBlock211, March 2015
    --------------------------------------------------------------------------------------------------------------------------
    INPUTS:
    -------

    fid:                  Matlab file identifier of a Rinex observations text file

    numSV:                total number of satellites with observations in
                          current observation block, integer

    numOfObsCodes:        column vector containing number of observation
                          types for each GNSS system. Order is the same as
                          GNSSsystems

    GNSSsystems:          cell array containing codes of GNSS systems included
                          in RINEX observationfile. Elements are strings.
                          ex. "G" or "E"

    obsCodeIndex:         cell with one cell element for each GNSS system.
                          Order is the same as GNSSsystems. Each cell element
                          contains an array of indices. These indices
                          indicate the observation types that should be
                          read for each GNSS system. ex. If one index for
                          GPS is 1 then the first observation type for GPS
                          should be read.

    readSS:                   Boolean, 0 or 1.
                              1 = read "Signal Strength" Indicators
                              0 = do not read "Signal Strength" Indicators

    readLLI:                  Boolean, 0 or 1.
                              1 = read "Loss-Of-Lock Indicators"
                              0 = do not read "Loss-Of-Lock Indicators"
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS:
    --------

    success:               Boolean. 1 if the function seems to be successful,
                          0 otherwise

    Obs:                  matrix [numSV x max_nObs] that stores all
                          observations of this observation block. max_nObs
                          is the highest number of observation codes that
                          any of the GNSS systems have. Which observation
                          types that are associated with what collumn will
                          vary between GNSS systems. SVlist will give
                          overview of what GNSS system each row is connected
                          to

    SVlist:               column cell [numSV x 1] that conatins the
                          identification code of each line of observation
                          block. ex. "G21". numSV is total number of
                          satellites minus amount of satellites removed.

    numSV:                numSV, unlike the input of same name, is the total
                          number of satellites minus amount of satellites
                          removed.

    LLI:                  matrix [numSV x max_nObs] that stores all
                          "loss-of-lock" indicators of this observation block.
                          max_nObs is the highest number of observation codes
                          that any of the GNSS systems have. Which observation
                          types that are associated with what collumn will
                          vary between GNSS systems. SVlist will give
                          overview of what GNSS system each row is connected
                          to

    SS:                   matrix [numSV x max_nObs] that stores all
                          "signal strength" indicators of this observation block.
                          max_nObs is the highest number of observation codes
                          that any of the GNSS systems have. Which observation
                          types that are associated with what collumn will
                          vary between GNSS systems. SVlist will give
                          overview of what GNSS system each row is connected
                          to

    eof:                  end-of-file flag; 1 if end-of-file was reached,
                          0 otherwise
    --------------------------------------------------------------------------------------------------------------------------
    """
    ## Initialize variables in case of input error
    success                     = np.nan
    eof                         = np.nan
    max_n_obs_Types             = np.nan
    Obs                         = np.nan
    LLI                         = np.nan
    SS                          = np.nan
    SVlist                      = np.nan
    removed_sat                 = np.nan
    desiredGNSSsystems          = np.nan

    ## -- Testing input arguments
    if type(numSV) != int:
        print(f'INPUT ERROR(rinexReadObsBlock304): The input argument numSV is of type {type(numSV)}.\n Must be of type double')
        success = 0
        return success

    nObsCodes = [int(x) for x in nObsCodes]
    ## Test type of numOfObsCodes
    if type(nObsCodes[0]) != int:
        print('INPUT ERROR(rinexReadObsBlock304): The input argument numOfObsTypes is of type %s.\n Must be of type double' % (type(nObsCodes)))
        success = 0
        return success


    ## Test size of numOfObsCodes
    if len(nObsCodes) != len(GNSSsystems):
        print('INPUT ERROR(rinexReadObsBlock304): The input argument numOfObsTypes must have same length as GNSSsystems')
        success = 0
        return success

    success = 1
    eof     = 0

    # Highest number of obs codes of any GNSS system
    max_n_obs_Types = max(nObsCodes)

    # Initialize variables
    Obs = np.empty([numSV, max_n_obs_Types])
    SVlist = [np.nan]*numSV
    if readLLI:
        LLI = np.empty([numSV, max_n_obs_Types])

    if readSS:
        SS  = np.empty([numSV, max_n_obs_Types])

    # number of satellites excluded so far
    removed_sat = 0
    desiredGNSSsystems = list(GNSSsystems.values())
    # Gobble up observation block
    for sat in np.arange(0,numSV):
        line = fid.readline().rstrip()
        if not line:
            return

        SV = line[0:3].strip() # Satellite code, ex. 'G11' or 'E03'
        if SV[0] not in desiredGNSSsystems:
            removed_sat +=1
        else:
            ## Index of current GNSS system
            GNSSsystemIndex = [i for i in GNSSsystems if GNSSsystems[i]==SV[0]][0]
            SVlist[sat - removed_sat] = SV # Store SV of current row
            n_obs_current_system = nObsCodes[GNSSsystemIndex-1]
            for obs_num in np.arange(0, n_obs_current_system):
                obsIndex = obsCodeIndex[GNSSsystemIndex][obs_num]
                charPos = 4+(obsIndex)*16
                ## check that the current observation of the current GNSS system
                ## is not on the list of obs types to be excluded
                ## stringlength of next obs.
                obsLen = min(14, len(line) - charPos)
                # read next obs
                newObs = line[charPos:charPos+obsLen].strip()
                # If observation missing, set to 0
                if newObs != '':
                    newObs = float(newObs)
                else:
                    newObs = 0
                # Store new obs
                Obs[sat - removed_sat, obs_num] = newObs
                # read LLI of current obs (if present)
                if readLLI:
                    if charPos+13<len(line):
                        newLLI = line[charPos+13]
                    else:
                        newLLI = ' '

                    if newLLI.isspace():
                        newLLI = -999
                    else:
                        newLLI = int(newLLI)
                    LLI[sat - removed_sat, obs_num] = newLLI


                if readSS:
                    # read SS of current obs (if present)
                    if charPos+14<len(line):
                        newSS = line[charPos+14]
                    else:
                        newSS = ' '

                    # if no SS set to -999
                    if newSS.isspace():
                        newSS = -999
                    else:
                        newSS = int(newSS)
                    SS[sat - removed_sat, obs_num]  = newSS



    ## -- Update number og satellites after satellites have been excluded
    numSV = numSV - removed_sat
    ## --Remove empty arrays
    SVlist = list(filter(None,SVlist))
    idx_keep = len(Obs) -1 -removed_sat + 1 # removing sats
    Obs = Obs[:idx_keep,:]
    return success, Obs,SVlist, numSV, LLI, SS, eof






def rinexReadObsBlockHead304(fid):
    """
    Reads the metadata in the head of a RINEX 3.xx observations block, NOT
    the header of the file.

    ATTENTION: Ignores all data in blocks with event flags with numbers
    greater than 1.

    Positioned in a RINEX 3.04 GNSS observations text file at the beginning
    of an observation block. In rinex 3.xx the line starts with '> '

    --------------------------------------------------------------------------------------------------------------------------
     INPUTS

     fid:              Python identifier of an open RINEX 3.04 GNSS
                       observations text file positioned at the beginning
                       of an observation block.
    --------------------------------------------------------------------------------------------------------------------------
     OUTPUTS

     success:          1 if function performs successfully, 0 otherwise

     epochflag:        Rinex observations epoch flag, as follows:
                           0: OK
                           1: power failure between previous and current epoch
                       From now on the "event flags":
                           2: start moving antenna
                           3: new site occupation
                           4: header information follows
                           5: external event (epoch is significant)

     clockOffset:          value of the receiver clock offset. If not present
                           in the metadata of the observations block
                           (it's optional RINEX 3.04 data)it is assumed to be
                           zero. If not zero implies that epoch, code, and
                           phase data have been corrected by applying
                           realtime-derived receiver clock offset

     date:                 time stamp of the observations block. Six-elements column-vector
                           as follows:
                               year: four-digits year (eg: 1959)
                               month: integers 1..12
                               day: integers 1..31
                               hour: integers 0..24
                               minute: integers 0..60
                               second: reals 0..60

     numSV:                number of satellites with observations in with
                           observations. This will include all satellite
                           systems.
    --------------------------------------------------------------------------------------------------------------------------

    """
    # Initialize variables
    success = 1
    eof     = 0
    date    = [0,0,0,0,0,0]
    numSV   = 0
    epochflag = 0
    clockOffset = 0
    noFlag = 1

    line = fid.readline().rstrip()
    if not line:
        eof = 1
        return success, epochflag, clockOffset, date, numSV, eof
    epochflag   = line[31]
    # skip to next block if event flag is more than 1
    while int(epochflag) > 1:
        noFlag = 0
        linejump = int(line[32:35])
        msg = f'WARNING(rinexReadsObsBlockHead304): Observations event flag encountered. Flag = {str(epochflag)}, hence {str(linejump)} lines were ignored. '
        for count in np.arange(0,linejump+1):
            line = fid.readline().rstrip()
        epochflag = int(line[31]) # changed from 30 to 31 29.08.2023

    numSV = int(line[32:35])
    clockOffset = 0
    if len(line) == 56:
        clockOffset = float(line[41:56])

    # Reads the time stamp of the observations block (6 numerical values)
    date = line[1::]
    date = [float(el) for el in line[1::].split(" ") if el != ""]
    date = date[:6]

    if noFlag == 0:
        msg2 = msg + 'Epoch date: %.4d %.2d %.2d %.2d:%.2d:%6.4f' % (date[0],date[1],date[2],date[3],date[4],date[5])
        print(msg2)


    return success, epochflag, clockOffset, date, numSV, eof


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------- RINEX 2.11 ---------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------

def readRinexObs211(filename, readSS=None, readLLI=None, includeAllGNSSsystems=None,includeAllObsCodes=None, \
                    desiredGNSSsystems=None, desiredObsCodes=None, desiredObsBands=None):
    """
    Program/function to read GNSS observations in RINEX V.2 observation files
    The main core of the program is 4 functions:
                                  rinexReadObsFileHeader211
                                  rinexReadObsBlockHead211
                                  rinexReadObsBlock211
                                  rinexFindNEpochs211


    To export every parameter use this code:

    GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
         obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
         rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success = \
         readRinexObs211(rinObsFilename)

    Tip: Unwanted outputs can be ignored using _

    --------------------------------------------------------------------------------------------------------------------------
    INPUTS

    filename:                 path and name of RINEX V.2 observation file,
                              string

    readSS:                   Boolean, 0 or 1.
                              1 = read "Signal Strength" Indicators
                              0 = do not read "Signal Strength" Indicators

    readLLI:                  Boolean, 0 or 1.
                              1 = read "Loss-Of-Lock Indicators"
                              0 = do not read "Loss-Of-Lock Indicators"

    includeAllGNSSsystems:    Boolean, 0 or 1.
                              1 = include alle GNSS systems(GPS, GLONASS, Galieo, BeiDou)
                              0 = include only GNSSsystems specified in desiredGNSSsystems

    includeAllObsTypes:       Boolean, 0 or 1.
                              1 = include all valid ObsTypes
                              0 = include only ObsTypes specified in desiredObsTypes

    desiredGNSSsystems:       array og strings containing  codes of desired
                              GNSSsystems to be included,
                              ex. ["G", "E"]
                              OBS: Must be string array, NOT char vector

    desiredObsTypes:          array of strings containing desired ObsTypes to be
                              included, ex. ["C", "L", "S", "D"]
                              OBS: Must be string array, NOT char vector

    desiredObsBands:          array of desired obs Bands to be included,
                              ex [1, 5]
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS

    GNSS_obs:                 cell containing a matrix for each GNSS system.
                              Each matrix is a 3D matrix containing all
                              observation of current GNSS system for all epochs.
                              Order of obsType index is same order as in
                              obsCodes cell

                              GNSS_obs{GNSSsystemIndex}(PRN, obsType, epoch)
                                              GNSSsystemIndex: double,
                                              1,2,...,numGNSSsystems
                                              PRN: double
                                              ObsType: double: 1,2,...,numObsTypes
                                              epoch: double

    GNSS_LLI:                 cell containing a matrix for each GNSS system
                              Each matrix stores loss of lock indicators for
                              each epoch for that GNSS system

    GNSS_SS:                  cell containing a matrix for each GNSS system
                              Each matrix stores signal strength indicators
                              for each epoch for that GNSS system

    GNSS_SVs:                 cell containing a matrix for each GNSS system.
                              Each matrix contains number of satellites with
                              obsevations for each epoch, and PRN for those
                              satellites

                              GNSS_SVs{GNSSsystemIndex}(epoch, j)
                                              j=1: number of observed satellites
                                              j>1: PRN of observed satellites

    time_epochs:              matrix conatining gps-week and "time of week"
                              for each epoch
                              time_epochs(epoch,i),   i=1: week
                                                      i=2: time-of-week in seconds (tow)

    nepochs:                  number of epochs with observations in rinex observation file.

    GNSSsystems:              cell array containing codes of GNSS systems included
                              in RINEX observationfile. Elements are strings.
                              ex. "G" or "E"

    obsCodes:                 Cell that defines the observation
                              codes available for all GNSS system. Each cell
                              element is another cell containing the codes for
                              that GNSS system. Each element in this cell is a
                              string with three-characters. The first
                              character (a capital letter) is an observation code
                              ex. "L" or "C". The second character (a digit)
                              is a frequency code. The third character(a Capital letter)
                              is the attribute, ex. "P" or "X"

    approxPosition:           array containing approximate position from rinex
                              observation file header. [X, Y, Z]

    max_sat:                  array conataining max PRN number for each GNSS
                              system. Follows same order as GNSSsystems

    tInterval:                observations interval; seconds.

    markerName:               name of the antenna marker; '' if not specified

    rinexVersion:             string. rinex observation file version

    recType:                  receiver type, char vector

    timeSystem:               three-character code string of the time system
                              used for expressing tfirstObs;
                              can be GPS, GLO or GAL;

    leapSec:                  number of leap seconds since 6-Jan-1980.
                              UTC=GPST-leapSec. NaN by default.
                              THIS IS RINEX V.2 OPTIONAL DATA

                              rinexHeader: cell column-vector containing the
                              following data:
                              rinexVersion:   RINEX version number; string.
                                              '' if not specified
                              rinexType:      RINEX file type; char

    gnssType:                 GNSS system of the satellites observed; can be
                              'G', 'R', 'E', 'C' or 'M' that stand for
                              GPS, GLONASS, GALILEO, BeiDou or Mixed ; char

    rinexProgr:               name of the software used to produce de RINEX
                              GPS obs file; '' if not specified


    rinexDate:                date/time of the RINEX file creation; '' if not
                              specified; char


    antDelta:                 column vector ot the three components of the
                              distance from the marker to the antenna,
                              in the following order - up, east and north;
                              null vector by default

    tFirstObs:                time stamp of the first observation record in the RINEX
                              observations file; column vector
                              [YYYY; MM; DD; hh; mm; ss.sssssss];
                              THIS IS CRITICAL DATA

    tLastObs:                 time stamp of the last observation record in the RINEX
                              observations file; column vector
                              [YYYY, MM, DD, hh, mm,ss.sssssss]. NaN by default.
                              THIS IS RINEX V.2 OPTIONAL DATA

    clockOffsetsON:           receiver clock offsets flag. O if no realtime-derived
                              receiver clock offset was applied to epoch,
                              code and phase data (in other words, if the
                              file only has raw data), 1 otherwise. 0 by default.
                              THIS IS RINEX V.2 OPTIONAL DATA

    GLO_Slot2ChannelMap:      map container that maps GLONASS slot numbers to
                              their respective channel number.
                              GLO_Slot2ChannelMap(slotnumber)

    success:                  Boolean. 1 if the reading of the RINEX
                              observations file seems to be successful,
                              0 otherwise
    --------------------------------------------------------------------------------------------------------------------------

    ADVICE: The function rinexFindNEpochs() calculates the amount of observation epochs in
    advance. It is recommended to manually add this information to the header if it is
    not included by default.
    --------------------------------------------------------------------------------------------------------------------------

    According to RINEX V.2 the observation type codes are:

    Observation code
            C: Pseudorange
    GPS: C/A, L2C
    Glonass: C/A
    Galileo: All
    L: Carrier phase
    D: Doppler frequency
    S: Raw signal strengths or SNR values as given by the receiver for the
    respective phase observations
    I: Ionosphere phase delay
    X: Receiver channel numbers

    Frequency code
    GPS Glonass Galileo SBAS
    1: L1 G1 E1 B1    (GPS,QZSS,SBAS,BDS)
    2: L2 G2 B1-2     (GLONASS)
    4: G1a            (Galileo)
    5: L5 E5a B2/B2a  (GPS, QZSS, SBAS, IRNSS)
    6: L6 E6 B3 G2a   (Galileo, QZSS, BDS, GLONASS)
    7: E5b B2/B2b     (Galileo)
    8: E5a+b E5a+b    (Galileo, BDS)
    9: S              (IRNSS)
    0: for type X     (all)

    Attribute:
    A = A channel     (Galileo,IRNSS,GLONASS)
    B = B channel     (Galileo,IRNSS,GLONASS)
    C = C channel     (Galiloe, IRNSS)
    C code-based  (SBAS, GPS, GLONASS, QZSS)
    D = Semi-codelss  (GPS)

    I = I channel     (GPS, Galileo, QZSS, BDS)
    L = L channel     (L2C GPS, QZSS)
    P channel     (GPS. QZSS)
    M = M code-based  (GPS)
    N = Codeless      (GPS)
    P = P code-based  (GPS, GLONASS)
    Pilot channel (BDS)

    Q = Q channel     (GPS, Galileo, QZSS, BDS)
    S = D channel     (GPS, Galileo, QZSS, BDS)
    M channel     (L2C GPS, QZSS)

    W = Based on Z-tracking (GPS)
    X = B+C channels  (Galileo, IRNSS)
    I+Q channels  (GPS, IRNSS)
    M+L channels  (GPS, QZSS)
    D+P channels  (GPS, QZSS, BDS)

    Y = Y code based  (GPS)
    Z = A+B+C channels(Galileo)
    D+P channels  (BDS)
    --------------------------------------------------------------------------------------------------------------------------
    """
    ## -- Setting None arguments
    if readSS is None:
        readSS = 1
    if readLLI is None:
        readLLI = 1
    if includeAllGNSSsystems is None and desiredGNSSsystems is None:
        includeAllGNSSsystems = 1
    if includeAllObsCodes is None and desiredObsCodes is None:
        includeAllObsCodes = 1
    if desiredGNSSsystems is None:
        desiredGNSSsystems = ['G','R','E','C']
    if desiredObsCodes is None:
        desiredObsCodes = ['C','L','S','D']
    if desiredObsBands is None:
        desiredObsBands = list(np.arange(1,10))

    ## Get the start time
    t = time.process_time()

    ## - Initialize variables in case of input error
    GNSS_obs       = np.nan
    GNSS_LLI       = np.nan
    GNSS_SS        = np.nan
    GNSS_SVs       = np.nan
    time_epochs    = np.nan
    nepochs        = np.nan
    GNSSsystems    = np.nan
    obsCodes       = np.nan
    approxPosition = np.nan
    max_sat        = np.nan
    tInterval      = np.nan
    markerName     = np.nan
    rinexVersion   = np.nan
    recType        = np.nan
    timeSystem     = np.nan
    leapSec        = np.nan
    gnssType       = np.nan
    rinexProgr     = np.nan
    rinexDate      = np.nan
    antDelta       = np.nan
    tFirstObs      = np.nan
    tLastObs       = np.nan
    clockOffsetsON = np.nan
    GLO_Slot2ChannelMap = np.nan

    ### --- Dict for storing data
    GNSS_obs = {}
    GPS = {}
    GLONASS = {}
    Galileo = {}
    BeiDou = {}
    GPS_LLI = {}
    GLONASS_LLI = {}
    Galileo_LLI = {}
    BeiDou_LLI = {}
    GPS_SS = {}
    GLONASS_SS = {}
    Galileo_SS = {}
    BeiDou_SS = {}

    ## -- Test if readSS is boolean
    if readSS!=1 and readSS!=0:
        print('INPUT ERROR(readRinexObs211): The input argument readSS must be either 1 or 0')
        success = 0
        return


    ## -- Test if readLLI is boolean
    if readLLI!=1 and readLLI!=0:
        print('INPUT ERROR(readRinexObs211): The input argument readLLI must be either 1 or 0')
        success = 0
        return

    max_GPS_PRN     = 36 # Max number of GPS PRN in constellation
    max_GLONASS_PRN = 36 # Max number of GLONASS PRN in constellation
    max_Galileo_PRN = 36 # Max number of Galileo PRN in constellation
    max_Beidou_PRN  = 100 # Max number of BeiDou PRN in constellation

    #  Read header of observation file
    [success, rinexVersion, gnssType, markerName, recType, antDelta,\
    GNSSsystems,numOfObsCodes, obsCodes, obsCodeIndex,tFirstObs, tLastObs, tInterval, \
    timeSystem, _, clockOffsetsON, rinexProgr, rinexDate,leapSec, approxPosition, GLO_Slot2ChannelMap, _, fid] = \
    rinexReadObsFileHeader211(filename, includeAllGNSSsystems, includeAllObsCodes,desiredGNSSsystems, desiredObsCodes, desiredObsBands)

    if success==0:
        return

    # Compute number of epochs with observations
    nepochs, tLastObs, tInterval, success = rinexFindNEpochs211(filename, tFirstObs, tLastObs, tInterval) #computes number of epochs in observation file

    if success==0:
        return

    # Number of GNSS systems
    nGNSSsystems = len(GNSSsystems)

    # Declare data cells, arrays and matrices
    GNSS_SVs = {}
    max_sat  =  np.zeros([nGNSSsystems,1])
    t_week = []
    t_tow = []

    GNSSsystems_full_names =  [""]*nGNSSsystems
    # Making dict for storin LLI and SS
    GNSS_LLI = {}
    GNSS_SS = {}


    # Create array for max_sat. Initialize cell elements in cell arrays
    for k in np.arange(0, nGNSSsystems):
        # Always extract scalar from max_sat[k] for shape arguments
        max_sat_k = int(max_sat[k].item()) if hasattr(max_sat[k], "item") else int(max_sat[k])
        if GNSSsystems[k+1] == 'G':
            GNSS_SVs['G'] = np.zeros([nepochs, max_sat_k + 1])
            GNSSsystems_full_names[k] = "GPS"
        elif GNSSsystems[k+1] == 'R':
            GNSS_SVs['R'] = np.zeros([nepochs, max_sat_k + 1])
            GNSSsystems_full_names[k] = "GLONASS"
        elif GNSSsystems[k+1] == 'E':
            GNSS_SVs['E'] = np.zeros([nepochs, max_sat_k + 1])
            GNSSsystems_full_names[k] = "Galileo"
        elif GNSSsystems[k+1] == 'C':
            GNSS_SVs['C'] = np.zeros([nepochs, max_sat_k + 1])
            GNSSsystems_full_names[k] = "BeiDou"
        else:
            print(f'ERROR(readRinexObs211): Only following GNSS systems are compatible with this program: GPS, GLONASS, Galileo, Beidou. {GNSSsystems[k]} is not valid')
            return

        curr_sys = GNSSsystems[k+1]
        GNSS_obs[curr_sys] = np.zeros([max_sat_k, numOfObsCodes[k], nepochs])


        curr_sys = GNSSsystems[k+1]
        GNSS_obs[curr_sys] = np.zeros([int(max_sat[k]), numOfObsCodes[k], nepochs])

        # Preallocation LLI and SS
        if readLLI:
            GNSS_LLI[curr_sys] = np.zeros([nGNSSsystems,1])
        else:
            GNSS_LLI[curr_sys] = np.nan

        if readSS:
            GNSS_SS[curr_sys] = np.zeros([nGNSSsystems,1])
        else:
            GNSS_SS[curr_sys] = np.nan


    GNSS_names = dict(zip(['G', 'R', 'E', 'C'],['GPS','GLONASS','Galileo','Beidou']))
    current_epoch      = 0

    # Initialize progress bar
    n_update_break = int(np.floor(nepochs/10)) #number of epoch before updating progressbar
    bar_format = '{desc}: {percentage:3.0f}%|{bar}| ({n_fmt}/{total_fmt})'

    with tqdm(total=100,desc ="Rinex observations are being read" , position=0, leave=True, bar_format=bar_format) as pbar:
        while 1:
            #  Read Obs Block Header
            success, _, _, date, numSV,SVlist_, eof = rinexReadObsBlockHead211(fid)
            if success==0 or eof==1:
                break
            # Read current block of observations
            success, Obs,SVlist, numSV, LLI, SS, eof = rinexReadObsBlock211(fid, numSV, numOfObsCodes, GNSSsystems, obsCodeIndex, readSS, readLLI, SVlist_)
            if success ==0 or eof==1:
                break

            current_epoch = current_epoch + 1
            ## -- Update progress bar every n_update_break epochs
            if np.mod(current_epoch, n_update_break) == 0:
                pbar.update(10)

            ## Convert date to GPS-week and "time-of-week"
            date[0] = float(str(tFirstObs[0][0])[0:2] + str(int(date[0])))  # change from 20 to 2020 to get tow, week correct
            week, tow = date2gpstime(int(date[0]), int(date[1]), int(date[2]), int(date[3]), int(date[4]), int(date[5]))

            ## Store GPS-week and "time-of-week" of current epoch
            t_week.append(week)
            t_tow.append(tow)
            time_epochs = np.column_stack((t_week,t_tow))

            ## Number of satellites with observations in this epoch, for each GNSS system
            nGNSS_sat_current_epoch = np.zeros([nGNSSsystems,1])

            ## Initialize dummy variables
            GNSS_obs_dum = {}
            GNSS_LLI_dum = {}
            GNSS_SS_dum  = {}
            for k in np.arange(0,nGNSSsystems):
                ## -- Initialize arrays of dummy variables
                GNSS_obs_dum[k+1] = np.zeros([int(max_sat[k]) +1, numOfObsCodes[k]]) # added +1 12.11.2022 to get PRN36 sats
                GNSS_LLI_dum[k+1] = np.zeros([int(max_sat[k]) +1, numOfObsCodes[k]])  # added +1 12.11.2022 to get PRN36 sats
                GNSS_SS_dum[k+1]  = np.zeros([int(max_sat[k]) +1, numOfObsCodes[k]]) # added +1 12.11.2022 to get PRN36 sats

            ## -- Iterate through satellites of epoch and store obs, LLI and SS
            for sat in np.arange(0,numSV):
                ## -- Get index of current GNSS system
                curr_sys = SVlist[sat][0]
                GNSSsystemIndex = [i for i in GNSSsystems if GNSSsystems[i]==curr_sys][0]
                ## --Increment amount of satellites this epoch for this GNSS system
                nGNSS_sat_current_epoch[GNSSsystemIndex-1] = nGNSS_sat_current_epoch[GNSSsystemIndex-1] + 1
                SV = int(SVlist[sat][1:3]) # Get just PRN number
                nObsTypes_current_sat = numOfObsCodes[GNSSsystemIndex-1] # Number of obs types for current satellite
                ## -- Store observations, LLI, and SS of current satellite this epoch
                GNSS_obs_dum[GNSSsystemIndex][SV][0:nObsTypes_current_sat] = Obs[sat,0:nObsTypes_current_sat] # removed -1 due to lack of C5X obs
                if readLLI:
                    GNSS_LLI_dum[GNSSsystemIndex][SV][0:nObsTypes_current_sat] = LLI[sat, 0:nObsTypes_current_sat] # fjernet "-1" 13.11.2022 siste koloonnen ble ikke med pga -1
                if readSS:
                    GNSS_SS_dum[GNSSsystemIndex][SV][0:nObsTypes_current_sat] = SS[sat, 0:nObsTypes_current_sat] # fjernet "-1" 13.11.2022

                # Store PRN number of current sat to PRNs of this epoch
                # GNSS_SVs[curr_sys][current_epoch-1, int(nGNSS_sat_current_epoch[GNSSsystemIndex-1])] = SV
                GNSS_SVs[curr_sys][current_epoch-1, int(nGNSS_sat_current_epoch[GNSSsystemIndex-1].item()) if hasattr(nGNSS_sat_current_epoch[GNSSsystemIndex-1], "item") else int(nGNSS_sat_current_epoch[GNSSsystemIndex-1])] = SV

            for k in np.arange(0,nGNSSsystems):
                curr_sys = GNSSsystems[k+1]
                ## --Set number of satellites with obs for each GNSS system this epoch
                # GNSS_SVs[curr_sys][current_epoch-1, 0]  = nGNSS_sat_current_epoch[k]
                GNSS_SVs[curr_sys][current_epoch-1, 0] = float(nGNSS_sat_current_epoch[k].item()) if hasattr(nGNSS_sat_current_epoch[k], "item") else float(nGNSS_sat_current_epoch[k])
                if curr_sys == 'G':
                    GPS[current_epoch] = GNSS_obs_dum[k+1]
                elif curr_sys == 'R':
                    GLONASS[current_epoch] = GNSS_obs_dum[k+1]
                elif curr_sys == 'E':
                    Galileo[current_epoch] = GNSS_obs_dum[k+1]
                    Galileo_LLI[current_epoch]  =GNSS_LLI_dum[k+1]
                elif curr_sys == 'C':
                    BeiDou[current_epoch] = GNSS_obs_dum[k+1]
                    BeiDou_LLI[current_epoch]  =GNSS_LLI_dum[k+1]


                if readLLI and curr_sys == 'G':
                    GPS_LLI[current_epoch]  =GNSS_LLI_dum[k+1]
                elif readLLI and curr_sys == 'R':
                    GLONASS_LLI[current_epoch] = GNSS_LLI_dum[k+1]
                elif readLLI and curr_sys == 'E':
                    Galileo_LLI[current_epoch]  = GNSS_LLI_dum[k+1]
                elif readLLI and curr_sys == 'C':
                    BeiDou_LLI[current_epoch]  = GNSS_LLI_dum[k+1]

                if readSS and curr_sys =='G':
                    GPS_SS[current_epoch] = GNSS_SS_dum[k+1]
                elif readSS and curr_sys =='R':
                    GLONASS_SS[current_epoch] = GNSS_SS_dum[k+1]
                elif readSS and curr_sys =='E':
                    Galileo_SS[current_epoch] = GNSS_SS_dum[k+1]
                elif readSS and curr_sys =='C':
                    BeiDou_SS[current_epoch]  = GNSS_SS_dum[k+1]



        ## -- Storing observation in dictionary
        GNSS_obs['G'] = GPS
        GNSS_obs['R'] = GLONASS
        GNSS_obs['E'] = Galileo
        GNSS_obs['C'] = BeiDou
        ## -- Storing "loss of lock indicaors"  in dict
        GNSS_LLI['G'] = GPS_LLI
        GNSS_LLI['R'] = GLONASS_LLI
        GNSS_LLI['E'] = Galileo_LLI
        GNSS_LLI['C'] = BeiDou_LLI
        ## -- Storing  SS in dict
        GNSS_SS['G'] = GPS_SS
        GNSS_SS['R'] = GLONASS_SS
        GNSS_SS['E'] = Galileo_SS
        GNSS_SS['C'] = BeiDou_SS

        del_sys = list(GNSS_obs.keys())
        for sys in del_sys: # Deleting systems with no observations
            if not GNSS_obs[sys]:
                del GNSS_obs[sys]

        ## -- Removing system not in desiredGNSSsystems
        del_GNSS_LLI = [k for k,v in GNSS_LLI.items() if k not in desiredGNSSsystems]
        del_GNSS_SS  = [k for k,v in GNSS_SS.items() if k not in desiredGNSSsystems]
        del_GNSS_obs = [k for k,v in GNSS_obs.items() if k not in desiredGNSSsystems]
        del_GNSS_names = [k for k,v in GNSS_names.items() if k not in desiredGNSSsystems]
        del_GNSSsystems = [k for k,v in GNSSsystems.items() if v not in desiredGNSSsystems]
        del_GNSS_SVs = [k for k,v in GNSS_SVs.items() if k not in desiredGNSSsystems]

        keep_idx = []
        for k, v in GNSSsystems.items():
            if v in desiredGNSSsystems:
                keep_idx.append(k-1)

        max_sat = max_sat[keep_idx]

        for k in del_GNSS_LLI:
            del GNSS_LLI[k]
        for k in del_GNSS_SS:
            del GNSS_SS[k]
        for k in del_GNSS_obs:
            del GNSS_obs[k]
        for k in del_GNSS_names:
            del GNSS_names[k]
        for k in del_GNSSsystems:
            del GNSSsystems[k]
        for k in del_GNSS_SVs:
            del GNSS_SVs[k]


        GNSSsystems = {i+1: v for i, v in enumerate(GNSSsystems.values())} #update keys

        if current_epoch!= nepochs and success == 1:
            print('ERROR(readRinexObs211): The amount of epochs calculated in advance(nepochs = %d) does not equal number og epochs prossesed(current_epoch = %d).\nCheck that header information concerning TIME OF FIRST OBS and TIME OF LAST OBS is correct.\n' %(nepochs, current_epoch))


        messages = {}
        if success == 1:
            messages[0]= 'INFO(readRinexObs211): The following GNSS systems have been read into the data:'
            for curr_sys in list(GNSSsystems.values()):
                k = list(GNSSsystems.values()).index(curr_sys)
                messages[k+1]= 'INFO(readRinexObs211): The following %s observation types have been registered:' % (GNSS_names[curr_sys])
                for obs in np.arange(0, len(obsCodes[k+1][curr_sys])):
                    if obs == 0:
                        messages[k+1]= messages[k+1] + ' %s' % (obsCodes[k+1][curr_sys][obs])
                    else:
                        messages[k+1]= messages[k+1] + ', %s' % (obsCodes[k+1][curr_sys][obs])

                if k == 0:
                    messages[0]= messages[0] + ' %s' % GNSS_names[GNSSsystems[1]]
                else:
                    messages[0]= messages[0] + ', %s' % GNSS_names[GNSSsystems[k+1]]

            for msg in np.arange(0,len(messages)):
                print(messages[msg])



        if readLLI:
            print('INFO(readRinexObs211): LLI have been read (if present in observation file)')
        else:
            print('INFO(readRinexObs211): LLI have not been read')


        if readSS:
            print('INFO(readRinexObs211): SS have been read (if present in observation file)')
        else:
            print('INFO(readRinexObs211): SS have not been read')


        ## --  Finding processing time
        et = time.process_time()  # get the end time
        e = et - t                # get execution time

        if e >= 3600:
            hours = np.floor(e/3600)
            minutes = np.floor((e-hours*3600)/60)
            seconds = e-hours*3600-minutes*60
            print('INFO(readRinexObs211): Total processing time: %d hours, %d minutes, %f seconds\n' % (hours, minutes, seconds))
        elif e>60:
            minutes = np.floor(e/60)
            seconds = e-minutes*60
            print('INFO(readRinexObs211): Total processing time: %d minutes, %f seconds\n' % (minutes, seconds))
        else:
            print('INFO(readRinexObs211): Total processing time: %f seconds\n\n' % (e))


    return GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
        obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
        rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success



def rinexFindNEpochs211(filename, tFirstObs, tLastObs, tInterval):
    """
    Function that computes number of epochs in Rinex 3.xx observation file.
    --------------------------------------------------------------------------------------------------------------------------
    INPUTS

    filename:         RINEX observation filename

    tFirstObs:        time stamp of the first observation record in the RINEX
                      observations file; column vector
                      [YYYY; MM; DD; hh; mm; ss.sssssss];

    tLastObs:         time stamp of the last observation record in the RINEX
                      observations file; column vector
                      [YYYY; MM; DD; hh; mm; ss.sssssss]. If this information
                      was not available in rinex observation header the
                      default value is Nan. In this case the variable is
                      determined in this function

    tInterval:        observations interval; seconds. If this information
                      was not available in rinex observation header the
                      default value is Nan. In this case the variable is
                      determined in this function.
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS:
    -------

    nepochs:          number of epochs in Rinex observation file with
                      observations

    tLastObs:         time stamp of the last observation record in the RINEX
                      observations file; column vector
                      [YYYY, MM, DD, hh, mm, ss.sssssss]. If this information
                      was not available in rinex observation header the
                      default value is Nan. In this case the variable is
                      determined in this function

    tInterval:        observations interval; seconds. If this information
                      was not available in rinex observation header the
                      default value is Nan.

    success:                  Boolean. 1 if the function seems to be successful,
                              0 otherwise
    --------------------------------------------------------------------------------------------------------------------------

    ADVICE: The function rinexFindNEpochs() calculates the amount of observation epochs in
    advance. This calculation will be incredibly more effective if TIME OF
    LAST OBS is included in the header of the observation file. It is
    strongly advized to manually add this information to the header if it is
    not included by default.
    --------------------------------------------------------------------------------------------------------------------------
    """
    # #  Testing input arguments
    success = 1
    nepochs = 0

    ## --Test if filename is valid format
    if type(filename) is not str:
        raise TypeError('INPUT ERROR(rinexFindNEpoch): The input argument filename'\
            'is of type %s. Must be of type string' % type(filename))


    ## --Open observation file
    fid = open(filename, 'rt')
    seconds_in_a_week = 604800
    #  tLastObs is in header
    if ~np.all(np.isnan(tLastObs)):  # endret 07.12.2022
        #  tInterval is not in header
        if np.isnan(tInterval):
            tInterval_found = 0
            first_epoch_found = 0
            #  calculate tInterval
            while not tInterval_found: #  calculate tInterval
                line = fid.readline().rstrip()
                #  start of new epoch
                pattern = r"\s*(\d{2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2}\.\d+)\s*"
                if re.match(pattern, line):
                # if '>' in line:
                    if not first_epoch_found: #  first epoch
                        first_epoch_time = [el for el in line.split(" ") if el != ""][0:6]
                        first_epoch_time = [float(el) for el in first_epoch_time]
                        first_epoch_found = 1
                    else: #  seconds epoch
                        second_epoch_time = [el for el in line.split(" ") if el != ""][0:6]
                        second_epoch_time = [float(el) for el in second_epoch_time]
                        tInterval = second_epoch_time[5]-first_epoch_time[5]
                        tInterval_found = 1

        # fid.close(); fid = open(filename, 'rt')
        file = open(filename, 'rt')
        tFirstObs = tFirstObs.astype(int)
        tLastObs = tLastObs.astype(int)
        rinex_lines = file.readlines()
        # idx_start = rinex_lines.index('                                                            END OF HEADER       \n')
        idx_start = [i for i, line in enumerate(rinex_lines) if line.strip() == "END OF HEADER"][0]
        rinex_lines = rinex_lines[idx_start::]
        pattern = r"\s*(\d{2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2}\.\d+)\s*"
        epoch_line = [line for line in rinex_lines if re.search(pattern, line)] # list with all the line thats defines a epoch
        nepochs = len(epoch_line)
        file.close()
    #  if tLastObs is not in header. Function counts number of epochs manually
    else:
    ## New code for finding last
        print('INFO(rinexFindEpochs211): The header of the rinex observation file does not contain TIME OF LAST OBS.\n' \
            'This will be calculated, but consider editing rinex header to include TIME OF LAST HEADER')

        ## -- Find tLastObs
        pattern = r"\s*(\d{2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2}\.\d+)\s*"
        tLastObs = find_match_in_file(filename, pattern)
        tLastObs = np.array([list(map(float, tLastObs))]).T
        try:
            tLastObs[0] = tFirstObs[0]
        except:
            tLastObs[0] = float('20' + str(tLastObs[0][0]))
            print("\n Neither tFirstObs or tLastObs exist in header. The year is therefor a guess for after\
                  year 2000")

        first_ep,second_ep = find_first_two_epochs(filename, pattern)

        ## Computing the tInterval if not present in the header
        if np.isnan(tInterval):
            first_epoch_time  = [float(el) for el in first_ep]
            second_epoch_time = [float(el) for el in second_ep]
            tInterval = second_epoch_time[5]-first_epoch_time[5]

        # nepochs = time_difference(tFirstObs, tLastObs)/tInterval
        # fid.close(); fid = open(filename, 'rt')
        file = open(filename, 'rt')
        rinex_lines = file.readlines()
        idx_start = [i for i, line in enumerate(rinex_lines) if line.strip() == "END OF HEADER"][0]
        rinex_lines = rinex_lines[idx_start::]
        pattern = r"\s*(\d{2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2}\.\d+)\s*"
        epoch_line = [line for line in rinex_lines if re.search(pattern, line)] # list with all the line thats defines a epoch
        nepochs = len(epoch_line)
        print('INFO(rinexFindNEpochs211): TIME OF LAST OBS has been found and amount of epochs have been computed')
        file.close()
    return int(nepochs), tLastObs, tInterval, success


def find_match_in_file(file_path, pattern):
    """
    Function that reads in a file backwards and looks for pattern from the
    bottom and up. Is used for finding tLastObs when not defined in header.
    Uses binary mode to save processing time.
    """
    with open(file_path, "rb") as f:
        # Go to the end of the file
        f.seek(0, 2)
        file_size = f.tell()
        # Read the file backwards
        while file_size > 0:
            # Read a chunk of the file
            chunk_size = min(1024, file_size)
            f.seek(file_size - chunk_size)
            chunk = f.read(chunk_size)
            # Split the chunk into lines
            lines = chunk.split(b"\n")
            # Process the lines in reverse order
            for line in reversed(lines):
                line = line.decode("utf-8").rstrip()
                match = re.match(pattern, line)
                if match:
                    return match.groups()
            # Repeat until the entire file has been processed
            file_size -= chunk_size
    # No match was found
    return None



def find_match_in_file_FWD(file_path, pattern):
    """
    Function that reads in a file and looks for pattern from the bottom and up.
    Is used for finding tLastObs when not defined in header. Uses binary mode
    to save processing time. Returns current line and next 5 lines after match.
    """
    with open(file_path, "rb") as f:
        # Go to the end of the file
        f.seek(0, 2)
        file_size = f.tell()
        # Read the file backwards
        while file_size > 0:
            # Read a chunk of the file
            chunk_size = min(1024, file_size)
            f.seek(file_size - chunk_size)
            chunk = f.read(chunk_size)
            # Split the chunk into lines
            lines = chunk.split(b"\n")
            # Process the lines in reverse order
            for i in range(len(lines)-1, -1, -1):
                line = lines[i].decode("utf-8").rstrip()
                match = re.match(pattern, line)
                if match:
                    return [line] + [lines[j].decode("utf-8").rstrip() for j in range(i+1, min(i+6, len(lines)))]
            # Repeat until the entire file has been processed
            file_size -= chunk_size
    # No match was found
    return None


def find_first_two_epochs(file_path, pattern):
    """
    Function that extracts the two first epochs for
    RINEX obs file to compute tInterval when not defined.
    """
    with open(file_path, 'r') as file:
        found_header = False
        count = 0
        matches = []
        for line in file:
            if not found_header:
                if line.strip() == "END OF HEADER":
                    found_header = True
                continue
            match = re.search(pattern, line)
            if match:
                count += 1
                matches.append(match.groups())
                if count >= 2:
                    break
        ep1, ep2  = matches
        ep1 = list(ep1)
        ep2 = list(ep2)
        return ep1,ep2


def find_nepochs(file_path, pattern):
    """Find number of epochs"""
    with open(file_path, 'r') as file:
        contents = file.read()
        return len(re.findall(pattern, contents))


def time_difference(arr1, arr2):
    """Find time difference between two arrays"""
    time1 = datetime(year=int(arr1[0]), month=int(arr1[1]), day=int(arr1[2]), hour=int(arr1[3]), minute=int(arr1[4]), second=int(arr1[5]))
    time2 = datetime(year=int(arr2[0]), month=int(arr2[1]), day=int(arr2[2]), hour=int(arr2[3]), minute=int(arr2[4]), second=int(arr2[5]))
    time_delta = (time2-time1).total_seconds()
    return time_delta

def rinexReadObsFileHeader211(filename, includeAllGNSSsystems, includeAllObsCodes,desiredGNSSsystems,desiredObsCodes, desiredObsBands):
    """
    Extracts relevant data from the header of a RINEX 3.xx GNSS observations
    file. Excludes undesired GNSS systems, obsevation codes and/or frequency
    bands.

    --------------------------------------------------------------------------------------------------------------------------
    INPUTS:
    ------

    filename:                     RINEX observation filename and path

    includeAllGNSSsystems:        Boolean, 0 or 1.
                                      1 = include alle GNSS systems
                                          (GPS, GLONASS, Galieo, BeiDou)
                                      0 = include only GNSSsystems
                                          specified in desiredGNSSsystems

    includeAllobsCodes:           Boolean, 0 or 1.
                                      1 = include all valid obsCodes
                                      0 = include only obsCodes
                                          specified in desiredobsCodes

    desiredGNSSsystems:           string array containing desired GNSSsystems
                                  to be included, ex. ["G", "E", "C"]

    desiredobsCodes:              string array containing desired obsCodes to
                                  be included, ex. ["C", "L", "S", "D"]

    desiredObsBands:              array of desired obs Bands to be included,
                                  ex [1, 5]

    NOTE: If both includeAllGNSSsystems and includeAllobsCodes Boolean are 1
          then the last three input arguments are optional to include and may
          be left blank without en error.
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS:
    -------

    success:                      1 if the reading of the RINEX observations
                                  file seems to be successful, 0 otherwise

    rinexVersion:                 string. rinex observation file version

    gnssType:                     GNSS system of the satellites observed; can
                                  be 'G', 'R', 'E', 'C' or 'M' that stand for
                                  GPS, GLONASS, GALILEO, BeiDou or Mixed; char

    markerName:                   name of the antenna marker; '' if not
                                  specified

    recType:                      Receiver type, char vector

    antDelta:                     column vector ot the three components of
                                  the distance from the marker to the antenna,
                                  in the following order - up, east and north;
                                  null vector by default

    GNSSsystems:                  cell array containing codes of GNSS systems
                                  included in RINEX observationfile. Elements
                                  are strings. ex. "G" or "E"

    numOfObsCodes:                column vector containing number of observation
                                  types for each GNSS system. Order is the same
                                  as GNSSsystems

    obsCodes:                     Cell that defines the observation
                                  codes available for all GNSS system. Each
                                  cell element is another cell containing the
                                  codes for that GNSS system. Each element in
                                  this cell is a string with three-characters.
                                  The first character (a capital letter) is
                                  an observation code ex. "L" or "C". The
                                  second character (a digit) is a frequency
                                  code. The third character(a Capital letter)
                                  is the attribute, ex. "P" or "X"

    obsCodeIndex:                 cell with one cell element for each GNSS
                                  system. Order is the same as GNSSsystems.
                                  Each cell element contains an array of
                                  indices. These indices indicate the
                                  observation types that should be read
                                  for each GNSS system. ex. If one index for
                                  GPS is 1 then the first observation type
                                  for GPS should  be read.

    tFirstObs:                    time stamp of the first observation record
                                  in the RINEX observations file; column vector
                                  [YYYY; MM; DD; hh; mm; ss.sssssss];
                                  THIS IS CRITICAL DATA

    tLastObs:                     time stamp of the last observation record
                                  in the RINEX observations file; column vector
                                  [YYYY; MM; DD; hh; mm;ss.sssssss].
                                  NaN by default.
                                  THIS IS RINEX V.2 OPTIONAL DATA

    tInterval:                    observations interval; seconds.

    timeSystem:                   three-character code string of the time
                                  system used for expressing tfirstObs;
                                  can be GPS, GLO or GAL;

    numHeaderLines:               number of lines in header

    rinexProgr:                   name of the software used to produce de
                                  RINEX GPS obs file; '' if not specified

    rinexDate:                    date/time of the RINEX file creation; ''
                                  if not specified; char

    leapSec:                      number of leap seconds since 6-Jan-1980.
                                  UTC=GPST-leapSec. NaN by default.
                                  THIS IS RINEX V.2 OPTIONAL DATA

    approxPosition:               array containing approximate position from
                                  rinex observation file header. [X, Y, Z]

    GLO_Slot2ChannelMap:          map container that maps GLONASS slot
                                  numbers to their respective channel number.
                                  GLO_Slot2ChannelMap(slotnumber)

    eof:                          end-of-file flag; 1 if end-of-file was reached,
                                  0 otherwise

    fid:                          Matlab file identifier of a Rinex
                                  observations text file
    --------------------------------------------------------------------------------------------------------------------------

    According to RINEX V.2 these codes are:

       Observation codes:
       ------------------
       C: Pseudorange
          GPS: C/A, L2C
          Glonass: C/A
          Galileo: All
       L: Carrier phase
       D: Doppler frequency
       S: Raw signal strengths or SNR values as given by the receiver for the
          respective phase observations
       I: Ionosphere phase delay
       X: Receiver channel numbers

       Frequency code:
       ---------------
       GPS Glonass Galileo SBAS
       1: L1 G1 E1 B1    (GPS,QZSS,SBAS,BDS)
       2: L2 G2 B1-2     (GLONASS)
       4: G1a            (Galileo)
       5: L5 E5a B2/B2a  (GPS, QZSS, SBAS, IRNSS)
       6: L6 E6 B3 G2a   (Galileo, QZSS, BDS, GLONASS)
       7: E5b B2/B2b     (Galileo)
       8: E5a+b E5a+b    (Galileo, BDS)
       9: S              (IRNSS)
       0: for type X     (all)

       Attribute:
       ----------
       A = A channel     (Galileo,IRNSS,GLONASS)
       B = B channel     (Galileo,IRNSS,GLONASS)
       C = C channel     (Galiloe, IRNSS)
           C code-based  (SBAS, GPS, GLONASS, QZSS)
       D = Semi-codelss  (GPS)

       I = I channel     (GPS, Galileo, QZSS, BDS)
       L = L channel     (L2C GPS, QZSS)
           P channel     (GPS. QZSS)
       M = M code-based  (GPS)
       N = Codeless      (GPS)
       P = P code-based  (GPS, GLONASS)
           Pilot channel (BDS)

       Q = Q channel     (GPS, Galileo, QZSS, BDS)
       S = D channel     (GPS, Galileo, QZSS, BDS)
           M channel     (L2C GPS, QZSS)

       W = Based on Z-tracking (GPS)
       X = B+C channels  (Galileo, IRNSS)
           I+Q channels  (GPS, IRNSS)
           M+L channels  (GPS, QZSS)
           D+P channels  (GPS, QZSS, BDS)

       Y = Y code based  (GPS)
       Z = A+B+C channels(Galileo)
           D+P channels  (BDS)
    -------------------------------------------------------------------------------------------------------------------------
    """
    eof         = 0
    success     = 1
    antDelta    = []
    timeSystem  = ''
    tFirstObs   = []
    tLastObs    = np.nan
    tInterval   = np.nan
    rinexProgr  = np.nan
    rinexDate   = np.nan
    obsCodes    = {}
    GNSSsystems = {}
    gnssType    = ""
    markerName  = ""
    numHeaderLines  = 0
    clockOffsetsON  = 0
    numGNSSsystems  = 0
    leapSec         = np.nan
    numOfObsCodes   = []
    rinexHeader     = {}
    approxPosition  = [0, 0, 0]
    obsCodeIndex = {}
    rinexVersion = np.nan
    recType = np.nan
    GLO_Slot2ChannelMap = np.nan

    ## -------Testing input arguments
    # Test if filename is valid format
    if type(filename) != str:
        print('INPUT ERROR(rinexReadsObsHeader211): The input argument filename is of type %s.\n Must be of type string or char' %(type(filename)))
        success = 0
        fid     = 0
        return success


    ## -- Open rinex observation file
    fid = open(filename,'r')
    if os.stat(filename).st_size == 0:
        raise ValueError('ERROR: This file seems to be empty')

    while 1: # Gobbling the header
        numHeaderLines = numHeaderLines + 1
        line = fid.readline().rstrip()

        if 'END OF HEADER' in line:
            break

        if numHeaderLines == 1: # if first line of header
            # store rinex version
            rinexVersion = line[0:9].strip()
            # store rinex type, ex. "N" or "O"
            rinexType = line[20]
            # if rinex file is not an observation file
            if rinexType != 'O':  # Rinex file is oservation file
                print('ERROR(rinexReadObsFileHeader211): the file is not a RINEX observations data file!')
                success = 0
                fid.close()
                return

            ## -- Check gnss type  ## Changend indent here 09.12.2022 (was apart of the if test above earlier, and thats wrong)
            gnssType = line[40] # reads the GNSS system type
            if gnssType not in [' ', 'G', 'R', 'C', 'E', 'M' ]:
                if gnssType in ['J', 'I', 'S']:
                    print('ERROR(rinexReadObsFileHeader211): This software is meant for reading GNSS data only.\
                           %s is an invalid satellite system type.' %(gnssType))
                else:
                    print('ERROR(rinexReadObsFileHeader211): %s is an unrecognized satellite system type.' %(gnssType))

                success = 0
                fid.close()
            ## -- If no system type, set G
            if gnssType == ' ':
                gnssType = 'G'


        if 'PGM / RUN BY / DATE' in line:
            rinexProgr = line[0:20] # rinex program
            rinexDate = line[40:60] # rinex date

        if 'MARKER NAME' in line:
            markerName = line.strip() # markername

        ## if no marker name, "MARKER" is read, so set to blank
        if 'Marker' in markerName:
            markerName = ''

        if 'ANTENNA: DELTA H/E/N' in line:
            for k in np.arange(0,3):
                line_ = [el for el in line.split(" ") if el != ""]
                antDelta = [line_[0],line_[1],line_[2]]



        if '# / TYPES OF OBSERV' in line:
            pattern_ep1 = r"\s*(\d{2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2}\.\d+)\s*"
            pattern_ep2 = r'\s+([A-Z]\d{2})' #pattern for next line in block header
            first_block_lines = find_match_in_file_FWD(filename, pattern_ep1)
            sat_lines = [line for line in first_block_lines if re.match(pattern_ep1, line) or re.match(pattern_ep2, line)]
            sat_lines = ''.join([line.strip() for line in sat_lines if re.match(pattern_ep1, line) or re.match(pattern_ep2, line)])[29:] #[29:] is remove the first part of the string
            GNSSsystems = list(set(re.findall(r"[a-zA-Z]", sat_lines)))
            GNSSsystems = {i+1: GNSSsystems[i] for i in range(len(GNSSsystems))} # make dict there index in list becomes key in dict
            numGNSSsystems = len(GNSSsystems)

            line_dum = []
            while '# / TYPES OF OBSERV' in line:
                line = line[0:60]     # deletes 'SYS / # / OBS TYPES'
                line_ = [el for el in line.split(" ") if el != ""]
                line_dum.extend(line_)
                line = fid.readline().rstrip()
            line_ = line_dum
            nObs = int(line_.pop(0)) # assingning nObs to variable and removing it from the list
            undesiredobsCodeIndex = []
            desiredObsCodeIndex = []
            ## is Sys amoung desired GNSS systems
            GNSSSystemObsCodes = {}  # Reset cell of obsCodes for this GNSS system
            obsCode_list = []
            for k in np.arange(0,nObs):
                obsCode = line_.pop(0)
                # Checking if obsCode is valid
                if len(obsCode) != 2 or obsCode[0] not in ['C', 'L', 'D','S', 'I', 'X','P'] or  \
                          obsCode[1] not in ['0', '1', '2','3', '4', '5', '6', '7', '8', '9']:
                    print('ERROR (rinexReadsObsHeader211):  obsCode %s is a not a standard RINEX 2.11 observation type!' %(obsCode))

                ## is obsCode amoung desired obscodes and frequency bands
                if includeAllObsCodes or obsCode[0] in desiredObsCodes and int(obsCode[1]) in desiredObsBands:
                     ## store obsCode if amoung desire obsCodes
                    obsCode_list.append(obsCode)
                    for sys in GNSSsystems.values():
                        GNSSSystemObsCodes[sys] =  obsCode_list
                        GNSSSystemObsCodes[sys] =  obsCode_list

                    desiredObsCodeIndex.append(k)
                else:
                    # store index of discareded obsCode
                    undesiredobsCodeIndex.append(k)

            for sys_indx in GNSSsystems.keys():
                sys = GNSSsystems[sys_indx]
                numOfObsCodes.append(len(GNSSSystemObsCodes[sys]))
                obsCodes[sys_indx] = GNSSSystemObsCodes
                obsCodes[sys_indx] = GNSSSystemObsCodes

            obsCodeIndex[numGNSSsystems] = desiredObsCodeIndex # Store indices of desired obsCodes

        if 'PRN / # OF OBS' in line:
            system_list = []
            while 'PRN / # OF OBS' in line:
                line = line[0:60]     # deletes 'SYS / # / OBS TYPES'
                line_ = [el for el in line.split(" ") if el != ""]
                Sys = line_.pop(0)[0] # assingning nObs to variable and removing it from the list
                if Sys not in ["G","R","E","C"]: # added this line 29.01.2023 to fix bug where Only one system and several lines with Obscodes in rinex file
                    continue
                if (includeAllGNSSsystems and Sys in ["G", "R", "E", "C"] or Sys in desiredGNSSsystems):
                    # numGNSSsystems  = numGNSSsystems + 1 # increment number of GNSS systems
                    system_list  = [Sys] + [s for s in system_list if s not in [Sys]]
                    numGNSSsystems  = len(system_list)

                    if Sys not in GNSSsystems.values():
                        GNSSsystems[numGNSSsystems] = str(Sys) # Store current GNSS system
                    else:
                        pass
                    # GNSSsystems[numGNSSsystems] = str(Sys) # Store current GNSS system
                    GNSSSystemObsCodes[Sys] =  obsCode_list
                    numOfObsCodes.append(len(GNSSSystemObsCodes[Sys]))
                    obsCodes[numGNSSsystems] = GNSSSystemObsCodes
                    obsCodeIndex[numGNSSsystems] = desiredObsCodeIndex # Store indices of desired obsCodes
                    GNSSSystemObsCodes[Sys] =  obsCode_list
                numHeaderLines = numHeaderLines + 1
                line = fid.readline().rstrip()
                line_ = [el for el in line.split(" ") if el != ""]



        if 'TIME OF FIRST OBS' in line:
            line = line[0:60]     #  deletes 'TIME OF FIRST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            for k in np.arange(0,6):
                tok = line_.pop(0)  # finds the substrings containing the components of the time of the first observation
                                      #(YYYY; MM; DD; hh; mm; ss.sssssss) and specifies
                                      # the Time System used in the
                                      # observations file (GPST, GLOT or
                                      # GALT)
                if k ==0:
                    yyyy = int(tok)
                elif k ==1:
                    mm = int(tok)
                elif k ==2:
                    dd = int(tok)
                elif k ==3:
                    hh = int(tok)
                elif k ==4:
                    mnt = int(tok)
                elif k ==5:
                    ss = float(tok)


            tFirstObs = np.array([[yyyy],[mm],[dd],[hh],[mnt],[ss]])

            # Get Time system
            aux = line_.pop(0)
            # aux = strtok(line);
            if aux == 'GPS':
                timeSystem = 'GPS'
            elif aux == 'GLO':
                timeSystem = 'GLO'
            elif aux == 'GAL':
                timeSystem = 'GAL'
            elif aux == 'BDT':
                timeSystem = 'BDT'

            else:
                if gnssType == 'G':
                    timeSystem = 'GPST'
                elif gnssType == 'R':
                    timeSystem = 'GLOT'
                elif gnssType == 'E':
                    timeSystem = 'GALT'
                elif gnssType == 'C':
                    timeSystem = 'BDT'
                else:
                    print('CRITICAL ERROR (rinexReadsObsHeader211):\n' \
                                       'The Time System of the RINEX observations file '\
                                       'isn''t correctly specified!\n')
                    success = 0
                    fid.close()


        if 'TIME OF LAST OBS' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            for k in np.arange(0,6):
                tok = line_.pop(0)
                if k ==0:
                    yyyy = int(tok)
                elif k ==1:
                    mm = int(tok)
                elif k ==2:
                    dd = int(tok)
                elif k ==3:
                    hh = int(tok)
                elif k ==4:
                    mnt = int(tok)
                elif k ==5:
                    ss = float(tok)

            tLastObs = np.array([[yyyy],[mm],[dd],[hh],[mnt],[ss]])

        if 'INTERVAL' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            tInterval = float(line_.pop(0))



          ## -- This is an optional record!
          # if 'RCV CLOCK OFFS APPL' in line:
          #     if (strtok(line)=='0'):
          #         clockOffsetsON = 0;
          #     elif (strtok(line)=='1'):
          #         clockOffsetsON = 1;
          #     else:
          #         success = 0;
          #         print('ERROR (rinexReadsObsHeader211): unrecognized receiver clock offsets flag!')
          #         fid.close()


           ## This is an optional record
        if 'LEAP SECONDS' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            leapSec = int(line_.pop(0))



           ## -- store approximate receiver position
        if 'APPROX POSITION XYZ' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            approxPosition = np.array([[float(line_[0])],[float(line_[1])],[float(line_[2])]])


         ## GLOANSS SLOTS
        if 'GLONASS SLOT / FRQ #' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            nGLOSat = int(line_.pop(0))
            slotNumbers = np.array([])
            channels = np.array([])
            for k in np.arange(0,nGLOSat):
                slotNumber = line_.pop(0)[1::]
                channel = int(line_.pop(0))
                slotNumbers = np.append(slotNumbers,slotNumber)
                channels = np.append(channels,channel)

                if np.mod(k+1, 8) == 0 and k+1 != 24:
                    # line = fgetl(fid); # end of line is reached so read next line
                    line = fid.readline().rstrip()
                    numHeaderLines = numHeaderLines + 1
                    line = line[0:60]     #  deletes 'TIME OF LAST OBS'
                    line_ = [el for el in line.split(" ") if el != ""]
                elif np.mod(k+1, 8) == 0 and k+1 == 24:
                    break

            GLO_Slot2ChannelMap = dict(zip(slotNumbers.astype(int),channels.astype(int)))

        if 'REC # / TYPE / VERS' in line:
            recType = line[20:40]



     # End of Gobbling Header Loop
    for k in np.arange(1,numGNSSsystems+1):
        # Give info if any of GNSS systems had zero of desired obscodes.
        if numOfObsCodes[k-1] == 0 or sum(tFirstObs) == 0:
            if GNSSsystems[k] == 'G':
                print('INFO: (rinexReadsObsHeader211)\nNone of the GPS satellites had any of the desired obsCodes\n\n')
            elif GNSSsystems[k] == 'R':
                print('INFO: (rinexReadsObsHeader211)\nNone of the GLONASS satellites had any of the desired obsCodes\n\n')

    ## store rinex header info
    rinexHeader['rinexVersion'] =rinexVersion
    rinexHeader['rinexType'] = rinexType
    rinexHeader['gnssType'] =gnssType
    rinexHeader['rinexProgr'] =rinexProgr
    rinexHeader['rinexDate'] =rinexDate


    print('INFO(rinexReadObsFileHeader211): Rinex header has been read')

    return success, rinexVersion, gnssType, markerName, recType, antDelta, GNSSsystems, numOfObsCodes, \
    obsCodes, obsCodeIndex,tFirstObs, tLastObs, tInterval,timeSystem, numHeaderLines, clockOffsetsON, \
    rinexProgr, rinexDate,leapSec, approxPosition, GLO_Slot2ChannelMap, eof, fid



def rinexReadObsBlock211(fid, numSV, nObsCodes, GNSSsystems, obsCodeIndex, readSS, readLLI, SVlist):
    """
    Reads all the observations from a RINEX observation block.

    Positioned at the beginning of the line immediately after the header of the
    observations block, reads all the observations in this block of a RINEX
    observations file. This function is meant to be used after using function
    rinexReadObsFileHeader211

    Based in the work of Antonio Pestana, rinexReadObsBlock211, March 2015
    --------------------------------------------------------------------------------------------------------------------------
    INPUTS:
    -------

    fid:                  Matlab file identifier of a Rinex observations text file

    numSV:                total number of satellites with observations in
                          current observation block, integer

    numOfObsCodes:        column vector containing number of observation
                          types for each GNSS system. Order is the same as
                          GNSSsystems

    GNSSsystems:          cell array containing codes of GNSS systems included
                          in RINEX observationfile. Elements are strings.
                          ex. "G" or "E"

    obsCodeIndex:         cell with one cell element for each GNSS system.
                          Order is the same as GNSSsystems. Each cell element
                          contains an array of indices. These indices
                          indicate the observation types that should be
                          read for each GNSS system. ex. If one index for
                          GPS is 1 then the first observation type for GPS
                          should be read.

    readSS:                   Boolean, 0 or 1.
                              1 = read "Signal Strength" Indicators
                              0 = do not read "Signal Strength" Indicators

    readLLI:                  Boolean, 0 or 1.
                              1 = read "Loss-Of-Lock Indicators"
                              0 = do not read "Loss-Of-Lock Indicators"
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS:
    --------

    success:               Boolean. 1 if the function seems to be successful,
                          0 otherwise

    Obs:                  matrix [numSV x max_nObs] that stores all
                          observations of this observation block. max_nObs
                          is the highest number of observation codes that
                          any of the GNSS systems have. Which observation
                          types that are associated with what collumn will
                          vary between GNSS systems. SVlist will give
                          overview of what GNSS system each row is connected
                          to

    SVlist:               column cell [numSV x 1] that conatins the
                          identification code of each line of observation
                          block. ex. "G21". numSV is total number of
                          satellites minus amount of satellites removed.

    numSV:                numSV, unlike the input of same name, is the total
                          number of satellites minus amount of satellites
                          removed.

    LLI:                  matrix [numSV x max_nObs] that stores all
                          "loss-of-lock" indicators of this observation block.
                          max_nObs is the highest number of observation codes
                          that any of the GNSS systems have. Which observation
                          types that are associated with what collumn will
                          vary between GNSS systems. SVlist will give
                          overview of what GNSS system each row is connected
                          to

    SS:                   matrix [numSV x max_nObs] that stores all
                          "signal strength" indicators of this observation block.
                          max_nObs is the highest number of observation codes
                          that any of the GNSS systems have. Which observation
                          types that are associated with what collumn will
                          vary between GNSS systems. SVlist will give
                          overview of what GNSS system each row is connected
                          to

    eof:                  end-of-file flag; 1 if end-of-file was reached,
                          0 otherwise
    --------------------------------------------------------------------------------------------------------------------------
    """
    ## Initialize variables in case of input error
    success                     = np.nan
    eof                         = np.nan
    max_n_obs_Types             = np.nan
    Obs                         = np.nan
    LLI                         = np.nan
    SS                          = np.nan
    removed_sat                 = np.nan
    desiredGNSSsystems          = np.nan
    obsCodeIndex = obsCodeIndex[list(obsCodeIndex.keys())[0]]

    ## Test type of numSV
    if type(numSV) != int:
        print(f'INPUT ERROR(rinexReadObsBlock211): The input argument numSV is of type {type(numSV)}.\n Must be of type double')
        success = 0
        return success

    nObsCodes = int(nObsCodes[0])
    success = 1
    eof     = 0

    # Highest number of obs codes of any GNSS system
    max_n_obs_Types = nObsCodes
    # Initialize variables
    Obs = np.empty([numSV, max_n_obs_Types])
    if nObsCodes > 5:
        factor = 2
        nLines = factor*numSV
    else:
        factor = 1
        nLines = factor*numSV

    if readLLI:
        LLI = np.empty([numSV, max_n_obs_Types])

    if readSS:
        SS  = np.empty([numSV, max_n_obs_Types])

    # number of satellites excluded so far
    removed_sat = 0
    desiredGNSSsystems = list(GNSSsystems.values())
    pattern = r"^\s*[0-9]+\.[0-9]+\s*"
    pattern2 = r'\s*\d+\.\d+'
    pattern4 = r'^\s*[A-Za-z]+[A-Za-z0-9]*(?:[\s]+[A-Za-z]+[A-Za-z0-9]*)*\s*$'
    ## -- Gobble up observation block
    for sat in np.arange(0,numSV):
        line = fid.readline().rstrip()
        while (not re.match(pattern,line) or not re.match(pattern2,line) or re.match(pattern4,line)) and line !="":
            sat_overview = line.split(" ")[-1]
            pattern3 = re.compile(r'[A-Z][0-9]{2}')
            sat_list = re.findall(pattern3, sat_overview)
            for s in sat_list:
                sys = s[0]
                if sys not in desiredGNSSsystems:
                    continue
                SVlist.append(s)
            line = fid.readline().rstrip()

        SV = SVlist[sat]
        if SV[0] not in desiredGNSSsystems:
            removed_sat +=1
        else:
            ## Index of current GNSS system
            GNSSsystemIndex = [i for i in GNSSsystems if GNSSsystems[i]==SV[0]][0]
            n_obs_current_system = nObsCodes
            nNew_line = 1 # counter to keep track in new line of same satellite
            for obs_num in np.arange(0, n_obs_current_system):
                # obsIndex = obsCodeIndex[GNSSsystemIndex][obs_num]
                obsIndex = obsCodeIndex[obs_num]
                # charPos = 4+(obsIndex)*16
                if obsIndex < 5:
                    charPos = 1+(obsIndex)*16
                else:
                    charPos = 1+(obsIndex-5)*16

                ## check that the current observation of the current GNSS system
                ## is not on the list of obs types to be excluded
                ## stringlength of next obs.
                obsLen = min(14, len(line) - charPos)
                # read next obs
                newObs = line[charPos:charPos+obsLen].strip()
                # If observation missing, set to 0
                if newObs != '':
                    newObs = float(newObs)
                else:
                    newObs = 0
               # Store new obs
                Obs[sat - removed_sat, obs_num] = newObs

                if readLLI:
                # read LLI of current obs (if present)
                    if charPos+13<len(line): ## endret til < (kun mindre enn) 13.11
                        newLLI = line[charPos+13] # loss of lock indicator ### endret fra 14 til 13 den 13.11.2022
                    else:
                        newLLI = ' '
                    if newLLI.isspace():
                        newLLI = -999
                    else:
                        newLLI = int(newLLI)
                     # Store LLI
                    LLI[sat - removed_sat, obs_num] = newLLI

                if readSS:
                    # read SS of current obs (if present)
                    if charPos+14<len(line): ## endret til < (kun mindre enn) 13.11
                        newSS = line[charPos+14] # signal strength endret fra 15 til 14 den 13.11.2022
                    else:
                        newSS = ' '
                    # if no SS set to -999
                    if newSS.isspace():
                        newSS = -999
                    else:
                        newSS = int(newSS)

                    ## Store SS
                    SS[sat - removed_sat, obs_num]  = newSS

                    # if np.mod(obs_num+1, 5) == 0 and nObsCodes > 5 and factor*sat < factor*numSV:
                    if np.mod(obs_num+1, 5) == 0 and nObsCodes>5 and nNew_line < factor:
                        nNew_line += 1
                        line = fid.readline().rstrip()


    ## -- Update number og satellites after satellites have been excluded
    numSV = numSV - removed_sat
    ## --Remove empty arrays
    idx_keep = len(Obs) -1 -removed_sat + 1 # removing sats
    Obs = Obs[:idx_keep,:]
    return success, Obs, SVlist, numSV, LLI, SS, eof





def rinexReadObsBlockHead211(fid):
    """
    Reads the metadata in the head of a RINEX 3.xx observations block, NOT
    the header of the file.

    ATTENTION: Ignores all data in blocks with event flags with numbers
    greater than 1!!!

    Positioned in a RINEX V.2 GNSS observations text file at the beginning
    of an observation block. In rinex 3.xx the line starts with '> '

    --------------------------------------------------------------------------------------------------------------------------
     INPUTS

     fid:              Matlab identifier of an open RINEX V.2 GNSS
                       observations text file positioned at the beginning
                       of an observation block.
    --------------------------------------------------------------------------------------------------------------------------
     OUTPUTS

     success:          1 if function performs successfully, 0 otherwise

     epochflag:        Rinex observations epoch flag, as follows:
                           0: OK
                           1: power failure between previous and current epoch
                       From now on the "event flags":
                           2: start moving antenna
                           3: new site occupation
                           4: header information follows
                           5: external event (epoch is significant)

     clockOffset:          value of the receiver clock offset. If not present
                           in the metadata of the observations block
                           (it's optional RINEX V.2 data)it is assumed to be
                           zero. If not zero implies that epoch, code, and
                           phase data have been corrected by applying
                           realtime-derived receiver clock offset

     date:                 time stamp of the observations block. Six-elements column-vector
                           as follows:
                               year: four-digits year (eg: 1959)
                               month: integers 1..12
                               day: integers 1..31
                               hour: integers 0..24
                               minute: integers 0..60
                               second: reals 0..60

     SVlist:               Dictionary with the system and PRN number for the current epoch in obsfile. Ex {G:[3,4,20,24,28]
                                                                                                           R:[1,3,4,9,12]}

     numSV:                number of satellites with observations in with
                           observations. This will include all satellite
                           systems.
    --------------------------------------------------------------------------------------------------------------------------

    """
    ## -- Initialize variables
    success = 1
    eof     = 0
    date    = [0,0,0,0,0,0]
    numSV   = 0
    epochflag = 0
    clockOffset = 0
    noFlag = 1
    SVlist = np.nan
    line = fid.readline().rstrip()
    if not line:
        eof = 1
        print('\nINFO(rinexReadObsBlockHead211): End of observations text file reached')
        return success, epochflag, clockOffset, date, numSV,SVlist, eof

    epochflag   = line[28]
    # skip to next block if event flag is more than 1
    while int(epochflag) > 1:
        noFlag = 0
        linejump = 0
        pattern6 = r'\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+(\d+)\s+([A-Z\dG]+)'
        while line:
            match = re.search(pattern6, line)
            if match:
                break
            line = fid.readline()
            linejump += 1
        epochflag = int(line[28])
        msg = 'WARNING(rinexReadsObsBlockHead211): Observations event flag encountered. Flag = %s. %s lines were ignored.' % (str(epochflag), str(linejump))
        print(msg)

    # Gets the number of used satellites in obs epoch
    numSV = int(line[30:32])
    # Gets the receiver clock offset. This is optional data!
    # clockOffset = 0
    # if len(line) == 56:
    #     clockOffset = float(line[41:56])

    ## -- Reads the time stamp of the observations block (6 numerical values)
    date = [el for el in line[1::].split(" ") if el != ""]
    date = date[:6]
    date = [float(el) for el in date]

    if noFlag == 0:
        msg2 = msg + '\nEpoch date = %.4d %.2d %.2d %.2d:%.2d:%6.4f' % (date[0],date[1],date[2],date[3],date[4],date[5])
        print(msg2)


    # SVlist = []
    # sat_overview = line.split(" ")[-1][len(str(numSV))::]
    # pattern = re.compile(r'[A-Z][0-9]{2}')
    # sat_list = re.findall(pattern, sat_overview)
    # for sat in sat_list: ## droppe denne forloopen ?? SVlist og sat_list inneholder det samme??
    #     SVlist.append(sat)

    SVlist = re.findall(r'[A-Z][0-9]{2}', line)

    return success, epochflag, clockOffset, date, numSV, SVlist, eof
