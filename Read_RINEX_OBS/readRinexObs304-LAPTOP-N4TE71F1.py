def readRinexObs304(filename, readSS, readLLI, includeAllGNSSsystems,includeAllObsCodes, \
                    desiredGNSSsystems, desiredObsCodes, desiredObsBands):
    
    #Program/function to read GNSS observations in RINEX 3.04 observation files
    #The main core of the program is 4 functions:
    #                              rinexReadObsFileHeader304
    #                              rinexReadObsBlockHead304
    #                              rinexReadObsBlock304
    #                              rinexFindNEpochs304
    #
    #The first three of these functions are based on three functions
    #written by António Pestana:
    #                              rinexReadsObsFileHeader211
    #                              rinexReadsObsBlockHead211
    #                              rinexReadsObsBlock211
    #%
    #The structure of the data storage is also based on simular work by 
    #Ola Øvstedal
    #%--------------------------------------------------------------------------------------------------------------------------
    #INPUTS 
    
    #filename:                 path and name of RINEX 3.04 observation file,
    #                          string
    
    #readSS:                   Boolean, 0 or 1. 
    #                          1 = read "Signal Strength" Indicators
    #                          0 = do not read "Signal Strength" Indicators
    
    #readLLI:                  Boolean, 0 or 1. 
    #                          1 = read "Loss-Of-Lock Indicators"
    #                          0 = do not read "Loss-Of-Lock Indicators"
    
    #includeAllGNSSsystems:    Boolean, 0 or 1. 
    #                          1 = include alle GNSS systems(GPS, GLONASS, Galieo, BeiDou)
    #                          0 = include only GNSSsystems specified in desiredGNSSsystems
    
    #includeAllObsTypes:       Boolean, 0 or 1. 
    #                          1 = include all valid ObsTypes
    #                          0 = include only ObsTypes specified in desiredObsTypes
    
    #desiredGNSSsystems:       array og strings containing  codes of desired 
    #                          GNSSsystems to be included, 
    #                          ex. ["G", "E"]
    #                          OBS: Must be string array, NOT char vector
    
    #desiredObsTypes:          array of strings containing desired ObsTypes to be
    #                          included, ex. ["C", "L", "S", "D"]
    #                          OBS: Must be string array, NOT char vector
    
    #desiredObsBands:          array of desired obs Bands to be included,
    #                          ex [1, 5]
    #%--------------------------------------------------------------------------------------------------------------------------
    #OUTPUTS
    
    #GNSS_obs:                 cell containing a matrix for each GNSS system.
    #                          Each matrix is a 3D matrix containing all 
    #                          observation of current GNSS system for all epochs. 
    #                          Order of obsType index is same order as in 
    #                          obsCodes cell
    
    #                          GNSS_obs{GNSSsystemIndex}(PRN, obsType, epoch)
    #                                          GNSSsystemIndex: double,
    #                                          1,2,...,numGNSSsystems
    #                                          PRN: double
    #                                          ObsType: double: 1,2,...,numObsTypes
    #                                          epoch: double
    
    #GNSS_LLI:                 cell containing a matrix for each GNSS system
    #                          Each matrix stores loss of lock indicators for
    #                          each epoch for that GNSS system
    
    #GNSS_SS:                  cell containing a matrix for each GNSS system
    #                          Each matrix stores signal strength indicators 
    #                          for each epoch for that GNSS system
    
    #GNSS_SVs:                 cell containing a matrix for each GNSS system.
    #                          Each matrix contains number of satellites with 
    #                          obsevations for each epoch, and PRN for those 
    #                          satellites
    
    #                          GNSS_SVs{GNSSsystemIndex}(epoch, j)  
    #                                          j=1: number of observed satellites
    #                                          j>1: PRN of observed satellites
    
    #time_epochs:              matrix conatining gps-week and "time of week" 
    #                          for each epoch
    #                          time_epochs(epoch,i),   i=1: week
    #                                                  i=2: time-of-week in seconds (tow)
    
    #nepochs:                  number of epochs with observations in rinex observation file.
    
    #GNSSsystems:              cell array containing codes of GNSS systems included 
    #                          in RINEX observationfile. Elements are strings.
    #                          ex. "G" or "E"
        
    #obsCodes:                 Cell that defines the observation
    #                          codes available for all GNSS system. Each cell
    #                          element is another cell containing the codes for
    #                          that GNSS system. Each element in this cell is a
    #                          string with three-characters. The first 
    #                          character (a capital letter) is an observation code 
    #                          ex. "L" or "C". The second character (a digit)
    #                          is a frequency code. The third character(a Capital letter)  
    #                          is the attribute, ex. "P" or "X"
    
    #approxPosition:           array containing approximate position from rinex
    #                          observation file header. [X, Y, Z]
    
    #max_sat:                  array conataining max PRN number for each GNSS
    #                          system. Follows same order as GNSSsystems
    
    #tInterval:                observations interval; seconds. 
    
    #markerName:               name of the antenna marker; '' if not specified
    
    #rinexVersion:             string. rinex observation file version                  
    
    #recType:                  receiver type, char vector
    
    #timeSystem:               three-character code string of the time system 
    #                          used for expressing tfirstObs; 
    #                          can be GPS, GLO or GAL; 
    
    #leapSec:                  number of leap seconds since 6-Jan-1980. 
    #                          UTC=GPST-leapSec. NaN by default. 
    #                          THIS IS RINEX 3.04 OPTIONAL DATA
    
    #                          rinexHeader: cell column-vector containing the 
    #                          following data:
    #                          rinexVersion:   RINEX version number; string. 
    #                                          '' if not specified
    #                          rinexType:      RINEX file type; char
    
    #gnssType:                 GNSS system of the satellites observed; can be 
    #                          'G', 'R', 'E', 'C' or 'M' that stand for 
    #                          GPS, GLONASS, GALILEO, BeiDou or Mixed ; char
    
    #rinexProgr:               name of the software used to produce de RINEX 
    #                          GPS obs file; '' if not specified
    #
    
    #rinexDate:                date/time of the RINEX file creation; '' if not
    #                          specified; char
    
    
    #antDelta:                 column vector ot the three components of the 
    #                          distance from the marker to the antenna, 
    #                          in the following order - up, east and north;
    #                          null vector by default
    
    #tFirstObs:                time stamp of the first observation record in the RINEX
    #                          observations file; column vector
    #                          [YYYY; MM; DD; hh; mm; ss.sssssss];
    #                          THIS IS CRITICAL DATA
    
    #tLastObs:                 time stamp of the last observation record in the RINEX
    #                          observations file; column vector 
    #                          [YYYY; MM; DD; hh; mm;ss.sssssss]. NaN by default. 
    #                          THIS IS RINEX 3.04 OPTIONAL DATA
    
    #clockOffsetsON:           receiver clock offsets flag. O if no realtime-derived
    #                          receiver clock offset was applied to epoch, 
    #                          code and phase data (in other words, if the 
    #                          file only has raw data), 1 otherwise. 0 by default. 
    #                          THIS IS RINEX 3.04 OPTIONAL DATA
    
    #GLO_Slot2ChannelMap:      map container that maps GLONASS slot numbers to 
    #                          their respective channel number.
    #                          GLO_Slot2ChannelMap(slotnumber)      
    
    #success:                  Boolean. 1 if the reading of the RINEX 
    #                          observations file seems to be successful, 
    #                          0 otherwise
    #%--------------------------------------------------------------------------------------------------------------------------
    
    #ADVICE: The function rinexFindNEpochs() calculates the amount of observation epochs in
    #advance. This calculation will be incredibly more effective if TIME OF
    #LAST OBS is included in the header of the observation file. It is
    #strongly advized to manually add this information to the header if it is 
    #not included by default.  
    #%--------------------------------------------------------------------------------------------------------------------------
    #%
    #  According to RINEX 3.04 the observation type codes are:
    #%
    #  Observation code
    #  C: Pseudorange 
    #     GPS: C/A, L2C
    #     Glonass: C/A
    #     Galileo: All
    #  L: Carrier phase
    #  D: Doppler frequency
    #  S: Raw signal strengths or SNR values as given by the receiver for the
    #     respective phase observations 
    #  I: Ionosphere phase delay
    #  X: Receiver channel numbers
    #%
    #  Frequency code
    #  GPS Glonass Galileo SBAS
    #  1: L1 G1 E1 B1    (GPS,QZSS,SBAS,BDS)
    #  2: L2 G2 B1-2     (GLONASS)
    #  4: G1a            (Galileo)
    #  5: L5 E5a B2/B2a  (GPS, QZSS, SBAS, IRNSS) 
    #  6: L6 E6 B3 G2a   (Galileo, QZSS, BDS, GLONASS)
    #  7: E5b B2/B2b     (Galileo)
    #  8: E5a+b E5a+b    (Galileo, BDS)
    #  9: S              (IRNSS)
    #  0: for type X     (all)
    #%
    #  Attribute:
    #  A = A channel     (Galileo,IRNSS,GLONASS)
    #  B = B channel     (Galileo,IRNSS,GLONASS)
    #  C = C channel     (Galiloe, IRNSS)
    #      C code-based  (SBAS, GPS, GLONASS, QZSS)
    #  D = Semi-codelss  (GPS)
    #  
    #  I = I channel     (GPS, Galileo, QZSS, BDS)  
    #  L = L channel     (L2C GPS, QZSS)
    #      P channel     (GPS. QZSS)
    #  M = M code-based  (GPS)
    #  N = Codeless      (GPS) 
    #  P = P code-based  (GPS, GLONASS)
    #      Pilot channel (BDS)
    #  
    #  Q = Q channel     (GPS, Galileo, QZSS, BDS)
    #  S = D channel     (GPS, Galileo, QZSS, BDS)
    #      M channel     (L2C GPS, QZSS)
    #  
    #  W = Based on Z-tracking (GPS)
    #  X = B+C channels  (Galileo, IRNSS)
    #      I+Q channels  (GPS, IRNSS)
    #      M+L channels  (GPS, QZSS)
    #      D+P channels  (GPS, QZSS, BDS)
    #%
    #  Y = Y code based  (GPS)
    #  Z = A+B+C channels(Galileo)
    #      D+P channels  (BDS)
    #%--------------------------------------------------------------------------------------------------------------------------
    #%%
    import numpy as np
    from rinexFindNEpochs304 import rinexFindNEpochs304
    from rinexReadObsFileHeader304 import rinexReadObsFileHeader304
    from tqdm import tqdm
    from rinexReadObsBlock304 import rinexReadObsBlock304
    from kepler2ecef import date2gpstime
    from rinexReadObsBlockHead304 import rinexReadObsBlockHead304
    import time
    

    ## Get the start time
    t = time.process_time()

    
    ##Initialize variables in case of input error
    
    GNSS_obs       = np.nan;
    GNSS_LLI       = np.nan;
    GNSS_SS        = np.nan;
    GNSS_SVs       = np.nan;
    time_epochs    = np.nan;
    nepochs        = np.nan;
    GNSSsystems    = np.nan;
    obsCodes       = np.nan;
    approxPosition = np.nan;
    max_sat        = np.nan;
    tInterval      = np.nan;
    markerName     = np.nan;
    rinexVersion   = np.nan;
    recType        = np.nan;
    timeSystem     = np.nan;
    leapSec        = np.nan;
    gnssType       = np.nan;
    rinexProgr     = np.nan;
    rinexDate      = np.nan;
    antDelta       = np.nan;
    tFirstObs      = np.nan;
    tLastObs       = np.nan;
    clockOffsetsON = np.nan;
    GLO_Slot2ChannelMap = np.nan;
    
    ### --- Temporary defining input variables
    # readLLI = 0
    # readSS = 0
    # includeAllGNSSsystems = 1
    # includeAllObsTypes = 1
    # includeAllObsCodes = 1
    # desiredGNSSsystems = ["G","R","E","C"]
    # desiredObsCodes = []
    # desiredObsBands = [1,2, 5]
    # desiredObsTypes = ["C", "L", "S", "D"]
    #                          OBS: Must be string array, NOT char vector
    
    test = []
    GNSS_obs = {}
    GPS = {} # dict for storing GPS obs
    GLONASS = {}
    Galileo = {}
    BeiDou = {}
    E = {}
    T = {}
    ### ----------------
    # filename = 'opec0020_3.04_kort.10o'
    
    
    ## Test if readSS is boolean
    if readSS!=1 and readSS!=0:
        print('INPUT ERROR(readRinexObs304): The input argument readSS must be either 1 or 0')
        success = 0
        # return
    
    
    #Test if readLLI is boolean
    if readLLI!=1 and readLLI!=0:
        print('INPUT ERROR(readRinexObs304): The input argument readLLI must be either 1 or 0')
        success = 0
        # return
    
    max_GPS_PRN     = 36 #Max number of GPS PRN in constellation
    max_GLONASS_PRN = 36 #Max number of GLONASS PRN in constellation
    max_Galileo_PRN = 36 #Max number of Galileo PRN in constellation
    max_Beidou_PRN  = 60 #Max number of BeiDou PRN in constellation
    
    ## -- Read header of observation file
    [success, rinexVersion, gnssType, markerName, recType, antDelta,\
    GNSSsystems,numOfObsCodes, obsCodes, obsCodeIndex,tFirstObs, tLastObs, tInterval, \
    timeSystem, _, clockOffsetsON, rinexProgr, rinexDate,leapSec, approxPosition, GLO_Slot2ChannelMap, _, fid] = \
    rinexReadObsFileHeader304(filename, includeAllGNSSsystems, includeAllObsCodes,desiredGNSSsystems, desiredObsCodes, desiredObsBands)
    
    if success==0:
        pass 
       # return
    
    #Compute number of epochs with observations
    nepochs, tLastObs, tInterval, success =\
        rinexFindNEpochs304(filename, tFirstObs, tLastObs, tInterval) #computes number of epochs in observation file
    
    if success==0:
        pass
       # return     # HUSK Å FJERN KOMMENTERINGEN ETTER DU HAR LAGET FUNKSJON AV SKRIPTET!!!!!!!!
    
    ## --Number of GNSS systems
    nGNSSsystems = len(GNSSsystems)
    
    ## Declare data cells, arrays and matrices
    # GNSS_obs = np.zeros([nGNSSsystems,1])
    
    
    # GNSS_obs = {}   ###### CONNETA UT DETTE 11.11
    
    # GNSS_SVs =  np.zeros([nGNSSsystems,1])
    GNSS_SVs = {}
    max_sat  =  np.zeros([nGNSSsystems,1])
    # time_epochs =  np.zeros([nGNSSsystems,2])
    t_week = []
    t_tow = []
    
    # time_epochs =  np.array([])
    # GNSSsystems_full_names =  np.zeros([nGNSSsystems,1])
    GNSSsystems_full_names =  [""]*nGNSSsystems
    # GNSS_obs = cell(nGNSSsystems,1);
    # GNSS_SVs = cell(nGNSSsystems,1);
    # max_sat = zeros(nGNSSsystems,1);
    # time_epochs = zeros(nepochs, 2);
    # GNSSsystems_full_names = cell(nGNSSsystems,1);
    
    # if readLLI:
    #    GNSS_LLI  = cell(nGNSSsystems,1)
    # else:
    #    GNSS_LLI = np.nan;
    
    # if readSS:
    #    GNSS_SS  = cell(nGNSSsystems,1);
    # else:
    #    GNSS_SS = np.nan;
    
    
    
    ## -- Create array for max_sat. Initialize cell elements in cell arrays
    # for k = 1:nGNSSsystems
    for k in range(0,nGNSSsystems):
       if GNSSsystems[k+1] == 'G':
           max_sat[k] = max_GPS_PRN
           GNSS_SVs['G'] = np.zeros([nepochs,int(max_sat[k] + 1)])
           GNSSsystems_full_names[k] = "GPS"
       elif GNSSsystems[k+1] == 'R':
           max_sat[k] = max_GLONASS_PRN
           GNSS_SVs['R'] = np.zeros([nepochs,int(max_sat[k] + 1)])
           GNSSsystems_full_names[k] = "GLONASS"
        
       elif GNSSsystems[k+1] == 'E':
           max_sat[k] = max_Galileo_PRN
           GNSS_SVs['E'] = np.zeros([nepochs,max_sat[k] + 1])
           GNSSsystems_full_names[k] = "Galileo"
    
       elif GNSSsystems[k+1] == 'C':
           max_sat[k] = max_Beidou_PRN
           GNSS_SVs['C'] = np.zeros([nepochs,max_sat[k] + 1])
           GNSSsystems_full_names[k] = "BeiDou"
       else:
           print('ERROR(readRinexObs304): Only following GNSS systems are compatible with this program: GPS, GLONASS, Galileo, Beidou. %s is not valid' % GNSSsystems[k])
           # return
           
       
       # GNSS_obs[k] = np.zeros([max_sat[k], numOfObsCodes[k], nepochs])
       curr_sys = GNSSsystems[k+1]
       GNSS_obs[curr_sys] = np.zeros([int(max_sat[k]), numOfObsCodes[k], nepochs])
    
           
    
    
    # GNSS_names = containers.Map({'G', 'R', 'E', 'C'}, {'GPS', 'GLONASS', 'Galileo', 'Beidou'});
    GNSS_names = dict(zip(['G', 'R', 'E', 'C'],['GPS','GLONASS','Galileo','Beidou']))
    
    
    current_epoch      = 0
    
    #Initialize progress bar
    n_update_break = np.floor(nepochs/10); #number of epoch before updating progressbar
    # wbar = waitbar(0, 'INFO(readRinexObs304): Rinex observations are being read. Please wait');
    # waitbar = tqdm(0, 'INFO(readRinexObs304): Rinex observations are being read. Please wait')
    waitbar = tqdm(total=50)
    
    while 1:
       #" Read Obs Block Header
       success, _, _, date, numSV, eof = rinexReadObsBlockHead304(fid)
       
       if success==0 or eof==1:
          break
       
    
       ## -- Read current block of observations
       success, Obs,SVlist, numSV, LLI, SS, eof = rinexReadObsBlock304(fid, numSV, numOfObsCodes, GNSSsystems, obsCodeIndex, readSS, readLLI)
    
       if success ==0 or eof==1:
          break
       
       current_epoch = current_epoch + 1
    
       ## -- Update progress bar every n_update_break epochs
       if np.mod(current_epoch, n_update_break) == 0:
           # msg = print('INFO(readRinexObs304): Rinex observations are being read. Please wait\n %2.0f %%' % (current_epoch/nepochs*100))
           waitbar.update((current_epoch/nepochs)*10)
           #waitbar.close()
            
    
       ## Convert date to GPS-week and "time-of-week"
       week, tow = date2gpstime(int(date[0]), int(date[1]), int(date[2]), int(date[3]), int(date[4]), int(date[5]))
      
       ## Store GPS-week and "time-of-week" of current epoch
       t_week.append(int(week))
       t_tow.append(int(tow))
       time_epochs = np.column_stack((t_week,t_tow))
       ## Number of satellites with observations in this epoch, for each GNSS system
       nGNSS_sat_current_epoch = np.zeros([nGNSSsystems,1])
      
       ## Initialize dummy variables
       # GNSS_obs_dum = np.zeros([nGNSSsystems,1])
       # GNSS_LLI_dum = np.zeros([nGNSSsystems,1])
       # GNSS_SS_dum  = np.zeros([nGNSSsystems,1])
       GNSS_obs_dum = {}
       GNSS_LLI_dum = {}
       GNSS_SS_dum  = {}
       for k in range(0,nGNSSsystems):
           ## Initialize cell elements of dummy variables
           # GNSS_obs_dum{k} = zeros(max_sat(k), numOfObsCodes(k)); 
           # GNSS_LLI_dum{k} = zeros(max_sat(k), numOfObsCodes(k)); 
           # GNSS_SS_dum{k}  = zeros(max_sat(k), numOfObsCodes(k)); 
           GNSS_obs_dum[k+1] = np.zeros([int(max_sat[k]), numOfObsCodes[k]])
           GNSS_LLI_dum[k+1] = np.zeros([int(max_sat[k]), numOfObsCodes[k]]) 
           GNSS_SS_dum[k+1]  = np.zeros([int(max_sat[k]), numOfObsCodes[k]]) 
    
       ## -- Iterate through satellites of epoch and store obs, LLI and SS
       for sat in range(0,numSV):
           #Get index of current GNSS system
           # GNSSsystemIndex = find([GNSSsystems{:}] == SVlist{sat}(1));
           curr_sys = SVlist[sat][0]
           GNSSsystemIndex = [i for i in GNSSsystems if GNSSsystems[i]==curr_sys][0]
          
           #Increment amount of satellites this epoch for this GNSS system
           nGNSS_sat_current_epoch[GNSSsystemIndex-1] = nGNSS_sat_current_epoch[GNSSsystemIndex-1] + 1
          
           #Get just PRN number
           SV = int(SVlist[sat][1:3])
          
           # Number of obs types for current satellite
           nObsTypes_current_sat = numOfObsCodes[GNSSsystemIndex-1]
          
           #store observations, LLI, and SS of current satellite this epoch
    
           # GNSS_obs_dum{GNSSsystemIndex}(SV,1:nObsTypes_current_sat) = Obs(sat, 1:nObsTypes_current_sat);
          
           # GNSS_obs_dum[GNSSsystemIndex][SV][list(np.arange(1,nObsTypes_current_sat+2))] = Obs[sat, list(np.arange(1,nObsTypes_current_sat+2))]
           # GNSS_obs_dum[GNSSsystemIndex][SV][nObsTypes_current_sat-1] = Obs[sat,nObsTypes_current_sat-1]
           GNSS_obs_dum[GNSSsystemIndex][SV][0:nObsTypes_current_sat-1] = Obs[sat,0:nObsTypes_current_sat-1]
    
    
          
           if readLLI:
              # GNSS_LLI_dum{GNSSsystemIndex}(SV,1:nObsTypes_current_sat) = LLI(sat, 1:nObsTypes_current_sat);
              GNSS_LLI_dum[GNSSsystemIndex][SV][nObsTypes_current_sat] = LLI[sat, nObsTypes_current_sat]
         
           if readSS:
              # GNSS_SS_dum{GNSSsystemIndex}(SV,1:nObsTypes_current_sat)  = SS(sat, 1:nObsTypes_current_sat)
              GNSS_SS_dum[GNSSsystemIndex][SV][nObsTypes_current_sat] = SS[sat, nObsTypes_current_sat]
          
    
           #Store PRN number of current sat to PRNs of this epoch
           # GNSS_SVs{GNSSsystemIndex}(current_epoch, nGNSS_sat_current_epoch(GNSSsystemIndex) +1) = SV; 
           # GNSS_SVs[curr_sys][current_epoch, int(nGNSS_sat_current_epoch[GNSSsystemIndex-1]) +1] = SV 
           GNSS_SVs[curr_sys][current_epoch-1, int(nGNSS_sat_current_epoch[GNSSsystemIndex-1])] = SV 
    
       
       # E = {}
       # T = {}
    
       for k in range(0,nGNSSsystems):
           curr_sys = GNSSsystems[k+1]
           #Set number of satellites with obs for each GNSS system this epoch
           # GNSS_SVs{k}(current_epoch, 1)  = nGNSS_sat_current_epoch(k) ;
           # GNSS_SVs[curr_sys][current_epoch, 1]  = nGNSS_sat_current_epoch[k]
           GNSS_SVs[curr_sys][current_epoch-1, 0]  = nGNSS_sat_current_epoch[k]
           #store dummy matrices of this epoch
           # GNSS_obs{k}(:,:,current_epoch) = GNSS_obs_dum{k}  
           # GNSS_obs[curr_sys][:,:,current_epoch] = GNSS_obs_dum[k+1]  
                 
           # GNSS_obs[curr_sys][current_epoch][:,:] = GNSS_obs_dum[k+1]   
           
           # OBS = {curr_sys: {current_epoch-1 : int(current_epoch)}}
           # test = {curr_sys: {current_epoch-1 : GNSS_obs_dum[k+1]}}
           
           # E =  GNSS_obs_dum[k+1]
           # test.append(GNSS_obs_dum[k+1])
           if curr_sys == 'G':
               GPS[current_epoch] = GNSS_obs_dum[k+1]
           elif curr_sys == 'R':
               GLONASS[current_epoch] = GNSS_obs_dum[k+1]
           elif curr_sys == 'E':
               Galileo[current_epoch] = GNSS_obs_dum[k+1]
           elif curr_sys == 'C':
               BeiDou[current_epoch] = GNSS_obs_dum[k+1]
               
           
           # E[current_epoch] =  GNSS_obs_dum[k+1]
           # test.append(GNSS_obs_dum[k+1])
           
           # if readLLI:
           #    GNSS_LLI{k}(:,:,current_epoch) = GNSS_LLI_dum{k}            
           # if readSS:
           #    GNSS_SS{k}(:,:,current_epoch)  = GNSS_SS_dum{k}          
           # T[curr_sys] = E
           # GNSS_obs[curr_sys] = E
           
    GNSS_obs['G'] = GPS 
    GNSS_obs['R'] = GLONASS
    GNSS_obs['E'] = Galileo
    GNSS_obs['C'] = BeiDou
    
    del_sys = list(GNSS_obs.keys())
    for sys in del_sys: # Deleting systems with no observations
        if not GNSS_obs[sys]:
            del GNSS_obs[sys]
    
    if current_epoch!= nepochs and success == 1:
        print('ERROR(readRinexObs304): The amount of epochs calculated in \
              advance(nepochs = %d) does not equal number og epochs prossesed(current_epoch = %d). \n\tCheck that header information concerning TIME OF FIRST OBS and TIME OF LAST OBS is correct.\n' %(nepochs, current_epoch))
         
    
    
    # close(wbar)
    # waitbar.close()
    
    if success == 1:
       # messages = cell(nGNSSsystems+1);
       messages = []
       messages.append('INFO(readRinexObs304): The following GNSS systems have been read into the data:')
       # for k = 1:nGNSSsystems
       for k in range(0,nGNSSsystems):
          # messages{k+1} = sprintf('INFO(readRinexObs304): The following %s observation types have been registered:', GNSS_names(GNSSsystems{k}));
           messages.append(('INFO(readRinexObs304): The following %s observation types have been registered:' % (GNSS_names[GNSSsystems[k+1]])))
           curr_sys = GNSSsystems[k+1]
          # for obs =1:length(obsCodes{k})
           for obs in range(0, len(obsCodes[k+1])):
               if obs == 1:
                   # messages{k+1} = append(messages{k+1}, sprintf(' %s', obsCodes{k}{obs}));
                   # messages[k+1] = messages[k+1] + " " + print(' %s' % (obsCodes[k][obs]))
                   messages.append((' %s' % (obsCodes[k+1][curr_sys])))
               else:
                   # messages{k+1} = append(messages{k+1}, sprintf(', %s', obsCodes{k}{obs}));
                    messages.append((', %s', obsCodes[k+1][curr_sys]))
            
           if k == 1:
               # messages{1} = append(messages{1}, sprintf(' %s', GNSS_names(GNSSsystems{k})));
               messages.append(' %s' % GNSS_names[GNSSsystems[k+1]])
           else:
               # messages{1} = append(messages{1}, sprintf(', %s', GNSS_names(GNSSsystems{k})));
               messages.append((', %s' % GNSS_names[GNSSsystems[k+1]]))
         
       
       # for msg = 1:length(messages)
       for msg in range(0,len(messages)):
           print(messages[msg])
    
    if readLLI:
       print('INFO(readRinexObs304): LLI have been read (if present in observation file)')
    else:
       print('INFO(readRinexObs304): LLI have not been read')
    
    
    if readSS:
       print('INFO(readRinexObs304): SS have been read (if present in observation file)')
    else:
       print('INFO(readRinexObs304): SS have not been read')
    
    
    # get the end time
    et = time.process_time()
    
    # get execution time
    e = et - t
    
    # e = cputime-t;
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