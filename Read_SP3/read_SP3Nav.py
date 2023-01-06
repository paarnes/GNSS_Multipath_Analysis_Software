def readSP3Nav(filename, desiredGNSSsystems=None):
        
   import numpy as np
   import copy

    #Function that reads the GNSS satellite position data from a SP3 position
    #file. The function has been tested with sp3c and sp3d. NOTE: It is
    #advised that any use of this function is made through the parent function
    #"read_multiple_SP3Nav.m", as it has more functionality. 
    #--------------------------------------------------------------------------------------------------------------------------
    #INPUTS
    
    #filename:             path and filename of sp3 position file, string
    
    #desiredGNSSsystems:   array. Contains string. Each string is a code for a
    #                      GNSS system that should have its position data stored 
    #                      in sat_positions. Must be one of: "G", "R", "E",
    #                      "C". If left undefined, it is automatically set to
    #                      ["G", "R", "E", "C"]
    #--------------------------------------------------------------------------------------------------------------------------
    #OUTPUTS
    
    #sat_positions:    cell. Each cell elements contains position data for a
    #                  specific GNSS system. Order is defined by order of 
    #                  navGNSSsystems. Each cell element is another cell that 
    #                  stores position data of specific satellites of that 
    #                  GNSS system. Each of these cell elements is a matrix 
    #                  with [X, Y, Z] position of a epoch in each row.
    
    #                  sat_positions{GNSSsystemIndex}{PRN}(epoch, :) = [X, Y, Z]
    
    #epoch_dates:      matrix. Each row contains date of one of the epochs. 
    #                  [nEpochs x 6]
    
    #navGNSSsystems:   array. Contains string. Each string is a code for a
    #                  GNSS system with position data stored in sat_positions.
    #                  Must be one of: "G", "R", "E", "C"
    
    #nEpochs:          number of position epochs, integer
    
    #epochInterval:    interval of position epochs, seconds
    
    #success:          boolean, 1 if no error occurs, 0 otherwise
    #--------------------------------------------------------------------------------------------------------------------------

   #Function that reads the GNSS satellite position data from a SP3 position
   #file. The function has been tested with sp3c and sp3d. NOTE: It is
   #advised that any use of this function is made through the parent function
   #"read_multiple_SP3Nav.m", as it has more functionality. 
   #--------------------------------------------------------------------------------------------------------------------------
   #INPUTS

   #filename:             path and filename of sp3 position file, string

   #desiredGNSSsystems:   array. Contains string. Each string is a code for a
   #                      GNSS system that should have its position data stored 
   #                      in sat_positions. Must be one of: "G", "R", "E",
   #                      "C". If left undefined, it is automatically set to
   #                      ["G", "R", "E", "C"]
   #--------------------------------------------------------------------------------------------------------------------------
   #OUTPUTS

   #sat_positions:    cell. Each cell elements contains position data for a
   #                  specific GNSS system. Order is defined by order of 
   #                  navGNSSsystems. Each cell element is another cell that 
   #                  stores position data of specific satellites of that 
   #                  GNSS system. Each of these cell elements is a matrix 
   #                  with [X, Y, Z] position of a epoch in each row.

   #                  sat_positions{GNSSsystemIndex}{PRN}(epoch, :) = [X, Y, Z]

   #epoch_dates:      matrix. Each row contains date of one of the epochs. 
   #                  [nEpochs x 6]

   #navGNSSsystems:   array. Contains string. Each string is a code for a
   #                  GNSS system with position data stored in sat_positions.
   #                  Must be one of: "G", "R", "E", "C"

   #nEpochs:          number of position epochs, integer

   #epochInterval:    interval of position epochs, seconds

   #success:          boolean, 1 if no error occurs, 0 otherwise
   #--------------------------------------------------------------------------------------------------------------------------



   #
   max_GNSSsystems = 4

   max_GPS_PRN     = 36 #Max number of GPS PRN in constellation
   max_GLONASS_PRN = 36 #Max number of GLONASS PRN in constellation
   max_Galileo_PRN = 36 #Max number of Galileo PRN in constellation
   max_Beidou_PRN  = 60 #Max number of BeiDou PRN in constellation
   max_sat = [max_GPS_PRN, max_GLONASS_PRN, max_Galileo_PRN, max_Beidou_PRN]

   #Initialize variables
   success = 1

   ## --- Open nav file

   try:
       fid = open(filename,'r')
   except:
       success = 0
       raise ValueError('No file selected!')
       
   if desiredGNSSsystems is None:
       desiredGNSSsystems = ["G", "R", "E", "C"]

   #GNSS system order
   navGNSSsystems = ["G", "R", "E", "C"];
   #Map mapping GNSS system code to GNSS system index
   # GNSSsystem_map = containers.Map(navGNSSsystems, [1, 2, 3, 4]);
   GNSSsystem_map = dict(zip(navGNSSsystems,[1, 2, 3, 4]))

   sat_pos = {}

   # Read header
   headerLine = 0
   line = fid.readline().rstrip()

   # All header lines begin with '*'
   while '*' not in line[0]:
       headerLine = headerLine + 1
       
       if headerLine == 1:
          sp3Version = line[0:2]
          
          # Control sp3 version    
          if '#c' not in sp3Version and '#d' not in sp3Version:
              print('ERROR(readSP3Nav): SP3 Navigation file is version %s, must be version c or d!' % (sp3Version))
              # [sat_positions, epoch_dates, navGNSSsystems, nEpochs, epochInterval] = deal(NaN)
              success = 0
              # return
        
          
          # Control that sp3 file is a position file and not a velocity file
          Pos_Vel_Flag = line[2]

          if 'P' not in Pos_Vel_Flag:
              print('ERROR(readSP3Nav): SP3 Navigation file is has velocity flag, should have position flag!')
              # [sat_positions, epoch_dates, navGNSSsystems, nEpochs, epochInterval] = deal(NaN);
              success = 0
              # return
          
          #Store coordinate system and amount of epochs
          CoordSys = line[46:51]
          nEpochs = int(line[32:39])
       
       
       if headerLine == 2:
          # Store GPS-week, "time-of-week" and epoch interval[seconds]
          GPS_Week = int(line[3:7])
          tow      = float(line[8:23])
          epochInterval = float(line[24:38])

       
       if headerLine == 3:
          
          #initialize list for storing indices of satellites to be excluded
          RemovedSatIndex = []
          
          
          if '#c' in sp3Version:
           nSat = int(line[4:6])
          else:
           nSat = int(line[3:6])
          
          ## -- Remove beginning of line
          line = line[9:60]
          
          #Initialize array for storing the order of satellites in the SP3
          #file(ie. what PRN and GNSS system index)
          
          GNSSsystemIndexOrder = []
          PRNOrder = []
           
          
          # Keep reading lines until all satellite IDs have been read
          for k in range(0,nSat):          
             # Control that current satellite is amoung desired systems
             if np.in1d(line[0], desiredGNSSsystems):
                
                 ## -- Get GNSSsystemIndex from map container
                 GNSSsystemIndex = GNSSsystem_map[line[0]]
                 #Get PRN number/slot number
                 PRN = int(line[1:3])
                 #remove satellite that has been read from line
                 line = line[3::]
                 #Store GNSSsystemIndex and PRN in satellite order arrays
                 GNSSsystemIndexOrder.append(GNSSsystemIndex)
                 PRNOrder.append(PRN)
                 
                 #if current satellite ID was last of a line, read next line
                 #and increment number of headerlines
                 if np.mod(k+1,17)==0 and k != 0: # FEILEN ER HER. DERFOR BLIR IKKE PR5 satellitten med
                     line = fid.readline().rstrip()
                     line = line[9:60]
                     headerLine = headerLine + 1
             #If current satellite ID is not amoung desired GNSS systems,
             #append its index to array of undesired satellites
             else:
                 RemovedSatIndex.append(k)
                 GNSSsystemIndexOrder.append(np.nan)
                 PRNOrder.append(np.nan)
                 
                 #if current satellite ID was last of a line, read next line
                 #and increment number of headerlines
                 if np.mod(k+1,17)==0 and k != 0: 
                     line = fid.readline().rstrip()
                     line = line[9:60]
                     headerLine = headerLine + 1;
       # Read next line
       line = fid.readline().rstrip()

   # Initialize matrix for epoch dates
   epoch_dates = []
   
   sys_dict = {}
   PRN_dict = {}

   PRN_dict_GPS = {}
   PRN_dict_Glonass = {}
   PRN_dict_Galileo = {}
   PRN_dict_BeiDou = {}
   
   # Read satellite positions of every epoch    
   ini_sys = list(GNSSsystem_map.keys())[0]
   for k in range(0,nEpochs):
       #Store date of current epoch
       epochs = line[3:31].split(" ")
       epochs = [x for x in epochs if x != "" ] # removing ''
       epoch_dates.append(epochs)
       
       # Store positions of all satellites for current epoch
       obs_dict = {}
       obs_dict_GPS = {}
       obs_dict_Glonass = {}
       obs_dict_Galileo = {}
       obs_dict_BeiDou = {}
       
       for i in range(0,nSat):
           #read next line
           line = fid.readline().rstrip()
           #if current satellite is amoung desired systems, store positions
           if np.in1d(i, RemovedSatIndex,invert = True):
               #Get PRN and GNSSsystemIndex of current satellite for
               #previously stored order
               PRN = PRNOrder[i]
               GNSSsystemIndex = GNSSsystemIndexOrder[i]
               # Store position of current satellite in correct location in
               sys_keys = list(GNSSsystem_map.keys())
               sys_values = list(GNSSsystem_map.values())
               sys_inx = sys_values.index(GNSSsystemIndex)
               sys = sys_keys[sys_inx]
               obs = line[5:46].split(" ")
               obs = [x for x in obs if x != "" ]

                   
               if sys != ini_sys:
                   ini_sys = sys

               obs_dict[str(PRN)]  = obs[:]  
               PRN_dict[int(k)] = obs_dict
               
               if sys == 'G':
                   obs_G = [x for x in obs if x != "" ]
                   obs_dict_GPS[PRN]  = obs_G.copy()
                   PRN_dict_GPS[k] = obs_dict.copy()
               elif sys =='R':
                   obs_R = [x for x in obs if x != "" ]
                   obs_dict_Glonass[PRN]  = obs_R.copy()
                   PRN_dict_Glonass[k] = obs_dict_Glonass.copy()
               elif sys =='E':
                   obs_E = [x for x in obs if x != "" ]
                   obs_dict_Galileo[PRN]  = obs_E.copy()
                   PRN_dict_Galileo[k] = obs_dict_Galileo.copy()
               elif sys =='C':
                   obs_C = [x for x in obs if x != "" ]
                   obs_dict_BeiDou[PRN]  = obs_C.copy()
                   PRN_dict_BeiDou[k] = obs_dict_BeiDou.copy()
                   
                   
                   
           sys_dict['G'] = PRN_dict_GPS
           sys_dict['R'] = PRN_dict_Glonass
           sys_dict['E'] = PRN_dict_Galileo
           sys_dict['C'] = PRN_dict_BeiDou
           sat_pos[sys] = sys_dict[sys]
                   
       #Get next line
       line = fid.readline().rstrip()
       
   #the next line should be eof. If not, raise warning
   try:
       line = fid.readline().rstrip()
   except:
       print('ERROR(readSP3Nav): End of file was not reached when expected!!')
       success = 0
        
   #remove NaN values
   GNSSsystemIndexOrder = [x for x in GNSSsystemIndexOrder if x != 'nan']
   PRNOrder = [x for x in GNSSsystemIndexOrder if x != 'nan']

   print('SP3 Navigation file "%s" has been read successfully.' %(filename))
   ## Remove GNSS systems not present in navigation file
   # sat_pos['GNSS_systems']  = sat_pos['GNSS_systems'][np.unique(GNSSsystemIndexOrder)]
   # navGNSSsystems = navGNSSsystems(np.unique(GNSSsystemIndexOrder))
   return sat_pos, epoch_dates, navGNSSsystems, nEpochs, epochInterval


# sat_pos, epoch_dates, navGNSSsystems, nEpochs, epochInterval =readSP3Nav('testfile4.SP3')