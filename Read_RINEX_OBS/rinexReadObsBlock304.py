def rinexReadObsBlock304(fid, numSV, nObsCodes, GNSSsystems, obsCodeIndex, readSS, readLLI):
    """
    #Reads all the observations from a RINEX observation block.
    
    #Positioned at the beginning of the line immediately after the header of the
    #observations block, reads all the observations in this block of a RINEX
    #observations file. This function is meant to be used after using function
    #rinexReadObsFileHeader304
    
    #Based in the work of Antonio Pestana, rinexReadObsBlock211, March 2015
    #--------------------------------------------------------------------------------------------------------------------------
    #INPUTS
    
    #fid:                  Matlab file identifier of a Rinex observations text file
    
    #numSV:                total number of satellites with observations in
    #                      current observation block, integer
    
    #numOfObsCodes:        column vector containing number of observation
    #                      types for each GNSS system. Order is the same as
    #                      GNSSsystems
    
    #GNSSsystems:          cell array containing codes of GNSS systems included 
    #                      in RINEX observationfile. Elements are strings.
    #                      ex. "G" or "E"
    
    #obsCodeIndex:         cell with one cell element for each GNSS system.
    #                      Order is the same as GNSSsystems. Each cell element
    #                      contains an array of indices. These indices
    #                      indicate the observation types that should be
    #                      read for each GNSS system. ex. If one index for
    #                      GPS is 1 then the first observation type for GPS
    #                      should be read.
    
    #readSS:                   Boolean, 0 or 1. 
    #                          1 = read "Signal Strength" Indicators
    #                          0 = do not read "Signal Strength" Indicators
    
    #readLLI:                  Boolean, 0 or 1. 
    #                          1 = read "Loss-Of-Lock Indicators"
    #                          0 = do not read "Loss-Of-Lock Indicators"
    #--------------------------------------------------------------------------------------------------------------------------
    #OUTPUTS
    
    #success:               Boolean. 1 if the function seems to be successful, 
    #                      0 otherwise
    
    #Obs:                  matrix [numSV x max_nObs] that stores all 
    #                      observations of this observation block. max_nObs 
    #                      is the highest number of observation codes that 
    #                      any of the GNSS systems have. Which observation
    #                      types that are associated with what collumn will
    #                      vary between GNSS systems. SVlist will give
    #                      overview of what GNSS system each row is connected
    #                      to
    
    #SVlist:               column cell [numSV x 1] that conatins the 
    #                      identification code of each line of observation 
    #                      block. ex. "G21". numSV is total number of 
    #                      satellites minus amount of satellites removed.
    
    #numSV:                numSV, unlike the input of same name, is the total 
    #                      number of satellites minus amount of satellites 
    #                      removed.
    
    #LLI:                  matrix [numSV x max_nObs] that stores all 
    #                      "loss-of-lock" indicators of this observation block. 
    #                      max_nObs is the highest number of observation codes 
    #                      that any of the GNSS systems have. Which observation
    #                      types that are associated with what collumn will
    #                      vary between GNSS systems. SVlist will give
    #                      overview of what GNSS system each row is connected
    #                      to
    
    #SS:                   matrix [numSV x max_nObs] that stores all 
    #                      "signal strength" indicators of this observation block. 
    #                      max_nObs is the highest number of observation codes 
    #                      that any of the GNSS systems have. Which observation
    #                      types that are associated with what collumn will
    #                      vary between GNSS systems. SVlist will give
    #                      overview of what GNSS system each row is connected
    #                      to
    
    #eof:                  end-of-file flag; 1 if end-of-file was reached, 
    #                      0 otherwise
    #--------------------------------------------------------------------------------------------------------------------------
    """
    import numpy as np
    import os
    
    
    ## Initialize variables in case of input error
    success                     = np.nan;
    eof                         = np.nan;
    max_n_obs_Types             = np.nan;
    Obs                         = np.nan;
    LLI                         = np.nan;
    SS                          = np.nan;
    SVlist                      = np.nan;
    removed_sat                 = np.nan;
    desiredGNSSsystems          = np.nan;
    
    ## Testing input arguments
    
    ## Test type of GNSS systems
    # if ~isa(GNSSsystems, 'cell')
    #    sprintf(['INPUT ERROR(rinexReadObsBlock304): The input argument GNSSsystems',...
    #         ' is of type #s.\n Must be of type cell'],class(GNSSsystems))
    #     success = 0;
    #     return;  
    # end
    
    ## Test type of numSV
    if type(numSV) != int: 
       print('INPUT ERROR(rinexReadObsBlock304): The input argument numSV is of type %s.\n Must be of type double' % (type(numSV)))
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
    
    
    ## Test type of obsCodeIndex
    # if ~isa(obsCodeIndex, 'cell')
    #    sprintf(['INPUT ERROR(rinexReadObsBlock304): The input argument numOfObsTypes',...
    #         ' is of type #s.\n Must be of type cell'],class(obsCodeIndex))
    #     success = 0;
    #     return;  
    # end
    
    ## Test type of cell elements of removed ObsCodeIndex
    # for k = 1:length(GNSSsystems)
    #     if ~isa(obsCodeIndex{k}, 'double')
    #         sprintf(['INPUT ERROR(rinexReadObsBlock304): Element with index #d in',...
    #             ' removedObsTypesIndex is not of type double as it should be. \nVariable has type #s'],...
    #             k, class(obsCodeIndex{k}))
    #         success = 0;
    #         return;
    #     end
    # end
    
    ##
    
    success = 1;
    eof     = 0;
    
    # Highest number of obs codes of any GNSS system
    max_n_obs_Types = max(nObsCodes)
    
    # Initialize variables
    # Obs = zeros(numSV, max_n_obs_Types) ;
    Obs = np.empty([numSV, max_n_obs_Types]) 
    # SVlist = cell(numSV,1);
    # SVlist = dict(numSV)
    # SVlist = np.arange(numSV).reshape(numSV,1)   
    # SVlist.fill(numSV)
    # SVlist = np.chararray(numSV).reshape(numSV,1)   
    # SVlist.fill(numSV)
    SVlist = [np.nan]*numSV
    if readLLI:
       # LLI = zeros(numSV, max_n_obs_Types) ;
       LLI = np.empty([numSV, max_n_obs_Types]) 
    
    if readSS:
       # SS  = zeros(numSV, max_n_obs_Types) ;
       SS  = np.empty([numSV, max_n_obs_Types]) 
    
                  
    
    # number of satellites excluded so far
    removed_sat = 0                                 
    # desiredGNSSsystems = str(GNSSsystems)      # DETTE MÅ TROLIG ENDRES!!   
    desiredGNSSsystems = list(GNSSsystems.values())
    # Gobble up observation block
    # for sat = 1:numSV
    for sat in range(0,numSV):
        
       line = fid.readline().rstrip()   
    
           
      # if line == -1
      #   eof = 1;
      #   disp(['ERROR (rinexReadsObsBlock211): the end of the '...
      #         'observations text file was reached unexpectedly'])
      #   success = 0;
      #   return
      # end
    
      # SV = line(1:3); # Satellite code, ex. 'G11' or 'E03'
      
       if not line:
           return
      
       SV = line[0:3].strip() # Satellite code, ex. 'G11' or 'E03'
       if SV[0] not in desiredGNSSsystems:
           removed_sat +=1
       else:
           ## Index of current GNSS system
           # GNSSsystemIndex = find([GNSSsystems{:}] == SV(1)); 
            # GNSSsystemIndex = {i for i in GNSSsystems if GNSSsystems[i]==SV[0]} # set virkelig?
           GNSSsystemIndex = [i for i in GNSSsystems if GNSSsystems[i]==SV[0]][0] 
    #      SVlist{sat - removed_sat,1} = SV; # Store SV of current row
           SVlist[sat - removed_sat] = SV # Store SV of current row
           
    #        n_obs_current_system = nObsCodes(GNSSsystemIndex);
           n_obs_current_system = nObsCodes[GNSSsystemIndex-1]
           
            # for obs_num = 1:n_obs_current_system
           for obs_num in range(0, n_obs_current_system):
    #           obsIndex = obsCodeIndex{GNSSsystemIndex}(obs_num);
               obsIndex = obsCodeIndex[GNSSsystemIndex][obs_num]
               # charPos = 4+(obsIndex-1)*16;
               charPos = 4+(obsIndex)*16
               ## check that the current observation of the current GNSS system
    #          ## is not on the list of obs types to be excluded
             
               ## stringlength of next obs. 
               # obsLen = min(14, length(line) - charPos); 
               obsLen = min(14, len(line) - charPos)
               # read next obs
               # newObs = strtrim(line(charPos:charPos+obsLen)); 
               newObs = line[charPos:charPos+obsLen].strip()
               # If observation missing, set to 0
               # if ~isempty(newObs):
                   
               if newObs != '':
                   newObs = float(newObs)
               else:
                   newObs = 0;
               # Store new obs
               Obs[sat - removed_sat, obs_num] = newObs 
                
               if readLLI:
               # read LLI of current obs (if present)
                   if charPos+13<len(line): ## endret til < (kun mindre enn) 13.11
                       newLLI = line[charPos+13] # loss of lock indicator ### endret fra 14 til 13 den 13.11.2022
                   else:
                       newLLI = ' '
                
    
                # if no LLI set to -999
                   # if isspace(newLLI):
                   #     newLLI = -999
                   # else:
                   #     newLLI = int(newLLI)
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
                       newSS = ' ';
                    
                
                    # if no SS set to -999
                   if newSS.isspace():
                       newSS = -999;
                   else:
                       newSS = int(newSS)
    
                    #Store SS
                   SS[sat - removed_sat, obs_num]  = newSS         
     
    
    
    # # update number og satellites after satellites have been excluded
    numSV = numSV - removed_sat
    
    #remove empty cell elements
    # SVlist(cellfun('isempty',SVlist))   = [];
    SVlist = list(filter(None,SVlist))
    # Obs[end-removed_sat+1:end, :]       = [];
    idx_keep = len(Obs) -1 -removed_sat + 1 # removing sats
    Obs = Obs[:idx_keep,:]
    # Obs[-1-removed_sat+1::, :]       = [];
    return success, Obs,SVlist, numSV, LLI, SS, eof

# end
