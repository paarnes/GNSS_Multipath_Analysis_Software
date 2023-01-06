from rinexReadObsFileHeader304 import rinexReadObsFileHeader304
from rinexReadObsBlockHead304 import rinexReadObsBlockHead304
from rinexReadObsBlock304 import rinexReadObsBlock304
from kepler2ecef import date2gpstime
import numpy as np
import os

filename = 'opec0020_3.04_kort.10o'

includeAllGNSSsystems = 1
includeAllObsCodes = 1
desiredGNSSsystems= ["G","R", "E", "C"]
desiredObsCodes = ["C", "L", "S", "D"]
desiredObsBands = [1,2,5]


success, rinexVersion, gnssType, markerName, recType, antDelta,GNSSsystems,nObsCodes, \
obsCodes, obsCodeIndex,tFirstObs, tLastObs, tInterval,timeSystem, numHeaderLines, clockOffsetsON, \
rinexProgr, rinexDate,leapSec, approxPosition, GLO_Slot2ChannelMap, eof, fid = rinexReadObsFileHeader304(filename, includeAllGNSSsystems, includeAllObsCodes,desiredGNSSsystems,
                                  desiredObsCodes, desiredObsBands)


success, epochflag, clockOffset, date, numSV, eof= rinexReadObsBlockHead304(fid)




## ----
#--------------------------------------------------------------------------------------------------------------------------

## -- Imports 
# % readSS:                   Boolean, 0 or 1. 
# %                           1 = read "Signal Strength" Indicators
# %                           0 = do not read "Signal Strength" Indicators

# % readLLI:                  Boolean, 0 or 1. 
# %                           1 = read "Loss-Of-Lock Indicators"
# %                           0 = do not read "Loss-Of-Lock Indicators"




# filename = 'opec0020_3.04_kort.10o'
# fid = open(filename,'r') 
readSS = 0
readLLI =0 

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





success, Obs,SVlist, numSV, LLI, SS, eof = rinexReadObsBlock304(fid, numSV, nObsCodes, GNSSsystems, obsCodeIndex, readSS, readLLI)




rinexFindNEpochs304(filename, tFirstObs, tLastObs, tInterval)

# --------------------------------------------------------------------------------------------------------------------------
# #  Testing input arguments
success = 1
nepochs = 0

# #  Test if filename is valid format
# if ~isa(filename,'char')&&~isa(filename,'string')
#     sprintf(['INPUT ERROR(rinexFindNEpoch): The input argument filename',...
#         'is of type # s.\n Must be of type string or char'],class(filename))
#     success = 0;
#     return
# end

# #  Test if filename is valid format
# if ~isa(tFirstObs,'double')||length(tFirstObs)~=6
#     sprintf(['INPUT ERROR(rinexFindNEpoch): The input argument tFirstObs',...
#         ' is of type # s and has length # d.\n Must be of type double with length 6'],...
#         class(tFirstObs), length(tFirstObs))
#     success = 0;
#     return
# end

# #  Test if filename is valid format
# if (~isa(tLastObs,'double') && length(tLastObs)~=6) && ~any(isnan(tLastObs))
#     sprintf(['INPUT ERROR(rinexFindNEpoch): The input argument tLastObs',...
#         ' is of type # s and has length # d.\n Must be of type double with length 6, or Nan'],...
#         class(tLastObs), length(tLastObs))
#     success = 0;
#     return
# end
    
# filename = 'opec0020_3.04_kort.10o'

# #  open observation file
# fid = open(filename, 'rt')
# seconds_in_a_week = 604800
# #  tLastObs is in header
# #%%
# # if ~isnan(tLastObs):
# if not True in np.isnan(tLastObs):
#     #  tInterval is not in header
#    if np.isnan(tInterval):
#        tInterval_found = 0
#        first_epoch_found = 0
#        #  calculate tInterval
#        while not tInterval_found: #  calculate tInterval
#            line = fid.readline().rstrip()
#            #  start of new epoch
#            # answer = strfind(line,'>');
#            if '>' in line:
#            # if ~isempty(answer) 
#                if not first_epoch_found: #  first epoch
#                    # first_epoch_time = cell2mat(textscan(line(2:end),'# f', 6));
#                    first_epoch_time = line[1::]
#                    first_epoch_time = [float(el) for el in line[1::].split(" ") if el != ""]
#                    first_epoch_time = first_epoch_time[:6]
#                    first_epoch_found = 1;
#                else: #  seconds epoch
#                    # second_epoch_time = cell2mat(textscan(line(2:end),'# f', 6));
#                    second_epoch_time = line[1::]
#                    second_epoch_time = [float(el) for el in line[1::].split(" ") if el != ""]
#                    second_epoch_time = second_epoch_time[:6]
                
#                    # tInterval = second_epoch_time(6)-first_epoch_time(6);
#                    tInterval = second_epoch_time[5]-first_epoch_time[5]
#                    tInterval_found = 1;

#    tFirstObs = tFirstObs.astype(int) 
#    tLastObs = tLastObs.astype(int) 
#    week1,tow1 = date2gpstime(tFirstObs[0][0],tFirstObs[1][0],tFirstObs[2][0],tFirstObs[3][0],tFirstObs[4][0],tFirstObs[5][0])
#    week2,tow2 = date2gpstime(tLastObs[0][0],tLastObs[1][0],tLastObs[2][0],tLastObs[3][0],tLastObs[4][0],tLastObs[5][0])
#    n = ((week2*seconds_in_a_week+tow2) - (week1*seconds_in_a_week+tow1)) / tInterval
#    nepochs = np.floor(n) + 1;  

# #  if tLastObs is not in header. Function counts number of epochs manually 
# #  OBS: This can be VERY inefficent
# else:
#    print('INFO(rinexFindEpochs304): The header of the rinex observation file does not contain TIME OF LAST OBS.\n' \
#         'Calculating nepochs will will require considerably more time. Consider editing rinex header to include\n' \
#         'TIME OF LAST HEADER')
    
#    line = fid.readline().rstrip()
#    end_of_header_reached = 0

#    while  not end_of_header_reached: #  Read past header
#        # answer = strfind(line,'END OF HEADER'); #  [] if string not found 
#        if 'END OF HEADER' in line: 
#        # if not np.isempty(answer):
#            line = fid.readline().rstrip()

    
#    while line != '':
#        # answer = strfind(line,'>'); #  [] if string not found 
#        # if ~isempty(answer)
#        if '>' in line:
#           nepochs = nepochs + 1; 
#        line = fid.readline().rstrip()
   
# print('INFO(rinexFindNEpochs304): Amount of epochs have been computed')
# fid.close()

# return nepochs, tLastObs, tInterval, success





# ## Testing input arguments

# ## Test type of GNSS systems
# # if ~isa(GNSSsystems, 'cell')
# #    sprintf(['INPUT ERROR(rinexReadObsBlock304): The input argument GNSSsystems',...
# #         ' is of type #s.\n Must be of type cell'],class(GNSSsystems))
# #     success = 0;
# #     return;  
# # end

# ## Test type of numSV
# if type(numSV) != int: 
# # if ~isa(numSV, 'double')
#    print('INPUT ERROR(rinexReadObsBlock304): The input argument numSV is of type %s.\n Must be of type double' % (type(numSV)))
#    success = 0
    

# nObsCodes = [int(x) for x in nObsCodes]
# ## Test type of numOfObsCodes
# if type(nObsCodes[0]) != int: 
# # if ~isa(numSV, 'double')
#    print('INPUT ERROR(rinexReadObsBlock304): The input argument numOfObsTypes is of type %s.\n Must be of type double' % (type(nObsCodes)))
#    success = 0
# # if ~isa(nObsCodes, 'double')
# #    sprintf(['INPUT ERROR(rinexReadObsBlock304): The input argument numOfObsTypes',...
# #         ' is of type #s.\n Must be of type double'],class(nObsCodes))
# #     success = 0;
# #     return;  
# # end


# ## Test size of numOfObsCodes
# if len(nObsCodes) != len(GNSSsystems):
#     print('INPUT ERROR(rinexReadObsBlock304): The input argument numOfObsTypes must have same length as GNSSsystems')
#     success = 0



# ## Test type of obsCodeIndex
# # if ~isa(obsCodeIndex, 'cell')
# #    sprintf(['INPUT ERROR(rinexReadObsBlock304): The input argument numOfObsTypes',...
# #         ' is of type #s.\n Must be of type cell'],class(obsCodeIndex))
# #     success = 0;
# #     return;  
# # end

# ## Test type of cell elements of removed ObsCodeIndex
# # for k = 1:length(GNSSsystems)
# #     if ~isa(obsCodeIndex{k}, 'double')
# #         sprintf(['INPUT ERROR(rinexReadObsBlock304): Element with index #d in',...
# #             ' removedObsTypesIndex is not of type double as it should be. \nVariable has type #s'],...
# #             k, class(obsCodeIndex{k}))
# #         success = 0;
# #         return;
# #     end
# # end

# ##

# success = 1;
# eof     = 0;

# # Highest number of obs codes of any GNSS system
# max_n_obs_Types = max(nObsCodes)

# # Initialize variables
# # Obs = zeros(numSV, max_n_obs_Types) ;
# Obs = np.empty([numSV, max_n_obs_Types]) 
# # SVlist = cell(numSV,1);
# # SVlist = dict(numSV)
# # SVlist = np.arange(numSV).reshape(numSV,1)   
# # SVlist.fill(numSV)
# # SVlist = np.chararray(numSV).reshape(numSV,1)   
# # SVlist.fill(numSV)
# SVlist = [np.nan]*numSV
# if readLLI:
#    # LLI = zeros(numSV, max_n_obs_Types) ;
#    LLI = np.empty([numSV, max_n_obs_Types]) 

# if readSS:
#    # SS  = zeros(numSV, max_n_obs_Types) ;
#    SS  = np.empty([numSV, max_n_obs_Types]) 

              

# # number of satellites excluded so far
# removed_sat = 0                                 
# # desiredGNSSsystems = str(GNSSsystems)      # DETTE MÅ TROLIG ENDRES!!   
# desiredGNSSsystems = list(GNSSsystems.values())
# # Gobble up observation block
# # for sat = 1:numSV
# for sat in range(0,numSV):
    
#    line = fid.readline().rstrip()   

       
#   # if line == -1
#   #   eof = 1;
#   #   disp(['ERROR (rinexReadsObsBlock211): the end of the '...
#   #         'observations text file was reached unexpectedly'])
#   #   success = 0;
#   #   return
#   # end

#   # SV = line(1:3); # Satellite code, ex. 'G11' or 'E03'
#    SV = line[0:3] # Satellite code, ex. 'G11' or 'E03'
#    if SV[0] not in desiredGNSSsystems:
#        removed_sat +=1
#    else:
#        ## Index of current GNSS system
#        # GNSSsystemIndex = find([GNSSsystems{:}] == SV(1)); 
#         # GNSSsystemIndex = {i for i in GNSSsystems if GNSSsystems[i]==SV[0]} # set virkelig?
#        GNSSsystemIndex = [i for i in GNSSsystems if GNSSsystems[i]==SV[0]][0] 
# #      SVlist{sat - removed_sat,1} = SV; # Store SV of current row
#        SVlist[sat - removed_sat] = SV # Store SV of current row
       
# #        n_obs_current_system = nObsCodes(GNSSsystemIndex);
#        n_obs_current_system = nObsCodes[GNSSsystemIndex-1]
       
#         # for obs_num = 1:n_obs_current_system
#        for obs_num in range(0, n_obs_current_system):
# #           obsIndex = obsCodeIndex{GNSSsystemIndex}(obs_num);
#            obsIndex = obsCodeIndex[GNSSsystemIndex][obs_num]
#            # charPos = 4+(obsIndex-1)*16;
#            charPos = 4+(obsIndex)*16
#            ## check that the current observation of the current GNSS system
# #          ## is not on the list of obs types to be excluded
         
#            ## stringlength of next obs. 
#            # obsLen = min(14, length(line) - charPos); 
#            obsLen = min(14, len(line) - charPos)
#            # read next obs
#            # newObs = strtrim(line(charPos:charPos+obsLen)); 
#            newObs = line[charPos:charPos+obsLen].strip()
#            # If observation missing, set to 0
#            # if ~isempty(newObs):
               
#            if newObs != '':
#                newObs = float(newObs)
#            else:
#                newObs = 0;
#            # Store new obs
#            Obs[sat - removed_sat, obs_num] = newObs 
            
#            if readLLI:
#            # read LLI of current obs (if present)
#                if charPos+14<=len(line):
#                    newLLI = line[charPos+14] # loss of lock indicator
#                else:
#                    newLLI = ' '
            

#             # if no LLI set to -999
#                # if isspace(newLLI):
#                #     newLLI = -999
#                # else:
#                #     newLLI = int(newLLI)
#                if newLLI.isspace():
#                    newLLI = -999
#                else:
#                    newLLI = int(newLLI)
#                 # Store LLI
#                LLI[sat - removed_sat, obs_num] = newLLI        
            
    
#            if readSS:
#                # read SS of current obs (if present)
#                if charPos+15<=len(line):
#                    newSS = line[charPos+15] # signal strength
#                else:
#                    newSS = ' ';
                
            
#                 # if no SS set to -999
#                if newSS.isspace():
#                    newSS = -999;
#                else:
#                    newSS = int(newSS)

#                 #Store SS
#                SS[sat - removed_sat, obs_num]  = newSS         
 


# # # update number og satellites after satellites have been excluded
# numSV = numSV - removed_sat

# #remove empty cell elements
# # SVlist(cellfun('isempty',SVlist))   = [];
# SVlist = list(filter(None,SVlist))
# # Obs[end-removed_sat+1:end, :]       = [];
# Obs[-1-removed_sat+1::, :]       = [];

# # end
