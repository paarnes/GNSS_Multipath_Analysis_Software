# function [success, rinexVersion, gnssType, markerName, recType, antDelta,...
# GNSSsystems,numOfObsCodes, obsCodes, obsCodeIndex,tFirstObs, tLastObs, tInterval, ...
# timeSystem, numHeaderLines, clockOffsetsON, rinexProgr, rinexDate,leapSec, approxPosition, GLO_Slot2ChannelMap, eof, fid] = ...
# rinexReadObsFileHeader304(filename, includeAllGNSSsystems, includeAllObsCodes,...
# desiredGNSSsystems, desiredObsCodes, desiredObsBands)
import numpy as np
import math, os 
filename = 'opec0020_3.04_kort.10o'

includeAllGNSSsystems = 1
includeAllObsCodes = 1
desiredGNSSsystems= ["G", "E", "C"]
desiredObsCodes = ["C", "L", "S", "D"]
desiredObsBands = [1,2,5]

eof         = 0                     
success     = 1                     
warnings    = 0                     
# antDelta    = [0; 0; 0]  
antDelta    = []             
timeSystem  = ''                    
# tFirstObs   = [0; 0; 0; 0; 0; 0]  
tFirstObs   = []    
# tLastObs    = NaN                   
# tInterval   = NaN                   
# rinexProgr  = NaN                   
# rinexDate   = NaN                   
obsCodes    = {}                    
GNSSsystems = {}                    
gnssType    = ""                    
markerName  = ""                    
numHeaderLines  = 0                 
clockOffsetsON  = 0                 
numGNSSsystems  = 0                 
# leapSec         = NaN               
numOfObsCodes   = []                
rinexHeader     = {}                
approxPosition  = [0, 0, 0]         
obsCodeIndex = {}                   
# rinexVersion = NaN                  
# recType = NaN                       
# GLO_Slot2ChannelMap = NaN           
## Testing input arguments


# Test if filename is valid format
if type(filename) != str:
    typ = type(filename)
    print('INPUT ERROR(rinexReadsObsHeader304): The input argument filename is of type %s.\n Must be of type string or char' %(typ))
    success = 0
    fid     = 0
    


# # Test if includeAllGNSSsystems is boolean
# if includeAllGNSSsystems~=1 && includeAllGNSSsystems~=0
#     sprintf(['INPUT ERROR(rinexReadsObsHeader304): The input argument includeAllGNSSsystems',...
#         'must be either 1 or 0'])
#     success = 0;
#     fid     = 0;
#     return
# end

# # Test if desiredGNSSsystems is valid type, given that includeAllGNSSsystems is 0 
# if ~isa(desiredGNSSsystems, 'string') && includeAllGNSSsystems == 0
#     sprintf(['INPUT ERROR(rinexReadsObsHeader304): The input argument desiredGNSSsystems',...
#         ' must be of type string(string array) as long as \nincludeAllGNSSsystems is 1. Variable is of type #s.'],...
#     class(desiredGNSSsystems))
#     success = 0;
#     fid     = 0;
#     return 
# end

# # Test if desiredGNSSsystems is empty, given that includeAllGNSSsystems is 0
# if isempty(desiredGNSSsystems) && includeAllGNSSsystems == 0
#     sprintf(['INPUT ERROR(rinexReadsObsHeader304): The input argument desiredGNSSsystems',...
#         ' can not be empty as long as includeAllGNSSsystems is 0'])
#     success = 0;
#     fid     = 0;
#     return
# end 

# # Test if includeAllobsCodes is boolean
# if includeAllObsCodes~=1 && includeAllObsCodes~=0
#     sprintf(['INPUT ERROR(rinexReadsObsHeader304): The input argument includeAllobsCodes',...
#         'must be either 1 or 0'])
#     success = 0;
#     fid     = 0;
#     return
# end

# # Test if desiredobsCodes is valid type, given that includeAllobsCodes is 0
# if ~isa(desiredObsCodes, 'string') && includeAllObsCodes == 0
#     sprintf(['INPUT ERROR(rinexReadsObsHeader304): The input argument desiredobsCodes',...
#         ' must be of type string(string array) as long as \nincludeAllObsTpes is 1. Variable is of type #s'],...
#         class(desiredObsCodes))
#     success = 0;
#     fid     = 0;
#     return 
# end

# # Test if desiredGNSSsystems is empty, given that includeAllGNSSsystems is 0
# if isempty(desiredObsCodes) && includeAllObsCodes == 0
#     sprintf(['INPUT ERROR(rinexReadsObsHeader304): The input argument desiredobsCodes',...
#         ' can not have length 0 as long as includeAllobsCodes is 0'])
#     success = 0;
#     fid     = 0;
#     return
# end



##

# Open rinex observation file
fid = open(filename,'r') 

if os.stat(filename).st_size == 0:
    raise ValueError('ERROR: This file seems to be empty')
    

# if fid == -1 # 
#     success = 0;
#     eof = 0;
#     disp('ERROR(rinexReadObsFileHeader304): RINEX observation file not found!')
#     return
# end

while 1: # Gobbling the header
     numHeaderLines = numHeaderLines + 1;
     # line = fgetl(fid); # returns -1 if only reads EOF
     line = fid.readline().rstrip()
     
     # if line == -1 # eof reached
     #   success = 0;
     #   eof = 1;
     #   disp('ERROR(rinexReadObsFileHeader304): End of file reached unexpectedly. No observations were read.')
     #   break
     
     
     ##
     if 'END OF HEADER' in line:
         break
     
     ##
     if numHeaderLines == 1: # if first line of header
       # store rinex version
       # rinexVersion = strtrim(line(1:9));
       rinexVersion = line[0:9]
       # store rinex type, ex. "N" or "O"
       # rinexType = line(21);
       rinexType = line[20]
       
       # if rinex file is not an observation file
       if rinexType != 'O':  # Rinex file is oservation file
         print('ERROR(rinexReadObsFileHeader304): the file is not a RINEX observations data file!')
         success = 0;
         fid.close()
         

       
       # Check gnss type
       gnssType = line[40] # reads the GNSS system type
       if gnssType not in [' ', 'G', 'R', 'C', 'E', 'M' ]:
           if gnssType in ['J', 'I', 'S']:
               print('ERROR(rinexReadObsFileHeader304): This software is meant for reading GNSS data only.\
                     %s is an invalid satellite system type.' %(gnssType)) 
           else:
               print('ERROR(rinexReadObsFileHeader304): %s is an unrecognized satellite system type.' %(gnssType))
          
           success = 0
           fid.close()
           
       
       
       # if no system type, set G
       if gnssType == ' ':
           gnssType = 'G'
           
     
      ##
    # answer = strfind(line,'PGM / RUN BY / DATE');
     if 'PGM / RUN BY / DATE' in line:
         rinexProgr = line[0:20] # rinex program
         rinexDate = line[40:60] # rinex date

    

    # answer = strfind(line,'MARKER NAME');
     if 'MARKER NAME' in line:
       # markerName = strtok(line); # markername
        markerName = line.strip() # markername
     
      # if no marker name, "MARKER" is read, so set to blank
     if 'Marker' in markerName:
         markerName = ''


     
      
     # answer = strfind(line,'ANTENNA: DELTA H/E/N');
     if 'ANTENNA: DELTA H/E/N' in line:
         for k in range(0,3):
             line_ = [el for el in line.split(" ") if el != ""]
             antDelta = [line_[0],line_[1],line_[2]]



      # Section describing what GNSS systems are present, and their obs types
      # answer = strfind(line,'SYS / # / OBS TYPES');
     
     if 'SYS / # / OBS TYPES' in line:
        # line = strtrim(line(1:60));     # deletes 'SYS / # / OBS TYPES'
        line = line[0:60]     # deletes 'SYS / # / OBS TYPES'
        line_ = [el for el in line.split(" ") if el != ""]
        Sys = line_.pop(0) # assingning system to variable and removing it from the list
        nObs = int(line_.pop(0))
        
        # [Sys, line] = strtok(line);      # reads GNSS system of this line
        # [nObs, line]    = strtok(line)      ; # Num ObsCodes of this GNSS system
        # nObs = str2double(nObs);
           
        # # array for storing indeces of undesired ObsCodes for this GNSS
        # # system
        undesiredobsCodeIndex = []
        desiredObsCodeIndex = []
           
            # is Sys amoung desired GNSS systems
        if (includeAllGNSSsystems and Sys in ["G", "R", "E", "C"] or Sys in desiredGNSSsystems):
      
              numGNSSsystems  = numGNSSsystems + 1; # increment number of GNSS systems
              GNSSsystems[numGNSSsystems] = str(Sys) # Store current GNSS system
              # GNSSsystems{numGNSSsystems} = string(Sys); # Store current GNSS system
       
              GNSSSystemObsCodes = {}; # Reset cell of obsCodes for this GNSS system
              
              for k in range(0,nObs):
              # for k = 1:nObs 
                  # [obsCode, line] = strtok(line)      ; # read obsCode
                  # line_ = [el for el in line.split(" ") if el != ""]
                  obsCode = line_.pop(0)
                  # Checking if obsCode is valid
                  if len(obsCode) != 3 or obsCode[0] not in ['C', 'L', 'D','S', 'I', 'X'] or  \
                  obsCode[1] not in ['0', '1', '2', '4', '5', '6', '7', '8', '9'] or \
                      obsCode[2] not in ['A', 'B', 'C', 'D', 'I', 'L', 'M', 'N', 'P', 'Q', 'S', 'W', 'X', 'Y', 'Z']:
                    print('ERROR (rinexReadsObsHeader304):  obsCode %s is a not a standard RINEX 3.04 observation type!' %(obsCode))     
                    success = 0
                    fid.close()
                    
      
                  # is obsCode amoung desired obscodes and frequency bands
                  if includeAllObsCodes or obsCode[0] in desiredObsCodes and int(obsCode[1]) in desiredObsBands:
                      ## store obsCode if amoung desire obsCodes
                      # GNSSSystemObsCodes = [GNSSSystemObsCodes, obsCode];
                      # desiredObsCodeIndex = [desiredObsCodeIndex, k];
                       GNSSSystemObsCodes[Sys] =  obsCode
                       desiredObsCodeIndex.append(k)
                  else:
                      # store index of discareded obsCode
                      # undesiredobsCodeIndex = [undesiredobsCodeIndex, k];
                      undesiredobsCodeIndex.append(k)
      
                  # Every 13 obsCodes is at end of line. In this case read
                  # next line and continue
                  if np.mod(k+1, 13) == 0:
                      numHeaderLines = numHeaderLines + 1;
                      line = fid.readline().rstrip()
                      line = line[0:60]     # deletes 'SYS / # / OBS TYPES'
                      line_ = [el for el in line.split(" ") if el != ""]
                      # line = fgetl(fid); # returns -1 if only reads EOF
                      # line = strtrim(line(1:60));     # deletes 'SYS / # / OBS TYPES'
                  
                  
              numOfObsCodes.append(len(GNSSSystemObsCodes))
              obsCodes[numGNSSsystems] = GNSSSystemObsCodes
              obsCodeIndex[numGNSSsystems] = desiredObsCodeIndex # Store indices of desired obsCodes
              # numOfObsCodes   = [numOfObsCodes; length(GNSSSystemObsCodes)]; # Store number of obsCodes for this GNSS system 
              # obsCodes{numGNSSsystems,1}              = GNSSSystemObsCodes; # Store obsCodes for this GNSS system 
              # obsCodeIndex{numGNSSsystems,1}          = desiredObsCodeIndex; # Store indices of desired obsCodes
        else:
            # If GNSS system is not desired, skip to last line connected or it
            lines2Skip = math.floor(nObs/13)
            for i in range(0,lines2Skip):
            # for i = 1:lines2Skip
                numHeaderLines = numHeaderLines + 1;
                line = fid.readline().rstrip()
                # line = fgetl(fid); # returns -1 if only reads EOF

       
         
    
     if 'TIME OF FIRST OBS' in line:
         line = line[0:60]     #  deletes 'TIME OF FIRST OBS'
         line_ = [el for el in line.split(" ") if el != ""]
         
         for k in range(0,6):
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

     
         # tFirstObs = [yyyy; mm; dd; hh; mnt; ss];
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
                 print('CRITICAL ERROR (rinexReadsObsHeader304):\n' \
                                  'The Time System of the RINEX observations file '\
                                  'isn''t correctly specified!\n')
                 success = 0
                 fid.close()


     if 'TIME OF LAST OBS' in line:  
         line = line[0:60]     #  deletes 'TIME OF LAST OBS'
         line_ = [el for el in line.split(" ") if el != ""]
         for k in range(0,6):

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


     
      ##
      # answer = strfind(line,'INTERVAL'); # This is an optional record
     if 'INTERVAL' in line:  
         tInterval = int(line_.pop(0));
      
     

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
         line = line[0:60]     #  deletes 'TIME OF LAST OBS'
         line_ = [el for el in line.split(" ") if el != ""]
         nGLOSat = int(line_.pop(0))
         slotNumbers = np.array([])
         channels = np.array([])
         for k in range(0,nGLOSat):
             # slotNumber = int(line_.pop(0)[1::])
             slotNumber = line_.pop(0)[1::]
             channel = int(line_.pop(0))
             slotNumbers = np.append(slotNumbers,slotNumber)
             channels = np.append(channels,channel)

             if np.mod(k+1, 8) == 0:
                 # line = fgetl(fid); # end of line is reached so read next line
                 line = fid.readline().rstrip()
                 numHeaderLines = numHeaderLines + 1
                 line = line[0:60]     #  deletes 'TIME OF LAST OBS'
                 line_ = [el for el in line.split(" ") if el != ""]

   

         # GLO_Slot2ChannelMap = containers.Map(slotNumbers, channels);
         GLO_Slot2ChannelMap = dict(zip(slotNumbers,channels.astype(int)))
       
    
      
         if 'REC # / TYPE / VERS' in line: 
             recType = line[20:40]

      
  
# End of Gobbling Header Loop

for k in range(0,numGNSSsystems):    
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
            



# store rinex header info
# rinexHeader = {rinexVersion; rinexType; gnssType; rinexProgr; rinexDate};
rinexHeader['rinexVersion'] =rinexVersion 
rinexHeader['rinexType'] = rinexType 
rinexHeader['gnssType'] =gnssType 
rinexHeader['rinexProgr'] =rinexProgr 
rinexHeader['rinexDate'] =rinexDate 


print('INFO(rinexReadObsFileHeader304): Rinex header has been read')


######### end rinexReadObsFileHeader304 #########