function [success, rinexVersion, gnssType, markerName, recType, antDelta,...
GNSSsystems,numOfObsCodes, obsCodes, obsCodeIndex,tFirstObs, tLastObs, tInterval, ...
timeSystem, numHeaderLines, clockOffsetsON, rinexProgr, rinexDate,leapSec, approxPosition, GLO_Slot2ChannelMap, eof, fid] = ...
rinexReadObsFileHeader304(filename, includeAllGNSSsystems, includeAllObsCodes,...
desiredGNSSsystems, desiredObsCodes, desiredObsBands)


eof         = 0                     ;
success     = 1                     ;
warnings    = 0                     ;
antDelta    = [0; 0; 0]             ;
timeSystem  = ''                    ;
tFirstObs   = [0; 0; 0; 0; 0; 0]    ;
tLastObs    = NaN                   ;
tInterval   = NaN                   ;
rinexProgr  = NaN                   ;
rinexDate   = NaN                   ;
obsCodes    = {}                    ;
GNSSsystems = {}                    ;
gnssType    = ""                    ;
markerName  = ""                    ;
numHeaderLines  = 0                 ;
clockOffsetsON  = 0                 ;
numGNSSsystems  = 0                 ;
leapSec         = NaN               ;
numOfObsCodes   = []                ;
rinexHeader     = {}                ;
approxPosition  = [0, 0, 0]         ;
obsCodeIndex = {}                   ;
rinexVersion = NaN                  ;
recType = NaN                       ;
GLO_Slot2ChannelMap = NaN           ;
## Testing input arguments


# Test if filename is valid format
if ~isa(filename,'char') && ~isa(filename,'string')
    sprintf(['INPUT ERROR(rinexReadsObsHeader304): The input argument filename',...
        ' is of type #s.\n Must be of type string or char'],class(filename))
    success = 0;
    fid     = 0;
    return;
end

# Test if includeAllGNSSsystems is boolean
if includeAllGNSSsystems~=1 && includeAllGNSSsystems~=0
    sprintf(['INPUT ERROR(rinexReadsObsHeader304): The input argument includeAllGNSSsystems',...
        'must be either 1 or 0'])
    success = 0;
    fid     = 0;
    return
end

# Test if desiredGNSSsystems is valid type, given that includeAllGNSSsystems is 0 
if ~isa(desiredGNSSsystems, 'string') && includeAllGNSSsystems == 0
    sprintf(['INPUT ERROR(rinexReadsObsHeader304): The input argument desiredGNSSsystems',...
        ' must be of type string(string array) as long as \nincludeAllGNSSsystems is 1. Variable is of type #s.'],...
    class(desiredGNSSsystems))
    success = 0;
    fid     = 0;
    return 
end

# Test if desiredGNSSsystems is empty, given that includeAllGNSSsystems is 0
if isempty(desiredGNSSsystems) && includeAllGNSSsystems == 0
    sprintf(['INPUT ERROR(rinexReadsObsHeader304): The input argument desiredGNSSsystems',...
        ' can not be empty as long as includeAllGNSSsystems is 0'])
    success = 0;
    fid     = 0;
    return
end 

# Test if includeAllobsCodes is boolean
if includeAllObsCodes~=1 && includeAllObsCodes~=0
    sprintf(['INPUT ERROR(rinexReadsObsHeader304): The input argument includeAllobsCodes',...
        'must be either 1 or 0'])
    success = 0;
    fid     = 0;
    return
end

# Test if desiredobsCodes is valid type, given that includeAllobsCodes is 0
if ~isa(desiredObsCodes, 'string') && includeAllObsCodes == 0
    sprintf(['INPUT ERROR(rinexReadsObsHeader304): The input argument desiredobsCodes',...
        ' must be of type string(string array) as long as \nincludeAllObsTpes is 1. Variable is of type #s'],...
        class(desiredObsCodes))
    success = 0;
    fid     = 0;
    return 
end

# Test if desiredGNSSsystems is empty, given that includeAllGNSSsystems is 0
if isempty(desiredObsCodes) && includeAllObsCodes == 0
    sprintf(['INPUT ERROR(rinexReadsObsHeader304): The input argument desiredobsCodes',...
        ' can not have length 0 as long as includeAllobsCodes is 0'])
    success = 0;
    fid     = 0;
    return
end



##

# Open rinex observation file
fid         = fopen(filename,'rt')  ;

if fid == -1 # 
    success = 0;
    eof = 0;
    disp('ERROR(rinexReadObsFileHeader304): RINEX observation file not found!')
    return
end

while 1 # Gobbling the header
    
  numHeaderLines = numHeaderLines + 1;
  line = fgetl(fid); # returns -1 if only reads EOF
  
  if line == -1 # eof reached
    success = 0;
    eof = 1;
    disp('ERROR(rinexReadObsFileHeader304): End of file reached unexpectedly. No observations were read.')
    break
  end
  
  ##
  answer = strfind(line,'END OF HEADER'); # [] if the string isn't found
  if ~isempty(answer) # the end of the header was found
    break
  end
  
  ##
  if numHeaderLines == 1 # if first line of header
    # store rinex version
    rinexVersion = strtrim(line(1:9));
    # store rinex type, ex. "N" or "O"
    rinexType = line(21);
    
    # if rinex file is not an observation file
    if rinexType ~= 'O'  # Rinex file is oservation file
      disp('ERROR(rinexReadObsFileHeader304): the file is not a RINEX observations data file!')
      success = 0;
      fclose(fid);
      return
    end
    
    # Check gnss type
    gnssType = line(41); # reads the GNSS system type
    if ~ismember(gnssType, [' ' 'G' 'R' 'C' 'E' 'M' ])
        if ismember(gnssType, ['J' 'I' 'S'])
            disp(['ERROR(rinexReadObsFileHeader304): This software is meant for reading GNSS data only"' gnssType '"' ' is an invalid satellite '...
                                       'system type.']) 
        else
            disp(['ERROR(rinexReadObsFileHeader304): "' gnssType '"' ' is an unrecognized satellite '...
                                       'system type.'])
        end
      success = 0;
      fclose(fid);
      return
    end
    
    # if no system type, set G
    if strcmp(gnssType,' ') 
      gnssType = 'G';
    end
  end
  
  ##
  answer = strfind(line,'PGM / RUN BY / DATE');
    if ~isempty(answer)
      rinexProgr = strtrim(line(1:20)); # rinex program
      rinexDate = strtrim(line(41:60)); # rinex date
    end
    
  ##
  answer = strfind(line,'MARKER NAME');
  if ~isempty(answer)
    markerName = strtok(line); # markername
    
    # if no marker name, "MARKER" is read, so set to blank
    if strcmp(markerName,'MARKER') # if no marker name, "MARKER" is read, so set to blank
        markerName = '';
    end
  end
  
  ##
  answer = strfind(line,'ANTENNA: DELTA H/E/N');
  if ~isempty(answer)
    for k = 1:3
      [number, line] = strtok(line); # finds the substring containing the
                                     # deltas of the antenna relative to the marker
      antDelta (k,1) = str2double(number);
    end
  end
  
  ##
  # Section describing what GNSS systems are present, and their obs types
  answer = strfind(line,'SYS / # / OBS TYPES');
  
  if ~isempty(answer)
    line = strtrim(line(1:60));     # deletes 'SYS / # / OBS TYPES'

    [Sys, line] = strtok(line);     # reads GNSS system of this line
    [nObs, line]    = strtok(line)      ; # Num ObsCodes of this GNSS system
    nObs = str2double(nObs);
        
    # array for storing indeces of undesired ObsCodes for this GNSS
    # system
    undesiredobsCodeIndex = [];
    desiredObsCodeIndex = [];
        
        # is Sys amoung desired GNSS systems
        if (includeAllGNSSsystems && ismember(Sys, ["G", "R", "E", "C"])) || ismember(Sys, desiredGNSSsystems)

            numGNSSsystems  = numGNSSsystems + 1; # increment number of GNSS systems
            GNSSsystems{numGNSSsystems} = string(Sys); # Store current GNSS system



            GNSSSystemObsCodes = {}; # Reset cell of obsCodes for this GNSS system

            for k = 1:nObs 

                [obsCode, line] = strtok(line)      ; # read obsCode

                # Checking if obsCode is valid
                if size(obsCode,2) ~= 3 || ~ismember(obsCode(1),['C' 'L' 'D' 'S' 'I' 'X']) || ...
                        ~ismember(obsCode(2),['0' '1' '2' '4' '5' '6' '7' '8' '9']) || ...
                        ~ismember(obsCode(3),['A' 'B' 'C' 'D' 'I' 'L' 'M' 'N' 'P' 'Q' 'S' 'W' 'X' 'Y' 'Z']) 
                  sprintf(['ERROR (rinexReadsObsHeader304):  obsCode #s'...
                        ' is a not a standard RINEX 3.04 observation type!'], obsCode);     
                  success = 0;
                  fclose(fid);
                  return
                end

                # is obsCode amoung desired obscodes and frequency bands
                if includeAllObsCodes || ...
                        (ismember(obsCode(1), desiredObsCodes) && ismember(str2double(obsCode(2)), desiredObsBands))
                    
                    # store obsCode if amoung desire obsCodes
                    GNSSSystemObsCodes = [GNSSSystemObsCodes, obsCode];
                    desiredObsCodeIndex = [desiredObsCodeIndex, k];
                else
                    
                    # store index of discareded obsCode
                    undesiredobsCodeIndex = [undesiredobsCodeIndex, k];
                end

                # Every 13 obsCodes is at end of line. In this case read
                # next line and continue
                if mod(k, 13) == 0
                    numHeaderLines = numHeaderLines + 1;
                    line = fgetl(fid); # returns -1 if only reads EOF
                    line = strtrim(line(1:60));     # deletes 'SYS / # / OBS TYPES'
                end
            end    

            numOfObsCodes   = [numOfObsCodes; length(GNSSSystemObsCodes)]; # Store number of obsCodes for this GNSS system 
            obsCodes{numGNSSsystems,1}              = GNSSSystemObsCodes; # Store obsCodes for this GNSS system 
            obsCodeIndex{numGNSSsystems,1}          = desiredObsCodeIndex; # Store indices of desired obsCodes
        else
            # If GNSS system is not desired, skip to last line connected ot
            # it
            lines2Skip = floor(nObs/13);
            for i = 1:lines2Skip
                numHeaderLines = numHeaderLines + 1;
                line = fgetl(fid); # returns -1 if only reads EOF
            end
        end    
  end
    
      
    
 ## 
  answer = strfind(line,'TIME OF FIRST OBS');
  if ~isempty(answer)
  line = strtrim(line(1:60)); # deletes 'TIME OF FIRST OBS'
  
  for k = 1:6
    [tok, line] = strtok(line); # finds the substrings containing
                                # the components of the time of the
                                # first observation (YYYY; MM; DD;
                                # hh; mm; ss.sssssss) and specifies
                                # the Time System used in the
                                # observations file (GPST, GLOT or
                                # GALT)
    switch k
      case 1
        yyyy    = str2num(tok);
      case 2
        mm      = str2num(tok);
      case 3
       dd       = str2num(tok);
     case 4
       hh       = str2num(tok);
     case 5
      mnt       = str2num(tok);
     otherwise
      ss        = str2num(tok);
    end
  end
  
  tFirstObs = [yyyy; mm; dd; hh; mnt; ss];
  
  
  # Get Time system
  aux = strtok(line); 
  switch aux
    case 'GPS'
      timeSystem = 'GPS';
    case 'GLO'
      timeSystem = 'GLO';
    case 'GAL'
      timeSystem = 'GAL';
    case 'BDT'
      timeSystem = 'BDT';   
    
    # If no Time System defined, use info from gnssType
    otherwise
      switch gnssType 
      case 'G'
        timeSystem = 'GPST';
      case 'R'
        timeSystem = 'GLOT';
      case 'E'
        timeSystem = 'GALT';
      case 'C'
        timeSystem = 'BDT';
      otherwise
        fprintf(['CRITICAL ERROR (rinexReadsObsHeader304):\n'...
            'The Time System of the RINEX observations file '...
            'isn''t correctly specified!\n'])
        success = 0;
        fclose(fid);
        return
      end
   end
  end
  
  ##
  answer = strfind(line,'TIME OF LAST OBS'); # This is an optional record
  if ~isempty(answer)
      
    for k = 1:6
      [tok, line] = strtok(line); # finds the substrings containing
                                  # the components of the time of the
                                  # last observation (YYYY; MM; DD;
                                  # hh; mm; ss.sssssss)
      switch k
        case 1
          yyyy  = str2num(tok);
        case 2
          mm    = str2num(tok);
        case 3
          dd    = str2num(tok);
        case 4
          hh    = str2num(tok);
        case 5
          mnt   = str2num(tok);
        otherwise
          ss    = str2num(tok);
      end
    end
    tLastObs = [yyyy; mm; dd; hh; mnt; ss];
  end
  
  ##
  answer = strfind(line,'INTERVAL'); # This is an optional record
  if ~isempty(answer)
    tInterval = str2double(strtok(line));
  end
  
  ##
  answer = strfind(line,'RCV CLOCK OFFS APPL'); # This is an optional record!
  if ~isempty(answer)
    if (strtok(line)=='0') 
      clockOffsetsON = 0;
    elseif (strtok(line)=='1') 
      clockOffsetsON = 1;
    else
      success = 0;
      disp(['ERROR (rinexReadsObsHeader304): unrecognized '...
            'receiver clock offsets flag!'])
      fclose(fid);
      return
    end
  end
  
  ##
  answer = strfind(line,'LEAP SECONDS');# This is an optional record
  if ~isempty(answer)
    leapSec = str2double(strtok(line));
  end
  
  ##
  
  # store approximate receiver position
  answer = strfind(line, 'APPROX POSITION XYZ');
  if ~isempty(answer)
      for k = 1:3
          [pos_k, line] = strtok(line);
          approxPosition(k) = str2double(pos_k);
      end
  end
  
  ##
  answer = strfind(line, 'GLONASS SLOT / FRQ #');
  if ~isempty(answer)
      [nGLOSat, line] = strtok(line);
      nGLOSat = str2double(nGLOSat);
      slotNumbers = zeros(1, nGLOSat);
      channels = zeros(1, nGLOSat);
      for k = 1:nGLOSat
          
          
          [slotNumber, line] = strtok(line);
          slotNumber = str2double(slotNumber(2:end));
          
          [channel, line] = strtok(line);
          channel = str2double(channel);
          
          slotNumbers(k) = slotNumber;
          channels(k) = channel;
          
          if mod(k, 8) == 0
              line = fgetl(fid); # end of line is reached so read next line
              numHeaderLines = numHeaderLines + 1;
          end
      end
      
      GLO_Slot2ChannelMap = containers.Map(slotNumbers, channels);
  end

  ##
  answer = strfind(line, 'REC # / TYPE / VERS');
  if ~isempty(answer)
      recType = line(21:40);
  end
  
  
end # End of Gobbling Header Loop

for k = 1:numGNSSsystems
    
    # Give info if any of GNSS systems had zero of desired obscodes.
    if numOfObsCodes(k) == 0||sum(tFirstObs) == 0
        switch GNSSsystems{k}
            case "G"
                fprintf(['INFO: (rinexReadsObsHeader304)\nNone of the GPS satellites '...
                        'had any of the desired obsCodes\n\n']);
            case "R"
                fprintf(['INFO: (rinexReadsObsHeader304)\nNone of the GLONASS satellites '...
                        'had any of the desired obsCodes\n\n']);
            case "E"
                fprintf(['INFO: (rinexReadsObsHeader304)\nNone of the Galileo satellites '...
                        'had any of the desired obsCodes\n\n']);
            case "C"
                fprintf(['INFO: (rinexReadsObsHeader304)\nNone of the BeiDou satellites '...
                        'had any of the desired obsCodes\n\n']);
        end
    end
end



# store rinex header info
rinexHeader = {rinexVersion; rinexType; gnssType; rinexProgr; rinexDate};
disp('INFO(rinexReadObsFileHeader304): Rinex header has been read')
end

######### end rinexReadObsFileHeader304 #########