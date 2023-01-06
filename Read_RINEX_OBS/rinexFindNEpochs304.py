import numpy as np
from datetime import date
from numpy import fix
def rinexFindNEpochs304(filename, tFirstObs, tLastObs, tInterval):

    #  Function that computes number of epochs in Rinex 3.xx observation file.
    # --------------------------------------------------------------------------------------------------------------------------
    #  INPUTS
    
    #  filename:         RINEX observation filename
    
    #  tFirstObs:        time stamp of the first observation record in the RINEX
    #                    observations file; column vector 
    #                    [YYYY; MM; DD; hh; mm; ss.sssssss]; 
    
    #  tLastObs:         time stamp of the last observation record in the RINEX
    #                    observations file; column vector
    #                    [YYYY; MM; DD; hh; mm; ss.sssssss]. If this information
    #                    was not available in rinex observation header the
    #                    default value is Nan. In this case the variable is
    #                    determined in this function
    
    #  tInterval:        observations interval; seconds. If this information
    #                    was not available in rinex observation header the
    #                    default value is Nan. In this case the variable is
    #                    determined in this function.
    # --------------------------------------------------------------------------------------------------------------------------
    #  OUTPUTS
    
    #  nepochs:          number of epochs in Rinex observation file with
    #                    observations
    
    #  tLastObs:         time stamp of the last observation record in the RINEX
    #                    observations file; column vector
    #                    [YYYY; MM; DD; hh; mm; ss.sssssss]. If this information
    #                    was not available in rinex observation header the
    #                    default value is Nan. In this case the variable is
    #                    determined in this function
    
    #  tInterval:        observations interval; seconds. If this information
    #                    was not available in rinex observation header the
    #                    default value is Nan.
    
    #  success:                  Boolean. 1 if the function seems to be successful, 
    #                            0 otherwise
    # --------------------------------------------------------------------------------------------------------------------------
    
    #  ADVICE: The function rinexFindNEpochs() calculates the amount of observation epochs in
    #  advance. This calculation will be incredibly more effective if TIME OF
    #  LAST OBS is included in the header of the observation file. It is
    #  strongly advized to manually add this information to the header if it is 
    #  not included by default. 
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
        
    
    
    #  open observation file
    fid = open(filename, 'rt')
    seconds_in_a_week = 604800
    #  tLastObs is in header
    # if ~isnan(tLastObs):
    # if not True in np.isnan(tLastObs):
    if ~np.all(np.isnan(tLastObs)):  # endret 07.12.2022 
        #  tInterval is not in header
        if np.isnan(tInterval):
            tInterval_found = 0
            first_epoch_found = 0
            #  calculate tInterval
            while not tInterval_found: #  calculate tInterval
                line = fid.readline().rstrip()
                #  start of new epoch
                # answer = strfind(line,'>');
                if '>' in line:
                # if ~isempty(answer) 
                    if not first_epoch_found: #  first epoch
                        # first_epoch_time = cell2mat(textscan(line(2:end),'# f', 6));
                        first_epoch_time = line[1::]
                        first_epoch_time = [float(el) for el in line[1::].split(" ") if el != ""]
                        first_epoch_time = first_epoch_time[:6]
                        first_epoch_found = 1;
                    else: #  seconds epoch
                        # second_epoch_time = cell2mat(textscan(line(2:end),'# f', 6));
                        second_epoch_time = line[1::]
                        second_epoch_time = [float(el) for el in line[1::].split(" ") if el != ""]
                        second_epoch_time = second_epoch_time[:6]
                    
                        # tInterval = second_epoch_time(6)-first_epoch_time(6);
                        tInterval = second_epoch_time[5]-first_epoch_time[5]
                        tInterval_found = 1;
       
        fid.close(); fid = open(filename, 'rt')
        tFirstObs = tFirstObs.astype(int) 
        tLastObs = tLastObs.astype(int) 
        rinex_lines = fid.readlines()
        epoch_line = [line for line in rinex_lines if '>' in line] # list with all the line thats defines a epoch
        nepochs = len(epoch_line)
        # week1,tow1 = date2gpstime(tFirstObs[0][0],tFirstObs[1][0],tFirstObs[2][0],tFirstObs[3][0],tFirstObs[4][0],tFirstObs[5][0])
        # week2,tow2 = date2gpstime(tLastObs[0][0],tLastObs[1][0],tLastObs[2][0],tLastObs[3][0],tLastObs[4][0],tLastObs[5][0])
        # n = ((week2*seconds_in_a_week+tow2) - (week1*seconds_in_a_week+tow1)) / tInterval
        # nepochs = np.floor(n) + 1;  
    
    #  if tLastObs is not in header. Function counts number of epochs manually 
    #  OBS: This can be VERY inefficent
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
    
        tLastObs = np.array([[yyyy],[mm],[dd],[hh],[mnt],[ss]]).astype(int) 
        print('INFO(rinexFindNEpochs304): TIME OF LAST OBS has been found and amount of epochs have been computed')
        fid.close()
    return int(nepochs), tLastObs, tInterval, success


def date2gpstime(year,month,day,hour,minute,seconds):
    """
    Computing GPS-week nr.(integer) and "time-of-week" from year,month,day,hour,min,sec
    Origin for GPS-time is 06.01.1980 00:00:00 UTC
    """
    
    t0=date.toordinal(date(1980,1,6))+366
    t1=date.toordinal(date(year,month,day))+366 
    week_flt = (t1-t0)/7;
    week = fix(week_flt);
    tow_0 = (week_flt-week)*604800;
    tow = tow_0 + hour*3600 + minute*60 + seconds;
    
    return week, tow