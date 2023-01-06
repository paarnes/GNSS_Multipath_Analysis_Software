# rinexReadObsBlockHead304(fid)
# Reads the metadata in the head of a RINEX 3.xx observations block, NOT 
# the header of the file.
#
# ATTENTION: Ignores all data in blocks with event flags with numbers
# greater than 1!!!
#
# Positioned in a RINEX 3.04 GNSS observations text file at the beginning
# of an observation block. In rinex 3.xx the line starts with '> '

# Based on the work of Ant�nio Pestana, rinexReadObsBlockHead211.m,
#--------------------------------------------------------------------------------------------------------------------------
# INPUTS

# fid:              Matlab identifier of an open RINEX 3.04 GNSS 
#                   observations text file positioned at the beginning
#                   of an observation block.
#--------------------------------------------------------------------------------------------------------------------------
# OUTPUTS

# success:          1 if function performs successfully, 0 otherwise
# 
# epochflag:        Rinex observations epoch flag, as follows:
#                       0: OK
#                       1: power failure between previous and current epoch
#                   From now on the "event flags":
#                       2: start moving antenna
#                       3: new site occupation
#                       4: header information follows
#                       5: external event (epoch is significant)

# clockOffset:          value of the receiver clock offset. If not present 
#                       in the metadata of the observations block 
#                       (it's optional RINEX 3.04 data)it is assumed to be 
#                       zero. If not zero implies that epoch, code, and 
#                       phase data have been corrected by applying 
#                       realtime-derived receiver clock offset

# date:                 time stamp of the observations block. Six-elements column-vector
#                       as follows:
#                           year: four-digits year (eg: 1959)
#                           month: integers 1..12
#                           day: integers 1..31
#                           hour: integers 0..24
#                           minute: integers 0..60
#                           second: reals 0..60

# numSV:                number of satellites with observations in with
#                       observations. This will include all satellite
#                       systems.
#--------------------------------------------------------------------------------------------------------------------------
##

filename = 'opec0020_3.04_kort.10o'
fid = open(filename,'r') 


# Initialize variables
success = 1     ;
eof     = 0     ;
date    = []    ;
numSV   = 0     ;
epochflag   = [];
clockOffset = [];
noFlag = 1;



line = fid.readline().rstrip()  # DENNE KODE KAN SIKKER FJERNES NÅR TING SETTES I SAMMEN NÅR FID GÅR INN
while 'END OF HEADER' not in line:
    line = fid.readline().rstrip()
    
line = fid.readline().rstrip()    

# end of observation file is reached at an expected point
# if line == -1
#   eof = 1;
#   disp(['INFO(rinexReadObsBlockHead304): End of observations '...
#         'text file reached'])
#   return
# end

# The first thing to do: the reading of the epoch flag
epochflag   = line[31]

# skip to next block if event flag is more than 1
while int(epochflag) > 1:
    noFlag = 0 
    linejump = int(line[32:35])
    msg = 'WARNING(rinexReadsObsBlockHead304): Observations event flag encountered. Flag = %s. %s lines were ignored.' % (str(epochflag), str(linejump))
    print(msg)
    for count in range(0,linejump+1):
    # for count=1:linejump + 1 # skip over current obs block and go to next one
        line = fid.readline().rstrip()
    epochflag = int(line[30])


# Gets the number of used satellites in obs epoch
numSV = int(line[32:35])

# Gets the receiver clock offset. This is optional data!
clockOffset = 0
# if size(line,2) == 56:
if len(line) == 56:
    clockOffset = int(line[41:56])


# # Reads the time stamp of the observations block (6 numerical values)
# date = lin;
date = line[1::]
date = [float(el) for el in line[1::].split(" ") if el != ""]
date = date[:6]

if noFlag == 0:
    msg2 = msg + '\nEpoch date = %.4d %.2d %.2d %.2d:%.2d:%6.4f' % (date[0],date[1],date[2],date[3],date[4],date[5])
    print(msg2)
    



# # return success, epochflag, clockOffset, date, numSV, eof