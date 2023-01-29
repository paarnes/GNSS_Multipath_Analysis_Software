from datetime import datetime
from scipy.interpolate import lagrange

def preciseOrbits2ECEF(sys, PRN, date_, dates, epochInterval, nEpochs, sat_positions, navGNSSsystems):
    """
     Function that finds positions of speficied satelite at nearest epochs.
     Then interpolates position at specified time from these epoch positions
    --------------------------------------------------------------------------------------------------------------------------
     INPUTS
    
     sys:            Satellite system, string, 
                       ex. "E" or "G"
    
     PRN:            Satellite identification number, integer
    
     date_:             array, date of specified epoch in form of [year, month, day, hour, min, sec]
    
     dates:            matrix. Each row contains date of one of the epochs in 
                       the SP3 orbit file. 
                       [nEpochs x 6]
    
    
     epochInterval:    interval of position epochs in SP3 file, seconds
    
     nEpochs:          number of position epochs in SP3 file, integer
    
     sat_positions:    Dict. Each cell elements contains position data for a
                       specific GNSS system. Order is defined by order of 
                       navGNSSsystems. Each cell element is another cell that 
                       stores position data of specific satellites of that 
                       GNSS system. Each of these cell elements is a matrix 
                       with [X, Y, Z] position of a epoch in each row.
    
                       sat_positions[GNSSsystemIndex][PRN][epoch, :] = [X, Y, Z]
    
     navGNSSsystems:   Dict. Contains char. Each string is a code for a
                       GNSS system with position data stored in sat_positions.
                       Must be one of: 'G', 'R', 'E', 'C'
    --------------------------------------------------------------------------------------------------------------------------
     OUTPUTS
    
     X, Y, Z:          ECEF coordinates of satellite at desired time computed
                       by interpolation 
    --------------------------------------------------------------------------------------------------------------------------
    """

    ## Formating the date 
    date_ = datetime(*date_)
    dates =  np.array(np.array(dates, dtype=float),dtype=int)
    datetime_array = create_datetime_array(dates)
    datetime_array = find_closest_dates(datetime_array, date_)
    
    GNSSsystemIndex = [idx for idx,val in enumerate(navGNSSsystems) if navGNSSsystems[idx]==sys][0] 
    curr_sys = navGNSSsystems[GNSSsystemIndex]
    sat_coord = []
    for idx in np.arange(0,len(sat_positions[curr_sys])):
        curr_ep = sat_positions[curr_sys][idx][PRN]
        sat_coord.append(curr_ep[0])

    sp3_data = []
    for idx in np.arange(0,len(datetime_array)):
        sp3_data.append((datetime_array[idx],sat_coord[idx][0],sat_coord[idx][1],sat_coord[idx][2]))
        
    x,y,z = interpolate_ecef(sp3_data, date_)
    return x,y,z




def interpolate_ecef(sp3_data, time):
    times = [t[0] for t in sp3_data]
    xs = [t[1] for t in sp3_data]
    ys = [t[2] for t in sp3_data]
    zs = [t[3] for t in sp3_data]
    
    # compute the time differences between the input time and the times of each satellite position
    time_diffs = [(time - t).total_seconds() for t in times]
    
    x_interp = lagrange(time_diffs, xs)
    y_interp = lagrange(time_diffs, ys)
    z_interp = lagrange(time_diffs, zs)
    
    return x_interp(0), y_interp(0), z_interp(0)


def create_datetime_array(dates):
    # Create an empty list to store the datetime objects
    datetime_list = []
    
    # Iterate over the rows of the array
    for row in dates:
        # Convert the elements of the row to integers
        year, month, day, hour, minute, second = row
        # Create a datetime object
        date = datetime(year, month, day, hour, minute, second)
        # Append the datetime object to the list
        datetime_list.append(date)
    
    # Convert the list to numpy array
    datetime_array = np.array(datetime_list)
    return datetime_array


from datetime import timedelta

def find_closest_dates(datetime_array, target_date):
    diff = abs(datetime_array - target_date)
    sort_indices = np.argsort(diff)[:6]
    return datetime_array[sort_indices]












# def interpolate_ecef(sp3_data, time):
#     times = [t[0] for t in sp3_data]
#     xs = [t[1] for t in sp3_data]
#     ys = [t[2] for t in sp3_data]
#     zs = [t[3] for t in sp3_data]
    
#     # compute the time differences between the input time and the times of each satellite position
#     time_diffs = [(time - t).total_seconds() for t in times]
    
#     x_interp = lagrange(time_diffs, xs)
#     y_interp = lagrange(time_diffs, ys)
#     z_interp = lagrange(time_diffs, zs)
    
#     return x_interp(0), y_interp(0), z_interp(0)

# # Example usage:
# sp3_data = [
#     (datetime(2022, 1, 1, 0, 0, 0), 13882271.938,-21710006.216 ,5357125.462),
#     (datetime(2022, 1, 1, 0, 5, 0), 13848517.625,-21468350.061 ,6298554.892)
# ]


# time = datetime(2022, 1, 1, 0, 4, 59,999999)

# koord = interpolate_ecef(sp3_data, time)
# print(koord)