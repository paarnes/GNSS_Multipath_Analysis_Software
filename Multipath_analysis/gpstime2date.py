import numpy as np
from datetime import datetime, timedelta
def gpstime2date(week, tow):
    """
    Calculates date from GPS-week number and "time-of-week" to Gregorian calendar.
    
    
    Example:
    week = 2236
    tow = 35898
    date = gpstime2date(week,tow) --> 2022-11-13 09:58:00  (13 november 2022)
        
    Parameters
    ----------
    week : GPS-week  
    tow : "Time of week" 
    
    Returns
    -------
    date : The date given in the Gregorian calender ([year, month, day, hour, min, sec]) 
    
    """
    hour = np.floor(tow/3600)
    res = tow/3600 - hour
    min_ = np.floor(res*60)
    res = res*60-min_
    sec = res*60
    
    # if hours is more than 24, extract days built up from hours
    days_from_hours = np.floor(hour/24)
    # hours left over
    hour = hour - days_from_hours*24
    
    ## -- Computing number of days
    days_to_start_of_week = week*7
    
    # Origo of GPS-time: 06/01/1980 
    # t0 = date.toordinal(date(1980,1,6))+366
    t0 = datetime(1980,1,6)
    # t0 = t0.strftime("%Y %m %d")
    # t1 = t0 + days(days_to_start_of_week + days_from_hours); 
    t1 = t0 + timedelta(days=(days_to_start_of_week + days_from_hours))
    
    
    ## --  Formating the date to "year-month- day"
    t1 = t1.strftime("%Y %m %d")
    t1_ = [int(i) for i in t1.split(" ")]
    
    [year, month, day] = t1_
    
    date_ = [year, month, day, hour, min_, sec]
    return date_

# ## TESTING OUT NEW FUNTION
# from datetime import datetime, timedelta

# def gpstime2date(week, time_of_week):
#     # GPS epoch (January 6, 1980)
#     gps_epoch = datetime(1980, 1, 6)
    
#     # Calculate time in seconds
#     time_in_seconds = (week * 7 * 24 * 60 * 60) + time_of_week
#     # Round to nearest second
#     time_in_seconds = round(time_in_seconds, 0)
#     # Add seconds to GPS epoch
#     date = gps_epoch + timedelta(seconds=time_in_seconds)
#     # return date in the format of a list
#     return [date.year, date.month, date.day, date.hour, date.minute, date.second]

# week = 2190
# tow = 518399
# gpstime2date(week, tow)
# gps_to_date(week, tow)
