import datetime

import numpy as np
# from datetime import datetime
from math import sin, cos, sqrt, atan2

import os, sys,numpy as np
sys.path.append(r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Read_RINEX_OBS')
sys.path.append(r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Read_RINEX_nav')
sys.path.append(r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Kepler2ECEF')
from read_rinex3_nav import read_rinex3_nav
from kepler2ecef import Satkoord2
from kepler2ecef import date2gpstime
from kepler2ecef import extract_nav_message
from readRinexObs304 import readRinexObs304
from tqdm import trange,tqdm
from compute_azimut_elev import compute_azimut_elev
import pandas as pd



def compute_GLO_coord_from_nav(ephemerides, time_epochs):
    """
    Function that use broadcast ephemerides (from rinex nav file) and interpolates to current time using 4th order Runge-kutta.  
    
    Parameters: 
        
    ephemerides: A array with ephemerides for current satellite
    time_epochs: An n_obs X 2 sized array with weeks and tow read from the rinex observation file.
    
    Return:
        
        
    """
    ## Parameters (from GLONASS Interface Control Document 1988)
    GM         = 398600.44e9     # Gravitational constant [m3/s2]   (product of the mass of the earth and and gravity constant)
    omega_e    = 7.292115e-5     # Earth rotation rate    [rad/sek]
    c          = 299792458       # Speed of light         [m/s]
    a          = 6378136         # Semi major axis PZ-90   [m]
    f          = 1/298.257839303 # Inverse flattening
    C_20       = -1082.63e-6     # Second zonal coefficient of spherical harmonic expression.
    ephemerides[0] = 9999
    ephemerides = ephemerides.astype(float)

    ## Read in data:
    toc = [ephemerides[1],ephemerides[2] ,ephemerides[3] ,ephemerides[4],ephemerides[5],ephemerides[6]] # year,month,day,hour,minute,second
    week,toc = date2gpstime(int(ephemerides[1]),int(ephemerides[2]) ,int(ephemerides[3]) ,int(ephemerides[4]),int(ephemerides[5]),int(ephemerides[6]))
    tauN = ephemerides[6]   # SV clock bias (sec) (-TauN)
    gammaN = ephemerides[7] # SV relative frequency bias (+GammaN)

    # week_rec, tow_rec = time_epochs[0,1]  # extracting tow 
    week_rec, tow_rec = time_epochs  # extracting tow 
    x_te = ephemerides[10]  # X-coordinates at t_e in PZ-90 [km]
    y_te = ephemerides[14]  # Y-coordinates at t_e in PZ-90 [km]
    z_te = ephemerides[18]  # Z-coordinates at t_e in PZ-90 [km]
    
    vx_te = ephemerides[11] # Velocity component at t_e in PZ-90 (v_x) [km/s]
    vy_te = ephemerides[15] # Velocity component at t_e in PZ-90 (v_y) [km/s]
    vz_te = ephemerides[19] # Velocity component at t_e in PZ-90 (v_z) [km/s]
    
    J_x = ephemerides[12]   # Moon and sun acceleration at t_e [km/sec**2]
    J_y = ephemerides[16]   # Moon and sun acceleration at t_e [km/sec**2]
    J_z = ephemerides[20]   # Moon and sun acceleration at t_e [km/sec**2]

    ## -- Convert from UTC to GPST by adding leap seconds
    leap_sec = get_leap_seconds(week,toc) # Get correct leap sec based on date
    toc_gps_time = toc + leap_sec # convert from UTC to GPST
    
    ## -- Find time difference
    if week_rec == week:
        tdiff = tow_rec - toc_gps_time
    else:
        time_eph = format_date_string(week,toc_gps_time)
        time_rec = format_date_string(week_rec, tow_rec)
        tdiff = (time_rec - time_eph).total_seconds()

    ## -- Clock correction (except for general relativity which is applied later)
    clock_err = tauN + tdiff * (gammaN)
    clock_rate_err = gammaN

    init_state = np.empty(6)
    init_state[0] = x_te
    init_state[1] = y_te
    init_state[2] = z_te
    init_state[3] = vx_te
    init_state[4] = vy_te
    init_state[5] = vz_te
    init_state = 1000*init_state        # converting from km to meters
    acc = 1000*np.array([J_x, J_y,J_z]) # converting from km to meters
    state = init_state
    tstep = 90
    if tdiff < 0:
        tt = -tstep
    elif tdiff > 0:
        tt = tstep
    while abs(tdiff) > 1e-9:
        if abs(tdiff) < tstep:
            tt = tdiff
        k1 = glonass_diff_eq(state, acc)
        k2 = glonass_diff_eq(state + k1*tt/2, -acc)
        k3 = glonass_diff_eq(state + k2*tt/2, -acc)
        k4 = glonass_diff_eq(state + k3*tt, -acc)
        state += (k1 + 2*k2 + 2*k3 + k4)*tt/6.0
        tdiff -= tt

    pos = state[0:3]
    vel = state[3:6]
    return pos, vel, clock_err, clock_rate_err


def format_date_string(week,toc):
    """
    Function for formating dates to be able to subtract with seconds and float
    """
    year,month,day,hour,min_,sec = np.array(gpstime2date(week, toc))
    time = str(int(year)) + "/" + str(int(month)) + "/" + str(int(day)) + " " + str(int(hour)) \
        + ":" + str(int(min_)) + ":" + str(sec)[0:9]      
    sec_dum = time.split(':')[-1]
    time = time.replace(sec_dum, str(format(float(sec_dum), '.6f'))) # removing decimals if more than 3
    time = datetime.datetime.strptime(time, "%Y/%m/%d %H:%M:%S.%f")   
    return time
        


def glonass_diff_eq(state, acc):
    """
    State is a vector containing x,y,z,vx,vy,vz from navigation message
    """
    J2 = 1.0826257e-3       # Second zonal coefficient of spherical harmonic expression.
    mu = 3.9860044e14       # Gravitational constant [m3/s2]   (product of the mass of the earth and and gravity constant)
    omega = 7.292115e-5     # Earth rotation rate    [rad/sek]
    ae = 6378136.0          # Semi major axis PZ-90   [m]
    r = np.sqrt(state[0]**2 + state[1]**2 + state[2]**2)
    ders = np.zeros(6)
    if r**2 < 0:
        return ders
    a = 1.5 * J2 * mu * (ae**2)/ (r**5)
    b = 5 * (state[2]**2) / (r**2)
    c = -mu/(r**3) - a*(1-b)
    
    ders[0:3] = state[3:6]
    ders[3] = (c + omega**2)*state[0] + 2*omega*state[4] + acc[0]
    ders[4] = (c + omega**2)*state[1] - 2*omega*state[3] + acc[1]
    ders[5] = (c - 2*a)*state[2] + acc[2]
    return ders




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
    
    import numpy as np
    from datetime import datetime, timedelta

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
 

def utc_to_gpst(t_utc):
    """
    Convert from UTC to GPST by adding leap seconds.
    """
    t_gpst = t_utc + get_leap_seconds(t_utc)
    return t_gpst


def get_leap_seconds(week,tow):
    """
    Add leap seconds based on date. Input is week and tow for current obs.
    """
    year,month,day,_,_,_ = gpstime2date(week, tow) # convert to gregorian date
    time = (year,month,day)
    if time <= (2006, 1, 1):
        raise ValueError("Have no history on leap seconds before 2006")
    elif time <= (2009, 1, 1):
        return 14
    elif time <= (2012, 7, 1):
        return 15
    elif time <= (2015, 7, 1):
        return 16
    elif time <= (2017, 1, 1):
        return 17
    else:
        return 18


def date2gpstime(year,month,day,hour,minute,seconds):
    """
    Computing GPS-week nr.(integer) and "time-of-week" from year,month,day,hour,min,sec
    Origin for GPS-time is 06.01.1980 00:00:00 UTC
    """
    from datetime import date
    from numpy import fix
    
    t0=date.toordinal(date(1980,1,6))+366
    t1=date.toordinal(date(year,month,day))+366 
    week_flt = (t1-t0)/7;
    week = fix(week_flt);
    tow_0 = (week_flt-week)*604800;
    tow = tow_0 + hour*3600 + minute*60 + seconds;
    
    return week, tow





