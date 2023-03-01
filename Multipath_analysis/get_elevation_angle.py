from Geodetic_functions import *
import numpy as np 
from numpy import fix,log,fmod,arctan,arctan2,sqrt
from get_elevation_angle import ECEF2enu
from gpstime2date import gpstime2date
from preciseOrbits2ECEF import preciseOrbits2ECEF
import datetime
import warnings
warnings.filterwarnings(action='ignore', message='invalid value encountered in fmod')

def get_elevation_angle(sys, PRN, week, tow, sat_positions, nEpochs, epoch_dates, epochInterval, navGNSSsystems, x_e):
    """
    Calculates elevation angle of a satelite with specified PRN at specified
    epoch, viewed from defined receiver position
    
    --------------------------------------------------------------------------------------------------------------------------
    INPUTS:
    ------
    
    sys:              Satellite system, string. ex. "E" or "G"
    
    PRN:              Satellite identification number, integer
     
    week:             GPS-week number, float
    
    tow:              "time-of-week", float
  

    sat_positions:    dictionary containing satellite navigation ephemeris of all 
                      satellites of observation period, of each GNSS system. The
                      structure is like this sat_positions[systemcode][epoch][PRN]. 
                      Then you get X,Y and Z coordinates. Ex: sat_positions['G'][100][24]
                      will extract GPS position at epoch 100 for PRN 24.
                      
     
    nEpochs:          integer, The number of navigation ephemeris epochs for this satellite.
    
    epoch_dates:      The gregorian date for each epock in SP3 file
    
    
    epochInterval:    float. The epoch interval in seconds. The time difference between epochs. 
    
                     
   
    navGNSSsystems:   list, conatins codes of GNSS systems with navigation data. 
                      ex: ['G', 'R']
                      
    x_e:             Coordinates, in ECEF reference frame, of receiver station ex. [X,Y,Z]
    
    
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS:
    --------
    
    elevation_angle: Elevation angle of specified satelite at specified
                      epoch, view from specified receiver station. Unit: Degrees 
   
    missing_nav_data: Boolean, 1 if orbit data for current satellite is
                      missing from sp3 file, o otherwise
                      
    Xs:               float. The computed X-coordinate
    Ys:               float. The computed Y-coordinate
    Zs:               float. The computed Z-coordinate
    
    --------------------------------------------------------------------------------------------------------------------------
    """
    
    ## -- Define GRS80 ellipsoid parameters
    a       = 6378137
    f       = 1/298.257222100882711243
    b       = a*(1-f)
    
    missing_nav_data = 0
    ## Get date in form of [year, month, week, day, min, sec] from GPS-week and tow
    date_ = gpstime2date(week, round(tow,1)) ## added round to prevent 59.99999 seconds   
    Xs, Ys, Zs = preciseOrbits2ECEF(sys, PRN, date_, epoch_dates, epochInterval, nEpochs, sat_positions, navGNSSsystems)
    if all([Xs,Ys,Zs]) == 0:
         missing_nav_data = 1
         elevation_angle = 0
         azimut_angle = 0 
    else:
        ##-- Define vector from receiver to satellite
        dx_e = np.array([Xs,Ys,Zs]) -  x_e

        ## -- Get geodetic coordinates of receiver
        lat,lon,h = ECEF2geodb(a,b,x_e[0][0],x_e[1][0],x_e[2][0])
            
        ##-- Transform dx_e vector to local reference frame
        e,n,u = ECEF2enu(lat,lon, dx_e[0][0],dx_e[1][0],dx_e[2][0])    
        ## -- Calculate elevation and azimut angle from receiver to satellite
        if not np.isnan(u) and not np.isnan(e) and not np.isnan(n):
            elevation_angle = np.rad2deg(atanc(u, sqrt(e**2 + n**2)))
            
            # Azimut computation and quadrant correction 
            if (e> 0 and n< 0) or (e < 0 and n < 0):
                azimut_angle = np.rad2deg(arctan(e/n)) + 180
            elif e < 0 and n > 0:
                azimut_angle = np.rad2deg(arctan(e/n)) + 360
            else:
                azimut_angle = np.rad2deg(arctan(e/n))

        else:
            elevation_angle = np.nan
            azimut_angle = np.nan

    return elevation_angle,azimut_angle, missing_nav_data, float(Xs), float(Ys), float(Zs)
