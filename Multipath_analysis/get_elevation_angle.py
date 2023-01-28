from Geodetic_functions import *

def get_elevation_angle(sys, PRN, week, tow, sat_positions, nEpochs, epoch_dates, epochInterval, navGNSSsystems, x_e):
    """
    Calculates elevation angle of a satelite with specified PRN at specified
    epoch, viewed from defined receiver position
    
    --------------------------------------------------------------------------------------------------------------------------
    INPUTS
    sys:            Satellite system, string, 
                     ex. "E" or "G"
    
    PRN:            Satellite identification number, integer
     
    week:             GPS-week number, integer
    
    tow:              "time-of-week", integer
  
    x_e:              Coordinates, in ECEF reference frame, of receiver station,
                      ex. [X;Y;Z]
    
    ephemerisCell:    cell containing satellite navigation ephemeris of all 
                      satellites of observation period, of each GNSS system. 
                      Each cell element contains ephemeris of one GNSS system.
                      Order of cells is defined by navGNSSsystems.
                      Each cell elements at the next level is a matrix 
                      containining ephemeris of one satellite. The ephemeris of 
                      one navigation block is conatined in one line.
                      ephemerisCell{GNSSsystemIndex}{PRN}
    
                      ephemerisCell{GNSSsystemIndex}{PRN}(navBlock, ephemeris)
     
    n_nav_blocks:     cell, one cell element for each GNSS system. Each cell
                      element is another cell with one element for each
                      satellite. Each of these elements contains the number
                      of navigation ephemeris block for this satellite.
                      n_nav_blocks{GNSSsystemIndex}{PRN}
   
    navGNSSsystems:   cell, conatins codes of GNSS systems with navigation 
                      ephemeris stored in ephemerisCell. Elements are char.
                      ex. 'G' or 'E'
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS:
   
    elevation_angle: Elevation angle of specified satelite at specified
                      epoch, view from specified receiver station. Unit: Degrees 
   
    missing_nav_data: Boolean, 1 if orbit data for current satellite is
                      missing from sp3 file, o otherwise
    --------------------------------------------------------------------------------------------------------------------------
    """
    import numpy as np 
    # from math import sqrt,sin,cos,tan,pi,atan,atan2,asin
    from numpy import fix,log,fmod,arctan,arctan2,sqrt
    from get_elevation_angle import ECEF2enu
    from gpstime2date import gpstime2date
    from preciseOrbits2ECEF import preciseOrbits2ECEF
    import datetime
    import warnings
    warnings.filterwarnings(action='ignore', message='invalid value encountered in fmod')
    
    ## -- Define GRS80 ellipsoid parameters
    a       = 6378137
    f       = 1/298.257222100882711243
    b       = a*(1-f)
    
    
    missing_nav_data = 0
    
    # get date in form of [year, month, week, day, min, sec] from GPS-week
    #  "time-of-week"
    # date_ = gpstime2date(week, tow)
    date_ = gpstime2date(week, round(tow,1)) ## added round to prevent 59.99999 seconds
    # date_ = datetime.datetime(int(date_[0]),int(date_[1]),int(date_[2]),int(date_[3]),int(date_[4]),int(date_[5]))
    
    
    Xs, Ys, Zs = preciseOrbits2ECEF(sys, PRN, date_, epoch_dates, epochInterval, nEpochs, sat_positions, navGNSSsystems)
    # interpol_coord = np.array([Xs,Ys,Zs])
    if all([Xs,Ys,Zs]) == 0:
         missing_nav_data = 1
         elevation_angle = 0
         azimut_angle = 0 # added 05.01.2023
    else:
        # Define vector from receiver to satellite
        dx_e = np.array([Xs,Ys,Zs]).reshape(len(np.array([Xs,Ys,Zs])),1) -  x_e.reshape(1,len(x_e))
    
        # Get geodetic coordinates of receiver
        station_x_g = ECEF2geodb(a,b,x_e[0][0],x_e[1][0],x_e[2][0])
    
        # transform dx_e vector to local reference frame
        # dx_l = ECEF2local(station_x_g, dx_e);
        lat = station_x_g[0]
        lon = station_x_g[1]
        dx_l = ECEF2enu(lat,lon, dx_e[0][0],dx_e[1][0],dx_e[2][0])
        # Local vector components
        e = dx_l[0]
        n = dx_l[1]
        u = dx_l[2]
    
        # Calculate elevation and azimut angle from receiver to satellite
        # if u != np.nan and e != np.nan and n != np.nan:
        if not np.isnan(u) and not np.isnan(e) and not np.isnan(n):
            elevation_angle = atanc(u, sqrt(e**2 + n**2))*180/pi
            
            # Azimut computation and quadrant correction 
            if (e> 0 and n< 0) or (e < 0 and n < 0):
                azimut_angle = (arctan(e/n)*(180/pi) + 180)
            elif e < 0 and n > 0:
                azimut_angle = (arctan(e/n)*(180/pi) + 360)
            else:
                azimut_angle = (arctan(e/n)*(180/pi))

        else:
            elevation_angle = np.nan
            azimut_angle = np.nan

    return elevation_angle,azimut_angle, missing_nav_data, float(Xs), float(Ys), float(Zs)
