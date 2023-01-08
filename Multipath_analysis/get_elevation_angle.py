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
    from math import sqrt,sin,cos,tan,pi,atan,atan2,asin
    from numpy import fix,log,fmod,arctan,arctan2
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
    date_ = gpstime2date(week, tow)
    # date_ = datetime.datetime(int(date_[0]),int(date_[1]),int(date_[2]),int(date_[3]),int(date_[4]),int(date_[5]))
    
    
    Xs, Ys, Zs = preciseOrbits2ECEF(sys, PRN, date_, epoch_dates, epochInterval, nEpochs, sat_positions, navGNSSsystems)
    
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
        if u != np.nan and e != np.nan and n != np.nan:
            elevation_angle = atanc(u, sqrt(e**2 + n**2))*180/pi
            
            # Azimut computation and quadrant correction 
            if (e> 0 and n< 0) or (e < 0 and n < 0):
                azimut_angle = (atan(e/n)*(180/pi) + 180)
            elif e < 0 and n > 0:
                azimut_angle = (atan(e/n)*(180/pi) + 360)
            else:
                azimut_angle = (atan(e/n)*(180/pi))

        else:
            elevation_angle = np.nan
            azimut_angle = np.nan

        

    return elevation_angle,azimut_angle, missing_nav_data



def ECEF2geodb(a,b,X,Y,Z):
    '''
    Konverter fra kartesiske ECEF-koordinater til geodetiske koordinater vha Bowrings metode.

    Parameters
    ----------
    a : Store halvakse
    b : Lille halvakse
    X : X-koordinat
    Y : Y-koordinat
    Z : Z-koordinat

    Returns
    -------
    lat : Breddegrad
    lon : Lengdegrad
    h :   Høyde

    '''
    import numpy as np 
    from math import sqrt,sin,cos,tan,pi,atan,atan2,asin
    from numpy import fix,log,fmod,arctan,arctan2
    from get_elevation_angle import Nrad
    
    e2m = (a**2 - b**2)/b**2
    e2  = (a**2 - b**2)/a**2
    rho = sqrt(X**2 +Y**2)
    my  = atan((Z*a)/(rho*b))
    lat = atan(( Z +e2m*b*(sin(my))**3)/(rho - e2*a*(cos(my))**3))
    lon = atan(Y/X)
    N   = Nrad(a,b,lat)
    h   = rho*cos(lat) + Z*sin(lat) - N*( 1 - e2*(sin(lat))**2)
    return lat, lon, h


def ECEF2enu(lat,lon,dX,dY,dZ):
    """
    Konverterer fra ECEF til lokaltoposentrisk koordinatsystem ENU.
    """
    import numpy as np 
    from math import sqrt,sin,cos,tan,pi,atan,atan2,asin
    from numpy import fix,log,fmod,arctan,arctan2,array
    dP_ECEF = array([dX, dY, dZ]).reshape((3,1))
    
    M = array([[-sin(lon), cos(lon), 0], 
        [-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat)], 
        [cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)]])
    
    dP_ENU = M @ dP_ECEF
    
    e = float(dP_ENU[0]) 
    n = float(dP_ENU[1])
    u = float(dP_ENU[2])
    return dP_ENU


def atanc(y,x):
    import numpy as np 
    from math import sqrt,sin,cos,tan,pi,atan,atan2,asin
    from numpy import fix,log,fmod,arctan,arctan2
    z=atan2(y,x)
    atanc=fmod(2*pi + z, 2*pi)
    return atanc

def Nrad(a,b,lat):
    '''
    Funksjonen beregner Normalkrumningsradiusen for den gitte breddegraden. På engelsk
    "Earth's prime-vertical radius of curvature", eller "The Earth's transverse radius of curvature".
    Den står ortogonalt på M (meridiankrumningsradiusen) for den gitte breddegraden. Dvs øst-vest. 

    Parameters
    ----------
    a : Store halvakse
    b : Lille halvakse
    lat : Breddegrad

    Returns
    -------
    N : Normalkrumningsradiusen

    '''
    import numpy as np 
    from math import sqrt,sin,cos,tan,pi,atan,atan2,asin
    from numpy import fix,log,fmod,arctan,arctan2
    e2 = (a**2 - b**2)/a**2
    N = a/(1 - e2*sin(lat)**2)**(1/2)
    return N