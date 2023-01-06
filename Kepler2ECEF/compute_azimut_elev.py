from Geodetic_functions import *

def compute_azimut_elev(X,Y,Z,xm,ym,zm):
    """
    Computes the satellites azimute and elevation angel based on satellitte and
    reciever ECEF-coordinates. Can take both single coordinates (float) and list(arrays). 
    
    Unit: Degree.

    Parameters
    ----------
    X : Satellite X-coordinate (float or array)
    Y : Satellite Y-coordinate (float or array)
    Z : Satellite Z-coordinate (float or array)
    xm : Reciever X-coordinate
    ym : Reciever Y-coordinate
    zm : Reciever X-coordinate

    Returns
    -------
    az: Azimut in degrees
    elev: Elevation angel in degrees
    """
    from math import sin,cos,tan,asin,acos,atan,pi
    import numpy as np

    ## -- WGS 84 datumsparametre:
    a   =  6378137.0         # store halvakse
    b   =  6356752.314245   # lille halvakse
    f   =  (a -b)/a         # flattrykning
    e2  = (a**2 - b**2)/a**2   # eksentrisitet
    
    ## -- Beregner bredde og lengdegrad til mottakeren:
    lat,lon,h = ECEF2geodb(a,b,xm,ym,zm)
    
    ## --Finner vektordifferansen:
    dX = (X - xm)
    dY = (Y - ym)
    dZ = (Z - zm)
    
    ## -- Transformerer koordinatene over til lokalttoposentrisk system:
    # east = []; north = []; up = []
    if X.shape == (): # if only float put in, not list or array
        east,north,up = ECEF2enu(lat,lon,dX,dY,dZ)
        ## -- Computes the azimut angle and elevation angel for current coordinates (in degrees)
        if (east > 0 and north < 0) or (east < 0 and north < 0):
            az = (atan(east/north)*(180/pi) + 180)
        elif east < 0 and north > 0:
           az = atan(east/north)*(180/pi) + 360
        else:
            az = atan(east/north)*(180/pi)
        elev = asin(up/(sqrt(east**2 + north**2 + up**2)))*(180/pi)
    else:
        east = np.array([]); north = np.array([]); up = np.array([])
        for i in range(0,len(dX)):    
            east_,north_,up_ = ECEF2enu(lat,lon,dX[i],dY[i],dZ[i])
            east = np.append(east,east_)
            north = np.append(north,north_)
            up = np.append(up,up_)
      
        ## -- Computes the azimut angle and elevation angel for list  coordinates (in degrees)
        az = []; elev = []
        for p in range(0,len(dX)):
            # Kvadrantkorreksjon 
            if (east[p]> 0 and north[p]< 0) or (east[p] < 0 and north[p] < 0):
                az.append(atan(east[p]/north[p])*(180/pi) + 180)
            elif east[p] < 0 and north[p] > 0:
               az.append(atan(east[p]/north[p])*(180/pi) + 360)
            else:
                az.append(atan(east[p]/north[p])*(180/pi))
            elev.append(asin(up[p]/(sqrt(east[p]**2 + north[p]**2 + up[p]**2)))*(180/pi)) 

    return az,elev

