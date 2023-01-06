def Satkoord2_GLO(ephemerides,t_rec,xm,ym,zm):
    """
    Funksjonen beregner satellittkordinater og korrigerer for jordrotasjonen. From state vector to ECEF.
    See: https://gssc.esa.int/navipedia/index.php/GLONASS_Satellite_Coordinates_Computation
    """
    
    from numpy import sqrt,pi, fmod,cos, sin, tan
    from math import atan

    ## Parameters (from GLONASS Interface Control Document 1988)
    GM         = 398600.44e9     # Gravitational constant [m3/s2]   (product of the mass of the earth and and gravity constant)
    omega_e    = 7.292115e-5     # Earth rotation rate    [rad/sek]
    c          = 299792458       # Speed of light         [m/s]
    a          = 6378136         # Semi major axis PZ-90   [m]
    f          = 1/298.257839303 # Inverse flattening
    C_20       = -1082.63e-6     # Second zonal coefficient of spherical harmonic expression.


    ## Read in data:
    toc = [ephemerides[1],ephemerides[2] ,ephemerides[3] ,ephemerides[4],ephemerides[5],ephemerides[6]] # year,month,day,hour,minute,second     
        
    toe  = ephemerides[9]   # Time of ephemerides (i think??). Time of week in seconds (UTC)    
    t_k  = t_rec - toe      # Difference between reciever time and ephemerides time.     
    
    x_te = ephemerides[10]  # X-coordinates at t_e in PZ-90 [km]
    y_te = ephemerides[14]  # Y-coordinates at t_e in PZ-90 [km]
    z_te = ephemerides[18]  # Z-coordinates at t_e in PZ-90 [km]
    
    vx_te = ephemerides[11] # Velocity component at t_e in PZ-90 (v_x) [km/s]
    vy_te = ephemerides[15] # Velocity component at t_e in PZ-90 (v_y) [km/s]
    vz_te = ephemerides[19] # Velocity component at t_e in PZ-90 (v_z) [km/s]
    
    J_x = ephemerides[12]   # Moon and sun acceleration at t_e [km/sec**2]
    J_y = ephemerides[16]   # Moon and sun acceleration at t_e [km/sec**2]
    J_z = ephemerides[20]   # Moon and sun acceleration at t_e [km/sec**2]
    
    k_num = ephemerides[17] # Frequency slot/number [-7..+13]
    ## -- Transformation to an intertial reference frame
    theta_g0 = siderial_time(toc[0], toc[1], toc[2], utc=0, long=0) ### RIktig?? MÅ gjøre om til sekunder!!!!
    theta_ge = theta_g0 + omega_e*(toe - 3) # - 3 hours
    ## -- Transformation to an inertial reference frame:
    ## Position:
    x_a = x_te*cos(theta_ge) - y_te*sin(theta_ge)
    y_a = x_te*sin(theta_ge) + y_te*cos(theta_ge)
    z_a = z_te 
    ## Velocity:
    v_xa = vx_te*cos(theta_ge) - vy_te*sin(theta_ge) - omega_e*y_te
    v_ya = vx_te*cos(theta_ge) + vy_te*cos(theta_ge) + omega_e*x_te
    v_za = vz_te
    
    ## -- Defining parameters (OL-"overline")
    r = sqrt(x_a**2 + y_a**2 + z_a**2)
    GM_OL = GM/r**2
    x_a_OL = x_a/r
    y_a_OL = y_a/r
    z_a_OL = z_a/r
    rho_OL = a/r
    r = sqrt(x_a**2 + y_a**2 + z_a**2)
    
    ## --Defining the differential equations
    dx = vx_te
    dy = vy_te
    dx = vz_te
    
    dvx = -GM/r**3*x_te + 3/2*C_20*GM*a**2/r**5*x_te*(1 - 5*z_te**2/r**2) + J_x + omega_e**2*x_te + 2*omega_e*vy_te
    dvy = -GM/r**3*y_te + 3/2*C_20*GM*a**2/r**5*y_te*(1 - 5*z_te**2/r**2) + J_y + omega_e**2*y_te - 2*omega_e*vx_te
    dvz = -GM/r**3*z_te + 3/2*C_20*GM*a**2/r**5*z_te*(3 - 5*z_te**2/r**2) + J_z    
    
    
    
    return 




def f1(vx_te):    
    return vx_te

def f2(vy_te):    
    return vy_te

def f3(vz_te):    
    return vz_te

def f4()


def julian_date(year, month, day, utc=0):
    """
    Returns the Julian date, number of days since 1 January 4713 BC 12:00.
    utc is UTC in decimal hours. If utc=0, returns the date at 12:00 UTC.
    """
    if month > 2:
        y = year
        m = month
    else:
        y = year - 1
        m = month + 12
    d = day
    h = utc/24
    if year <= 1582 and month <= 10 and day <= 4:
        # Julian calendar
        b = 0
    elif year == 1582 and month == 10 and day > 4 and day < 15:
        # Gregorian calendar reform: 10 days (5 to 14 October 1582) were skipped.
        # In 1582 after 4 October follows the 15 October.
        d = 15
        b = -10
    else:
        # Gregorian Calendar
        a = int(y/100)
        b = 2 - a + int(a/4)
    jd = int(365.25*(y+4716)) + int(30.6001*(m+1)) + d + h + b - 1524.5
    return(jd)


def siderial_time(year, month, day, utc=0, long=0):
    """
    Returns the siderial time in decimal hours. Longitude (long) is in 
    decimal degrees. If long=0, return value is Greenwich Mean Siderial Time 
    (GMST).
    
    https://www.nies.ch/doc/astro/sternzeit.en.php
    """
    jd = julian_date(year, month, day)
    t = (jd - 2451545.0)/36525
    # Greenwich siderial time at 0h UTC (hours)
    st = (24110.54841 + 8640184.812866 * t +
          0.093104 * t**2 - 0.0000062 * t**3) / 3600
    # Greenwich siderial time at given UTC
    st = st + 1.00273790935*utc
    # Local siderial time at given UTC (longitude in degrees)
    st = st + long/15
    st = st % 24
    return(st)

def GMST(year,month,day,hour,minute,second):
    """
    Funtions computes the sideral time at 0h (UT)
    http://www.astro.sunysb.edu/metchev/AST443/times.html
    """
    import pandas as pd
    ts = pd.Timestamp(year = year,  month = month, day = day, hour = hour, minute= minute, second = second, tz = 'Europe/Berlin')
    JD = ts.to_julian_date()
    # GMST (Greenwich Mean Sidereal Time) at 0h UT = 24110.54841 + 8640184.812866 TU + 0.093104TU2 - 6.2x10-6TU3
    T_U = (JD -2451545.0)/36525 #  is the number of Julian Centuries since J2000.0
    GMST_0 = 24110.54841 + 8640184.812866*T_U + 0.093104*T_U**2 - 6.2e10-6*T_U**3
    return GMST_0
    



#########################


