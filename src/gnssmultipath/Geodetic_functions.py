from numpy import fix,array,log,fmod,arctan,arctan2,arcsin,sqrt,sin,cos,pi,arange
import numpy as np
from datetime import datetime,timedelta
from datetime import date
from typing import List, Union, Literal
from numpy import ndarray
import warnings
warnings.filterwarnings("ignore")

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
    e2m = (a**2 - b**2)/b**2
    e2  = (a**2 - b**2)/a**2
    rho = sqrt(X**2 +Y**2)
    my  = arctan((Z*a)/(rho*b))
    lat = arctan(( Z +e2m*b*(sin(my))**3)/(rho - e2*a*(cos(my))**3))
    lon = arctan2(Y,X)
    N   = Nrad(a,b,lat)
    h   = rho*cos(lat) + Z*sin(lat) - N*( 1 - e2*(sin(lat))**2)
    return lat, lon, h



def ECEF2enu(lat,lon,dX,dY,dZ): ## added this new function 28.01.2023
    """
    Convert from ECEF to a local toposentric coordinate system (ENU)
    """
    ## -- Ensure that longitude is bewtween -180 and 180
    if -2*pi < lon < - pi:
        lon = lon + 2*pi
    elif pi < lon < 2*pi:
        lon = lon - 2*pi

    ## -- Compute sin and cos before putting in to matrix to gain speed
    sin_lon = np.sin(lon)
    cos_lon = np.cos(lon)
    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)

    dP_ECEF = np.array([[dX, dY, dZ]]).T
    ## -- Propagate through rotation matrix
    M = np.array([[-sin_lon, cos_lon, 0],
        [-sin_lat*cos_lon, -sin_lat*sin_lon, cos_lat],
        [cos_lat*cos_lon, cos_lat*sin_lon, sin_lat]])

    dP_ENU = np.dot(M, dP_ECEF)

    e = dP_ENU[0, 0]
    n = dP_ENU[1, 0]
    u = dP_ENU[2, 0]

    return e, n, u






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

    e2 = (a**2 - b**2)/a**2
    N = a/(1 - e2*sin(lat)**2)**(1/2)
    return N



def compute_satellite_azimut_and_elevation_angle(X, Y, Z, xm, ym, zm):
    """
    Computes the satellites azimute and elevation angle based on satellitte and
    receiver ECEF-coordinates. Utilizes vectorization (no for loops) for
    better performance.

    Unit: Degree.

    Parameters
    ----------
    X : Satellite X-coordinate (numpy array)
    Y : Satellite Y-coordinate (numpy array)
    Z : Satellite Z-coordinate (numpy array)
    xm : Reciever X-coordinate (numpy array)
    ym : Reciever Y-coordinate (numpy array)
    zm : Reciever X-coordinate (numpy array)

    Returns
    -------
    az : Azimut angle in degrees  (numpy array)
    elev: Elevation angle in degrees  (numpy array)


    NOT IN USE AT THE MOMENT
    """

    ## -- WGS 84 ellipsoid:
    a   =  6378137.0         # semi-major ax
    b   =  6356752.314245    # semi minor ax

    # Compute latitude and longitude for the receiver
    lat,lon,h = ECEF2geodb(a,b,xm,ym,zm)

    # Find coordinate difference between satellite and receiver
    dX = (X - xm)
    dY = (Y - ym)
    dZ = (Z - zm)

    # Convert from ECEF to ENU (east,north, up)
    east, north, up = np.vectorize(ECEF2enu)(lat,lon,dX,dY,dZ)

    # Calculate azimuth angle and correct for quadrants
    azimuth = np.rad2deg(np.arctan(east/north))
    azimuth = np.where((east > 0) & (north < 0) | ((east < 0) & (north < 0)), azimuth + 180, azimuth)
    azimuth = np.where((east < 0) & (north > 0), azimuth + 360, azimuth)

    # Calculate elevation angle
    elevation = np.rad2deg(np.arctan(up / np.sqrt(east**2 + north**2)))
    elevation = np.where((elevation <= 0) | (elevation >= 90), np.nan, elevation) # Set elevation angle to NaN if not in the range (0, 90)

    return azimuth, elevation





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

    ## -- WGS 84 datumsparametre:
    a   =  6378137.0         # store halvakse
    b   =  6356752.314245    # lille halvakse

    ## -- Beregner bredde og lengdegrad til mottakeren:
    lat,lon,h = ECEF2geodb(a,b,xm,ym,zm)

    ## --Finner vektordifferansen:
    dX = (X - xm)
    dY = (Y - ym)
    dZ = (Z - zm)

    ## -- Transformerer koordinatene over til lokalttoposentrisk system:
    if isinstance(X, float): # if only float put in, not list or array
        east,north,up = ECEF2enu(lat,lon,dX,dY,dZ)
        ## -- Computes the azimut angle and elevation angel for current coordinates (in degrees)
        if (east > 0 and north < 0) or (east < 0 and north < 0):
            az = np.rad2deg(arctan(east/north)) + 180
        elif east < 0 and north > 0:
            az = np.rad2deg(arctan(east/north)) + 360
        else:
            az = np.rad2deg(arctan(east/north))
        # elev = arcsin(up/(sqrt(east**2 + north**2 + up**2)))*(180/pi)
        elev = np.rad2deg(atanc(up, sqrt(east**2 + north**2)))
        if not 0 < elev < 90: # if the satellite is below the horizon (elevation angle is below zero)
            elevation_angle = np.nan
    else:
        east = np.array([]); north = np.array([]); up = np.array([])
        for i in np.arange(0,len(dX)):
            east_,north_,up_ = ECEF2enu(lat,lon,dX[i],dY[i],dZ[i])
            east = np.append(east,east_)
            north = np.append(north,north_)
            up = np.append(up,up_)

        ## -- Computes the azimut angle and elevation angel for list  coordinates (in degrees)
        az = []; elev = []
        for p in np.arange(0,len(dX)):
            # # Kvadrantkorreksjon
            if (east[p]> 0 and north[p]< 0) or (east[p] < 0 and north[p] < 0):
                az.append(np.rad2deg(arctan(east[p]/north[p])) + 180)
            elif east[p] < 0 and north[p] > 0:
                az.append(np.rad2deg(arctan(east[p]/north[p])) + 360)
            else:
                az.append(np.rad2deg(arctan(east[p]/north[p])))
            # elev.append(arcsin(up[p]/(sqrt(east[p]**2 + north[p]**2 + up[p]**2)))*(180/pi))
            elev.append(np.rad2deg(atanc(up, sqrt(east**2 + north**2))))

    return az,elev


def atanc(y,x):
    z=arctan2(y,x)
    atanc=fmod(2*pi + z, 2*pi)
    return atanc


def date2gpstime(year,month,day,hour,minute,seconds):
    """
    Computing GPS-week nr.(integer) and "time-of-week" from year,month,day,hour,min,sec
    Origin for GPS-time is 06.01.1980 00:00:00 UTC
    """
    t0=date.toordinal(date(1980,1,6))+366
    t1=date.toordinal(date(year,month,day))+366
    week_flt = (t1-t0)/7
    week = fix(week_flt)
    tow_0 = (week_flt-week)*604800
    tow = tow_0 + hour*3600 + minute*60 + seconds

    return week, tow


def filter_array_on_system(eph_data,system):
    """
    Filtering the array with ephemerids based on system code "G", "E" "C" etc.
    """
    mask = np.char.startswith(eph_data[:,0],system)
    return eph_data[mask]


def extract_nav_message(data,PRN,tidspunkt):
    """
    Funksjonen samler de satellittene som har de spesifiserte PRN nummrene,
    samt de som ligger nærmest spesifisert tidspunkt.
    """
    samla_data = gathering_sat_by_PRN(data,PRN);
    Sat_liste = find_message_closest_in_time(samla_data, tidspunkt)
    return Sat_liste


def gathering_sat_by_PRN(data,PRN):
    """
    Funksjonen bruker rinex-data i form av array ordnet ved funksjonen read_rinex2_nav eller read_rines3_nav
    """
    m = len(data)
    Sat_data  = np.zeros((1,36))
    for k in np.arange(0,m):
        PRN_ = np.array([])

        if len(data[k,0]) == 3: # RINEX v3 navfiles have letters in addition. EX: G01.
            sat = data[k,0][1::]
        else:
            sat = data[k,0]

        if int(sat) == PRN:
            PRN_ =np.append(PRN_,data[k,:])
            PRN_ = PRN_.reshape(1,len(PRN_))
        if np.size(PRN_) != 0:
            Sat_data  = np.concatenate([Sat_data , PRN_], axis=0)

    if all(Sat_data[0,:].astype(float)) == 0:
        Sat_data  = np.delete(Sat_data , (0), axis=0)

    return Sat_data


def find_message_closest_in_time(data,tow_mot):
    """
    Funksjonen plukker ut den linjen i et datasett som ligg nærmest tidspunktet vi ønsker i bestemme koordinatene for.
    """
    towSat = [] # Tom liste for å lagre tidspunktene
    length = len(data)
    towSat = np.array([])
    for i in np.arange(0,length):
        week, tow = date2gpstime(2000 + int(data[i,1]), int(data[i,2]), int(data[i,3]), int(data[i,4]), int(data[i,5]), int(data[i,6]))
        towSat = np.append(towSat,tow)
    # Finner verdien som ligger nærmest
    # index = int(min(abs(tow_mot - towSat)))
    index = np.abs(tow_mot - towSat).argmin()
    GNSS_linjer = data[index,:]

    return GNSS_linjer




def Satkoord2(efemerider,tow_mot,xm,ym,zm):
    """
    Funksjonen beregner satellittkordinater og korrigerer for jordrotasjonen. Fra keplerelemtenter til ECEF.
    """

    GM         = 3.986005e14      # Produktet av jordas masse og gravitasjonskonstanten
    omega_e    = 7.2921151467e-5  # [rad/sek]
    c          = 299792458        # Lyshastigheten [m/s]


    #Leser inn data:
    M0         = efemerider[13]
    delta_n    = efemerider[12]
    e          = efemerider[15]
    A          = efemerider[17]**2
    OMEGA      = efemerider[20]
    i0         = efemerider[22]
    omega      = efemerider[24]
    OMEGA_dot  = efemerider[25]
    i_dot      = efemerider[26]
    C_uc       = efemerider[14]
    C_us       = efemerider[16]
    C_rc       = efemerider[23]
    C_rs       = efemerider[11]
    C_ic       = efemerider[19]
    C_is       = efemerider[21]
    toe        = efemerider[18]

    n0  = sqrt(GM/A**3)   #(rad/s)
    t_k = tow_mot - toe

    n_k = n0 + delta_n   #Koorigert  midlere bevegelse
    M_k = M0 + n_k*t_k   #Midlere anomali (rad/s)


    #Beregner eksentrisk anomali
    E_old = M_k
    for i in arange(0,10):
        E = M_k+ e*sin(E_old)
        dE = fmod(E - E_old, 2*pi)
        if abs(dE) < 1.e-12:
           break
        E_old = E

    E = fmod(E+2*pi,2*pi)

    cosv = (cos(E) - e)/(1 - e*cos(E))
    sinv = (sqrt(1 - e**2)*sin(E))/(1-e*cos(E))
    tanv = sinv/cosv

    ## -- Kvadrantkorreksjon
    if sinv > 0 and cosv < 0 or sinv < 0 and cosv < 0:
        v = arctan(tanv) + pi
    elif sinv < 0 and cosv > 0:
        v = arctan(tanv) + 2*pi
    else:
        v = arctan(tanv)


    theta   = v + omega

    ## -- Korreksjon på baneparametere:
    du_k    = C_uc*cos(2*theta) + C_us*sin(2*theta)
    dr_k    = C_rc*cos(2*theta) + C_rs*sin(2*theta)
    di_k    = C_ic*cos(2*theta) + C_is*sin(2*theta)

    ## -- Korrigerte baneparametere:
    u_k     = theta + du_k
    r_k     = A*(1 - e*cos(E)) + dr_k
    i_k     = i0 + i_dot*t_k + di_k


    ## -- Korrigert lengdegrad for baneplanknuten:
    OMEGA_k = OMEGA + (OMEGA_dot - omega_e)*t_k - omega_e*toe


    ## -- Satelittens posisjon i banen:
    x = r_k*cos(u_k)
    y = r_k*sin(u_k)

    ## --ECEF-koordinater for satellitten:
    X = x*cos(OMEGA_k) - y*sin(OMEGA_k)*cos(i_k)
    Y = x*sin(OMEGA_k) + y*cos(OMEGA_k)*cos(i_k)
    Z = y*sin(i_k)

    ## Relativistisk klokkekorreksjon
    dT_rel = (-2/c**2)*sqrt(A*GM)*e*sin(E)

    ## --Tester om korreksjonen skal utføres.
    if (abs(xm)) > 1.0 and (abs(ym)) > 1.0 and (abs(zm)) > 1.0:
        TRANS  = 0
        TRANS0 = 0.075
        j = 0
        while(abs(TRANS0 - TRANS) > 1e-10):
            j = j +1
            if (j > 20):
                print('Feil, Gangtids-rotasjonen konvergerer ikke!')
                break

            TRANS = TRANS0
            OMEGA_k = OMEGA + (OMEGA_dot - omega_e)*t_k - omega_e*(toe + TRANS)
            X = x*cos(OMEGA_k) - y*sin(OMEGA_k)*cos(i_k)
            Y = x*sin(OMEGA_k) + y*cos(OMEGA_k)*cos(i_k)
            Z = y*sin(i_k)
            #Regner ut avstanden mottaker-satellitt
            dX = (X - xm)
            dY = (Y - ym)
            dZ = (Z - zm)
            DS = sqrt(dX**2 + dY**2 + dZ**2)
            #Regner ny verdi for signalets gangtid
            TRANS0 = DS/c

    else:
        #Jordrotasjonen skal ikke utføres
        OMEGA_k = OMEGA + (OMEGA_dot - omega_e)*t_k - omega_e*toe
        X = x*cos(OMEGA_k) - y*sin(OMEGA_k)*cos(i_k)
        Y = x*sin(OMEGA_k) + y*cos(OMEGA_k)*cos(i_k)
        Z = y*sin(i_k)

    return X,Y,Z,dT_rel



def compute_GLO_coord_from_nav(ephemerides, time_epochs):
    """
    Function that use broadcast ephemerides (from rinex nav file) and interpolates to current time using 4th order Runge-kutta.

    Parameters:

    ephemerides: A array with ephemerides for current satellite
    time_epochs: An n_obs X 2 sized array with weeks and tow read from the rinex observation file.

    Return:


    """
    ephemerides[0] = 9999 #reming char from string (ex G01 -> 9999)
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
    # time = datetime.datetime.strptime(time, "%Y/%m/%d %H:%M:%S.%f")
    time = datetime.strptime(time, "%Y/%m/%d %H:%M:%S.%f")
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
    der_state = np.zeros(6)
    if r**2 < 0:
        return der_state
    a = 1.5 * J2 * mu * (ae**2)/ (r**5)
    b = 5 * (state[2]**2) / (r**2)
    c = -mu/(r**3) - a*(1-b)

    der_state[0:3] = state[3:6]
    der_state[3] = (c + omega**2)*state[0] + 2*omega*state[4] + acc[0]
    der_state[4] = (c + omega**2)*state[1] - 2*omega*state[3] + acc[1]
    der_state[5] = (c - 2*a)*state[2] + acc[2]
    return der_state





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
    t0 = datetime(1980,1,6)
    t1 = t0 + timedelta(days=(days_to_start_of_week + days_from_hours))

    ## --  Formating the date to "year-month- day"
    t1 = t1.strftime("%Y %m %d")
    t1_ = [int(i) for i in t1.split(" ")]

    [year, month, day] = t1_

    date_ = [year, month, day, hour, min_, sec]
    return date_



def gpstime2date_arrays(week: Union[List[int], np.ndarray], tow: Union[List[float], np.ndarray]) -> np.ndarray:
    """
    Calculates date from GPS-week number and "time-of-week" to Gregorian calendar.

    Example:
    -------
    .. code-block:: python

    time_epochs = array([[  2190.        , 518399.99999988],
                         [  2190.        , 518429.99999988],
                         [  2190.        , 531539.99999988],
                         [  2190.        , 531569.99999988]])

    greg_time = gpstime2date_arrays(time_epochs[:,0], time_epochs[:,1])

    Parameters
    ----------
    week : array/list, GPS-week numbers.
    tow  : array/list, "Time of week" values.

    Returns
    -------
    dates : array, An array of dates given in the Gregorian calendar ([year, month, day, hour, min, sec]).
    """
    # Convert week and tow to regular Python lists
    week = week.tolist()
    tow = tow.tolist()
    total_seconds = [w * 7 * 24 * 3600 + t for w, t in zip(week, tow)] # Calculate the time since the GPS epoch (6th January 1980) for each week and tow
    delta = [timedelta(seconds=s) for s in total_seconds] # Calculate the timedelta from the GPS epoch for each week and tow
    gps_epoch = datetime(1980, 1, 6) # The GPS reference epoch
    dates = [gps_epoch + d for d in delta] # Calculate the final date by adding the timedelta to the GPS epoch
    date_array = np.array([[date.year, date.month, date.day, date.hour, date.minute, date.second] for date in dates]) # Create an array by extracting date components

    return date_array


def gpstime2date_arrays_with_microsec(week: Union[List[int], np.ndarray], tow: Union[List[float], np.ndarray]) -> np.ndarray:
    """
    Calculates date from GPS-week number and "time-of-week" to Gregorian calendar with microsecond precision.

    Example:
    -------
    .. code-block:: python

    time_epochs = array([[  2190.        , 518399.99999988],
                         [  2190.        , 518429.99999988],
                         [  2190.        , 531539.99999988],
                         [  2190.        , 531569.99999988]])

    greg_time = gpstime2date_arrays(time_epochs[:,0], time_epochs[:,1])

    Parameters
    ----------
    week : array/list, GPS-week numbers.
    tow  : array/list, "Time of week" values.

    Returns
    -------
    dates : array, An array of dates given in the Gregorian calendar ([year, month, day, hour, min, sec, microsec]).
    """
    # Convert week and tow to NumPy arrays for efficient processing
    week = np.asarray(week, dtype=int)
    tow = np.asarray(tow, dtype=float)

    # Calculate the total seconds since GPS epoch
    total_seconds = week * 7 * 24 * 3600 + tow

    # Create a timedelta object for each time in total_seconds
    gps_epoch = datetime(1980, 1, 6)  # The GPS reference epoch
    deltas = [timedelta(seconds=s) for s in total_seconds]

    # Calculate final dates and extract components
    dates = [gps_epoch + delta for delta in deltas]
    date_array = np.array([
        [date.year, date.month, date.day, date.hour, date.minute, date.second, date.microsecond]
        for date in dates
    ])

    return date_array


def convert_to_datetime_vectorized(time_datetime: np.ndarray) -> list:
    """
    Convert numpy array of datetime components to datetime strings with specified format.
    """
    # Convert datetime components to datetime objects
    datetimes = [datetime(*components) for components in time_datetime]
    # Format datetime objects to strings
    datetime_strings = np.array([dt.strftime("%Y-%m-%d %H:%M:%S") for dt in datetimes]).tolist()
    return datetime_strings


def gpstime_to_utc_datefmt(time_epochs_gpstime: np.ndarray) -> list:
    """
    Coverters form GPS time to UTC with formatting.
    Ex output: "2022-01-01 02:23:30".

    """
    time_datetime = gpstime2date_arrays(*time_epochs_gpstime.T)
    return convert_to_datetime_vectorized(time_datetime)



def date2gpstime_vectorized(gregorian_date_array):
    """
    Convert from gregorian dates to GPS time by using
    np.vectorization.
    """
    years, months, days, hours, minutes, seconds = gregorian_date_array.astype(float).astype(int).T
    weeks, tows = np.vectorize(date2gpstime)(years, months, days, hours, minutes, seconds)
    return np.round(weeks), np.round(tows)



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





if __name__=="__main__":
    time_epochs = array([[  2190.        , 518401.40],
                         [  2190.        , 531570.00]])


    date_array = gpstime2date_arrays_with_microsec(time_epochs[:,0], time_epochs[:,1])



