
'''
@author: Per Helge Aarnes
'''

from math import sqrt,sin,cos,tan,pi,atan,atan2,asin
from numpy import fix,array,log,fmod,arctan,arctan2

def Marc(a,b,lat):
    '''
    Funksjonen beregner buelengden langs meridianen. Meridian-buelengden er altså 
    avstanden mellom to punkter som har samme lengdegrad.

    Parameters
    ----------
    a : Store halvakse
    b : Lille halvakse
    lat : Breddegrad

    Returns
    -------
    B : Meridianbuelengden

    '''    
    f = (a - b)/a
    b0 = a*(1 - 1/2*f + 1/16*f**2 + 1/32*f**3)

    B = b0*(lat - (3/4*f + 3/8*f**2 + 15/128*f**3)*sin(2*lat)
            + (15/64*f**2 + 15/64*f**3)*sin(4*lat)
            - 35/384*f**3*sin(6*lat))
    return B




def Mrad(a,b,lat):
    '''
    Funksjonen beregner merdidiankrumningsradiusen for den gitte
    breddegraden. Altså retning nord-sør. På engelsk:  "Earth's meridional radius of curvature"

    Parameters
    ----------
    a : Store halvakse
    b : Lille halvakse
    lat : Breddegrad

    Returns
    -------
    M : Meridiankrumningsradius 

    '''
    e2 = (a**2 - b**2)/a**2
    M = a*(1 - e2)/(1 - e2*sin(lat)**2)**(3/2)
    
    return M




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



def arctanc(y, x):
    z = arctan2(y, x)
    return fmod(2*pi + z, 2*pi)



def geod2TMgrid(a,b,lat,lon,lat0,lon0,scale,fnorth,feast):
    '''
    Funksjonen konverterer geodetiske koordinater til kartplanet (TM) ved hjelp av taylor rekkeutvikling. 
    Tap av presisjon er noen tidels mm for hver konvertering (neglisjerbart i de fleste tilfeller)

    Parameters
    ----------
    a : store halvakse til ellipsoiden
    b : lille halvakse til ellipsoiden
    lat : breddegrad
    lon : lengdegrad
    lat0 : startbreddegrad, 0 gir ekvator (eks 58 grader for NGO 1948)
    lon0 : startlengdegrad,lengdegraden til tangeringsmeredianen til \n
    TMgriden. Eks (9 grader for UTM sone 32)
    scale : skaleringsfaktor, (0.9996 for UTM)
    fnorth : Falsk nord i meter (0 for UTM)
    feast : Falsk øst i meter (500 000 for UTM)

    Returns
    -------
    north : Nordkoordinaten gitt i en sone
    east : Østkoordinaten gitt i en sone

    '''
    
    B =Marc(a,b,lat) - Marc(a,b,lat0)
    e2 = (a**2 - b**2)/a**2
    eps2 = e2/(1 - e2)*cos(lat)**2
    N = Nrad(a,b,lat)
    l = lon - lon0
    
    x = B + 1/2*l**2*N*sin(lat)*cos(lat) \
        + 1/24*l**4*N*sin(lat)*cos(lat)**3*(5 - tan(lat)**2 + 9*eps2 + 4*eps2**2) \
        + 1/720*l**6*N*sin(lat)*cos(lat)**5*(61 - 58*tan(lat)**2 + tan(lat)**4)

    y = l*N*cos(lat) + 1/6*l**3*N*cos(lat)**3*(1 - tan(lat)**2 + eps2) \
        + 1/120*l**5*N*cos(lat)**5*(5 - 18*tan(lat)**2 + tan(lat)**4)
    
    north = x*scale
    east = y*scale
    
    north = north + fnorth
    east  = east + feast
    
    return north, east




def TMgrid2geod(a,b,north,east,lat0,lon0,scale,fnorth,feast):
    '''
    Funksjonen konverterer kartplan-koordinater til geodetisk bredde- og lengdegrad. 
    Tap av presisjon er noen tidels mm for hver konvertering (neglisjerbart i de fleste tilfeller)

    Parameters
    ----------
    a : store halvakse til ellipsoiden
    b : lille halvakse til ellipsoiden
    north : Nordkoordinaten gitt i en sone
    east : Østkoordinaten gitt i en sone
    lat0 : startbreddegrad, 0 gir ekvator (eks 58 grader for NGO 1948)
    lon0 : startlengdegrad,lengdegraden til tangeringsmeredianen til \n
    TMgriden. Eks (9 grader for UTM sone 32)
    scale : skaleringsfaktor, (0.9996 for UTM)
    fnorth : Falsk nord i meter (0 for UTM)
    feast : Falsk øst i meter (500 000 for UTM)

    Returns
    -------
    lat : Breddegrad
    lon : Lengdegrad

    '''

    x = (north-fnorth)/scale
    y = (east-feast)/scale
    latf = footlat(a,b,x,lat0)    
    M = Mrad(a,b,latf)    
    N = Nrad(a,b,latf)    
    e2 = (a**2-b**2)/a**2    
    E = e2/(1 - e2)*cos(latf)**2   
    
        
    lat = latf - 1/2*y**2*tan(latf)/(M*N) \
          + 1/24*y**4*tan(latf)/(M*N**3)*(5 + 3*tan(latf)**2 + E - 9*E*tan(latf)**2 - 4*E**2) \
          - 1/720*y**6*tan(latf)/(M*N**5)*(61 + 90*tan(latf)**2 + 45*tan(latf)**4)

    lon = y/(N*cos(latf)) \
        - 1/6*y**3/(N**3*cos(latf))*(1 + 2*tan(latf)**2 + E) \
        + 1/120*y**5/(N**5*cos(latf))*(5 + 28*tan(latf)**2 + 24*tan(latf)**4) + lon0
        
        
    return lat, lon
          


def deg2rad(deg):
    '''
    Konverterer fra grader til radianser

    Parameters
    ----------
    deg : Vinkel i grader

    Returns
    -------
    rad : Vinkel i radianser

    '''
    rad = deg*(pi/180)
    return rad



def footlat(a,b,x,lat0):
    '''
    Beregner fotpunktsbredden som er meridian-buelengden 
    mellom ekvator og målestasjonspunktet (kun en mellomregningsstørrelse)

    Parameters
    ----------
    a : Store halvakse
    b : Lillehalvakse
    x : X-koordinat (Nord)
    lat0 : Nullbredde

    Returns
    -------
    latf : Fotpunktsbredde

    '''
    f = (a - b)/a
    b0 = a*(1 - 1/2*f + 1/16*f**2 + 1/32*f**3)

    B = Marc(a, b, lat0) + x

    latf = B/b0 + (3/4*f + 3/8*f**2 + 21/256*f**3)*sin(2*B/b0) \
           + (21/64*f**2 + 21/64*f**3)*sin(4*B/b0) \
           + 151/768*f**3*sin(6*B/b0)
       
    return latf






def geod2ECEF(a,b,lat,lon,h):
    '''
    Funksjonen tar inn geodetisk bredde og lengde og konverterer til 
    kartesiske koordinater(ECEF-jordsentrisk og jordfast koordinatsystem)

    Parameters
    ----------
    a   : Store halvakse
    b   : lille halvakse
    lat : Breddegrad
    lon : Lengdegrad
    h   : Høyde

    Returns
    -------
    X : X-koordinat
    Y : Y-koordinat
    Z : Z-koordinat

    '''
    
    N = Nrad(a,b,lat)
    X = (N+h)*cos(lat)*cos(lon)
    Y = (N+h)*cos(lat)*sin(lon) 
    Z = ((b**2/a**2)*N + h)*sin(lat)
    return X, Y, Z



def rad2deg(rad):
    '''
    Konvererer fra radianer til grader

    Parameters
    ----------
    rad : Radianer

    Returns
    -------
    deg : Grader

    '''
    deg = rad*(180/pi)
    return deg
    

def dms2deg(d,m,s):
    '''
    Konverterer fra grader,minutter og sekunder til desimalgrader
    Parameters
    ----------
    d : Degrees (grader)
    m : Minutes (minutter)
    s : Seconds (sekunder)

    Returns
    -------
    deg : Desimalgrader

    '''
    deg = abs(d) + m/60 + s/3600
    return deg



def dms2rad(d,m,s):
    '''
    Konverterer fra grader,minutter og sekunder til desimalradianer
    
    Parameters
    ----------
    d : Degrees (grader)
    m : Minutes (minutter)
    s : Seconds (sekunder)

    Returns
    -------
    deg : Desimalradianer

    '''
    deg = dms2deg(d,m,s)
    rad = deg2rad(deg)
    return rad


def deg2dms(deg):
    '''
    Konverterer fra desimalradianer til grader,minutter og sekunder.
    
    Parameters
    ----------
    deg : Desimalradianer

    Returns
    -------
    d : Degrees (grader)
    m : Minutes (minutter)
    s : Seconds (sekunder)

    '''
    if deg > 0:
        d = fix(deg) #bruker fix-funksjonen for å få kun heltall.
        m = fix(abs((deg -d)*100)*0.6) # Bruker if- for at de negative tallene skal bli rett. 
        s = abs(((deg-d-(m/60))*3600))
    else:   
        d = fix(deg)
        m = fix(abs((deg - d)*100)*0.6)
        s = abs(((deg - d + (m/60))*3600))
    return d, m, s


def rad2dms(rad): 
    '''
    Konverterer fra desimalradianer til grader,minutter og sekunder.
    
    Parameters
    ----------
    deg : Desimalradianer

    Returns
    -------
    d : Degrees (grader)
    m : Minutes (minutter)
    s : Seconds (sekunder)

    '''
    [deg] = rad2deg(rad)
    [d,m,s] = deg2dms(deg)
    return



def ECEF2geod(a, b, X, Y, Z):
    '''
    Konverter fra kartesiske ECEF-koordinater til geodetiske koordinater vha interasjon.

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
    N = None
    e2 = (a**2 - b**2)/a**2
    p = sqrt(X**2 + Y**2)
    lat_new = arctan(Z/p)
    epsilon = 1e-10
    lat = 0
    while abs(lat_new - lat) > epsilon:
        lat = lat_new
        N = Nrad(a, b, lat)
        lat_new = arctan(Z/p + N*e2*sin(lat)/p)

    lat = lat_new
    lon = arctan(Y/X)
    h = p*cos(lat) + Z*sin(lat) - N*(1 - e2*sin(lat)**2)

    return lat, lon, h


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
    my  = atan((Z*a)/(rho*b))
    lat = atan(( Z +e2m*b*(sin(my))**3)/(rho - e2*a*(cos(my))**3))
    lon = atan(Y/X)
    N   = Nrad(a,b,lat)
    h   = rho*cos(lat) + Z*sin(lat) - N*( 1 - e2*(sin(lat))**2)
    return lat, lon, h





def ECEF2geodv(a, b, X, Y, Z):
    """
    Konverter fra kartesiske ECEF-koordinater til geodetiske koordinater vha H. Vermeilles metode.
    
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
    """

    e2 = (a**2 - b**2)/a**2
    p = (X**2 + Y**2)/a**2
    q = (1 - e2)/a**2*Z**2
    r = (p + q - e2**2)/6
    s = e2**2*(p*q)/(4*r**3)
    t = (1 + s + sqrt(s*(2 + s)))**(1/3)
    u = r*(1 + t + 1/t)
    v = sqrt(u**2 + e2**2*q)
    w = e2*(u + v - q)/(2*v)
    k = sqrt(u + v + w**2) - w
    D = k*sqrt(X**2 + Y**2)/(k + e2)

    lat = 2*arctan(Z/(D + sqrt(D**2 + Z**2)))
    lon = arctan(Y/X)
    h = (k + e2 - 1)/k*sqrt(D**2 + Z**2)

    return lat, lon, h




def ECEF2enu(lat,lon,dX,dY,dZ):
    """
    Konverterer fra ECEF til lokaltoposentrisk koordinatsystem ENU.
    """

    dP_ECEF = array([dX, dY, dZ]).reshape((3,1))
    
    M = array([[-sin(lon), cos(lon), 0], 
        [-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat)], 
        [cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)]])
    
    dP_ENU = M @ dP_ECEF
    
    e = float(dP_ENU[0]) 
    n = float(dP_ENU[1])
    u = float(dP_ENU[2])
    return e, n, u


def TMcorr(a, b, x1, y1, x2, y2, lat0):
    '''
    Funksjonen beregner kartprojektsjonskorreksjoner for avstand og retning (ds og daz)
    Parameters
    ----------
    a : Store halvakse
    b : Lille halvakse
    x1 : X-koordinat punk 1
    y1 : Y-koordinat punk 1
    x2 : X-koordinat punk 2
    y2 : Y-koordinat punk 2
    lat0 : Nullbredde(startbredde)

    Returns
    -------
    daz : Korreksjon asimut
    ds : Korreksjon avstand
    '''
    latf = footlat(a, b, x1, lat0)

    az1 = (atan2((y2-y1), (x2-x1)))
    Ra = Rrad(a, b, latf, az1)

    s = sqrt((x2-x1)**2+(y2-y1)**2)

    daz = -((x2-x1)*(2*y1+y2))/(6*Ra**2)
    ds = (s/(6*Ra**2))*(y1**2+y1*y2+y2**2)
    return daz, ds


def TMconv(a, b, x, y, lat0):
    '''
    Funksjonen beregner meridiankonvergensen "gamma" i et punkt.
    Parameters
    ----------
    a : Store halvakse
    b : Lille halvakse
    x : x-koordinat
    y : y-koordinat
    lat0 : Nullbredde(startbredde) (eks 58 grader for NGO1948)
    Returns
    -------
    gamma : Meridiankonvergens (Grid-convergence) i punktet (x,y)

    '''
    e2 = (a**2 - b**2) / a**2
    latf = footlat(a, b, x, lat0)
    N = Nrad(a, b, latf)
    eps2 = (e2/(1-e2))*cos(latf)**2
    gamma = ((y*tan(latf)/N) - ((y**3*tan(latf))/(3*N**3)*(1 + (tan(latf))**2 - eps2 - 2*eps2**2)))
    return gamma



def TMscale(a, b, x, y, lat0):
    '''

    Parameters
    ----------
    a : Store halvakse
    b : Lille halvakse
    x : X-koordinat 
    y : Y-koordinat 

    Returns
    -------
    m : Skaleringsfaktor

    '''
    latf = footlat(a, b, x, lat0)
    M = Mrad(a, b, latf)
    N = Nrad(a, b, latf)
    m = 1 + (y**2)/(2*M*N)
    return m


def Rrad(a,b,lat,az1): 
    '''
    Funksjonen beregner krumningsradiusen i sikteretningen vha Eulers likning
    Parameters
    ----------
    a : Store halvakse
    b : Lille halvakse
    lat : Breddegrad
    az1 : Retningsvinkel (asimut)

    Returns
    -------
    Ra : Krumningsradiusen i den aktuelle retningen (asimuten)
    '''
    M = Mrad(a,b,lat)
    N = Nrad(a,b,lat)
    Ra = (M*N)/((M*sin(az1)**2)+(N*cos(az1)**2))
    return Ra

def red2ell2(d2,h1,h2):
    '''
    Reduserer avstanden fra korde til ellipsoiden
    
    Parameters
    ----------
    d2 : Avstand (korde)
    h1 : Høyde 1
    h2 : Høyde 2

    Returns
    -------
    d4 : Avstand ellipsoiden
    '''
    R  = 6371000
    d4 = 2*R*asin(sqrt((d2**2 -(h2 - h1)**2)/(4*(R+h1)*(R+h2))))
    return d4


def geod2Mgrid(a, b, lat, lon):
    '''
    Konverterer fra geodetiske koordinater til koordinater gitt i 
    Normal Mercator projeksjon (stående sylinder). Det som er kommentert ut skal i følge notatet
    til Jon Glenn være med. Dette samsvarer ikke med Wikipedia, og gir feil resultat 
    sammenlignet med en converter på internett. Hvis det leddet ikke er med stemmer det overrens.
    
    Parameters
    ----------
    a   : Store halvakse
    b   : lille halvakse
    lat : Breddegrad
    lon : Lengdegrad
    
    Returns
    -------
    north : Nordkoordinat gitt med Mercator projeksjon
    east : Østkoordinat gitt med Mercator projeksjon
    '''
    e = sqrt((a ** 2 - b ** 2)/(a ** 2))
    x = a*log(tan(pi/4 + lat/2)) # + ((a*e)/2)*log((1 - e*sin(lat))/(1 + e*sin(lat))) 
    y = a*lon
    north = x
    east = y
    return north, east




def red2ell(d1,z,R,k,H,N):
    '''
    Reduserer skråavstanden ned til ellipsoiden. Korrigerer senitvinkelen.

    Parameters
    ----------
    d1 : målt avstand
    z : målt zenit-vinkel
    R : gjennomsnittlig krumningsradius (6390000 i Norge)
    k : refraksjonskoeffisient
    H : snitthøyde for prosjekt (ortometrisk)
    N : geoide-snitthøyde

    Returns
    -------
    s : skråavstand redusert til ellipsoiden

    '''    
    # Korrigere zenit-vinkel (Skogseth s.152) 
    dz = -((d1*sin(z))/(2*R)) * (1-(k/(sin(z))**2))
    z = z + dz
    # Horisontering (Skogseth s.153)
    d_h = d1*sin(z)
    # Reduksjon til ellipsoide (Meyer s.124)
    s = d_h*(R/(R + N + H))
    
    return s



def geod2_indir(a, b, lat1, lon1, lat2, lon2):
    '''
    Indirekte geodetisk problem

    Parameters
    ----------
    a   : Store halvakse
    b   : lille halvakse
    lat1 : Breddegrad punkt 1
    lon1 : Lengdegrad punkt 1
    lat2 : Breddegrad punkt 2
    lon2 : Lengdegrad punkt 2

    Returns
    -------
    az1 : asimut 1
    az2 : asimut 2
    s   : 

    '''
    epsilon = 1e-10
    dlon_old = -1

    f = (a-b)/a
    e2 = (a ** 2-b ** 2)/b ** 2

    beta1 = atan(b/a*tan(lat1))
    beta2 = atan(b/a*tan(lat2))

    dlon = lon2-lon1
    dlon_new = dlon

    while (abs(dlon_new-dlon_old) > epsilon):

        dlon_old = dlon_new
    
        X = cos(beta1)*sin(beta2)-sin(beta1)*cos(beta2)*cos(dlon_old)
        Y = cos(beta2)*sin(dlon_old)
        Z = sin(beta1)*sin(beta2)+cos(beta1)*cos(beta2)*cos(dlon_old)
    
        sigma = atan2(sqrt(X ** 2+Y ** 2), Z)
        az1 = atan2(Y, X)
        az0 = asin(sin(az1)*cos(beta1))
    
        sigma1 = atan2(tan(beta1), cos(az1))
        sigma2 = sigma1+sigma
    
        K = (f+f ** 2)/4*cos(az0) ** 2-f ** 2/4*cos(az0) ** 4
    
        dlon_new = dlon+f*sin(az0)*((1-K-K ** 2)*sigma+K*sin(sigma) *
                                    cos(sigma1+sigma2)+K ** 2*sin(sigma)*cos(sigma)*cos(2*(sigma1+sigma2)))

    dlon = dlon_new

    az2 = atan2(cos(beta1)*sin(dlon), (cos(beta1)*sin(beta2)*cos(dlon)-sin(beta1)*cos(beta2)))

    if az2 < pi:
        az2 = az2+pi
    else:
        az2 = az2-pi
        
    g = e2*cos(az0) ** 2
    H = 1/8*g-1/16*g ** 2+37/1024*g ** 3
    b0 = b*(1+1/4*g-3/64*g ** 2+5/256*g ** 3)

    s = b0*(sigma-2*H*sin(sigma)*cos(sigma1+sigma2)-H ** 2/2*sin(2*sigma) *
            cos(2*(sigma1+sigma2))-H ** 3/3*sin(3*sigma)*cos(3*(sigma1+sigma2)))
    return az1, az2, s



def rotation_matrix(rx,ry,rz):
    """
    Build rotation matrix
    
    Parameters
    ----------
    rx : Rotation along x axis
    ry : Rotation along x axis
    rz : Rotation along x axis

    Returns
    -------
    R : Rotation matrix
    """
    Rx = array([[1, 0, 0],
              [0, cos(rx), -sin(rx)],
              [0, sin(rx),  cos(rx)]])
    
    Ry = array([[cos(ry), 0, sin(ry)],
                  [0, 1, 0],
                  [-sin(ry), 0, cos(ry)]])
    
    
    Rz = array([[cos(rz), -sin(rz), 0],
                      [sin(rz),  cos(rz), 0],
                      [0,        0,       1]])
    
    R = Rx@Ry@Rz
    return R


def EUREF89_2_NGO(lat,lon,h):
    # --- Params WGS84
    a      = 6378137                # Store halvakse WGS84
    f      = 1/298.257222100882711  # Flattrykning
    b      = a*(1-f)                # Lille halvakse WGS84

    # Convert from geodetic to ECEF
    X,Y,Z = geod2ECEF(a, b, lat, lon, h)
    P = array([[X],
               [Y],
               [Z]])
    
    # Parameters 7-parameter transformation (NMBU campus)
    T = array([[-313.368],
               [125.818],
               [-626.643]])
    
    m = (1 + 7.781959e-6)
    
    rx = dms2rad(0, 0, -2.336248)
    ry = dms2rad(0, 0, -1.712020)
    rz = dms2rad(0, 0,  1.169871)
    R = rotation_matrix(rx,ry,rz)

    P = T + m*R@P

    # Modified Bessel ellipsoid
    a = 6377492.0176
    f = 1/299.15281285
    b = a*(1 - f)
    
    # TM projection
    lat0 = deg2rad(58)
    lon0 = deg2rad(8.38958333333337)  # axis 2
    scale = 1
    fnorth = 0
    feast = 0
    
    # Convert from ECEF to geodetic
    X = P[0, 0]
    Y = P[1, 0]
    Z = P[2, 0]
    lat, lon, h = ECEF2geodb(a, b, X,Y,Z)
    
    # Convert from geodetic to projection (NGO48)
    N, E = geod2TMgrid(a, b, lat, lon, lat0, lon0, scale, fnorth, feast)

    return N,E,h

# N,E,h =EUREF89_2_NGO(deg2rad(63.03530998),deg2rad(8.91120223),300)
# print(N,E,h)

def ellip_and_proj_params(ellip,proj=None,zone=None):
    """
    Function that extract correct parameters for ellipsoids and projections. Proj can be UTM or GH (Gauss-Krüger).
    Zone can be from 1-8 (GK) or 31-33 (UTM)
    
    """
    
    if proj == None and zone == None:
        if ellip == 'WGS84':
            a      = 6378137            # Store halvakse WGS84
            f      = 1/298.257223563    # Flattrykning
            b      = a*(1-f)            # Lille halvakse WGS84
            lat0   = 0                  # Ekvator er nullbredde
            lon0   = deg2rad(9)         # Lengdegraden til sentralmeridianen til sone 32
            scale  = 0.9996             # UTM skalering
            fnorth = 0                  # Ingen falsk nord i UTM
            feast  = 500000             # Falsk øst
        elif ellip == 'GRS80':
            a      = 6378137                # Store halvakse WGS84
            f      = 1/298.257222100882711  # Flattrykning
            b      = a*(1-f)                # Lille halvakse WGS84
            lat0   = 0                      # Ekvator er nullbredde
            lon0   = deg2rad(9)             # Lengdegraden til sentralmeridianen til sone 32
            scale  = 0.9996                 # UTM skalering
            fnorth = 0                      # Ingen falsk nord i UTM
            feast  = 500000                 # Falsk øst
            
        elif ellip == "Mod.Bessel":
            a      = 6377492.0176              # Store halvakse modified Bessel
            f      = 1 /299.15281285           # Flattrykning modified Bessel
            b      = a*(1-f)                   # Lille halvakse modified Bessel
            lat0   = deg2rad(58)               # Ekvator er nullbredde
            lon0   = deg2rad(10.7229166666667) # Lengdegraden til sentralmeridianen til akse 3 (NGO)
            scale  = 1                         # UTM skalering
            fnorth = 0                         # Ingen falsk nord i UTM
            feast  = 0                         # Falsk øst
                
        elif ellip =="ED50":
            a      = 6378388.000    # Store halvakse ED50 (international 1924)
            f      = 1/297.000      # Falttrykning  ED50 (international 1924)
            b      = 6356752.3142   # Lille halvakse ED50 (international 1924)
            lat0   = 0              # Ekvator er nullbredde
            lon0   = deg2rad(9)     # Lengdegraden til sentralmeridianen til sone 32
            scale  = 0.9996         # UTM skalering
            fnorth = 0              # Ingen falsk nord i UTM
            feast  = 500000         # Falsk øst
    else:
        
        if proj == "UTM" and zone == "31":
            lon0 = deg2rad(3)         # Lengdegraden til sentralmeridianen til zone 31
        if proj == "UTM" and zone == "32":
            lon0 = deg2rad(9)         # Lengdegraden til sentralmeridianen til zone 32
        if proj == "UTM" and zone == "33":
            lon0 = deg2rad(15)         # Lengdegraden til sentralmeridianen til zone 32
            
            
        if proj == "GK" and zone == "1":
            lon0 = deg2rad(6.0562500000000306)   # Lengdegraden til sentralmeridianen til akse 1
        if proj == "GK" and zone == "2":
            lon0 = deg2rad(8.38958333333337)     # Lengdegraden til sentralmeridianen til zone 2
        if proj == "GK" and zone == "3":
            lon0 = deg2rad(10.7229166666667)     # Lengdegraden til sentralmeridianen til zone 3
        if proj == "GK" and zone == "4":
            lon0 = deg2rad(13.2229166666667)   # Lengdegraden til sentralmeridianen til akse 4
        if proj == "GK" and zone == "5":
            lon0 = deg2rad(16.88958333333337)     # Lengdegraden til sentralmeridianen til zone 5
        if proj == "GK" and zone == "6":
            lon0 = deg2rad(20.8895833333334)     # Lengdegraden til sentralmeridianen til zone 6
        if proj == "GK" and zone == "7":
            lon0 = deg2rad(24.8895833333334)   # Lengdegraden til sentralmeridianen til akse 7
        if proj == "GK" and zone == "8":
            lon0 = deg2rad(29.05625)   # Lengdegraden til sentralmeridianen til akse 8
    
        
        if ellip == 'WGS84':
            a      = 6378137            # Store halvakse WGS84
            f      = 1/298.257223563    # Flattrykning
            b      = a*(1-f)            # Lille halvakse WGS84
            lat0   = 0                  # Ekvator er nullbredde
            scale  = 0.9996             # UTM skalering
            fnorth = 0                  # Ingen falsk nord i UTM
            feast  = 500000             # Falsk øst
        elif ellip == 'GRS80':
            a      = 6378137                # Store halvakse WGS84
            f      = 1/298.257222100882711  # Flattrykning
            b      = a*(1-f)                # Lille halvakse WGS84
            lat0   = 0                      # Ekvator er nullbredde
            scale  = 0.9996                 # UTM skalering
            fnorth = 0                      # Ingen falsk nord i UTM
            feast  = 500000                 # Falsk øst
            
        elif ellip == "Mod.Bessel":
            a      = 6377492.0176              # Store halvakse modified Bessel
            f      = 1 /299.15281285           # Flattrykning modified Bessel
            b      = a*(1-f)                   # Lille halvakse modified Bessel
            lat0   = deg2rad(58)               # Ekvator er nullbredde
            scale  = 1                         # UTM skalering
            fnorth = 0                         # Ingen falsk nord i UTM
            feast  = 0                         # Falsk øst
                
        elif ellip =="ED50":
            a      = 6378388.000    # Store halvakse ED50 (international 1924)
            f      = 1/297.000      # Falttrykning  ED50 (international 1924)
            b      = 6356752.3142   # Lille halvakse ED50 (international 1924)
            lat0   = 0              # Ekvator er nullbredde
            scale  = 0.9996         # UTM skalering
            fnorth = 0              # Ingen falsk nord i UTM
            feast  = 500000         # Falsk øst

            
    return a,b,lat0,lon0,scale,fnorth,feast

# import pyproj
# from pyproj import CRS
# a,b,lat0,lon0,scale,fnorth,feast=ellip_and_proj_params("GRS80")
# x = 495507.6008
# y = 6989526.4622
# gamma = TMconv(a, b, x, y, lat0)
# gamma = rad2deg(gamma)
# skalering = TMscale(a, b, x, y, lat0)

# # Proj beregnet meridianconv
# crs_EUREF89 = CRS.from_user_input(25832)
# p = pyproj.Proj(crs_EUREF89)
# lat, lon = 63.03530998, 8.91120223
# facts = p.get_factors(lon, lat, radians=False)
# print(p(lon, lat), facts.meridian_convergence, facts.meridional_scale,  crs_EUREF89.name)



