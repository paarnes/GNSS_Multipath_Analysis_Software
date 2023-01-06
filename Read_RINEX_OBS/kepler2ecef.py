# from read_rinex2_nav import read_rinex2_nav
# from compute_azimut_elev import compute_azimut_elev

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

# week, tow = date2gpstime(2022,11,5,9,45,00)

def extract_nav_message(data,PRN,tidspunkt):
    """
    Funksjonen samler de satellittene som har de spesifiserte PRN nummrene, samt de som ligger nærmest spesifisert tidspunkt. 

    """
    samla_data = gathering_sat_by_PRN(data,PRN);
    Sat_liste = find_message_closest_in_time(samla_data, tidspunkt)
    return Sat_liste
    

def gathering_sat_by_PRN(data,PRN):
    """
    Funksjonen bruker rinex-data i form av array ordnet ved funksjonen read_rinex2_nav
    """
    import numpy as np
    j = 0;
    m = len(data) # kanskje m = len(data[0]) for å få ant col 
    # Sat_data = np.array(np.empty)
    Sat_data  = np.zeros((1,36))
    for k in range(0,m):
        PRN_ = np.array([])
        
        if len(data[k,0]) == 3: # RINEX v3 navfiles have letters in addition. EX: G01.
            sat = data[k,0][1::]
        else:
            sat = data[k,0]
        
        if int(sat) == PRN:
            # Sat_data[j,:] = data[k,:]
            PRN_ =np.append(PRN_,data[k,:])
            PRN_ = PRN_.reshape(1,len(PRN_))
        if np.size(PRN_) != 0:
            Sat_data  = np.concatenate([Sat_data , PRN_], axis=0)
            
    # check = str(Sat_data[0,:] == 0)
    # if 'False' not in check:
    #     Sat_data  = np.delete(Sat_data , (0), axis=0)
    if all(Sat_data[0,:].astype(float)) == 0:
        Sat_data  = np.delete(Sat_data , (0), axis=0)


    return Sat_data


def find_message_closest_in_time(data,tow_mot):
    """
    Funksjonen plukker ut den linjen i et datasett som ligg nærmest tidspunktet vi ønsker i bestemme koordinatene for. 
    """
    import numpy as np
    towSat = [] # Tom liste for å lagre tidspunktene
    length = len(data)
    towSat = np.array([])
    for i in range(0,length):
        week, tow = date2gpstime(2000 + int(data[i,1]), int(data[i,2]), int(data[i,3]), int(data[i,4]), int(data[i,5]), int(data[i,6]))
        towSat = np.append(towSat,tow)
    # Finner verdien som ligger nærmest
    # index = int(min(abs(tow_mot - towSat)))
    index = np.abs(tow_mot - towSat).argmin()
    GNSS_linjer = data[index,:]   
    
    return GNSS_linjer





def Satkoord2(efemerider,t,xm,ym,zm):
    """
    Funksjonen beregner satellittkordinater og korrigerer for jordrotasjonen. Fra keplerelemtenter til ECEF.
    """
    
    from numpy import sqrt,pi, fmod,cos, sin, tan
    from math import atan
    
    
    GM         = 3.986005e14      # Produktet av jordas masse og gravitasjonskonstanten
    omega_e    = 7.2921151467e-5  # [rad/sek]
    c          = 299792458        # Lyshastigheten [m/s]
    
    tow_mot  = t
    
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
    
    n0  = sqrt(GM/A**3);   #(rad/s)
    t_k = tow_mot - toe;
    n_k = n0 + delta_n;   #Koorigert  midlere bevegelse
    M_k = M0 + n_k*t_k;   #Midlere anomali (rad/s)
    
    
    #Beregner eksentrisk anomali
    E_old = M_k;
    for i in range(0,10):
        E = M_k+ e*sin(E_old);
        dE = fmod(E - E_old, 2*pi)
        if abs(dE) < 1.e-12:
           break
        E_old = E
        
    E      = fmod(E+2*pi,2*pi);
    
    cosv = (cos(E) - e)/(1 - e*cos(E));
    sinv = (sqrt(1 - e**2)*sin(E))/(1-e*cos(E))
    tanv = sinv/cosv
    
    ## -- Kvadrantkorreksjon
    if sinv > 0 and cosv < 0 or sinv < 0 and cosv < 0:
        v = atan(tanv) + pi
    elif sinv < 0 and cosv > 0:
        v = atan(tanv) + 2*pi
    else: 
        v = atan(tanv)

    
    theta   = v + omega;
    
    ## -- Korreksjon på baneparametere:
    du_k    = C_uc*cos(2*theta) + C_us*sin(2*theta);
    dr_k    = C_rc*cos(2*theta) + C_rs*sin(2*theta);
    di_k    = C_ic*cos(2*theta) + C_is*sin(2*theta);
    
    ## -- Korrigerte baneparametere:
    u_k     = theta + du_k;
    r_k     = A*(1 - e*cos(E)) + dr_k;
    i_k     = i0 + i_dot*t_k + di_k;
    
    
    ## -- Korrigert lengdegrad for baneplanknuten:
    OMEGA_k = OMEGA + (OMEGA_dot - omega_e)*t_k - omega_e*toe;
    
    
    ## -- Satelittens posisjon i banen:
    x = r_k*cos(u_k);
    y = r_k*sin(u_k); 
    
    ## --ECEF-koordinater for satellitten:
    X = x*cos(OMEGA_k) - y*sin(OMEGA_k)*cos(i_k);
    Y = x*sin(OMEGA_k) + y*cos(OMEGA_k)*cos(i_k);
    Z = y*sin(i_k);
    
    ## Relativistisk klokkekorreksjon
    dT_rel = (-2/c**2)*sqrt(A*GM)*e*sin(E);
    
    ## --Tester om korreksjonen skal utføres.
    if (abs(xm)) > 1.0 and (abs(ym)) > 1.0 and (abs(zm)) > 1.0:
        TRANS  = 0;
        TRANS0 = 0.075 ;
        j = 0;
        while(abs(TRANS0 - TRANS) > 1e-10):
            j = j +1;
            if(j > 20):
                print('Feil, Gangtids-rotasjonen konvergerer ikke!')
                break
            
            TRANS = TRANS0;
            OMEGA_k = OMEGA + (OMEGA_dot - omega_e)*t_k - omega_e*(toe + TRANS);
            X = x*cos(OMEGA_k) - y*sin(OMEGA_k)*cos(i_k);
            Y = x*sin(OMEGA_k) + y*cos(OMEGA_k)*cos(i_k);
            Z = y*sin(i_k);
            #Regner ut avstanden mottaker-satellitt
            dX = (X - xm);
            dY = (Y - ym);
            dZ = (Z - zm);
            DS = sqrt(dX**2 + dY**2 + dZ**2);
            #Regner ny verdi for signalets gangtid
            TRANS0 = DS/c;
        
    else:
        #Jordrotasjonen skal ikke utføres
        OMEGA_k = OMEGA + (OMEGA_dot - omega_e)*t_k - omega_e*toe;
        X = x*cos(OMEGA_k) - y*sin(OMEGA_k)*cos(i_k);
        Y = x*sin(OMEGA_k) + y*cos(OMEGA_k)*cos(i_k);
        Z = y*sin(i_k);
        
    return X,Y,Z,dT_rel



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

# data, header, n_eph = read_rinex2_nav('testfile.20n')
# tow_mot = 556559.999999882
# X0 =[];Y0 =[]; Z0 = []; dT_rel0 = []
# x = 3149785.9652
# y = 598260.8822
# z = 5495348.4927 
# for ep in range(0,len(data)):
#     X0_,Y0_,Z0_,dT_rel0_ = Satkoord2(data[ep,:],tow_mot,x,y,z)
#     X0.append(X0_)
#     Y0.append(Y0_)
#     Z0.append(Z0_)
#     dT_rel0.append(dT_rel0_)


# import pandas as pd
# df_satpos_ECEF = pd.DataFrame(list(zip(X0,Y0,Z0,dT_rel0)),columns=["X","Y","Z","dT_rel"])

# xm = 3149785.9652
# ym = 598260.8822
# zm = 5495348.4927 
# X = df_satpos_ECEF['X'].to_numpy()
# Y = df_satpos_ECEF['Y'].to_numpy()
# Z = df_satpos_ECEF['Z'].to_numpy()
# az, elev = compute_azimut_elev(X,Y,Z,xm,ym,zm)