def JD2GMST(JD): 
    #----------------------------- Begin Function -----------------------------
    #Purpose:
    #--------
    #Convert a specified Julian Date Vector to Greenwhich Mean Sidereal Time.
    #
    #Inputs:
    #-------
    #JD             [N x M x L]                         Julian Date Vector                      
    #       
    #
    #Outputs:
    #--------
    #GMST           [N x M x L]                         Greenwich Mean Sidereal
    #                                                   Time in degrees from
    #                                                   0-360
    #
    #References:
    #-----------
    #Approximate Sidereal Time, 
    #http://www.usno.navy.mil/USNO/astronomical-applications/...
    #astronomical-information-center/approx-sider-time
    #
    #Universal Sidereal Times, The Astronomical Almanac For The Year 2004
    #
    #Programed by: 
    #-------------
    #Darin Koblick 07-11-2010
    #--------------------------------------------------------------------------
    import pandas as pd
    ts = pd.Timestamp(year = 2011,  month = 11, day = 21, hour = 10, second = 49, tz = 'Europe/Berlin')
    ts.to_julian_date()
    ##Find the Julian Date of the previous midnight, JD0
    JD0 = NaN(size(JD));
    JDmin = floor(JD)-.5;
    JDmax = floor(JD)+.5;
    JD0(JD > JDmin) = JDmin(JD > JDmin);
    JD0(JD > JDmax) = JDmax(JD > JDmax);
    H = (JD-JD0)*24;       #Time in hours past previous midnight
    D = JD - 2451545.0;     #Compute the number of days since J2000
    D0 = JD0 - 2451545.0;   #Compute the number of days since J2000
    T = D/36525;           #Compute the number of centuries since J2000
    ##Calculate GMST in hours (0h to 24h) ... then convert to degrees
    GMST = mod(6.697374558 + 0.06570982441908.*D0  + 1.00273790935.*H + ...
        0.000026.*(T.^2),24).*15;
    return GMST