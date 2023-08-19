from read_SP3Nav import *
from readRinexNav import *
from gpstime2date import gpstime2date
from tqdm import tqdm
import numpy as np
from Geodetic_functions import ECEF2geodb, ECEF2enu, atanc
from numpy import arctan,sqrt
from preciseOrbits2ECEF import preciseOrbits2ECEF
import warnings
warnings.filterwarnings(action='ignore', message='invalid value encountered in fmod')



def computeSatElevations(GNSS_SVs, GNSSsystems, approxPosition,\
    nepochs, time_epochs, max_sat, sp3_nav_filename_1, sp3_nav_filename_2, sp3_nav_filename_3):
    """
    Function that computes the satellite elevation angles of all satellites
    of all GNSS systems at each epoch of the observations period.
    --------------------------------------------------------------------------------------------------------------------------
     INPUTS
    
     GNSS_SVs:                 Dict containing number of satellites with 
                               observations and PRN for each satellite, 
                               for each GNSS system
    
                               GNSS_SVs[GNSSsystemIndex][epoch,j]     
                                           j=0: (index null) number of observed GPS satellites 
                                           j>1: PRN of observed satellites
    
     GNSSsystems:              Dict containing different GNSS 
                               systems included in RINEX file. Elements are strings. 
                               Must be either "G","R","E" or "C".  
    
     approxPosition:           array containing approximate position from rinex
                               observation file header. [X, Y, Z]
    
     nepochs:                  number of epochs with observations in 
                               rinex observation file.
    
     time_epochs:              matrix conatining gps-week and "time of week" 
                               for each epoch
                               time_epochs(epoch,i),   i=1: week
                                                       i=2: time-of-week in seconds (tow)
    
     max_sat:                  array that stores max satellite PRN number for 
                               each of the GNSS systems. Follows same order as GNSSsystems
    
     nav_filename:             string, path and filename of rinex3.xx navigation file.
    
     almanac_nav_filename:     string, path and filename of sen almanac
                               navigation filename for GLONASS
    --------------------------------------------------------------------------------------------------------------------------
     OUTPUTS
    
     sat_elevation_angles:     Dict contaning satellite elevation angles at each
                               epoch, for each GNSS system. 
                               sat_elevation_angles[GNSSsystemIndex][epoch, PRN]
                               
     sat_azimut_angles:        Dict contaning satellite azimuth angles at each
                               epoch, for each GNSS system. 
                               sat_azimut_angles[GNSSsystemIndex][epoch, PRN]
                               
                               
     sat_coordinates:          Dict contaning satellite coordinates for each
                               epoch, for each GNSS system. 
                               sat_coordinates[GNSSsystemIndex][epoch, PRN]
    --------------------------------------------------------------------------------------------------------------------------
    """
    
    nGNSSsystems = len(GNSSsystems)
    sat_coordinates = {}
    sat_elevation_angles = {}
    sat_azimut_angles    = {}
    
    if all(approxPosition == 0):
        print('ERROR(computeSatElevations): The approximate receiver position is missing.\n',\
            'Please chack that APPROX POSITION XYZ is in header of Rinex file.\n',\
            'Elevation angles will not be computed\n\n.')
        return
    
    ## -- Testing whether two and tree sp3 files are defined
    two_sp3_files = 0
    three_sp3_files = 0
    
    if sp3_nav_filename_2 != "":
        two_sp3_files = 1
        if sp3_nav_filename_3 != "":
          three_sp3_files = 1
       
        
    
    ## ---  Read first SP3 file
    sat_positions_1, epoch_dates_1, navGNSSsystems_1, nEpochs_1, epochInterval_1, success = readSP3Nav(sp3_nav_filename_1)
    # if two SP3 files inputted by user, read second SP3 file
    if two_sp3_files:
        sat_positions_2, epoch_dates_2, navGNSSsystems_2, nEpochs_2, epochInterval_2, success = readSP3Nav(sp3_nav_filename_2)
    else:
        # sat_positions_2, epoch_dates_2, navGNSSsystems_2, nEpochs_2, epochInterval_2 = deal(NaN)
        sat_positions_2, epoch_dates_2, navGNSSsystems_2, nEpochs_2, epochInterval_2 = np.nan,np.nan,np.nan,np.nan,np.nan
    
    # if three SP3 files inputted by user, read third SP3 file
    if three_sp3_files:
        sat_positions_3, epoch_dates_3, navGNSSsystems_3, nEpochs_3, epochInterval_3, success = readSP3Nav(sp3_nav_filename_3)
    else:
        # sat_positions_3, epoch_dates_3, navGNSSsystems_3, nEpochs_3, epochInterval_3 = deal(NaN);
        sat_positions_3, epoch_dates_3, navGNSSsystems_3, nEpochs_3, epochInterval_3 = np.nan,np.nan,np.nan,np.nan,np.nan
    
    ## Combine data from different SP3 files 
    if two_sp3_files:
        sat_positions, epoch_dates, navGNSSsystems, nEpochs, epochInterval, success = combineSP3Nav(three_sp3_files,\
            sat_positions_1, epoch_dates_1, navGNSSsystems_1, nEpochs_1, epochInterval_1,\
            sat_positions_2, epoch_dates_2, navGNSSsystems_2, nEpochs_2, epochInterval_2,\
            sat_positions_3, epoch_dates_3, navGNSSsystems_3, nEpochs_3, epochInterval_3,GNSSsystems)
    else:
        sat_positions = sat_positions_1;
        epoch_dates = epoch_dates_1;
        navGNSSsystems = navGNSSsystems_1;
        nEpochs = nEpochs_1;
        epochInterval = epochInterval_1;
    
    satMissingData = []
    
    
    bar_format = '{desc}: {percentage:3.0f}%|{bar}| ({n_fmt}/{total_fmt})'
    for k in tqdm(range(0,nGNSSsystems),desc='Looping through the systems',position=0,leave=True,bar_format=bar_format):
        n_ep = int(nepochs)
        nSat = int(max_sat[k]) + 1
        # Initialize data matrix for current GNSSsystem
        sat_elevation_angles[k] = np.zeros([n_ep, nSat])
        sat_azimut_angles[k]    = np.zeros([n_ep, nSat]) 
        X = np.full((n_ep, nSat), np.nan)  # Array for storing X-coordinate
        Y = np.full((n_ep, nSat), np.nan)  # Array for storing Y-coordinate
        Z = np.full((n_ep, nSat), np.nan)  # Array for storing Z-coordinate
        
        sys = GNSSsystems[k+1]
        sat_coordinates[sys] = {} 
        if sys in navGNSSsystems: 
            curr_pos = {}  # dict for storing data   
            for epoch in tqdm(range(0,nepochs),desc='Satellite elevation angles are being calculated for system %s of %s' %(k+1, nGNSSsystems),position=0, leave=False, bar_format=bar_format):             
                ##-- GPS Week and time of week of current epoch
                week = time_epochs[epoch,0]
                tow  = time_epochs[epoch,1]
                ## -- Satellites in current epoch that should have elevation computed
                # SVs = np.nonzero(FirstLastObsEpochOverview[k][epoch, :])[0] # COMMENTED OUT 11.03.2023 to prevent computing elevation angles for satellites not visalble in the current epoc
                SVs  = GNSS_SVs[sys][epoch,1::][GNSS_SVs[sys][epoch,1::] != 0].astype(int) #extract only nonzero and convert to integer
                n_sat = len(SVs)
                for sat in np.arange(0,n_sat):
                    # Get satellite elevation angle of current sat at current epoch
                    PRN = int(SVs[sat])
                    elevation_angle, azimut_angle, missing_nav_data, X[epoch,PRN], Y[epoch,PRN], Z[epoch,PRN] = \
                        get_elevation_angle(GNSSsystems[k+1], SVs[sat], week, tow, sat_positions, nEpochs,\
                       epoch_dates, epochInterval, navGNSSsystems, approxPosition);
                    curr_pos[int(SVs[sat])] = np.array([X[:,PRN],Y[:,PRN],Z[:,PRN]]).T
                    sat_elevation_angles[k][epoch,SVs[sat]] = elevation_angle
                    sat_azimut_angles[k][epoch,SVs[sat]] = azimut_angle
                    if missing_nav_data:
                        ## Combine PRN number and GNSS system 
                        SVN = str(GNSSsystems[k]) + str(SVs(sat))
                        if not any(SVN in satMissingData):
                            satMissingData.append(SVN)
                sat_coordinates[sys]  = curr_pos
                           
  

        
    if satMissingData:
        print('INFO(computeSatElevations): The following satellites had missing orbit data in SP3 file.',\
        'Their elevation angles were set to 0 for these epochs') 
        print(satMissingData)
    
        
    print('INFO(computeSatElevations): Satellite elevation angles have been computed')

    return sat_elevation_angles, sat_azimut_angles, sat_coordinates



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
    
    elevation_angle:  Elevation angle of specified satelite at specified
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
            if not 0 < elevation_angle < 90: # if the satellite is below the horizon (elevation angle is below zero)
                elevation_angle = np.nan
            
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
