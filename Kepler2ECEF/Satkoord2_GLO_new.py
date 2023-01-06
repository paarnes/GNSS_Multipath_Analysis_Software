import datetime
import json
import warnings
from abc import ABC, abstractmethod
from collections import defaultdict
from enum import IntEnum
from typing import Dict, List, Optional

import numpy as np
import numpy.polynomial.polynomial as poly
# from datetime import datetime
from math import sin, cos, sqrt, fabs, atan2

import os, sys,numpy as np
sys.path.append(r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Read_RINEX_OBS')
sys.path.append(r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Read_RINEX_nav')
sys.path.append(r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Kepler2ECEF')
from read_rinex3_nav import read_rinex3_nav
# from read_SP3Nav import readSP3Nav
from kepler2ecef import Satkoord2
from kepler2ecef import date2gpstime
from kepler2ecef import extract_nav_message
from readRinexObs304 import readRinexObs304
from scipy.interpolate import lagrange
import matplotlib.pyplot as plt
from tqdm import trange,tqdm
from compute_azimut_elev import compute_azimut_elev
import pandas as pd


def get_leap_seconds(time):
  """
  Add leap seconds based on years
  """
  if time <= GPSTime.from_datetime(datetime.datetime(2006, 1, 1)):
    raise ValueError("Don't know how many leap seconds to use before 2006")
  elif time <= GPSTime.from_datetime(datetime.datetime(2009, 1, 1)):
    return 14
  elif time <= GPSTime.from_datetime(datetime.datetime(2012, 7, 1)):
    return 15
  elif time <= GPSTime.from_datetime(datetime.datetime(2015, 7, 1)):
    return 16
  elif time <= GPSTime.from_datetime(datetime.datetime(2017, 1, 1)):
    return 17
  else:
    return 18


def gpst_to_utc(t_gpst):
    t_utc = t_gpst - get_leap_seconds(t_gpst)
    if utc_to_gpst(t_utc) - t_gpst != 0:
      return t_utc + 1
    else:
      return t_utc


def utc_to_gpst(t_utc):
    t_gpst = t_utc + get_leap_seconds(t_utc)
    return t_gpst


def GLONASSEphemeris(ephemerides):
    ## Parameters (from GLONASS Interface Control Document 1988)
    GM         = 398600.44e9     # Gravitational constant [m3/s2]   (product of the mass of the earth and and gravity constant)
    omega_e    = 7.292115e-5     # Earth rotation rate    [rad/sek]
    c          = 299792458       # Speed of light         [m/s]
    a          = 6378136         # Semi major axis PZ-90   [m]
    f          = 1/298.257839303 # Inverse flattening
    C_20       = -1082.63e-6     # Second zonal coefficient of spherical harmonic expression.
    ephemerides[0] = 9999
    ephemerides = ephemerides.astype(float)

    ## Read in data:
    toc = [ephemerides[1],ephemerides[2] ,ephemerides[3] ,ephemerides[4],ephemerides[5],ephemerides[6]] # year,month,day,hour,minute,second
    week,toc = date2gpstime(int(ephemerides[1]),int(ephemerides[2]) ,int(ephemerides[3]) ,int(ephemerides[4]),int(ephemerides[5]),int(ephemerides[6]))
    tauN = ephemerides[6]   # SV clock bias (sec) (-TauN)
    gammaN = ephemerides[7] # SV relative frequency bias (+GammaN)
    
    # toe  = ephemerides[9]   # Time of ephemerides (i think??). Time of week in seconds (UTC)    
    # _,t_rec = date2gpstime(int(ephemerides[1]),int(ephemerides[2]) ,int(ephemerides[3]) ,int(ephemerides[4]),int(ephemerides[5]),int(ephemerides[6]))
    _,t_rec = date2gpstime(int(ephemerides[1]),int(ephemerides[2]) ,int(ephemerides[3]) ,int(23),int(15),int(ephemerides[6]))
    
    x_te = ephemerides[10]  # X-coordinates at t_e in PZ-90 [km]
    y_te = ephemerides[14]  # Y-coordinates at t_e in PZ-90 [km]
    z_te = ephemerides[18]  # Z-coordinates at t_e in PZ-90 [km]
    
    vx_te = ephemerides[11] # Velocity component at t_e in PZ-90 (v_x) [km/s]
    vy_te = ephemerides[15] # Velocity component at t_e in PZ-90 (v_y) [km/s]
    vz_te = ephemerides[19] # Velocity component at t_e in PZ-90 (v_z) [km/s]
    
    J_x = ephemerides[12]   # Moon and sun acceleration at t_e [km/sec**2]
    J_y = ephemerides[16]   # Moon and sun acceleration at t_e [km/sec**2]
    J_z = ephemerides[20]   # Moon and sun acceleration at t_e [km/sec**2]
    


    # TODO should handle leap seconds better
    # toc_gps_time = utc_to_gpst(toc)
    toc_gps_time = toc +  18 # CHACE THIS!!
    # tdiff = time - toc_gps_time
    # tdiff = toe - toc_gps_time # toe??
    tdiff = t_rec - toc_gps_time # toe??

    # Clock correction (except for general relativity which is applied later)
    clock_err = tauN + tdiff * (gammaN)
    clock_rate_err = gammaN

    init_state = np.empty(6)
    init_state[0] = x_te
    init_state[1] = y_te
    init_state[2] = z_te
    init_state[3] = vx_te
    init_state[4] = vy_te
    init_state[5] = vz_te
    init_state = 1000*init_state         # converting from km to meters
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

def glonass_diff_eq(state, acc):
    """
    State is a cevtor containing x,y,z,vx,vy,vz from navigation message
    """
    J2 = 1.0826257e-3
    mu = 3.9860044e14
    omega = 7.292115e-5
    ae = 6378136.0
    r = np.sqrt(state[0]**2 + state[1]**2 + state[2]**2)
    ders = np.zeros(6)
    if r**2 < 0:
        return ders
    a = 1.5 * J2 * mu * (ae**2)/ (r**5)
    b = 5 * (state[2]**2) / (r**2)
    c = -mu/(r**3) - a*(1-b)
    
    ders[0:3] = state[3:6]
    ders[3] = (c + omega**2)*state[0] + 2*omega*state[4] + acc[0]
    ders[4] = (c + omega**2)*state[1] - 2*omega*state[3] + acc[1]
    ders[5] = (c - 2*a)*state[2] + acc[2]
    return ders



class GPSTime:
  """
  GPS time class to add and subtract [week, tow]
  """
  def __init__(self, week, tow):
    self.week = week
    self.tow = tow
    self.seconds_in_week = 604800

  @classmethod
  def from_datetime(cls, datetime):
    week, tow = datetime_to_tow(datetime)
    return cls(week, tow)

  @classmethod
  def from_glonass(cls, cycle, days, tow):
    # https://en.wikipedia.org/wiki/GLONASS
    # Day number (1 to 1461) within a four-year interval
    # starting on 1 January of the last leap year
    t = datetime.datetime(1992, 1, 1, 0, 0, 0, 0, None)
    t += datetime.timedelta(days=cycle*(365*4+1)+(days-1))
    # according to Moscow decree time.
    t -= datetime.timedelta(hours=3)
    t += datetime.timedelta(seconds=tow)
    ret = cls.from_datetime(t)
    return utc_to_gpst(ret)

  @classmethod
  def from_meas(cls, meas):
    return cls(meas[1], meas[2])

  def __sub__(self, other):
    if isinstance(other, type(self)):
      return (self.week - other.week)*self.seconds_in_week + self.tow - other.tow
    elif isinstance(other, float) or isinstance(other, int):
      new_week = self.week
      new_tow = self.tow - other
      while new_tow < 0:
        new_tow += self.seconds_in_week
        new_week -= 1
      return GPSTime(new_week, new_tow)
    raise NotImplementedError(f"subtracting {other} from {self}")

  def __add__(self, other):
    if isinstance(other, float) or isinstance(other, int):
      new_week = self.week
      new_tow = self.tow + other
      while new_tow >= self.seconds_in_week:
        new_tow -= self.seconds_in_week
        new_week += 1
      return GPSTime(new_week, new_tow)
    raise NotImplementedError(f"adding {other} from {self}")

  def __lt__(self, other):
    return self - other < 0

  def __gt__(self, other):
    return self - other > 0

  def __le__(self, other):
    return self - other <= 0

  def __ge__(self, other):
    return self - other >= 0

  def __eq__(self, other):
    return self - other == 0

  def as_datetime(self):
    return tow_to_datetime(self.tow, self.week)

  def as_unix_timestamp(self):
    return (gpst_to_utc(self).as_datetime() - datetime.datetime(1970, 1, 1)).total_seconds()

  @property
  def day(self):
    return int(self.tow/(24*3600))

  def __repr__(self):
    return f"GPSTime(week={self.week}, tow={self.tow})"

def datetime_to_tow(t):
    """
    Convert a Python datetime object to GPS Week and Time Of Week.
    Does *not* convert from UTC to GPST.
    Fractional seconds are supported.
    Parameters
    ----------
    t : datetime
      A time to be converted, on the GPST timescale.
    mod1024 : bool, optional
      If True (default), the week number will be output in 10-bit form.
    Returns
    -------
    week, tow : tuple (int, float)
      The GPS week number and time-of-week.
    """
    # DateTime to GPS week and TOW
    wk_ref = datetime.datetime(2014, 2, 16, 0, 0, 0, 0, None)
    refwk = 1780
    wk = (t - wk_ref).days // 7 + refwk
    tow = ((t - wk_ref) - datetime.timedelta((wk - refwk) * 7.0)).total_seconds()
    return wk, tow



def tow_to_datetime(tow, week):
    """
    Convert a GPS Week and Time Of Week to Python datetime object.
    Does *not* convert from GPST to UTC.
    Fractional seconds are supported.
    Parameters
    ----------
    tow : time of week in seconds
    weeks : gps week
    Returns
    -------
    t : datetime
      Python datetime
    """
    #  GPS week and TOW to DateTime
    t = datetime.datetime(1980, 1, 6, 0, 0, 0, 0, None)
    t += datetime.timedelta(seconds=tow)
    t += datetime.timedelta(weeks=week)
    return t


def get_leap_seconds(time):
  # TODO use library for this
  import datetime
  if time <= GPSTime.from_datetime(datetime.datetime(2006, 1, 1)):
    raise ValueError("Don't know how many leap seconds to use before 2006")
  elif time <= GPSTime.from_datetime(datetime.datetime(2009, 1, 1)):
    return 14
  elif time <= GPSTime.from_datetime(datetime.datetime(2012, 7, 1)):
    return 15
  elif time <= GPSTime.from_datetime(datetime.datetime(2015, 7, 1)):
    return 16
  elif time <= GPSTime.from_datetime(datetime.datetime(2017, 1, 1)):
    return 17
  else:
    return 18


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

sp3 = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Multipath_analysis/NMBU_SAMSUNG_S20_ALL.20p'
data,header, n_eph = read_rinex3_nav(sp3,dataframe='no')
ephemerides = data[836]
pos_init = np.array([float(ephemerides[10])*1000,float(ephemerides[14])*1000,float(ephemerides[18])*1000])
fasit = np.array([-11472.653194, -22189.176009, -5320.765921])*1000


#%%
pos, vel, clock_err, clock_rate_err = GLONASSEphemeris(ephemerides)

