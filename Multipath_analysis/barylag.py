import numpy as np 
import numpy.matlib
import warnings
warnings.filterwarnings(action='ignore', message='invalid value encountered in fmod')

def barylag(data,x):

    """
    Interpolates the given data using the Barycentric
    Lagrange Interpolation formula. Vectorized to remove all loops
   
    data - a two column vector where column one contains the
           nodes and column two contains the function value 
           at the nodes
    p - interpolated data. Column one is just the 
        fine mesh x, and column two is interpolated data 
   
    Reference:
   
    (1) Jean-Paul Berrut & Lloyd N. Trefethen, "Barycentric Lagrange 
        Interpolation" 
        http://web.comlab.ox.ac.uk/oucl/work/nick.trefethen/berrut.ps.gz
    (2) Walter Gaustschi, "Numerical Analysis, An Introduction" (1997) pp. 94-95
   
   
    Written by: Greg von Winckel       03/07/04
    Contact:    gregvw@chtm.unm.edu 
   
    Copyright (c) 2009, Greg von Winckel
    All rights reserved.
    
    """
    
    M = len(data)
    
    if M ==0:
        p = np.array([np.nan]) ## added this 08.12.2022 to try to speed up the algoritm
        return p
    
    # N = len(x)
    N = 1 # length i always 1 ???????????
    
    
    ## -- Compute the barycentric weights
    try:
        X = np.matlib.repmat(data[:,0].reshape(len(data[:,0]),1),1,M)
    except:
        X = np.array([])
    
    ## -- Matrix of weights
    W = np.matlib.repmat(1/np.prod(np.transpose(X) - X  + np.eye(M),1),N,1);  #MATLAB USE TRANSPOSE HERE
    
    ## -- Get distances between nodes and interpolation points
    xdist = np.matlib.repmat(x,1,M) - np.matlib.repmat(data[:,0].reshape(len(data[:,0])),N,1)
    
    ## -- Find all of the elements where the interpolation point is on a node
    fixi,fixj = np.where(xdist==0)
    
    
    ## -- Use NaNs as a place-holder
    xdist[fixi,fixj] = np.nan # hvis feilmeld, gjør om xdist til float
    H = W/xdist
    
    ## -- Compute the interpolated polynomial
    p=(H@data[:,1])/np.sum(H,1) # added @ 20.11.2022
    
    ## -- Replace NaNs with the given exact values. 
    p[fixi] = data[fixj,1]
    return p



# from scipy.interpolate import lagrange
# from datetime import datetime, timedelta

# def interpolate_ecef_datetime(sp3_data, time):
#     """
#     Interpolates ECEF coordinates from SP3 data using the Barycentric Lagrange Interpolation formula.
    
#     Parameters:
#         - sp3_data (List[Tuple[datetime, float, float, float]]): List of tuples containing the datetime, x, y, z coordinates of the satellite
#         - time (datetime): The datetime for which the ECEF coordinates are to be interpolated
    
#     Returns:
#         - Tuple[float, float, float]: Interpolated ECEF coordinates (x, y, z) at the given datetime
#     """
#     times, xs, ys, zs = zip(*sp3_data)
#     times_seconds = [ (time-datetime(1970,1,1)).total_seconds() for time in times] 
#     time_seconds = (time-datetime(1970,1,1)).total_seconds()
#     poly_x = lagrange(times_seconds, xs)
#     poly_y = lagrange(times_seconds, ys)
#     poly_z = lagrange(times_seconds, zs)
#     return poly_x(time_seconds), poly_y(time_seconds), poly_z(time_seconds)


# sp3_data = [
#     (datetime(2022, 1, 1, 0, 0, 0), 13882271.938,-21710006.216 ,5357125.462),
#     (datetime(2022, 1, 1, 1, 0, 0), 13882271.938,-21710006.216 ,5357125.462),
#     (datetime(2022, 1, 1, 2, 0, 0), 13882271.938,-21710006.216 ,5357125.462),
#     (datetime(2022, 1, 1, 3, 0, 0), 13882271.938,-21710006.216 ,5357125.462),
#     (datetime(2022, 1, 1, 4, 0, 0), 13882271.938,-21710006.216 ,5357125.462),
# ]

# koord = interpolate_ecef(sp3_data, 300)