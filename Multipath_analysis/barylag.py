# import numpy as np 
# import numpy.matlib
def barylag(data,x):

    """
    #
    # barylag.m
    #
    # Interpolates the given data using the Barycentric
    # Lagrange Interpolation formula. Vectorized to remove all loops
    #
    # data - a two column vector where column one contains the
    #        nodes and column two contains the function value 
    #        at the nodes
    # p - interpolated data. Column one is just the 
    #     fine mesh x, and column two is interpolated data 
    #
    # Reference:
    #
    # (1) Jean-Paul Berrut & Lloyd N. Trefethen, "Barycentric Lagrange 
    #     Interpolation" 
    #     http://web.comlab.ox.ac.uk/oucl/work/nick.trefethen/berrut.ps.gz
    # (2) Walter Gaustschi, "Numerical Analysis, An Introduction" (1997) pp. 94-95
    #
    #
    # Written by: Greg von Winckel       03/07/04
    # Contact:    gregvw@chtm.unm.edu 
    #
    # Copyright (c) 2009, Greg von Winckel
    # All rights reserved.
    
    """
    import numpy as np 
    import numpy.matlib
    import warnings
    warnings.filterwarnings(action='ignore', message='invalid value encountered in fmod')
    
    M = len(data)
    
    if M ==0:
        p = np.array([np.nan]) ## added this 08.12.2022 to try to speed up the algoritm
        return p
    
    # N = len(x)
    N = 1 # length i always 1 ???????????
    
    
    # Compute the barycentric weights
    # if not data:
    #     X = np.array([])
    # else:
    #     X=np.matlib.repmat(data[:,0].reshape(len(data[:,0]),1),1,M)
    
    try:
        X = np.matlib.repmat(data[:,0].reshape(len(data[:,0]),1),1,M)
    except:
        X = np.array([])
    
    # matrix of weights
    W = np.matlib.repmat(1/np.prod(np.transpose(X) - X  + np.eye(M),1),N,1);  #MATLAB USE TRANSPOSE HERE
    
    # Get distances between nodes and interpolation points
    # xdist=np.matlib.repmat(x,1,M)-np.matlib.repmat(data[:,0].reshape(len),N,1)

    xdist = np.matlib.repmat(x,1,M) - np.matlib.repmat(data[:,0].reshape(len(data[:,0])),N,1)
    
    # Find all of the elements where the interpolation point is on a node
    fixi,fixj = np.where(xdist==0)
    
    
    ## -- Use NaNs as a place-holder
    xdist[fixi,fixj] = np.nan # hvis feilmeld, gjør om xdist til float
    H = W/xdist
    
    ## -- Compute the interpolated polynomial
    p=(H@data[:,1])/np.sum(H,1) # added @ 20.11.2022
    
    ## -- Replace NaNs with the given exact values. 
    p[fixi] = data[fixj,1]
    return p

# x = 29.999999046325684
# test = np.array([[0, 1.504235654600000009e+07],[3.000000000000000000e+02 ,1.552872639699999988e+07]])
# barylag(test,x)