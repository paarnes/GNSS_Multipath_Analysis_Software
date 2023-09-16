import numpy as np 
import numpy.matlib
import warnings
warnings.filterwarnings("ignore")

def barylag(data, time_diff):

    """
    Interpolates the given data using the Barycentric Lagrange Interpolation formula. 
    Vectorized to remove all loops. Based on work by Greg von Winckel.

    Parameters:
    ----------
    data      : numpy array. A two column vector where column one contains the
                nodes (time) and column two contains the satellite coordinate
                at the nodes (numpy array).

    time_diff : float. The time difference between first epoch and current
           
    Returns:
    -------        
    p - interpolated data. Column one is just the 
        fine mesh time_diff, and column two is interpolated data 
   

    """
    M = len(data)
    N = 1 # computing each  coordinate component seprately (X,Y,Z)
    if M ==0:
        p = np.array([np.nan]) # return np.nan if no data
        return p
    
    #  Compute the barycentric weights
    try:
        X = np.matlib.repmat(data[:,0].reshape(-1, 1), 1, M)
    except:
        X = np.array([])
    
    W = np.matlib.repmat(1/np.prod(np.transpose(X) - X  + np.eye(M), 1), N, 1) # Matrix of weights
    xdist = np.matlib.repmat(time_diff, 1, M) - np.matlib.repmat(data[:,0].reshape(len(data[:,0])), N, 1) # Get distances between nodes and interpolation points
    fixi,fixj = np.where(xdist==0) # Find all of the elements where the interpolation point is on a node
    xdist[fixi,fixj] = np.nan # Use NaNs as a place-holder
    H = W/xdist

    # Compute the interpolated polynomial
    p = H@data[:,1]/np.sum(H,1)
    # Replace NaNs with the given exact values. 
    p[fixi] = data[fixj,1]

    return p

