import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as _mptl
import matplotlib.cm as cm
from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os, sys

plt.rcParams['axes.axisbelow'] = True
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)
plt.rc('figure', figsize=(14, 9),dpi = 170)



def make_polarplot(multipath,sat_elevation,sat_azimut,GNSSsystemName,range1_code,graph_dir): 
    """
    The function makes a polar plot that shows the multipath effect as a function 
    of azimut and elevation angle."
    
    
    """
    
    ## -- Setting some arguments
    vmin     = 0          # Multipath minimum. Har med fargingen av multipathverdiene ift fargeskalaen
    # vmax     = 1.5        # Multipath max. Verdier på 0.8 meter og over får max fargemetning.
    vmax     = 2*np.nanmean(np.abs(multipath)) # Multipath max color by 2 times the means value
    cmap     = 'cividis'  # Fargeskalen
    # cmap = 'jet'
    dpi_fig  = 300        # Oppløsningen på figurene
    # dpi_fig  = 300        # Oppløsningen på figurene
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'},figsize=(16,12),dpi=170)
    fig.subplots_adjust(left=None, bottom=0.1, right=None, top=None, wspace=None, hspace=None)
    ax.set_rlim(bottom=90, top=0)
    _,num_sat = multipath.shape
    
    color = np.arange(1,0,-0.1)
    color = color.reshape(len(color),1)
    pc = ax.imshow(color,cmap = cmap,vmin = vmin, vmax = vmax,data = multipath, origin='upper', extent=[0,0,0,0])
    # cbar = fig.colorbar(pc, ax=ax, orientation='vertical',cmap = cmap,shrink=0.55,pad=.04,aspect=15)
    cbar = fig.colorbar(pc, ax=ax, orientation='vertical',shrink=0.55,pad=.04,aspect=15) #removed cmap due to warning 10.01.2023

    cbar.ax.set_title('MP[m]',fontsize=18,pad=15)
    cbar.ax.tick_params(labelsize=18) 
        
    ax.set_rticks([10 ,20 ,30, 40, 50, 60, 70, 80, 90])  # Less radial ticks
    # ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    ax.tick_params(axis='both',labelsize=18,pad=4)
    ax.grid(True)
    ax.set_theta_zero_location("N")  # theta=0 at the top
    ax.set_theta_direction(-1)  # theta increasing clockwise
    ax.set_title("Polarplott av flerveis-interferens som funksjon av asimut og elevasjonsvinkel \n Signal: %s (%s)" % (range1_code, GNSSsystemName), va='bottom',fontsize=28)
    # plt.xlabel('Dato: {}'.format(curr_date),fontsize=22,labelpad =14)
    for PRN in range(0,num_sat):
        mp_est = multipath[:,PRN]
        if all(np.isnan(mp_est)) == True:
            continue
        else:
            sat_el = sat_elevation[:,PRN]
            sat_az = sat_azimut[:,PRN]
            c = list(abs(mp_est))
            c = cm.jet((c-np.min(c))/(np.max(c)-np.min(c)))
            ax.scatter(np.radians(sat_az), sat_el, c = abs(mp_est), cmap = cmap ,marker = 'o', s = 100, edgecolor='none',vmin = vmin, vmax = vmax)
        
    # plt.show()
    # plt.ioff()
    filename = 'MP_GPS_L1' + "_" + GNSSsystemName + "_" + range1_code + '.png'
    fig.savefig(graph_dir + "/" + filename, dpi=dpi_fig, orientation='landscape')
    plt.close()


             
        


