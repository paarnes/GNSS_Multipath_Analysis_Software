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
    
    TEST UT Å LEGG fjern marking "o" og bruk heller sammentrukken linje. 
    
    
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


def make_skyplot(azimut_currentSys, elevation_currentSys, GNSSsystemName,graph_dir):
    """
    Generates a skyplot for GPS based on azimuth and elevation angles.
    azimuth: list of azimuth angles in degrees
    elevation: list of elevation angles in degrees
    title: title of the skyplot
    """
    GNSS_Name2Code =  dict(zip(['GPS', 'GLONASS', 'Galileo', 'BeiDou'], ['G', 'R', 'E', 'C']))
    sys_code = GNSS_Name2Code[GNSSsystemName]
    num_sat = azimut_currentSys.shape[1]
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'},figsize=(16,10),dpi=180)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rlim(bottom=90, top=0)
    ax.set_title("Skyplot for %s" % (GNSSsystemName), va='bottom',fontsize=28)
    for PRN in np.arange(0,num_sat):
        sat_el = elevation_currentSys[:,PRN]

        if all(np.isnan(sat_el)) == True or all(sat_el)==0:
            continue
        else:
            sat_az = azimut_currentSys[:,PRN]
        # Convert azimut angles to radians
        azimuth_rad = np.deg2rad(sat_az)        
        # Plot the satellite positions on the skyplot
        PRN_ = sys_code+str(PRN)
        PRN_ = sys_code+str(PRN).zfill(2)
        # ax.scatter(azimuth_rad, sat_el,label=PRN_)
        line = ax.plot(azimuth_rad, sat_el,label=PRN_,linewidth=5.5)

    ax.set_rticks([10 ,20 ,30, 40, 50, 60, 70, 80, 90])  # Less radial ticks
    # ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    ax.tick_params(axis='both',labelsize=18,pad=7)
    ax.grid(True)
    # ax.legend(fontsize=14,bbox_to_anchor=(1.40, 0.5),fancybox=True, shadow=True,ncol=2,loc='center right')
    legend = ax.legend(fontsize=14,bbox_to_anchor=(1.40, 0.5),fancybox=True, shadow=True,ncol=2,loc='center right')
    ## Set the linewidth of each legend object (then not dependent of linewith in plot)
    for legobj in legend.legendHandles:
        legobj.set_linewidth(3.5)
    filename = 'Skyplot_' + GNSSsystemName + '.png'
    filename2 = 'Skyplot_' + GNSSsystemName + '.pdf'
    # fig.savefig(graph_dir + "/" + filename, dpi=300, orientation='landscape')
    fig.savefig(graph_dir + "/" + filename2, orientation='landscape',bbox_inches='tight')
    
    return
        


