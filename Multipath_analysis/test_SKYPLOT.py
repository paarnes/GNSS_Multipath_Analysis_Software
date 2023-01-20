import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import pickle
import pandas as pd


plt.rcParams['axes.axisbelow'] = True
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)
plt.rc('figure', figsize=(14, 9),dpi = 170)


file_to_read = open(r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\Multipath_analysis/Output_Files/analysisResults.pkl', "rb")
file_to_read = open(r'C:\Users\perhe\OneDrive\Documents\Python_skript\Output_Files_08012023/analysisResults.pkl', "rb")
analysisResults = pickle.load(file_to_read)
  
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
    # plt.title(title)
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
        ax.scatter(azimuth_rad, sat_el,label=PRN_)

    ax.set_rticks([10 ,20 ,30, 40, 50, 60, 70, 80, 90])  # Less radial ticks
    # ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    ax.tick_params(axis='both',labelsize=18,pad=7)
    ax.grid(True)
    ax.legend(fontsize=14,bbox_to_anchor=(1.40, 0.5),fancybox=True, shadow=True,ncol=2,loc='center right')
    filename = 'Skyplot_' + GNSSsystemName + '.png'
    # filename2 = 'Skyplot_' + GNSSsystemName + '.pdf'
    # fig.savefig(graph_dir + "/" + filename, dpi=300, orientation='landscape')
    # fig.savefig(filename2, orientation='landscape',bbox_inches='tight')
    # plt.show()
    
    return

azimut_currentSys = analysisResults['Sat_position']['G']['Azimut']
elevation_currentSys = analysisResults['Sat_position']['G']['Elevation']  
GNSSsystemName ='GPS'
graphDir = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\Multipath_analysis'
make_skyplot(azimut_currentSys,elevation_currentSys,GNSSsystemName,graphDir)

#%%
GNSS_Name2Code =  dict(zip(['GPS', 'GLONASS', 'Galileo', 'BeiDou'], ['G', 'R', 'E', 'C']))
for sys in analysisResults['GNSSsystems']:
    curr_sys = GNSS_Name2Code[sys]
    azimut_currentSys = analysisResults['Sat_position'][curr_sys]['Azimut']
    elevation_currentSys = analysisResults['Sat_position'][curr_sys]['Elevation'] 
    make_skyplot(azimut_currentSys,elevation_currentSys,sys,graphDir)

