import warnings
import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
from matplotlib.ticker import MaxNLocator
from .plotResults import set_linewidt_for_each_object

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)

plt.rcParams['axes.axisbelow'] = True
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)
plt.rc('figure', figsize=(14, 9),dpi = 170)

def make_polarplot(analysisResults, graph_dir):
    """
    The function makes a polar plot that shows the multipath effect as a function
    of azimut and elevation angle."
    """
    GNSS_Name2Code =  dict(zip(['GPS', 'GLONASS', 'Galileo', 'BeiDou'], ['G', 'R', 'E', 'C']))

    for system in analysisResults['GNSSsystems']:
        curr_sys = GNSS_Name2Code[system]
        bands_curr_sys = analysisResults[system]['Bands']
        try:
            sat_elevation = analysisResults['Sat_position'][curr_sys]['elevation']
            sat_azimut = analysisResults['Sat_position'][curr_sys]['azimuth']
        except:
            logger.warning("INFO(GNSS_MultipathAnalysis): Polarplot of multipath is not possible for %s. Satellite azimuth and elevation angles are missing.", system)
            continue


        vmax_list = []
        ## -- Finding larges mean multipath value for scale on cbar (vmax)
        for band in bands_curr_sys:
            codes_curr_sys = analysisResults[system][band]['Codes']
            codes_curr_sys = [ele for ele in codes_curr_sys if ele != []] # removing empty list if exist
            for code in codes_curr_sys:
                if code not in list(analysisResults[system][band].keys()):
                    continue
                else:
                    multipath = analysisResults[system][band][code]['multipath_range1']
                    vmax_list.append(round(2*np.nanmean(np.abs(multipath)),1)) # Multipath max color by 2 times the means value
        ## --Do the plotting
        for band in bands_curr_sys:
            codes_curr_sys = analysisResults[system][band]['Codes']
            codes_curr_sys = [ele for ele in codes_curr_sys if ele != []] # removing empty list if exist
            for code in codes_curr_sys:
                if code not in list(analysisResults[system][band].keys()):
                    continue
                else:
                    multipath = analysisResults[system][band][code]['multipath_range1']
                    range1_code = code
                    ## -- Setting some arguments
                    vmin     = 0          # Multipath minimum. Har med fargingen av multipathverdiene ift fargeskalaen
                    # vmax     = 1.5        # Multipath max. Verdier på 0.8 meter og over får max fargemetning.
                    vmax     = np.max(vmax_list) # Multipath max color by 2 times the means value
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
                    kwargs = {'format': '%.1f'}
                    cbar = fig.colorbar(pc, ax=ax, orientation='vertical',shrink=0.55,pad=.04,aspect=15,**kwargs) #removed cmap due to warning 10.01.2023
                    # cbar = fig.colorbar(pc, ax=ax, orientation='vertical',shrink=0.55,pad=.04,aspect=15) #removed cmap due to warning 10.01.2023
                    cbar.ax.set_title('MP[m]',fontsize=18,pad=15)
                    cbar.ax.tick_params(labelsize=18)
                    ax.set_rticks([10 ,20 ,30, 40, 50, 60, 70, 80, 90])  # Less radial ticks
                    # ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
                    ax.tick_params(axis='both',labelsize=18,pad=4)
                    ax.grid(True)
                    ax.set_theta_zero_location("N")  # theta=0 at the top
                    ax.set_theta_direction(-1)  # theta increasing clockwise
                    ax.set_title("Polar plot of the multipath effect as funtion of azimut and elevation angle for \n Signal: %s (%s)" % (range1_code, system), va='bottom',fontsize=28)
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

                    filename = 'MP_' + system + "_" + range1_code + '.png'
                    fig.savefig(graph_dir + "/" + filename, dpi=dpi_fig, orientation='landscape', bbox_inches='tight')
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
        line = ax.plot(azimuth_rad, sat_el,label=PRN_,linewidth=5.5, solid_capstyle='round') #solid_capstyle='round' makes rounded edges on lines

    ax.set_rticks([10 ,20 ,30, 40, 50, 60, 70, 80, 90])  # Less radial ticks
    # ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    ax.tick_params(axis='both',labelsize=18,pad=7)
    ax.grid(True)
    # ax.legend(fontsize=14,bbox_to_anchor=(1.40, 0.5),fancybox=True, shadow=True,ncol=2,loc='center right')
    legend = ax.legend(fontsize=14,bbox_to_anchor=(1.40, 0.5),fancybox=True, shadow=True,ncol=2,loc='center right')
    ## Set the linewidth of each legend object (then not dependent of linewith in plot)
    set_linewidt_for_each_object(legend, 3.5)

    filename = 'Skyplot_' + GNSSsystemName + '.png'
    filename2 = 'Skyplot_' + GNSSsystemName + '.pdf'
    # fig.savefig(graph_dir + "/" + filename, dpi=300, orientation='landscape')
    fig.savefig(graph_dir + "/" + filename2, orientation='landscape',bbox_inches='tight')

    return


def make_polarplot_SNR(analysisResults, GNSS_obs,GNSSsystems, obsCodes, graphDir):
    """
    The function makes a polar plot that shows the Signal to noise ration (SNR) a function
    of azimut and elevation angle."
    """
    SNR_obs = {}
    for sys in GNSS_obs.keys():
        GNSSsystemIndex = [k for k in GNSSsystems if GNSSsystems[k] == sys][0]
        SNR_codes = [SNR_code for SNR_code in obsCodes[GNSSsystemIndex][sys] if 'S' in SNR_code[0]]
        SNR_obs[sys] = {}
        for SNR_code in SNR_codes:
            SNR_idx = obsCodes[GNSSsystemIndex][sys].index(SNR_code)
            SNR_data =  np.array([epoch[:, SNR_idx] for epoch in GNSS_obs[sys].values()])
            SNR_data[SNR_data ==0] = np.nan
            SNR_obs[sys][SNR_code] = SNR_data

    GNSS_Name2Code =  dict(zip(['GPS', 'GLONASS', 'Galileo', 'BeiDou'], ['G', 'R', 'E', 'C']))

    for system in analysisResults['GNSSsystems']:
        curr_sys = GNSS_Name2Code[system]
        try:
            sat_elevation = analysisResults['Sat_position'][curr_sys]['elevation']
            sat_azimut = analysisResults['Sat_position'][curr_sys]['azimuth']
        except:
            logger.warning(f"INFO(GNSS_MultipathAnalysis): Polarplot of SNR is not possible for {system}. Satellite azimuth and elevation angles are missing." )
            continue
        vmax_list = []
        ## -- Finding larges mean multipath value for scale on cbar (vmax)
        # for code in SNR_obs[curr_sys].keys():
        #     SNR = SNR_obs[curr_sys][code]
        #     vmax_list.append(round(1.5*np.nanmean(np.abs(SNR)),1)) # SNR max color by 2 times the means value
        ## --Do the plotting
        for code in SNR_obs[curr_sys].keys():
            SNR = SNR_obs[curr_sys][code]
            if np.all(np.isnan(SNR)):
                logger.warning(f"INFO(GNSS_MultipathAnalysis): Polarplot of SNR not possible for {code} for {system}. The RINEX file does not contain data for this code for this system." )
                continue
            range1_code = code
            ## -- Setting some arguments
            # vmin     = 0  # Multipath minimum. Har med fargingen av multipathverdiene ift fargeskalaen
            vmin     = round(np.nanmin(np.abs(SNR)),0)  # Multipath minimum. Har med fargingen av multipathverdiene ift fargeskalaen
            vmax     = round(1.05*np.nanmax(np.abs(SNR)),0)
            # cmap     = 'cividis'  # Fargeskalen
            cmap       = 'jet'
            dpi_fig  = 300        # Oppløsningen på figurene
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'},figsize=(16,12),dpi=170)
            fig.subplots_adjust(left=None, bottom=0.1, right=None, top=None, wspace=None, hspace=None)
            ax.set_rlim(bottom=90, top=0)
            _,num_sat = SNR.shape

            color = np.arange(1,0,-0.1)
            color = color.reshape(len(color),1)
            pc = ax.imshow(color,cmap = cmap,vmin = vmin, vmax = vmax,data = SNR, origin='upper', extent=[0,0,0,0])
            num_ticks = 8
            tick_locs = MaxNLocator(nbins=num_ticks).tick_values(vmin, vmax)
            cbar = fig.colorbar(pc, ax=ax, orientation='vertical', ticks=tick_locs,shrink=0.55,pad=.04,aspect=15) #removed cmap due to warning 10.01.2023
            cbar.ax.set_title('SNR[dB-Hz]',fontsize=18,pad=15)
            cbar.ax.tick_params(labelsize=18)

            ## -- set color
            ax.set_facecolor('#373B44')
            # ax.set_facecolor((200/255, 270/255, 200/255))
            # set ticks color
            ax.tick_params(colors="white", labelcolor="white")
            ax.tick_params(axis="y", colors="white")
            ax.tick_params(axis="x", colors="black") # want the ticks outside the facecolor to be black

            ax.set_rticks([10 ,20 ,30, 40, 50, 60, 70, 80, 90])  # Less radial ticks
            # # ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
            ax.tick_params(axis='both',labelsize=18,pad=4)


            ax.grid(True)
            ax.set_theta_zero_location("N")  # theta=0 at the top
            ax.set_theta_direction(-1)  # theta increasing clockwise


            ax.set_title("Polar plot of the SNR as funtion of azimut and elevation angle for \n Signal: %s (%s)" % (range1_code, system), va='bottom',fontsize=28)
            for PRN in range(0,num_sat):
                SNR_est = SNR[:,PRN]
                if all(np.isnan(SNR_est)) == True:
                    continue
                else:
                    sat_el = sat_elevation[:,PRN]
                    sat_az = sat_azimut[:,PRN]
                    c = list(abs(SNR_est))
                    c = cm.jet((c-np.min(c))/(np.max(c)-np.min(c)))
                    ax.scatter(np.radians(sat_az), sat_el, c = abs(SNR_est), cmap = cmap ,marker = 'o', s = 100, edgecolor='none',vmin = vmin, vmax = vmax)
            filename = 'SNR_Polar_' + system + "_" + range1_code + '.png'
            fig.savefig(graphDir + "/" + filename, dpi=dpi_fig, orientation='landscape',bbox_inches='tight')
            plt.close()




def plot_SNR_wrt_elev(analysisResults,GNSS_obs, GNSSsystems, obsCodes, graphDir,tInterval):
    """
    Funtion that makes a subplot of the Signal to nosie ration (SNR) wrt to
    time and the satellites elevation angle.
    """
    SNR_obs = {}
    for sys in GNSS_obs.keys():
        GNSSsystemIndex = [k for k in GNSSsystems if GNSSsystems[k] == sys][0]
        SNR_codes = [SNR_code for SNR_code in obsCodes[GNSSsystemIndex][sys] if 'S' in SNR_code[0]]
        SNR_obs[sys] = {}
        for SNR_code in SNR_codes:
            SNR_idx = obsCodes[GNSSsystemIndex][sys].index(SNR_code)
            SNR_data =  np.array([epoch[:, SNR_idx] for epoch in GNSS_obs[sys].values()])
            SNR_data[SNR_data ==0] = np.nan
            SNR_obs[sys][SNR_code] = SNR_data

    GNSS_Name2Code =  dict(zip(['GPS', 'GLONASS', 'Galileo', 'BeiDou'], ['G', 'R', 'E', 'C']))
    for system in analysisResults['GNSSsystems']:
        curr_sys = GNSS_Name2Code[system]
        try:
            sat_elevation = analysisResults['Sat_position'][curr_sys]['elevation']
            sat_elevation[(sat_elevation < 0) | (sat_elevation > 90)] = np.nan # removes elevation angels when the sat not visable (below the horizon)
        except:
            logger.warning(f"INFO(GNSS_MultipathAnalysis): Plot of SNR wrt elevation angle is not possible for {system}. Satellite azimuth and elevation angles are missing." )
            continue
        for code in SNR_obs[curr_sys].keys():
            SNR = SNR_obs[curr_sys][code]
            if np.all(np.isnan(SNR)):
                logger.warning(f"INFO(GNSS_MultipathAnalysis): Plot of SNR wrt elevation is not possible for {code} for {system}. The RINEX file does not contain data for this code for this system." )
                continue
            range1_code = code
            ## -- Setting some arguments
            dpi_fig  = 300        # Oppløsningen på figurene
            # fig, ax = plt.subplots(figsize=(16,12),dpi=170)
            fig, ax = plt.subplots(nrows=2, ncols=1,sharex=False, squeeze=True,figsize=(16,11),dpi=160)
            fig.subplots_adjust(left=0.07, bottom=0.1, right=0.78, top=0.91, wspace=None, hspace=0.45)
            num_ep,num_sat = SNR.shape
            ## -- Time stamps
            t = np.arange(1,num_ep+1)*tInterval/60**2 # Convert to hours

            ## -- Subplot 1
            for PRN in range(1,num_sat):
                if not np.isnan(SNR[:,PRN]).all():
                    ax[0].plot(t,SNR[:,PRN], label='PRN%s' % (PRN),linewidth=0.7)

            ax[0].grid(True,linewidth=0.3)
            ax[0].set_xlim([0,t[-1]])
            ax[0].set_ylim(0,np.nanmax(SNR)+10)
            ax[0].set_title("Signal to noise ratio (SNR) as funtion of time for \n Signal: %s (%s)" % (range1_code, system), va='bottom',fontsize=22)
            ax[0].set_xlabel('Time $[h]$',fontsize=18,labelpad=10)
            ax[0].set_ylabel('[dB-Hz]',fontsize=18,labelpad=10)
            ax[0].tick_params(axis='both', labelsize=16)
            legend = ax[0].legend(loc='center right',fontsize=12,bbox_to_anchor=(1.25, 0.5), fancybox=True, shadow=True,ncol=2) # frame = legend.get_frame(); frame.set_facecolor((0.89701,0.79902,0.68137)); frame.set_edgecolor('black') #legend
            set_linewidt_for_each_object(legend, 1.5)


            ## -- Subplot 2
            ax[1].grid(True,linewidth=0.3)
            ax[1].set_xlim([0,90])
            ax[1].set_title("Signal to noise ratio (SNR) as funtion of elevation angle for \n Signal: %s (%s)" % (range1_code, system), va='bottom',fontsize=22)
            ax[1].set_ylim(0,np.nanmax(SNR)+10)
            ax[1].set_xlabel('Elevation angle $[^{\circ}]$',fontsize=18,labelpad=10)
            ax[1].set_ylabel('[dB-Hz]',fontsize=18,labelpad=10)
            ax[1].tick_params(axis='both', labelsize=16)
            for PRN in range(0,num_sat):
                SNR_PRN = SNR[:,PRN]
                if all(np.isnan(SNR_PRN)) == True:
                    continue
                else:
                    sat_el = sat_elevation[:,PRN]
                    ax[1].plot(sat_el, SNR_PRN, label='PRN%s' % (PRN),linewidth=0.7)
            legend = ax[1].legend(loc='center right',fontsize=12,bbox_to_anchor=(1.25, 0.5), fancybox=True, shadow=True,ncol=2) # frame = legend.get_frame(); frame.set_facecolor((0.89701,0.79902,0.68137)); frame.set_edgecolor('black') #legend
            set_linewidt_for_each_object(legend, 1.5)


            filename = 'SNR_' + system + "_" + range1_code + '.pdf'
            fig.savefig(graphDir + "/" + filename, orientation='landscape',bbox_inches='tight')
            plt.close()


