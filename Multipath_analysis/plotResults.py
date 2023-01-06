def plotResults(ion_delay_phase1, multipath_range1, sat_elevation_angles,\
    tInterval, currentGNSSsystem, range1_Code, range2_Code, phase1_Code, phase2_Code, graphDir):
    """
    # Function that plots the results of estimates of delay on signals. The
    # following plots are made:
    #                           ionosphere delay on phase 1 signal vs time
    
    #                           zenit mapped ionosphere delay on phase 1 
    #                           signal vs time
    
    #                           multipath delay on range 1 signal vs time
    
    #                           multipath delay on range 1 signal vs elevation
    #                           angle
    
    #                           a combined plot of both multipath plots. This
    #                           plot is also saved
    #--------------------------------------------------------------------------------------------------------------------------
    # INPUTS
    
    # ion_delay_phase1:     matrix containing estimates of ionospheric delays 
    #                       on the first phase signal for each PRN, at each epoch.
    
    #                       ion_delay_phase1(epoch, PRN)
    
    # multipath_range1:     matrix containing estimates of multipath delays 
    #                       on the first range signal for each PRN, at each epoch.
    
    #                       multipath_range1(epoch, PRN)
    
    # sat_elevation_angles: Array contaning satellite elevation angles at each
    #                       epoch, for current GNSS system. 
    
    #                       sat_elevation_angles(epoch, PRN)
    
    # tInterval:            observations interval; seconds. 
    
    # currentGNSSsystem:    string containing code for current GNSS system,
    #                       ex. "G" or "E"
    
    # range1_Code:          string, defines the observation type that will be 
    #                       used as the first range observation. ex. "C1X", "C5X"    
    
    # range2_Code:          string, defines the observation type that will be 
    #                       used as the second range observation. ex. "C1X", "C5X"
    
    # phase1_Code:          string, defines the observation type that will be 
    #                       used as the first phase observation. ex. "L1X", "L5X"
    
    # phase2_Code:          string, defines the observation type that will be 
    #                       used as the second phase observation. ex. "L1X", "L5X"
    
    # graphDir:             string. Gives path to where figure the combined
    #                       multipath figure should be saved.
    #--------------------------------------------------------------------------------------------------------------------------
    """
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from matplotlib import rc
    matplotlib.use('Agg') # dont want the plots to be displayed.
    
    n,m = ion_delay_phase1.shape
    plt.rcParams['axes.axisbelow'] = True
    rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    rc('text', usetex=True)
    plt.rc('figure', figsize=(14, 9),dpi = 170)
    # plt.rcParams['axes.facecolor'] = (0.95294, 0.92157, 0.87058)
    # graphDir = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Multipath_analysis' #FJERN!!!
    
    ## -- Map mapping GNSS system code to full name
    GNSSsystemName_map = dict(zip(['G', 'R', 'E', 'C'], ['GPS', 'GLONASS', 'Galileo', 'BeiDou']))
    GNSSsystemName = GNSSsystemName_map[currentGNSSsystem]
    
    ## -- Time stamps
    # t = np.arange(1,n+1)*tInterval # seconds
    t = np.arange(1,n+1)*tInterval/60**2 # Convert to hours
    
    ## ---------- Plotting ionospheric delay vs time --------------------------
    fig1, ax1 = plt.subplots(nrows=1, ncols=1,sharex=True, squeeze=True,figsize=(16,9),dpi=170)
    fig1.subplots_adjust(left=0.07, bottom=0.1, right=0.78, top=None, wspace=None, hspace=None)

    for PRN in range(1,m):
        if not np.isnan(ion_delay_phase1[:,PRN]).all():
            ax1.plot(t, ion_delay_phase1[:,PRN], label='PRN%s' % (PRN)) 
    ax1.set_title('Ionospheric delay vs time for %s, \n Signal 1: %s  Signal 2: %s' % (GNSSsystemName, phase1_Code, phase2_Code),fontsize=22)
    ax1.set_xlabel('Time $[h]$',fontsize=16,labelpad=10)
    ax1.set_ylabel('$[m]$',fontsize=16,labelpad=10)
    ax1.tick_params(axis='both', labelsize=16)
    legend = ax1.legend(loc='center right',fontsize=12,bbox_to_anchor=(1.25, 0.5), fancybox=True, shadow=True,ncol=2) # frame = legend.get_frame(); frame.set_facecolor((0.89701,0.79902,0.68137)); frame.set_edgecolor('black') #legend
    ax1.grid(color='k', linestyle='-', linewidth=0.1)         
    ax1.axhline(y=0.0, color='k', linestyle='-',linewidth=1)
    # plt.show()
    # fig1_name  = GNSSsystemName + "_" + 'ionospheric_delay' + '.png'
    # fig1.savefig(graphDir + "/" +  fig1_name, dpi=300)
    fig1_name  = GNSSsystemName + "_" + 'ionospheric_delay' + '.pdf'
    fig1.savefig(graphDir + "/" +  fig1_name,bbox_inches='tight')
    plt.close()
    
    
    
    
    ## ---------- Plot zenit mapped ionosphere delay vs time ------------------
    fig2, ax2 = plt.subplots(nrows=1, ncols=1,sharex=True, squeeze=True,figsize=(16,9),dpi=150)
    fig2.subplots_adjust(left=0.05, bottom=0.1, right=0.8, top=None, wspace=None, hspace=None)

    ion_delay_phase1_zenit = ion_delay_phase1*np.cos(np.arcsin(6371000/(6371000 + 450000) * np.sin((90-sat_elevation_angles)*np.pi/180)))
    excluded_PRN = []
    for PRN in range(1,m):
        ## only plot for PRN that have any observations
        if not np.isnan(ion_delay_phase1[:,PRN]).all():
            if not np.isnan(sat_elevation_angles[:, PRN]).all():
                ax2.plot(t, ion_delay_phase1_zenit[:, PRN], label='PRN%s' % (str(PRN)),linewidth=2) 
            else:
                excluded_PRN.append(PRN)
                
    if len(excluded_PRN) == 0:
        ax2.set_title('Zenit mapped ionosphere delay for %s, \n Signal 1: %s  Signal 2: %s' % (GNSSsystemName, phase1_Code, phase2_Code),fontsize=22)
    else:
        ax2.set_title('Zenit mapped ionosphere delay for %s,\
                      %s \n  Signal 1: %s  Signal 2: %s \n  Ekskluderte PRN: %s' %(GNSSsystemName,range1_Code,\
                          phase1_Code, phase2_Code, excluded_PRN),fontsize=22)

    ax2.set_xlabel('Time $[h]$',fontsize=16,labelpad=10)
    ax2.set_ylabel('$[m]$',fontsize=16,labelpad=10)
    ax2.tick_params(axis='both', labelsize=16)
    legend = ax2.legend(loc='center right',fontsize=12,bbox_to_anchor=(1.23, 0.5), fancybox=True, shadow=True,ncol=2) # frame = legend.get_frame(); frame.set_facecolor((0.89701,0.79902,0.68137)); frame.set_edgecolor('black') #legend
    ax2.grid(color='k', linestyle='-', linewidth=0.1)         
    ax2.axhline(y=0.0, color='k', linestyle='-',linewidth=1)
    # plt.show()
    # fig2_name  = GNSSsystemName + "_" + 'ionospheric_delay_zenit_mapped' + '.png'
    # fig2.savefig(graphDir + "/" +  fig2_name, dpi=300)
    fig2_name  = GNSSsystemName + "_" + 'ionospheric_delay_zenit_mapped' + '.pdf'
    # fig2.savefig(graphDir + "/" +  fig2_name)
    fig2.savefig(graphDir + "/" +  fig2_name,bbox_inches='tight')
    plt.close()

    
    
    
    
    
    ## ------------------- Plot multipath delay on range 1 signal vs time --------------------
    fig3, ax3 = plt.subplots(nrows=1, ncols=1,sharex=True, squeeze=True,figsize=(16,9),dpi=170)
    fig3.subplots_adjust(left=0.07, bottom=0.1, right=0.78, top=None, wspace=None, hspace=None)

    for PRN in range(1,m):
        if not np.isnan(multipath_range1[:,PRN]).all():
            ax3.plot(t, multipath_range1[:,PRN], label='PRN%s' % (PRN),linewidth=0.7)


    ## -- Crop figure by seting  y lim to mean values pluss minus 7 std
    y_mean = np.nanmean(multipath_range1)
    y_std  = np.nanstd(multipath_range1)
    ax3.set_xlim(0,t[-1])
    ax3.set_ylim(y_mean - 7*y_std, y_mean+7*y_std)
    ax3.set_title('Multipath effect vs time for the signal %s,\
                  %s \n Signal combination: %s-%s-%s' % (range1_Code, GNSSsystemName,range1_Code, phase1_Code, phase2_Code),fontsize=28)
    ax3.set_xlabel('Time $[h]$',fontsize=20,labelpad=10)
    ax3.set_ylabel('$[m]$',fontsize=20,labelpad=10)
    ax3.tick_params(axis='both', labelsize=18)
    legend = ax3.legend(loc='center right',fontsize=14,bbox_to_anchor=(1.28, 0.5), fancybox=True, shadow=True,ncol=2) # frame = legend.get_frame(); frame.set_facecolor((0.89701,0.79902,0.68137)); frame.set_edgecolor('black') #legend
    ## Set the linewidth of each legend object (then not dependent of linewith in plot)
    for legobj in legend.legendHandles:
        legobj.set_linewidth(1.5)

    ax3.grid(color='k', linestyle='-', linewidth=0.1)         
    ax3.axhline(y=0.0, color='k', linestyle='-',linewidth=0.4)
    # plt.show()
    
    filename  = '%s_%s_%s_MP_time.pdf' % (GNSSsystemName, range1_Code, range2_Code)
    filename2 = '%s_%s_%s_MP_time.png' % (GNSSsystemName, range1_Code, range2_Code)
    full_filename = graphDir + '/' + filename
    full_filename2 = graphDir + '/' + filename2

    fig3.savefig(graphDir + "/" +  filename)
    fig3.savefig(graphDir + "/" +  filename2, dpi=300)
    plt.close()




    ## ----- Plot multipath delay on range 1 signal vs elevation angles -----
    """
    SJEKK UT HVORFOR FØRSTE OG SISTE OPS ER NAN FOR SATELITTE ELEVASJONVINKENELE
    """
    excluded_PRN = []
    fig4, ax4 = plt.subplots(nrows=1, ncols=1,sharex=True, squeeze=True,figsize=(16,9),dpi=170)
    fig4.subplots_adjust(left=0.07, bottom=0.1, right=0.78, top=None, wspace=None, hspace=None)
    
    # sat_elevation_angles = sat_elevation_angles[1:-1,:] ## DETTE MÅ FJERNES. FINN UT HVORFOR LENGDEN IKKE STEMMER! MANGLER TO EPOKER i sat_elevation_angles
    # multipath_range1 = multipath_range1[1:-1,:] ## DETTE MÅ FJERNES. FINN UT HVORFOR LENGDEN IKKE STEMMER! MANGLER TO EPOKER i sat_elevation_angles
    # t = t[1:-1] # Matlab er det null der hvor de ikke er satellitt.
   
    for PRN in range(1,m):
        epoch_missing_sat_ele = np.where(np.isnan(sat_elevation_angles[:,PRN]))
        multipath_range1[epoch_missing_sat_ele,PRN] = np.nan
        
   
    for PRN in range(1,m):
        # only plot for PRN that have any observations
        if not np.isnan(multipath_range1[:,PRN]).all(): 
            if not np.isnan(sat_elevation_angles[:, PRN]).all(): # change to all
                ax4.plot(sat_elevation_angles[:, PRN], multipath_range1[:,PRN],  label='PRN%s' % (PRN), linewidth= 0.7) 
            else:
                # excluded_PRN.append(str(PRN) + ", ")
                excluded_PRN.append(PRN)


    ## -- Crop figure by seting  y lim to mean values pluss minus 7 std
    y_mean = np.nanmean(multipath_range1)
    y_std  = np.nanstd(multipath_range1)
    ax4.set_xlim(0,90)
    ax4.set_ylim(y_mean - 7*y_std, y_mean+7*y_std)
    if len(excluded_PRN) == 0:
        ax4.set_title('Multipath effect vs satellite elevation angle for the signal %s,\
                      %s \n Signal combination: %s-%s-%s' % (range1_Code, GNSSsystemName,range1_Code, phase1_Code, phase2_Code),fontsize=28)
    else:
        ax4.set_title('Multipath effect vs satellite elevation angle for the signal %s,\
                      %s \n Signal combination: %s-%s-%s, \n  Ekskluderte PRN: %s' %(range1_Code, GNSSsystemName,range1_Code,\
                          phase1_Code, phase2_Code, excluded_PRN),fontsize=28)
    ax4.set_xlabel('Elevation angle $[^{\circ}]$',fontsize=20,labelpad=10)
    ax4.set_ylabel('$[m]$',fontsize=20,labelpad=10)
    ax4.tick_params(axis='both', labelsize=18)
    legend = ax4.legend(loc='center right',fontsize=14,bbox_to_anchor=(1.28, 0.5), fancybox=True, shadow=True,ncol=2) # frame = legend.get_frame(); frame.set_facecolor((0.89701,0.79902,0.68137)); frame.set_edgecolor('black') #legend
    ## Set the linewidth of each legend object (then not dependent of linewith in plot)
    for legobj in legend.legendHandles:
        legobj.set_linewidth(1.5)

    ax4.grid(color='k', linestyle='-', linewidth=0.1)         
    ax4.axhline(y=0.0, color='k', linestyle='-',linewidth=1)
    # plt.show()
    
    filename  = '%s_%s_%s_MP_elevation.pdf' % (GNSSsystemName, range1_Code, range2_Code)
    filename2 = '%s_%s_%s_MP_elevation.png' % (GNSSsystemName, range1_Code, range2_Code)
    full_filename = graphDir + '/' + filename
    full_filename2 = graphDir + '/' + filename2

    fig4.savefig(graphDir + "/" +  filename)
    fig4.savefig(graphDir + "/" +  filename2, dpi=300)
    plt.close()
    
  
    ## -- Combine multipath delay on range 1 signal plots together ---
    fig5, ax5 = plt.subplots(nrows=2, ncols=1,sharex=False, squeeze=True,figsize=(16,11),dpi=160)
    fig5.subplots_adjust(left=0.07, bottom=0.1, right=0.78, top=0.91, wspace=None, hspace=0.45)
    # fig5.tight_layout()
    ## Ax 1
    for PRN in range(1,m):
        if not np.isnan(multipath_range1[:,PRN]).all():
            ax5[0].plot(t, multipath_range1[:,PRN], label='PRN%s' % (PRN),linewidth=0.7)
            
    ## -- Crop figure by seting  y lim to mean values pluss minus 7 std
    y_mean = np.nanmean(multipath_range1)
    y_std  = np.nanstd(multipath_range1)
    ax5[0].set_xlim(0,t[-1])
    ax5[0].set_ylim(y_mean - 7*y_std, y_mean+7*y_std)
    ax5[0].set_title('Multipath effect vs time for the signal %s,\
                  %s \n Signal combination: %s-%s-%s' % (range1_Code, GNSSsystemName,range1_Code, phase1_Code, phase2_Code),fontsize=22)
    ax5[0].set_xlabel('Time $[h]$',fontsize=20,labelpad=10)
    ax5[0].set_ylabel('$[m]$',fontsize=20,labelpad=10)
    ax5[0].tick_params(axis='both', labelsize=18)
    legend = ax5[0].legend(loc='center right',fontsize=12,bbox_to_anchor=(1.25, 0.5), fancybox=True, shadow=True,ncol=2) # frame = legend.get_frame(); frame.set_facecolor((0.89701,0.79902,0.68137)); frame.set_edgecolor('black') #legend
    ## Set the linewidth of each legend object (then not dependent of linewith in plot)
    for legobj in legend.legendHandles:
        legobj.set_linewidth(1.5)

    ax5[0].grid(color='k', linestyle='-', linewidth=0.08)         
    ax5[0].axhline(y=0.0, color='k', linestyle='-',linewidth=1)
    
    ## Ax 2
    excluded_PRN = []
    for PRN in range(1,m):
        epoch_missing_sat_ele = np.where(np.isnan(sat_elevation_angles[:,PRN]))
        multipath_range1[epoch_missing_sat_ele,PRN] = np.nan
        
   
    for PRN in range(1,m):
        # only plot for PRN that have any observations
        if not np.isnan(multipath_range1[:,PRN]).all(): 
            if not np.isnan(sat_elevation_angles[:, PRN]).all(): # change to all
                ax5[1].plot(sat_elevation_angles[:, PRN], multipath_range1[:,PRN],  label='PRN%s' % (PRN), linewidth= 0.7) 
            else:
                # excluded_PRN.append(str(PRN) + ", ")
                excluded_PRN.append(PRN)


    ## -- Crop figure by seting  y lim to mean values pluss minus 7 std
    y_mean = np.nanmean(multipath_range1)
    y_std  = np.nanstd(multipath_range1)
    ax5[1].set_xlim(0,90)
    ax5[1].set_ylim(y_mean - 7*y_std, y_mean+7*y_std)
    if len(excluded_PRN) == 0:
        ax5[1].set_title('Multipath effect vs satellite elevation angle for the signal %s,\
                      %s \n Signal combination: %s-%s-%s' % (range1_Code, GNSSsystemName,range1_Code, phase1_Code, phase2_Code),fontsize=22)
    else:
        ax5[1].set_title('Multipath effect vs satellite elevation angle for the signal %s,\
                      %s \n Signal combination: %s-%s-%s, \n  Ekskluderte PRN: %s' %(range1_Code, GNSSsystemName,range1_Code,\
                          phase1_Code, phase2_Code, excluded_PRN),fontsize=22)
    ax5[1].set_xlabel('Elevation angle $[^{\circ}]$',fontsize=20,labelpad=10)
    ax5[1].set_ylabel('$[m]$',fontsize=20,labelpad=10)
    ax5[1].tick_params(axis='both', labelsize=18)
    legend = ax5[1].legend(loc='center right',fontsize=12,bbox_to_anchor=(1.25, 0.5), fancybox=True, shadow=True,ncol=2) # frame = legend.get_frame(); frame.set_facecolor((0.89701,0.79902,0.68137)); frame.set_edgecolor('black') #legend
    ## Set the linewidth of each legend object (then not dependent of linewith in plot)
    for legobj in legend.legendHandles:
        legobj.set_linewidth(1.5)

    ax5[1].grid(color='k', linestyle='-', linewidth=0.08)         
    ax5[1].axhline(y=0.0, color='k', linestyle='-',linewidth=1)
    # plt.show()
    
    filename  = '%s_%s_%s_MP_combined.pdf' % (GNSSsystemName, range1_Code, range2_Code)
    filename2 = '%s_%s_%s_MP_combined.png' % (GNSSsystemName, range1_Code, range2_Code)
    full_filename = graphDir + '/' + filename
    full_filename2 = graphDir + '/' + filename2
    fig5.savefig(graphDir + "/" +  filename)
    fig5.savefig(graphDir + "/" +  filename2, dpi=300)
    plt.close()
    
    
       
    

    
    ## ----------- Combine ionospheric delay plots together ------------------

    fig6, ax6 = plt.subplots(nrows=2, ncols=1,sharex=True, squeeze=True,figsize=(16,11),dpi=160)
    fig6.subplots_adjust(left=0.07, bottom=0.1, right=0.78, top=0.91, wspace=None, hspace=0.4)

    for PRN in range(1,m):
        if not np.isnan(ion_delay_phase1[:,PRN]).all():
            ax6[0].plot(t, ion_delay_phase1[:,PRN], label='PRN%s' % (PRN),linewidth=2)
            
    ax6[0].set_title('Ionospheric delay vs time for %s \n Signal 1: %s  Signal 2: %s' % (GNSSsystemName, phase1_Code, phase2_Code),fontsize=22)
    # ax6[0].set_xlabel('Time $[h]$',fontsize=16,labelpad=10)
    ax6[0].set_ylabel('$[m]$',fontsize=16,labelpad=10)
    ax6[0].tick_params(axis='both', labelsize=16)
    legend = ax6[0].legend(loc='center right',fontsize=12,bbox_to_anchor=(1.25, 0.5), fancybox=True, shadow=True,ncol=2) # frame = legend.get_frame(); frame.set_facecolor((0.89701,0.79902,0.68137)); frame.set_edgecolor('black') #legend
    ax6[0].grid(color='k', linestyle='-', linewidth=0.1)         
    ax6[0].axhline(y=0.0, color='k', linestyle='-',linewidth=1)
        
    
    excluded_PRN = []
    for PRN in range(1,m):
        ## only plot for PRN that have any observations
        if not np.isnan(ion_delay_phase1[:,PRN]).all():
            if not np.isnan(sat_elevation_angles[:, PRN]).all():
                ax6[1].plot(t, ion_delay_phase1_zenit[:, PRN], label='PRN%s' % (str(PRN)),linewidth=2) 
            else:
                excluded_PRN.append(PRN)
                
    if len(excluded_PRN) == 0:
        ax6[1].set_title('Zenit mapped ionospheric  delay for %s \n Signal 1: %s  Signal 2: %s' % (GNSSsystemName, phase1_Code, phase2_Code),fontsize=22)
    else:
        ax6[1].set_title('Zenit mapped ionospheric  delay for %s\
                      %s \n  Signal 1: %s  Signal 2: %s \n  Ekskluderte PRN: %s' %(GNSSsystemName,range1_Code,\
                          phase1_Code, phase2_Code, excluded_PRN),fontsize=22)

    ax6[1].set_xlabel('Time $[h]$',fontsize=16,labelpad=10)
    ax6[1].set_ylabel('$[m]$',fontsize=16,labelpad=10)
    ax6[1].tick_params(axis='both', labelsize=16)
    legend = ax6[1].legend(loc='center right',fontsize=12,bbox_to_anchor=(1.25, 0.5), fancybox=True, shadow=True,ncol=2) # frame = legend.get_frame(); frame.set_facecolor((0.89701,0.79902,0.68137)); frame.set_edgecolor('black') #legend
    ax6[1].grid(color='k', linestyle='-', linewidth=0.1)         
    ax6[1].axhline(y=0.0, color='k', linestyle='-',linewidth=1)
    # plt.show()
    
    # fig6_name  = GNSSsystemName + "_" + 'ionospheric_delay_zenit_mapped' + '.png'
    # fig6.savefig(graphDir + "/" +  fig6_name, dpi=300)
    fig6_name  = GNSSsystemName + "_" + 'ionospheric_delay_combined' + '.pdf'
    # fig6.savefig(graphDir + "/" +  fig6_name)
    fig6.savefig(graphDir + "/" +  fig6_name,bbox_inches='tight')
    plt.close()

    return