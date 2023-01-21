import os, sys, pickle,numpy as np
from readRinexObs304 import readRinexObs304
from GNSS_MultipathAnalysis import GNSS_MultipathAnalysis

abs_path = os.path.abspath("GNSS_MultipathAnalysis.py")
base_path = os.path.normpath(os.getcwd() + os.sep + os.pardir)
## Path to TestData
relpath_to_testdata = 'TestData'
full_path_testdata = os.path.join(base_path, relpath_to_testdata) 

#Path to outputdir 
relpath_to_outputdir = 'Results'
full_path_ouputdir = os.path.join(base_path, relpath_to_outputdir) 
## -----  Defining input data --------

## Rinex observation file
# rinObsFilename = full_path_testdata  + '/ObservationFiles/' + 'OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx'
# rinObsFilename = full_path_testdata  + '/ObservationFiles/' + 'OPEC00NOR_S_20220010000_01D_30S_CONV.rnx'
# rinObsFilename = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Multipath_analysis/NMBUS_SAMSUNG_S20.20o'
rinObsFilename = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles/NMBUS_SAMSUNG_S20.20o'

## SP3 files
# sp3NavFilename_1 = full_path_testdata  + '/SP3/' + 'samsung2.SP3'
sp3NavFilename_1 = full_path_testdata  + '/SP3/' + 'test1.eph'
sp3NavFilename_2 = full_path_testdata  + '/SP3/' + 'test2.SP3'
sp3NavFilename_3 = full_path_testdata  + '/SP3/' + 'test3.SP3'


broadcastNav1 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_GN.rnx'
broadcastNav2 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_RN.rnx'
broadcastNav3 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_EN.rnx'
broadcastNav4 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_CN.rnx'

# broadcastNav1 =  full_path_testdata  + '/NavigationFiles/' + 'BRDC00IGS_R_20220010000_01D_MN.rnx'

## Broadcast ephemerides
# broadcastNav1 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_GN.rnx'
# broadcastNav2 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_RN.rnx'
# broadcastNav3 =  full_path_testdata  + '/NavigationFiles/' + 'ENDRET_BCEmerge_30_10_2020.20p'
broadcastNav3 =  full_path_testdata  + '/NavigationFiles/' + 'BCEmerge_30_10_2020.20p'


## -- Simple example for running analysis (no userdefined settings)
# analysisResults = GNSS_MultipathAnalysis(rinObsFilename, 
#                                          broadcastNav1=broadcastNav1,
#                                          broadcastNav2=broadcastNav3
#                                          )


# analysisResults = GNSS_MultipathAnalysis(rinObsFilename,
#                                           # desiredGNSSsystems=['G'],
#                                           outputDir = full_path_ouputdir,
#                                           broadcastNav1=broadcastNav1
#                                           )



analysisResults = GNSS_MultipathAnalysis(rinObsFilename,
                                            # desiredGNSSsystems=['G'],
                                          outputDir = full_path_ouputdir,
                                          broadcastNav1=broadcastNav3
                                          )

# analysisResults = GNSS_MultipathAnalysis(rinObsFilename,
#                                           sp3NavFilename_1=sp3NavFilename_1
#                                           )

# analysisResults = GNSS_MultipathAnalysis(rinObsFilename, 
#                                          broadcastNav1=broadcastNav1,
#                                          broadcastNav2=broadcastNav2,
#                                          broadcastNav3=broadcastNav3,
#                                          broadcastNav4=broadcastNav4
#                                          )


# from readRinexNav import *
# data, header, n_eph = read_rinex3_nav(broadcastNav1)
#%% -- Advanced example (more user defined settings)

## Parameters
phaseCodeLimit              = 0
ionLimit                    = 0
cutoff_elevation_angle      = 0
outputDir                   = full_path_ouputdir
plotEstimates               = True
includeResultSummary        = True
includeCompactSummary       = True
includeObservationOverview  = True
includeLLIOverview          = True


GNSSsystems = ["R"] # run analysis in GLONASS only
analysisResults = GNSS_MultipathAnalysis(rinObsFilename,
                       desiredGNSSsystems=GNSSsystems,
                       sp3NavFilename_1 = sp3NavFilename_1,
                       sp3NavFilename_2 = sp3NavFilename_2,
                       sp3NavFilename_3 = sp3NavFilename_3,
                       phaseCodeLimit = phaseCodeLimit,
                       ionLimit = ionLimit,
                       cutoff_elevation_angle = cutoff_elevation_angle,
                       outputDir = outputDir, 
                       plotEstimates = plotEstimates , 
                       includeResultSummary = includeResultSummary , 
                       includeCompactSummary = includeCompactSummary , 
                       includeObservationOverview = includeObservationOverview ,
                       includeLLIOverview = includeLLIOverview
                       )
 
  
## --  Read in the result file from a analysis
# result_path = os.path.join(full_path_ouputdir, "analysisResults.pkl")
file_to_read = open(r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo/analysisResults.pkl', "rb")
loaded_dictionary = pickle.load(file_to_read)


#%% Read multiple nav
import os, sys, pickle,numpy as np
abs_path = os.path.abspath("GNSS_Receiver_QC_2020.py")
base_path = os.path.normpath(os.getcwd() + os.sep + os.pardir)
## Path to TestData
relpath_to_testdata = 'TestData'
full_path_testdata = os.path.join(base_path, relpath_to_testdata) 

from readRinexNav import *

broadcastNav1 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_GN.rnx'
broadcastNav2 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_RN.rnx'
broadcastNav3 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_EN.rnx'
broadcastNav4 =  full_path_testdata  + '/NavigationFiles/' + 'OPEC00NOR_S_20220010000_01D_CN.rnx'


data, header, n_eph = read_rinex3_nav(broadcastNav2)
nav_list = [broadcastNav1,broadcastNav2,None,None]
nav_list = [i for i in nav_list if i is not None]
sat_pos2 = {}
for file in nav_list:
    data, header, n_eph =read_rinex3_nav(file)
    curr_sys = data[0][0][0]
    sat_pos2[curr_sys] = data 



#%% Plotting RMS value in bar plot (not here, just for now)
import os,pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib import patheffects
from matplotlib.ticker import MaxNLocator
from matplotlib import rc
# from add_value_labels_to_bar import add_value_labels_to_bar
from matplotlib.ticker import FormatStrFormatter



# fig = plt.figure(figsize = (10, 5))

# Tot = len(current_systems)
# Cols = 2
# if Tot == 1:
#     Cols = 1
# # Compute Rows required
# Rows = Tot // Cols 
# #     EDIT for correct number of rows:If one additional row is necessary -> add one:
# if Tot % Cols != 0:
#     Rows += 1

# current_systems = analysisResults['GNSSsystems']
# # current_systems = ['GPS','GLONASS','BeiDou']
# if len(current_systems) == 4:
#     # fig, ax = plt.subplots()
#     max_MP = []
#     fig, ax = plt.subplots(nrows=2, ncols=2,sharex=False,figsize=(18,12),dpi = 100)
#     fig.subplots_adjust(left=0.082, bottom=0.08, right=0.887, top=0.93, wspace=None, hspace=0.2)
#     for idx,sys in enumerate(current_systems):
#         if idx == 0:
#             row_idx = 0
#             col_idx = 0
#         elif idx == 1:
#             row_idx = 0
#             col_idx = 1
#         elif idx == 2:
#             row_idx = 1
#             col_idx = 0
#         elif idx == 3:
#             row_idx = 1
#             col_idx = 1
#         data_elw_rms = []
#         data_rms = []
#         data_codes = []
#         bands_curr_sys = analysisResults[sys]['Bands']
#         for band in bands_curr_sys:
#             codes_curr_sys = analysisResults[sys][band]['Codes']
#             codes_curr_sys = [ele for ele in codes_curr_sys if ele != []] # removing empty list if exist
#             for code in codes_curr_sys:
#                 elweight_rms_MP = analysisResults[sys][band][code]['elevation_weighted_average_rms_multipath_range1']
#                 rms_MP = analysisResults[sys][band][code]['rms_multipath_range1_averaged']
#                 data_rms.append(rms_MP)
#                 data_elw_rms.append(elweight_rms_MP)
#                 data_codes.append(code)
#         # creating the bar plot
#         max_MP.append(max(data_rms + data_elw_rms))
#         width = 0.35  # the width of the bars
#         x = np.arange(len(data_codes))  # the label locations
#         rects1 = ax[row_idx,col_idx].bar(x - width/2, data_rms, width, label='RMS')
#         rects2 = ax[row_idx,col_idx].bar(x + width/2, data_elw_rms, width, label='RMS (weighted)')
#         # Add some text for labels, title and custom x-axis tick labels, etc.
#         ax[row_idx,col_idx].set_ylabel('RMS [m]',fontsize=18,labelpad=20)
#         ax[row_idx,col_idx].set_title('RMS values for the multipath effect (%s)' %(sys),fontsize=24)
#         ax[row_idx,col_idx].set_xticks(x)
#         ax[row_idx,col_idx].set_xticklabels(data_codes)
#         # ax[row_idx,col_idx].set_ylim([0, max(max_MP)+0.08])
#         # plt.locator_params(axis='y', nbins=12)
#         ax[row_idx,col_idx].locator_params(tight=True, nbins=12)
#         ax[row_idx,col_idx].legend(fontsize=15,fancybox=True, shadow=True)
#         ax[row_idx,col_idx].tick_params(axis='both', labelsize= 15)
#         ax[row_idx,col_idx].grid(color='grey', linestyle='-', linewidth=0.3)
#     plt.setp(ax,ylim=(0,max(max_MP)+0.08))
#     plt.show()
#     # fig.savefig('Barplot_RMS.png', dpi=300, orientation='landscape')
#     fig.savefig('Barplot_RMS.pdf', orientation='landscape',bbox_inches='tight')
# else:
#     ## first find max value of RMS
#     max_MP = []
#     for idx,sys in enumerate(current_systems):
#         data_elw_rms = []
#         data_rms = []
#         bands_curr_sys = analysisResults[sys]['Bands']
#         for band in bands_curr_sys:
#             codes_curr_sys = [ele for ele in analysisResults[sys][band]['Codes'] if ele != []] # removing empty list if exist
#             for code in codes_curr_sys:
#                 elweight_rms_MP = analysisResults[sys][band][code]['elevation_weighted_average_rms_multipath_range1']
#                 rms_MP = analysisResults[sys][band][code]['rms_multipath_range1_averaged']
#                 data_rms.append(rms_MP)
#         max_MP.append(max(data_rms + data_elw_rms))
#     ## then do the plotting
#     for idx,sys in enumerate(current_systems):
#         data_elw_rms = []
#         data_rms = []
#         data_codes = []
#         bands_curr_sys = analysisResults[sys]['Bands']
#         fig, ax = plt.subplots(nrows=1, ncols=1,sharex=False,figsize=(18,12),dpi = 150)
#         fig.subplots_adjust(left=0.082, bottom=0.08, right=0.887, top=0.93, wspace=None, hspace=0.2)
#         for band in bands_curr_sys:
#             codes_curr_sys = analysisResults[sys][band]['Codes']
#             codes_curr_sys = [ele for ele in codes_curr_sys if ele != []] # removing empty list if exist
#             for code in codes_curr_sys:
#                 elweight_rms_MP = analysisResults[sys][band][code]['elevation_weighted_average_rms_multipath_range1']
#                 rms_MP = analysisResults[sys][band][code]['rms_multipath_range1_averaged']
#                 data_rms.append(rms_MP)
#                 data_elw_rms.append(elweight_rms_MP)
#                 data_codes.append(code)
#         # creating the bar plot
#         width = 0.35  # the width of the bars
#         x = np.arange(len(data_codes))  # the label locations
#         rects1 = ax.bar(x - width/2, data_rms, width, label='RMS')
#         rects2 = ax.bar(x + width/2, data_elw_rms, width, label='RMS (weighted)')
#         # Add some text for labels, title and custom x-axis tick labels, etc.
#         ax.set_ylabel('RMS [m]',fontsize=18,labelpad=20)
#         ax.set_title('RMS values for the multipath effect (%s)' %(sys),fontsize=24)
#         ax.set_xticks(x); ax.set_xticklabels(data_codes)
#         ax.locator_params(tight=True, nbins=12)
#         ax.legend(fontsize=15,fancybox=True, shadow=True)
#         ax.tick_params(axis='both', labelsize= 15)
#         ax.grid(color='grey', linestyle='-', linewidth=0.3)
#         plt.setp(ax,ylim=(0,max(max_MP)+0.08))
#         plt.show()
#         # fig.savefig('Barplot_RMS.png', dpi=300, orientation='landscape')
#         fileName ='Barplot_RMS_%s.pdf' % (sys)
#         fig.savefig(fileName, orientation='landscape',bbox_inches='tight')
 

#%%

# phaseCodeLimit              = 0;
# ionLimit                    = 0;
# cutoff_elevation_angle      = 0; 
# # outputDIR                   = r'C:\Users\perhe\Desktop\Data\Test2/'
# outputDir                   = r'C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS\Multipath_analysis/Test2'
# plotEstimates               = 1;
# saveWorkspace               = 1;
# includeResultSummary        = 1;
# includeCompactSummary       = 1;
# includeObservationOverview  = 1;
# includeLLIOverview          = 1;
# readLLI=1
# readSS = 1
# includeAllGNSSsystems   = 0
# includeAllObsCodes      = 0
# desiredGNSSsystems = ["G"]
# desiredObsCodes = ["C", "L"];               # code and phase observations
# desiredObsBands = list(np.arange(1,10))     # all carrier bands. Tot 9, but arange stops at 8 -> 10

# ## --- Read RINEX 3.0x observation file
# [GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
#     obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
#     rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success] = \
#     readRinexObs304(rinObsFilename, readSS, readLLI, includeAllGNSSsystems,includeAllObsCodes, desiredGNSSsystems,\
#     desiredObsCodes, desiredObsBands)
