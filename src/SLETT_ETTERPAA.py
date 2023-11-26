# import os
# os.chdir(r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\src")
# from GNSS_MultipathAnalysis import GNSS_MultipathAnalysis
# from readRinexObs import readRinexObs
# # rinobs = r"C:\Users\perhe\Desktop\test_GNSS\Bahadir\ANK21550.23O"
# # rinobs = r"C:\Users\perhe\OneDrive\Skrivebord\Ny mappe\P15600USA_R_20230010000_01D_15S_MO.rnx"
# rinobs = r"C:\Users\perhe\OneDrive\Skrivebord\Ny mappe\P15600USA_R_20230010000_01D_60S_MO.obs"


# # rinobs = r"C:\Users\perhe\Desktop\test_GNSS\Bahadir\ANK21550_RINEX3.obs"
# # sp3 = r"C:\Users\perhe\Desktop\test_GNSS\Bahadir\IGS0OPSFIN_20231550000_01D_15M_ORB.SP3"
# sp3 = r"C:\Users\perhe\OneDrive\Skrivebord\Ny mappe\P15600USA_R_20230010000_01D_MN.rnx"
# # analysisResults = GNSS_MultipathAnalysis(rinobs,sp3NavFilename_1=sp3,outputDir=r"C:\Users\perhe\Desktop\test_GNSS\Bahadir\Results",desiredGNSSsystems=['G'],include_SNR=True)
# desiredGNSSsystems=['G']
# analysisResults = GNSS_MultipathAnalysis(rinobs,broadcastNav1=sp3,outputDir=r"C:\Users\perhe\OneDrive\Skrivebord\Ny mappe",desiredGNSSsystems=desiredGNSSsystems,include_SNR=True)

# # GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
# #       obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
# #       rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success = readRinexObs(rinobs)







#%% sy428 orginalfiler

# import os
# os.chdir(r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\src")
# from GNSS_MultipathAnalysis import GNSS_MultipathAnalysis
# from readRinexObs import readRinexObs

# # outputDir = r"C:\Users\perhe\Desktop\test_GNSS\sy428\Riktig res"
# outputDir = r"C:\Users\perhe\Desktop\test_GNSS\sy428\test"
# rinobs = r"C:\Users\perhe\Desktop\test_GNSS\sy428\Orignalfiler\abmf2440.22o"
# # rinobs = r"C:\Users\perhe\Desktop\test_GNSS\sy428\Orignalfiler\abmf2440_org_fjernet_fcn.22o"
# # rinobs = r"C:\Users\perhe\Desktop\test_GNSS\sy428\Orignalfiler\abmf2440_crop.obs"
# sp3 = r"C:\Users\perhe\Desktop\test_GNSS\sy428\Orignalfiler\BRDC00IGS_R_20222440000_01D_MN.rnx"

# desiredGNSSsystems=['G',"R","C"]
# analysisResults = GNSS_MultipathAnalysis(rinobs,broadcastNav1=sp3,
#                                           outputDir=outputDir,
#                                           include_SNR=True)


#%% sy428 orginalfiler
import os
os.chdir(r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\src")
from gnssmultipath import GNSS_MultipathAnalysis
from gnssmultipath import readRinexObs

# outputDir = r"C:\Users\perhe\Desktop\test_GNSS\sy428\Riktig res"
# outputDir = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\sy428\test_org_2\ajac2440"
# outputDir2 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\sy428\test_org_2\ajac2440"
# outputDir3 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\sy428\test_org_2\alic2440"
# outputDir4 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\sy428\test_org_2\anmg2440"
outputDir_test = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\sy428\test_org_2\test"

# rinobs1 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\sy428\Orginalfiler2\aggo2440.22o"
# rinobs2 = r"C:/Users/perhe/OneDrive/Documents/Python_skript/test_GNSS/sy428/Orginalfiler2/ajac2440.22o"
rinobs3 = r"C:/Users/perhe/OneDrive/Documents/Python_skript/test_GNSS/sy428/Orginalfiler2/alic2440.22o"
# rinobs4 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\sy428\Orginalfiler2\anmg2440.22o"
# broad_nav = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\sy428\Orginalfiler2\BRDC00IGS_R_20222440000_01D_MN.rnx"
sp3 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\sy428\Orginalfiler2\gfz22254.sp3"
# desiredGNSSsystems=['G',"R","C"]
desiredGNSSsystems=['G']
analysisResults = GNSS_MultipathAnalysis(rinobs3,
                                          # broadcastNav1=broad_nav,
                                          sp3NavFilename_1=sp3,
                                          outputDir=outputDir_test,
                                           desiredGNSSsystems = desiredGNSSsystems,
                                          include_SNR=True)
#%%

import os
os.chdir(r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\src")
from GNSS_MultipathAnalysis import GNSS_MultipathAnalysis
from readRinexObs import readRinexObs

outputDir_test = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\sy428\test_org_2\test"
rin_NMBUS = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\NMBUS_SAMSUNG_S20.20o"
sp3_NMBUS = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\SP3\NMBUS_2020 10 30.SP3"
analysisResults = GNSS_MultipathAnalysis(rin_NMBUS,
                                          # broadcastNav1=broad_nav,
                                          sp3NavFilename_1=sp3_NMBUS,
                                          outputDir=outputDir_test,
                                           # desiredGNSSsystems = desiredGNSSsystems,
                                          plotEstimates=True)


#%%
import os
os.chdir(r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\src")
from GNSS_MultipathAnalysis import GNSS_MultipathAnalysis
from readRinexObs import readRinexObs

outputDir_test = r"C:\Users\perhe\Desktop\_new"
rin_NMBUS = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.rnx"
broad_NMBUS = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\BRDC00IGS_R_20220010000_01D_MN.rnx"
# sp3_NMBUS = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\SP3\Testfile_20220101.eph"

desiredGNSSsystems = ["G"]
analysisResults = GNSS_MultipathAnalysis(rin_NMBUS,
                                          broadcastNav1=broad_NMBUS,
                                         # sp3NavFilename_1=sp3_NMBUS,
                                         outputDir=outputDir_test,
                                          desiredGNSSsystems=desiredGNSSsystems,
                                         include_SNR=True)


#%%
import os
os.chdir(r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\src")
from GNSS_MultipathAnalysis import GNSS_MultipathAnalysis
from readRinexObs import readRinexObs

outputDir_erly = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\Erly Caldas de Lima\Results"
rin_erly = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\Erly Caldas de Lima\P15600USA_R_20230010000_01D_15S_MO.rnx"
broad_erly = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\Erly Caldas de Lima\P15600USA_R_20230010000_01D_MN.rnx"

analysisResults = GNSS_MultipathAnalysis(rin_erly,
                                         broadcastNav1=broad_erly,
                                         outputDir=outputDir_erly,
                                         desiredGNSSsystems=["R"],
                                         include_SNR=True)

#%% io bakadir
import os
os.chdir(r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\src")
from GNSS_MultipathAnalysis import GNSS_MultipathAnalysis
from readRinexObs import readRinexObs

outputDir = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\IO\Resultat"
outputDir = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\test"
rin = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\IO\YLOV3350_croped.22O"
sp3 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\IO\WUM0MGXULA_20223350000_01D_05M_ORB.SP3"

analysisResults = GNSS_MultipathAnalysis(rin,
                                         sp3NavFilename_1=sp3,
                                         outputDir=outputDir,
                                          # desiredGNSSsystems=["G"],
                                         include_SNR=True)


#%%

import os
from GNSS_MultipathAnalysis import GNSS_MultipathAnalysis

base_path = os.getcwd()
parent_dir = os.path.abspath(os.path.join(base_path, os.pardir))

## Path to TestData
outputDir = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\IO\\Resultat_more_user_arg"

rin = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\IO\YLOV3350.22O"
sp3 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\IO\WUM0MGXULA_20223350000_01D_05M_ORB.SP3"

# analysisResults = GNSS_MultipathAnalysis(rin,
#                                          sp3NavFilename_1=sp3,
#                                           desiredGNSSsystems=["R"],
#                                           outputDir=outputDir,
#                                           include_SNR=True,
#                                          use_LaTex = False)

## Parameters
GNSSsystems                 = ["R"] # run analysis in GLONASS only
phaseCodeLimit              = 6.667 #
ionLimit                    = 0.0667 #
cutoff_elevation_angle      = 20 # 10 degree elevation cutoff #
outputDir                   = outputDir #
plotEstimates               = True
plot_polarplot              = True
include_SNR                 = True
tLim_R                      = None
tLim_GEC                    = None
includeResultSummary        = True
includeCompactSummary       = True
includeObservationOverview  = True
includeLLIOverview          = True
use_LaTex                   = False


# sp3NavFilename_1=None,

## Rinex observation file
rinObsFilename1 = rin
## SP3 files
sp3NavFilename_1_opec = sp3


analysisResults = GNSS_MultipathAnalysis(rinObsFilename1,
                        # desiredGNSSsystems=GNSSsystems,
                        sp3NavFilename_1 = sp3NavFilename_1_opec,
                        phaseCodeLimit = phaseCodeLimit,
                        ionLimit = ionLimit,
                        cutoff_elevation_angle = cutoff_elevation_angle,
                        outputDir = outputDir,
                        plotEstimates = plotEstimates,
                        plot_polarplot=plot_polarplot,
                        include_SNR = include_SNR,
#                        tLim_R = tLim_R,
#                        tLim_GEC = tLim_GEC,
                        includeResultSummary = includeResultSummary,
                        includeCompactSummary = includeCompactSummary,
                        includeObservationOverview = includeObservationOverview,
                        includeLLIOverview = includeLLIOverview,
                        # use_LaTex = use_LaTex
                        )


#%% black749 second folder
import os
os.chdir(r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\src")
from gnssmultipath import GNSS_MultipathAnalysis
from gnssmultipath import readRinexObs

rin1 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\apsl2440.22o"
# rin11 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\apsl2440_croped.22o"
# rin111 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\apsl2440_croped.obs"
# rin1111 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\apsl2440_croped2.22o"
# rin2 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\beri2440.22o"
rin3 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\krs12440.22o"
rin4 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\wuh22440.22o"
# rin5 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\wuhn2440.22o"

brod1 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\BRDC00IGS_R_20222440000_01D_MN.rnx"
sp3 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\COD0MGXFIN_20222440000_01D_05M_ORB.SP3"
# outputDir1 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\Results_apsl2440"
# outputDir11 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\Results_apsl2440_croped2"
# outputDir111 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\Results_apsl2440 2"

# outputDir2 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\Results_beri2440"
# outputDir3 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\Results_krs12440"
outputDir4 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\Results_wuh22440"
# outputDir5 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\Results_wuhn2440"

analysisResults = GNSS_MultipathAnalysis(rin3,
                                           # broadcastNav1=brod1,
                                           sp3NavFilename_1=sp3,
                                          outputDir=r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\test2",
                                         # outputDir=outputDir4,
                                           # desiredGNSSsystems=["G"],
                                           # plotEstimates=True,
                                          use_LaTex=True,
                                         include_SNR=True)



#%%
import os
os.chdir(r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\src")
from readRinexObs import readRinexObs
rin1 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\apsl2440.22o"
rin4 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\wuh22440.22o"
rin5 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\wuhn2440.22o"
GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
      obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
      rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success = \
      readRinexObs(rin5)

#%%
from readRinexNav import read_rinex3_nav
import os
os.chdir(r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\src")
brod1 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\BRDC00IGS_R_20222440000_01D_MN.rnx"
brod1 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\NMBUS_SAMSUNG_S20.20n"
# brod1 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\BCEmerge_30_10_2020.20p"
data = read_rinex3_nav(brod1)


#%%
import numpy as np
data_array,header, n_eph = data
# Get the PRN and FCN columns
prn_column = data_array[:, 0]

# Find GLONASS data
glo_rows = np.char.startswith(prn_column, 'R')
glo_data = data_array[glo_rows]
unique_prns_name, unique_glo_prns_idx = np.unique(glo_data[:, 0], return_index=True)
PRN_data = glo_data[unique_glo_prns_idx]

prn_fcn_dict2 = {str(prn): int(fcn) for prn, fcn in zip(PRN_data[:,0], PRN_data[:,17].astype(float).astype(int))}




#%%
import numpy as np

data_array,header, n_eph = data
# Get the PRN and FCN columns
prn_column = data_array[:, 0]
fcn_column = data_array[:, 17]

# Create an empty dictionary to store FCN numbers for each PRN
fcn_dict = {}

# Use NumPy's unique function to get unique PRN numbers
glo_rows = np.char.startswith(prn_column, 'R')
unique_prn = prn_column[glo_rows]
# unique_prn = np.unique(prn_column.startswith('R'))

glo_data = data_array[glo_rows]

# Iterate through unique PRN numbers and extract corresponding FCN numbers
for prn in unique_prn:
    fcn_numbers = fcn_column[prn_column == prn]
    fcn_dict[prn] = fcn_numbers
    
    
#     def export_results_as_comp_pickle(dict_data):
#     data_out = pickle.dumps(dict_data)
#     size_orig = len(data_out)
#     datra

# def put_object_redis(object, redis_key): 
#     data_out = pickle.dumps(object) 
#     size_orig = len(data_out) 
#     data_out_compress = zstd.compress(data_out, level=22) 
#     print("Object compression ratio", size_orig/len(data_out_compress))
#     redis_client.set(redis_key, data_out_compress)


#%%
from PickleHandler import PickleHandler
from GNSS_MultipathAnalysis import GNSS_MultipathAnalysis

old = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\tests\analysisResults_OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.pkl"
old_res = PickleHandler.read_pickle(old)

ofile = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.rnx"
nfile = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\BRDC00IGS_R_20220010000_01D_MN.rnx"
new_res = GNSS_MultipathAnalysis(ofile,broadcastNav1=nfile,nav_data_rate=120)



#%%
from GNSS_MultipathAnalysis import GNSS_MultipathAnalysis
rinObsFilename1 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx"
sp3 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\SP3\Testfile_20220101.eph"
output = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\test2"
analysisResults = GNSS_MultipathAnalysis(rinObsFilename1, 
                                         sp3NavFilename_1=sp3,
                                         outputDir=output)


#%%
from GNSS_MultipathAnalysis import GNSS_MultipathAnalysis
rinObsFilename1 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx"
broad = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\BRDC00IGS_R_20220010000_01D_MN.rnx"
output = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\test3"
analysisResults = GNSS_MultipathAnalysis(rinObsFilename1, 
                                         broadcastNav1=broad,
                                         outputDir=output,
                                         nav_data_rate=120)


#%%
from gnssmultipath import GNSS_MultipathAnalysis

rinObsFilename1 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx"
broad = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\BRDC00IGS_R_20220010000_01D_MN.rnx"
sp3 = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\SP3\Testfile_20220101.eph"
output = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\test5"
analysisResults = GNSS_MultipathAnalysis(rinObsFilename1, 
                                           # broadcastNav1=broad,
                                           sp3NavFilename_1=sp3,
                                          # desiredGNSSsystems = ["R"],
                                         outputDir=output,
                                         nav_data_rate=120)


#%%
from gnssmultipath import PickleHandler

path_to_picklefile_codespace = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\test3\analysisResults_codespace.pkl"
path_to_picklefile = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\black749\test3\analysisResults.pkl"

res = PickleHandler.read_zstd_pickle(path_to_picklefile)
res_codespace = PickleHandler.read_zstd_pickle(path_to_picklefile_codespace)



#%%
from gnssmultipath import GNSS_MultipathAnalysis

outputDir_erly = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\Erly Caldas de Lima\Results3"
outputDir_erly = r"C:\Users\perhe\Desktop\TEST VEC\test1"
rin_erly = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\Erly Caldas de Lima\P15600USA_R_20230010000_01D_15S_MO.rnx"
broad_erly = r"C:\Users\perhe\OneDrive\Documents\Python_skript\test_GNSS\Erly Caldas de Lima\P15600USA_R_20230010000_01D_MN.rnx"

analysisResults = GNSS_MultipathAnalysis(rin_erly,
                                         broadcastNav1=broad_erly,
                                         outputDir=outputDir_erly,
                                          desiredGNSSsystems=["R"],
                                         include_SNR=True)



#%%
import os
from gnssmultipath import GNSS_MultipathAnalysis

outputDir_test = r"C:\Users\perhe\Desktop\TEST VEC\test2"
# rinObs = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\NMBUS_SAMSUNG_S20.20o"
rinObs = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx"
sp3_NMBUS = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\SP3\NMBUS_2020 10 30.SP3"
# rinnav = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\BCEmerge_30_10_2020.20p"
rinnav = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\BRDC00IGS_R_20220010000_01D_MN.rnx"
analysisResults = GNSS_MultipathAnalysis(rinObs,
                                           broadcastNav1=rinnav,
                                          # sp3NavFilename_1=sp3_NMBUS,
                                          outputDir=outputDir_test,
                                           # desiredGNSSsystems = desiredGNSSsystems,
                                          plotEstimates=True,
                                          include_SNR=True,
                                          nav_data_rate=120)

#%%
import os
from gnssmultipath import GNSS_MultipathAnalysis

outputDir_test = r"C:\Users\perhe\Desktop\TEST VEC\test3"
rinObs = r"C:\Users\perhe\Desktop\TEST VEC\test3\OPEC00NOR_S_20180010000_01D_01S_MO.rnx"
# rinnav = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\BRDC00IGS_R_20220010000_01D_MN.rnx"
sp3 = r"C:\Users\perhe\Desktop\TEST VEC\test3\COM_20180010000_.eph"
analysisResults = GNSS_MultipathAnalysis(rinObs,
                                           # broadcastNav1=rinnav,
                                           sp3NavFilename_1=sp3,
                                          outputDir=outputDir_test,
                                          plotEstimates=True,
                                          include_SNR=True)



#%%
import os
from gnssmultipath import GNSS_MultipathAnalysis
nav1 = r"C:\Users\perhe\Desktop\TEST VEC\test3\OPEC00NOR_S_20180010000_01D_EN.rnx"
nav2 = r"C:\Users\perhe\Desktop\TEST VEC\test3\OPEC00NOR_S_20180010000_01D_GN.rnx"
nav3 = r"C:\Users\perhe\Desktop\TEST VEC\test3\OPEC00NOR_S_20180010000_01D_RN.rnx"
nav4 = r"C:\Users\perhe\Desktop\TEST VEC\test3\OPEC00NOR_S_20180010000_01D_CN.rnx"
rinnav = [nav1,nav2,nav3,nav4]
outputDir_test = r"C:\Users\perhe\Desktop\TEST VEC\test3_2"
rinObs = r"C:\Users\perhe\Desktop\TEST VEC\test3\OPEC00NOR_S_20180010000_01D_01S_MO.rnx"
sp3 = r"C:\Users\perhe\Desktop\TEST VEC\test3\COM_20180010000_.eph"
analysisResults = GNSS_MultipathAnalysis(rinObs,
                                            broadcastNav1=rinnav,
                                           # sp3NavFilename_1=sp3,
                                          outputDir=outputDir_test,
                                          plotEstimates=True,
                                          include_SNR=True,
                                          nav_data_rate=120)

#%%
import os
from gnssmultipath import GNSS_MultipathAnalysis

outputDir_test = r"C:\Users\perhe\Desktop\TEST VEC\test5_nytest"
rinObs = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.rnx"
rinnav = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\BRDC00IGS_R_20220010000_01D_MN.rnx"
analysisResults = GNSS_MultipathAnalysis(rinObs,
                                            broadcastNav1=rinnav,
                                           # sp3NavFilename_1=sp3,
                                          outputDir=outputDir_test,
                                          plotEstimates=True,
                                          include_SNR=True,
                                          nav_data_rate=120)

#%%
import os
from gnssmultipath import GNSS_MultipathAnalysis

outputDir_test = r"C:\Users\perhe\Desktop\TEST VEC\test6_nytest"
rinObs = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\NMBUS_SAMSUNG_S20.20o"
sp3_NMBUS = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\SP3\NMBUS_2020 10 30.SP3"
analysisResults = GNSS_MultipathAnalysis(rinObs,
                                          sp3NavFilename_1=sp3_NMBUS,
                                          outputDir=outputDir_test,
                                          plotEstimates=True,
                                          save_results_as_compressed_pickle=True)

#%%
from gnssmultipath import PickleHandler
import numpy as np

path_to_picklefile = r"C:\Users\perhe\Desktop\TEST VEC\test5_nytest\analysisResults_OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.pkl"
res = PickleHandler.read_zstd_pickle(path_to_picklefile)

# mp_c1c = res["GPS"]["Band_1"]["C1C"]["multipath_range1"]
# el_c1c = res["Sat_position"]["G"]["elevation"]
# az_c1c = res["Sat_position"]["G"]["azimuth"]
# mask = np.isnan(mp_c1c)
# el_c1c[mask] = np.nan
# az_c1c[mask] = np.nan

# mp_c1c = res["GLONASS"]["Band_1"]["C1C"]["multipath_range1"]
# el_c1c = res["Sat_position"]["R"]["elevation"]
# az_c1c = res["Sat_position"]["R"]["azimuth"]
# mask = np.isnan(mp_c1c)
# el_c1c[mask] = np.nan
# az_c1c[mask] = np.nan

mp_c2x = res["BeiDou"]["Band_2"]["C2X"]["multipath_range1"]
el_c2x = res["Sat_position"]["C"]["elevation"]
az_c2x = res["Sat_position"]["C"]["azimuth"]
mask = np.isnan(mp_c2x)
el_c2x[mask] = np.nan
az_c2x[mask] = np.nan


