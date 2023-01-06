# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 07:02:56 2022

@author: perhe
"""

import pandas as pd
import numpy as np
columns = list(range(1,37))
mp_python =pd.read_csv('multipath_range1.csv')
mp_matlab =pd.read_csv('multipath_range1_matlab.csv')

mp_python = mp_python.iloc[:,1::]
# diff = pd.concat([mp_matlab,mp_python]).drop_duplicates(keep=False)
# diff = set(mp_matlab.columns).symmetric_difference(mp_python.columns)
# for ep in range(0,len(mp_matlab)):
#     diff = mp_matlab.iloc[ep,:] - mp_python.iloc[ep,:]
# diff = np.array([len(mp_matlab),36])
# for i in range(0,37):
#     diff[:,i] = (mp_matlab.iloc[:,i] - mp_python.iloc[:,i])
# PRN = 5
# diff = np.array([])
# for PRN in range(0,36):
#   diff=  np.append(diff,np.array(mp_matlab.iloc[:,PRN] - mp_python.iloc[:,PRN]))
    
# avvik = np.where(diff > 0.2)[0]
PRN = 5
diff = np.zeros([len(mp_matlab),37])
for ep in range(0,len(mp_matlab)):
    for PRN in range(0,36):
        diff[ep,PRN]  =  mp_matlab.iloc[ep,PRN] - mp_python.iloc[ep,PRN]
    
avvik = np.where(abs(diff) > 0.2)[0]


#%%
import pandas as pd
import numpy as np
mp_python =pd.read_csv('ion_delay_phase1_pyth.csv')
mp_matlab =pd.read_csv('ion_delay_phase1_mat.csv')

mp_python = mp_python.iloc[:,1::]
# diff = pd.concat([mp_matlab,mp_python]).drop_duplicates(keep=False)
# diff = set(mp_matlab.columns).symmetric_difference(mp_python.columns)
# for ep in range(0,len(mp_matlab)):
#     diff = mp_matlab.iloc[ep,:] - mp_python.iloc[ep,:]
# diff = np.array([len(mp_matlab),36])
# for i in range(0,37):
#     diff[:,i] = (mp_matlab.iloc[:,i] - mp_python.iloc[:,i])
# PRN = 5
# diff = np.array([])
# for PRN in range(0,36):
#   diff=  np.append(diff,np.array(mp_matlab.iloc[:,PRN] - mp_python.iloc[:,PRN]))
    
# avvik = np.where(diff > 0.2)[0]
sat = []
diff = np.zeros([len(mp_matlab),37])
for ep in range(0,len(mp_matlab)):
    for PRN in range(0,36):
        diff[ep,PRN]  =  mp_matlab.iloc[ep,PRN] - mp_python.iloc[ep,PRN]
    
avvik = np.where(abs(diff) > 0.000000001)[0]




