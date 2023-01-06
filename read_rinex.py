# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 23:26:29 2022

@author: perhe
"""

import georinex as gr
import subprocess
obs = gr.load(r'C:\Users\perhe\OneDrive\Skrivebord/opec0020_3.04_kort.10o',use=['G','R','E'])


rin_data = obs.to_dict()
observ = rin_data['data_vars']
obs_codes = observ.keys()

for ocode in obs_codes:
    print(observ[ocode]['data'])

