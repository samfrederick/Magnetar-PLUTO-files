#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 07:50:16 2020

@author: samfrederick
"""
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

df_1 = pd.read_csv('201206_VolAvg_Total_Bfield.csv', index_col='t')
df_2 = pd.read_csv('201206_VolAvg_Poloidal_Bfield.csv', index_col='t')
df_3 = pd.read_csv('201206_VolAvg_Toroidal_Bfield.csv', index_col='t')

fig, ax = plt.subplots(1, 1, figsize=(8, 6))
ax.plot(df_1.index, df_1.BTotal_v_avg, label='Total')
ax.plot(df_2.index, df_2.BPoloidal_v_avg, label='Poloidal')
ax.plot(df_3.index, df_3.BToroidal_v_avg, label='Toroidal')
ax.legend()

today = datetime.now().strftime('%Y%m%d_%H%M%S')
plt.savefig('201206_Volume_Avg_Bfield_'+today+'.png', dpi=300)
plt.close()