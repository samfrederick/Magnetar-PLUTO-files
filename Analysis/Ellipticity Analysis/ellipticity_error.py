#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 19:16:48 2021

@author: samfrederick
"""
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime


sns.set_style('ticks')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

os.chdir('/Users/samfrederick/Documents/GitHub/Magnetar-PLUTO-files/Analysis'
         '/Ellipticity Analysis')
ellip_err = pd.read_csv('delta_Izx_vs_angular_res.csv')

fig, ax = plt.subplots(1,1, figsize=(8,5.85))
ax.plot(ellip_err.iloc[:, 0], ellip_err.iloc[:, 1], marker='o')
ax.set_yscale('log')
ax.set_xlabel('Angular resolution ($n_{\\theta} = \pi/\Delta_{\\theta,\phi}$)',
              fontsize=18)
ax.set_ylabel('$\delta I_{zx}/I_0$', fontsize=18)
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)

plt.tight_layout()
today = datetime.now().strftime('%Y%m%d_%H%M')
plt.savefig('ellip_err_'+today+'.png', dpi=300)


