#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 07:48:06 2019

@author: samfrederick
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

magnetar_dat = pd.read_csv('Mcgill_catalog_data.csv')

#magnetar_dat = magnetar_dat.drop(['psqr','1/d'],axis=1)
magnetar_dat = magnetar_dat.drop([16,21])

fig = plt.figure(figsize=(12,9))

ax = plt.gca()
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(direction='in') 

ax.plot(magnetar_dat.Fgw.values,
        magnetar_dat.Wavestrain.values,
        '*',
        c='crimson',
        markeredgecolor='none',
        markersize = 12, 
        alpha = .8
        )


plt.errorbar(magnetar_dat.Fgw.values, magnetar_dat.Wavestrain.values, \
             yerr=(magnetar_dat['w_plus'] ,magnetar_dat['w_minus']), \
             fmt='none', c= 'crimson', alpha=.4)

ax.set_yscale('log')

plt.grid(linestyle='dashed')
plt.xlabel('Gravitational Wave Strain Frequency (Hz)', fontsize = 22)
plt.xticks(fontsize = 17)
plt.ylabel(r'Strain Sensitivity (Hz$^{-1/2}$)',fontsize = 22)
plt.yticks(fontsize = 17)
plt.title("Gravitational Wave Strain Estimates for McGill Catalog Sources",fontsize = 24)
#plt.ylim(10**(-40), 10**(-1))

plt.savefig('McGill_Strain_Estimates.png',bbox_inches='tight', dpi=300)