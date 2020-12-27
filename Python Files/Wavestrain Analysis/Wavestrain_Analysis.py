#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 08:40:40 2020

@author: samfrederick
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from datetime import datetime
from Detector_Curves import DECIGO_curve, BBO_curve

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Load McGill catalog data
df = pd.read_csv('Mcgill_catalog_data.csv')
df = df.drop(index=16)

# Detector sensitivity curves
freqs = np.logspace(-3,2, 100)
DECIGO_curve_vals = DECIGO_curve(freqs)
BBO_curve_vals = BBO_curve(freqs)

# Create plot
fig, ax = plt.subplots(1,1, figsize=(12,9))

# Plot detector curves
ax.loglog(freqs, DECIGO_curve_vals,  color = 'r')
ax.loglog(freqs, BBO_curve_vals,color = 'b')
ax.text(x = 3e-1, y = 4e-24, s='DECIGO', color='r', fontsize=15)
ax.text(x = 3e-1, y = 3e-25, s='BBO', color='b', fontsize=15)

# Plot mcgill source wavestrains, error bars
ax.loglog(df.Fgw, df.Wavestrain, '*', c='crimson',
        markeredgecolor='none', markersize = 12, alpha = .8)
ax.errorbar(df.Fgw, df.Wavestrain, yerr=(df['w_plus'] ,df['w_minus']),
             fmt='none', c= 'crimson', alpha=.4)

# Set plot formatting
ax.grid(linestyle='dashed')
ax.set_xlabel('Gravitational Wave Strain Frequency (Hz)', fontsize = 20)
ax.set_ylabel(r'Strain Sensitivity (Hz$^{-1/2}$)', fontsize = 20)
ax.set_title("Gravitational Wave Strain Estimates for McGill Catalog Sources",
             fontsize = 20)
ax.set_yscale('log')
ax.set_xscale('log')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.set_ylim(2e-30, 1e-22)
ax.set_xlim(1e-2, 1e1)
ax.tick_params('both', labelsize=15)

today = datetime.now().strftime('%Y%m%d_%H%M%S')
plt.savefig('McGill_Strain_Estimates_'+ today + '.png',
            bbox_inches='tight', dpi=300)
plt.close()



