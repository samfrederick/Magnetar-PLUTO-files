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
import os
import seaborn as sns
os.chdir('/Users/samfrederick/Documents/GitHub/Magnetar-PLUTO-files/Analysis'
         '/Wavestrain Analysis')

from Detector_Curves import DECIGO_curve, BBO_curve
sns.set_style('ticks')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Load McGill catalog data
df = pd.read_csv('Mcgill_catalog_data.csv')
df = df.drop(index=16)

# cgs
G = 6.6743E-8
c = 2.998E10 
R = 1E6
rho_c = 2E15
ellip = 0.087
I = (8*(np.pi**2-6)*R**5*rho_c) / (3*np.pi**3)
const = (4*np.pi**2*G/c**4)*I*ellip

fgw_vals = df['Fgw']
d_vals = df['Distance']*3.086E21 # convert dist from kpc to cm 
df['wavestrain_est'] = (const*fgw_vals**2)/d_vals

d_plserr_vals = d_vals + df['d_err_plus']*3.086E21

d_mnserr_vals = d_vals - df['d_err_minus']*3.086E21

df['w_pluserr'] = (const*fgw_vals**2)/d_plserr_vals
df['w_minuserr'] = (const*fgw_vals**2)/d_mnserr_vals


# Detector sensitivity curves
freqs = np.logspace(-3,2, 100)
DECIGO_curve_vals = DECIGO_curve(freqs)
BBO_curve_vals = BBO_curve(freqs)

# Create plot
fig, ax = plt.subplots(1,1, figsize=(12,8.0))

# Plot detector curves
ax.loglog(freqs, DECIGO_curve_vals,  color = 'r')
ax.loglog(freqs, BBO_curve_vals,color = 'b')
ax.text(x = 3e-1, y = 4e-24, s='DECIGO', color='r', fontsize=15)
ax.text(x = 3e-1, y = 3e-25, s='BBO', color='b', fontsize=15)

# Plot mcgill source wavestrains, error bars
ax.loglog(df['Fgw'], df['wavestrain_est'], '*', c='crimson',
        markeredgecolor='none', markersize = 12, alpha = .8)
ax.errorbar(df.Fgw, df.Wavestrain, yerr=(df['w_minuserr'] ,df['w_pluserr']),
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



