#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 07:50:16 2020

@author: samfrederick
"""
import os
os.chdir('/Users/samfrederick/Documents/GitHub/Magnetar-PLUTO-files/Analysis')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from differentiate import TimeDerivative # Must set path to Analysis folder
import math as mt
sns.set_style('darkgrid')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def ImportBfieldData():
    total = pd.read_csv('Bfield Analysis/201206_VolAvg_Total_Bfield.csv',
                        index_col='t')
    poloidal = pd.read_csv('Bfield Analysis/201206_VolAvg_Poloidal_Bfield.csv',
                           index_col='t')
    toroidal = pd.read_csv('Bfield Analysis/201206_VolAvg_Toroidal_Bfield.csv',
                           index_col='t')
    
    bfield_df = total.join(poloidal).join(toroidal)
    
    # Compute first time derivative
    for col in bfield_df:
        bfield_df = TimeDerivative(bfield_df, col, 't')
    
    # Compute second time derivative
    for col in [col for col in bfield_df.columns if col.startswith('d_')]:
        bfield_df = TimeDerivative(bfield_df, col, 't')
    
    return bfield_df
    
def PlotBfield(df):    
    fsize = 12
    fig, axs = plt.subplots(3, 1, figsize=(6, 7))
    
    time_lims = -0.5, mt.ceil(bfield_df.index[-1])
    
    comp_names = ['BPoloidal_v_avg', 'BToroidal_v_avg']
    d_comp_names = ['d_BPoloidal_v_avg', 'd_BToroidal_v_avg']
    dd_comp_names = ['d_d_BPoloidal_v_avg', 'd_d_BToroidal_v_avg']
    
    axs[0].set_title('Volume-averaged Component Fields', fontsize=fsize)
    axs[0].plot(df.index,
                df[comp_names],
                linewidth=2)
    axs[0].legend(labels=['$\\bar{B}_{p}$', '$\\bar{B}_{t}$'],
                  loc='lower left', fontsize=10, framealpha=1)
    axs[0].set_xlim(-1, 22)
    axs[0].tick_params(axis='both', labelsize=11) 
    
    axs[1].set_title('First Derivative', fontsize=fsize)
    axs[1].plot(df.index,
                df[d_comp_names],
                linewidth=2)
    axs[1].axhline(y=0, color='#949494', linestyle='--')
    axs[1].legend(labels=['$\partial_t$ $\\bar{B}_{p}$',
                          '$\partial_t$ $\\bar{B}_{t}$'],
                  loc='lower left',fontsize=10)
    axs[1].set_ylim(-2e17, 2e17)
    axs[1].set_xlim(-1, 22)
    axs[1].tick_params(axis='both', labelsize=11) 
    
    axs[2].set_title('Second Derivative', fontsize=fsize)
    axs[2].plot(df.index, 
                df[dd_comp_names],
                linewidth=2)
    axs[2].axhline(y=0, color='#949494', linestyle='--')
    axs[2].legend(labels=['$\partial^2_t$ $\\bar{B}_{p}$',
                          '$\partial^2_t$ $\\bar{B}_{t}$'],
                  loc='lower left',fontsize=10)
    axs[2].set_ylim(-1e18, 1e18)
    axs[2].set_xlim(-1, 22)
    axs[2].set_xlabel('Time (s)', fontsize=fsize)
    axs[2].tick_params(axis='both', labelsize=11) 
    
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top = 0.9, hspace=.5)
    plt.tight_layout()

    today = datetime.now().strftime('%Y%m%d_%H%M%S')
    plt.savefig('Bfield Analysis/Volume_Avg_Bfield.png', dpi=300)
    #plt.close()
    
# ----------------------------------------------------------------------------

bfield_df = ImportBfieldData()
PlotBfield(bfield_df.dropna())