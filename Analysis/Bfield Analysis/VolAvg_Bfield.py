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
    fig, axs = plt.subplots(3, 1, figsize=(9, 7))
    plt.subplots_adjust(hspace=.35)
    
    time_lims = -0.5, mt.ceil(bfield_df.index[-1])
    
    #axs[0].set_title(r'Stellar Ellipticity ($\epsilon$)')
    axs[0].plot(df.index,
                df[['BTotal_v_avg', 'BPoloidal_v_avg', 'BToroidal_v_avg']])
    axs[0].legend(labels=[r'$B_{total}$', '$B_{p}$', '$B_{t}$'], loc='upper right',
                  fontsize=8, framealpha=.5)
    #axs[0].set_ylim(-2e17, 5.1e17)
    axs[0].set_xlim(time_lims)
    
    axs[1].set_title('First Derivative')
    axs[1].plot(df.index,
                df[['d_BTotal_v_avg', 'd_BPoloidal_v_avg', 'd_BToroidal_v_avg']])
    axs[1].axhline(y=0, color='#949494', linestyle='--')
    #axs[1].legend(labels=[r'$\partial_t$ $\epsilon$'], loc='lower right')
    axs[1].set_ylim(-2e17, 2e17)
    axs[1].set_xlim(time_lims)
    
    axs[2].set_title('Second Derivative')
    axs[2].plot(df.index, 
                df[['d_d_BTotal_v_avg', 'd_d_BPoloidal_v_avg', 'd_d_BToroidal_v_avg']])
    axs[2].axhline(y=0, color='#949494', linestyle='--')
    #axs[2].legend(labels=[r'$\partial^2_t$ $\epsilon$'], loc='lower right')
    axs[2].set_ylim(-1e18, 1e18)
    axs[2].set_xlim(time_lims)
    axs[2].set_xlabel('Time (s)')
    
    today = datetime.now().strftime('%Y%m%d_%H%M%S')
    plt.savefig('Bfield Analysis/201206_Volume_Avg_Bfield_'+today+'.png', dpi=300)
    #plt.close()
    
# ----------------------------------------------------------------------------

bfield_df = ImportBfieldData()
PlotBfield(bfield_df)