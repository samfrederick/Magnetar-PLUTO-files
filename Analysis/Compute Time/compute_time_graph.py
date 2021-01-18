#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  6 06:29:09 2019

@author: samfrederick
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import seaborn as sns

os.chdir('/Users/samfrederick/Documents/GitHub/Magnetar-PLUTO-files/Analysis'
         '/Compute Time')

sns.set_style('ticks')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

compute_data = pd.read_csv('simple_compute_time.csv')

def f(x):
    # Interpolating quadratic function for data
    yval = 1.1263 + (-0.1435)*x + (0.01643)*x**2
    return yval

xinterp = np.arange(10,71,1)
yinterp = f(xinterp)

plt.figure(figsize=(8, 5.85))

# Configure tick positioning on both sides of axes   
ax = plt.subplot(111)    
#ax.xaxis.set_ticks_position('both')
#ax.yaxis.set_ticks_position('both')
#ax.tick_params(direction='in')     

plt.scatter(
        compute_data.nthetaphi.values,
        compute_data.timetosol.values,
        s = .1)

plt.errorbar(compute_data.nthetaphi.values,compute_data.timetosol.values,
             yerr=compute_data.err.values,fmt='o')

plt.plot(xinterp,yinterp,'--',color='lightcoral')
plt.xlabel('Angular resolution ($n_{\\theta} = \pi/\Delta_{\\theta,\phi}$)',fontsize = 18)
plt.xticks(fontsize = 15)
plt.ylabel('Time to Solution (hours)',fontsize = 18)
plt.yticks(fontsize = 15)

plt.tight_layout()
plt.savefig('compute_time.png', dpi=300)