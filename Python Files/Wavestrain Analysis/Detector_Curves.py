#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 22:06:59 2020

@author: samfrederick
"""
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt


def DECIGO_curve(f):
    """
    Noise curve for DECIGO. Numerical model via 
    https://arxiv.org/pdf/1101.3940.pdf Equation 5.
    
    Term 1: Shot noise
    Term 2: Radiation pressure noise
    Term 3: Acceleration noise
    """
    fp = 7.36 # Hz
    term1 = 7.05e-48*(1 + (f/fp)**2)
    term2 = (4.8e-51*((f)**-4))*(1/(1 + (f/fp)**2))
    term3 = 5.33e-52*((f)**-4)
    
    S = term1 + term2 + term3
    
    return S**(1/2)

def BBO_curve(f):
    """
    Noise curve for BBO. Numerical model via 
    https://arxiv.org/pdf/1101.3940.pdf Equation 6.
    """
    term1 = 2.00e-49*(f)**2
    term2 = 4.58e-49
    term3 = 1.26e-51*(f)**-4

    S = term1 + term2  + term3

    return S**(1/2)

def Plot_Detector_Curves():
    freqs = np.logspace(-3,2, 100)

    DECIGO_curve_vals = DECIGO_curve(freqs)
    BBO_curve_vals = BBO_curve(freqs)

    fig, ax = plt.subplots(1,1,figsize=(8,6))
    ax.plot(freqs, DECIGO_curve_vals,  color = 'r')
    ax.plot(freqs, BBO_curve_vals,color = 'b')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.text(x = 3e-3, y = 2e-20, s='DECIGO', color='r')
    ax.text(x = 2e-3, y = 1e-21, s='BBO', color='b')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('S$_h^{1/2}$ (Hz$^{-1/2}$)')
