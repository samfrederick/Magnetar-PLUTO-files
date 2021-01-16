#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 10:02:27 2021

@author: samfrederick
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def Decaying_Sinusoid(t, offset, a, lam, omega):
    return offset + a*np.exp(lam*t)*np.cos(omega*t)

df_concat = df.loc[7.45:,'ellip']
mean_val = df_concat.loc[15.6:].mean()
df_centered = df_concat - mean_val
df_centered = df_centered.to_frame()
df_concat = df_concat.to_frame()


p_optimal, p_covariance = curve_fit(Decaying_Sinusoid,
                                    df_concat.index,
                                    df_concat.ellip,
                                    p0=(0.072, -.04, -0.05, 2))

func_label = '$\hat{y}(t) = A + B e^{\lambda t}\cos(\omega_0 t)$'

fig, ax = plt.subplots(1, 1)
ax.plot(df.index, df.ellip, alpha=0.7, color='gray', label='Ellipticity ($\epsilon$)')
ax.plot(df_concat.index, df_concat.ellip, '.', color='blue',
        label='Fitting $\epsilon$ data')
ax.plot(np.linspace(0,30,601),
        Decaying_Sinusoid(np.linspace(0,30,601),
                          p_optimal[0], p_optimal[1], p_optimal[2], p_optimal[3]),
        color='red', alpha=0.6, label=func_label, linestyle='--')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Ellipticity')

print(p_optimal[0], p_optimal[1], p_optimal[2], p_optimal[3])