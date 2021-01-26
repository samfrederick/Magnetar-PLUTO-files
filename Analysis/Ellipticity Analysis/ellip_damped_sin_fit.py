#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 10:02:27 2021

@author: samfrederick
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

sns.set_style('darkgrid')
def Decaying_Sinusoid(t, offset, a, lam, omega, c):
    return offset + a*np.exp(lam*t)*np.cos(omega*t) + c*(t**-1)

fit_start = 7.5

df_concat = df.loc[fit_start:,'ellip']
mean_val = df_concat.loc[15.6:].mean()
df_centered = df_concat - mean_val
df_centered = df_centered.to_frame()
df_concat = df_concat.to_frame()


p_optimal, p_covariance = curve_fit(Decaying_Sinusoid,
                                    df_concat.index,
                                    df_concat.ellip,
                                    p0=(0.07, -.04, -0.05, 2, 0.01))

func_label = ('$\hat{y}(t) = A + B e^{\lambda t}\cos(\omega_0 t) '
              '+ C\left(\\frac{1}{t}\\right)$')

fig, ax = plt.subplots(1, 1, figsize=(6, 3))
ax.plot(df.index, df.ellip, alpha=0.7, color='gray',
        label='Ellipticity ($\epsilon$)')
ax.plot(df_concat.index, df_concat.ellip, '.', color='#1d5bc4',
        label='Fitting $\epsilon$ data')
ax.plot(np.linspace(fit_start, 30, 601),
        Decaying_Sinusoid(np.linspace(fit_start,30,601),
                          p_optimal[0], p_optimal[1], p_optimal[2],
                          p_optimal[3], p_optimal[4]),
        color='red', alpha=0.6, label=func_label, linestyle='--')

ax.axhline(y=0, color='#949494', linestyle='--')

plt.legend()
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Ellipticity', fontsize=12)
plt.legend(loc='lower right', fontsize=11)

plt.ylim(-0.15, 0.15)
plt.xlim(-1, 27)
ax.tick_params(axis='both', labelsize=11) 

fmt_str = 'A={:4.3f}, B={:4.3f}, Lambda={:4.3f}, Omega={:4.3f}, C={:4.3f}'
print(fmt_str.format(p_optimal[0], p_optimal[1],
                     p_optimal[2], p_optimal[3], 
                     p_optimal[4]))

plt.tight_layout()

plt.savefig('Ellipticity Analysis/Ellip_fit.png', dpi=300)