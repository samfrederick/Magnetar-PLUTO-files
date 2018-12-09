"""
Author: Sam Frederick
Date: 10-12-18

This script graphs the graviational potential for the computational domain 
given equations specifying the potential in three main radial regions: r = 0, 
0 < r < 1, and r > 1. 
"""

import numpy as np
import pylab as py

RPOT = 1.857595e20 # CORRECTED VALUE FOR GM/R, was 6.67e19 prior to correction (likely
# due to temporary M_STAR mass of 1.0e33)
GPRSQ = 1.4674e20 
G_CONST = 6.67e-8 
M_STAR = 2.785e33 
RHO_C = 2.2e15          
R = 1.e6 
CONST_PI = 3.1415926

r = np.linspace(0,2,70)
phivals = list()

for i in r:
    
    if i <= 1 and i != 0:
        phi = (-4*G_CONST*RHO_C*R**2*np.sin(CONST_PI*i))/(CONST_PI*CONST_PI*i) - RPOT
        phivals.append(phi)
         
    if i > 1:
        phi = -(G_CONST*M_STAR) / (R*i)
        phivals.append(phi)
     
    if i == 0:
        phi = 4*G_CONST*RHO_C*(-(R*R)/(CONST_PI)-(M_STAR)/(4*R*RHO_C))
        phivals.append(phi)
 
#print phivals
#print len(phivals)
#print len(r)


py.plot(r,phivals,"bo-")
py.title("Gravitational Potential (Phi) vs. Radius")
py.ylabel("Gravitational Potential")
py.xlabel("Radius")

py.show()