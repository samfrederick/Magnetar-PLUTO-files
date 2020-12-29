"""
Author: Sam Frederick
Date: 10-15-18

This script graphs density and pressure as functions of radius over the computational
domain for (0<r<=2.) given equations specifying these variables for r = 0, r < 1.0,
and r>= 1.0.
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
VACUUM = 5e10
K = 4.25e4

UNIT_DENSITY = 1.e8
UNIT_LENGTH = 1.e10
UNIT_VELOCITY = 1.e10

r = np.linspace(0,2,70)
pressure = list()
rho = list()

for x1 in r:
    if x1 < 1.0 and x1 != 0.0:
        rhoval = (RHO_C*np.sin(CONST_PI*x1))/(x1*CONST_PI) +VACUUM
        pressureval = K*rhoval*rhoval

    elif x1 == 0:
        rhoval = RHO_C + VACUUM
        pressureval = K*RHO_C*RHO_C

    elif x1 >= 1.0:
        rhoval = VACUUM;
        pressureval = (K*VACUUM*VACUUM);

    rhoval = rhoval / UNIT_DENSITY
    pressureval = pressureval / (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY)

    rho.append(rhoval)
    pressure.append(pressureval)


py.plot(r,rho,"bo-")
#py.plot(r,pressure,"ro-")
py.yscale("log")
py.title("Log(Density) vs. Radius ")
py.ylabel("Log(Density)")
py.xlabel("Radius")
py.show()
