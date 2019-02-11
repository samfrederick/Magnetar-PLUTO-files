"""
Mixed_Field_Visualization.py
Author: Sam Frederick


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
VACUUM = 1e8
K = 4.25e4
BMAX = 1e15
RMAX = 2.0

h  = 1e-4
Lambda  = 2.362

UNIT_DENSITY = 1.e8
UNIT_LENGTH = 1.e10
UNIT_VELOCITY = 1.e10

rin = np.linspace(0,1,99,endpoint=False)
rout = np.linspace(1,2,30)
r = np.concatenate((rin,rout),axis=0)
rghost = [2.01]
r = np.concatenate((r,rghost),axis=0)
print r
Alist = list()
bx1list = list()
bx2list = list()
bx3list = list()

x2 = 1.7*CONST_PI # Theta set to pi/4

Bx1 = 0.0
Bx2 = 0.0
Bx3 = 0.0

def A(x1):
    if x1 != 0:
      Aval = ((BMAX*R*R)/((Lambda*Lambda-1)*(Lambda*Lambda-1)*CONST_PI*x1))*  \
      (2*CONST_PI*((Lambda*CONST_PI*x1*np.cos(Lambda*CONST_PI*x1)-np.sin(Lambda*CONST_PI*x1))/  \
      (CONST_PI*Lambda*np.cos(CONST_PI*Lambda)-np.sin(CONST_PI*Lambda))) +  \
      ((1-Lambda*Lambda)*(CONST_PI*x1)*(CONST_PI*x1)-2)*np.sin(CONST_PI*x1) +  \
      2*CONST_PI*x1*np.cos(CONST_PI*x1))
    else:
        Aval = 0

    return (Aval)

def dA(x1):
    if x1 != 0 and x1 != 2.0:
      D_Aval = (A(x1+h) - A(x1-h))/(2*h)  /* Central Difference Approx */
    else
      D_Aval = 0     /* Forward Differnce Approx */



  return D_Aval


for x1 in r:
    Aval = A(x1)
    Alist.append(Aval)


"""
    Bx1 = Bx1 / (np.sqrt(UNIT_DENSITY)*UNIT_VELOCITY)
    Bx2 = Bx2 / (np.sqrt(UNIT_DENSITY)*UNIT_VELOCITY)
    Bx3 = Bx3 / (np.sqrt(UNIT_DENSITY)*UNIT_VELOCITY)

    bx1list.append(Bx1)
    bx2list.append(Bx2)
    bx3list.append(Bx3)


py.plot(r,bx1list,"bo-")
py.plot(r,bx2list,"ro-")
py.plot(r,bx3list,"go-")
"""
py.plot(r,Alist)
#py.yscale("log")
py.xlim(0,2.1)
#py.ylim(-1,1)
#py.ylabel("Log(Density)")
#py.xlabel("Radius")
py.show()
