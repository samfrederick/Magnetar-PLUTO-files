"""
Mixed_Field_Visualization.py
Author: Sam Frederick
2-11-19

This script graphs mixed magnetic field component equations from 
Haskell et al. (2008). The function A(r) is called the stream 
function, and its derivative is the function dA(r). This function 
and its derivative appear in the field component equations: Br, 
Btheta, and Bphi. 

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
#print r
Alist = list()
DAlist = list()
bx1list = list()
bx2list = list()
bx3list = list()

x2 = .25*CONST_PI # Theta set to pi/4


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
    if x1 != 0:
        DAval = (A(x1+h) - A(x1-h))/(2*h)
    else:
        DAval = 0


    return (DAval)


for x1 in r:
    Aval = A(x1)
    DAvalue = dA(x1)
    Alist.append(Aval)
    DAlist.append(DAvalue)

    if x1 != 0:
        Bx1 = (2*A(x1)*np.cos(x2))/((x1*R)*(x1*R))
        Bx2 = (-dA(x1)*np.sin(x2))/(x1*R*R)
        Bx3 = (Lambda*CONST_PI*A(x1)*np.sin(x2))/(x1*R*R)

    else:
        Bx1 = (2*A(0.0001)*np.cos(x2))/((0.0001*R)*(0.0001*R))
        Bx2 = (-dA(0.0001)*np.sin(x2))/(0.0001*R*R) # 
        Bx3 = 0

#    Bx1 = Bx1 / (np.sqrt(UNIT_DENSITY)*UNIT_VELOCITY)
#    Bx2 = Bx2 / (np.sqrt(UNIT_DENSITY)*UNIT_VELOCITY)
#    Bx3 = Bx3 / (np.sqrt(UNIT_DENSITY)*UNIT_VELOCITY)

    bx1list.append(Bx1)
    bx2list.append(Bx2)
    bx3list.append(Bx3)

print(bx3list)

py.plot(r,bx1list,"bo-",label="Br")
py.plot(r,bx2list,"ro-",label="Btheta")
py.plot(r,bx3list,"go-",label="Bphi")

#py.plot(r,Alist,"o-",label="A(r)")
#py.plot(r,DAlist,"o-",label="A'(r)")

#py.yscale("log")
py.xlim(0,2.1)
#py.ylim(-1,1)
py.ylabel("Strength of Field")
py.xlabel("Radius")
py.legend(loc="upper right")
py.show()
