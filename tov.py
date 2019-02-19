# -*- coding: utf-8 -*-
"""
tov.py

Use: This script solves the Tolman-Oppenheimer-Volkov (TOV) equations using 
a two component polytrope. 

Credit: 
This script is the work of Jolien Creighton and is detailed in a
presentation given on July 16th, 2012 entitled Lecture 2: Relativistic Stars

The code was transposed from the presentation slides by Sam Frederick on 
February 18th, 2019 for potential use in research. 

"""
import pylab
from scipy import integrate
from scipy.constants import pi, G, c, hbar, m_n
Msun = 1.98892e30

# piecewise polytrope equation of state
Gamma0 = 5.0/3.0 # low densities: soft non-relativistic degeneracy pressure
K0 = (3.0*pi**2)**(2.0/3.0)*hbar**2/(5.0*m_n**(8.0/3.0))
Gamma1 = 3.0 # high densities: stiffer equation of state
rho1 = 5e17
P1 = K0*rho1**Gamma0
K1 = P1/rho1**Gamma1

dr = 100.0
r = pylab.arange(10.0,20000.0,dr)
m = pylab.zeros_like(r)
P = pylab.zeros_like(r)


def eos(rho):
    if rho < rho1: 
        return K0*rho**Gamma0
    else: 
        return K1*rho**Gamma1
    
def inveos(P):
    if P < P1:
        return (P/K0)**(1.0/Gamma0)
    else:
        return (P/K1)**(1.0/Gamma1)

def tov(y,r):
    """ Tolman-Oppenheimer-Volkov Equations """
    Pval, mval = y[0],y[1]
    rho = inveos(Pval)
    dPdr = -G*(rho + Pval/c**2)*(mval+4.0*pi*r**3*Pval/c**2)
    dPdr = dPdr/(r*(r-2.0*G*mval/c**2))
    dmdr = 4.0*pi*r**2*rho
    return pylab.array([dPdr,dmdr])

def tovsolve(rhoc):
    """ Solves TOV equations given a central density """

    y = pylab.array(P[0],m[0])
    m[0] = 4.0*pi*r[0]**3*rhoc
    P[0] = eos(rhoc)
    i = 0 # integrate until density drops below zero 
    while P[i] > 0.0 and i < len(r) - 1:
        y= integrate.odeint(tov,y,r[i])
        print (y)
        P[i+1] = y[0]
        m[i+1] = y[1]
        i = i + 1
    # return mass and radius of star (units of stellar mass, km respectively)  
    return m[i-1]/Msun,r[i-1]/1000.0

# plot mass-radius curve
rhoc = pylab.logspace(17.5,20) # logspace range of central densities
M = pylab.zeros_like(rhoc)
R = pylab.zeros_like(rhoc)
for i in range(len(rhoc)):
    M[i], R[i] = tovsolve(rhoc[i])

pylab.plot(R,M)
pylab.xlabel('Radius (km)')
pylab.ylabel('Mass (solar)')
pylab.grid()
pylab.show()        