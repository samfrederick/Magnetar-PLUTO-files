# -*- coding: utf-8 -*-
"""
SLyEOS.py
Author: Sam Frederick
Created on Fri Feb 22 11:28:26 2019

Use: This script builds on tov_solver by utilizing the SLy EOS in place of 
piecewise polytropes. The analytic representation of SLy EOS was obtained via
Haensel, Potekhin (2004). The basis for the TOV-Solver code is based off work 
by Jolien Creighton and is detailed in a presentation given on July 16th, 2012
entitled Lecture 2: Relativistic Stars. The presentation can be accessed at 
http://www.ictp-saifr.org/schoolgr/Lecture2Creighton.pdf. 

"""
import math
import numpy as np
from bisect import bisect_left
import pylab
from scipy import integrate
from astropy import constants as const 

""" Physical Constants """
pi = np.pi
G = const.G.cgs.value
c = const.c.cgs.value
hbar = const.hbar.cgs.value
m_n = const.m_n.cgs.value
Msun = const.M_sun.cgs.value

""" Constants for zeta function """ 
a1 = 6.22
a2 = 6.121
a3 = 0.005925
a4 = 0.16326
a5 = 6.48
a6 =11.4971
a7 = 19.105
a8 = 0.8938
a9 = 6.54
a10 = 11.4950
a11 = -22.775
a12 = 1.5707
a13 = 4.3
a14 = 14.08
a15 = 27.80
a16 = -1.653
a17 = 1.50
a18 = 14.67

""" Constants and variables for TOV Solver """
Upper  = 2.0e6 #cm 
Lower = 1e3 #cm 
dr = (Upper - Lower)/19990.0
r = pylab.arange(Lower,Upper,dr)
m = pylab.zeros_like(r)
Pvals = pylab.zeros_like(r)
Plist = list()
Rlist = list()
Rholist = list()

""" Variables for HD search list"""
HDPlist = list()
HDrholist = np.logspace(1,16,50000)

"""
TOV-Solver Functions
"""
def tov(y,r):
    """ Tolman-Oppenheimer-Volkov Equations """
    Pval, mval = y[0],y[1]
    #print (Pval)
    rho = zetainverse(Pval)
    dPdr = -G*(rho + Pval/c**2)*(mval+4.0*pi*r**3*Pval/c**2)
    dPdr = dPdr/(r*(r-2.0*G*mval/c**2))
    dmdr = 4.0*pi*r**2*rho
    return pylab.array([dPdr,dmdr])

def tovsolve(rhoc):
    """ Solves TOV equations given a central density """
    m[0] = 4.0*pi*r[0]**3*rhoc
    Pvals[0] = expzeta(rhoc)
    y = pylab.array([Pvals[0],m[0]])
    i = 0 # integrate until density drops below zero 

    while Pvals[i] > 0.0 and i < len(r) - 1:
        y= integrate.odeint(tov,y,[r[i],r[i+1]])
        y = y[-1]       
        Pvals[i+1] = y[0]
        Plist.append(y[0])
        Rlist.append(r[i])
        rhoi = zetainverse(y[0])
        Rholist.append(rhoi)
        m[i+1] = y[1]
        i = i + 1
        if bool==1:
            newfile.write("%e  %e  %e\n" % (r[i], rhoi,y[0]))
    # return mass and radius of star (units of stellar mass, km respectively)  
    return m[i-1]/Msun,r[i-1]/1e6
    #print "Rho_c: %.3e," % rhoc, "Stellar Mass: %.3f," % (m[i-1]/Msun),"Stellar Radius: %.3f" % (r[i-1]/1e6)


    
"""
Functions for expressing analytic version of SLy EOS
"""    
""" Helper function """
def f0(x):
    value = (math.exp(x)+1)**(-1)
    return value 

""" Inputs log(rho) and outputs log(pressure) """
def zeta(xi):
    value = (a1+a2*xi+a3*xi**3)/(1+a4*xi)
    value = value*f0(a5*(xi-a6))
    value = value + (a7+a8*xi)*f0(a9*(a10-xi))
    value = value + (a11+a12*xi)*f0(a13*(a14-xi))
    value = value + (a15+a16*xi)*f0(a17*(a18-xi))
    return value 

""" 
Inputs exponential value for rho, converts to log10 and outputs exponential
value for pressure.   
"""
def expzeta(rho): # corresponds to eos() in tov_solver.py
    xi = np.log10(rho)
    value = 10**(zeta(xi))
    return value

""" 
Search a large list of pressure values and find closest value corresponding to
requested value, record its index, and return corresponding rho at corresponding
index. Purpose is to go the opposite direction of expzeta() by inputting pressure
and outputting density. 
"""
def zetainverse(Pval): # corresponds to inveos() in tov_solver.py
    position = bisect_left(HDPlist,Pval)
    rho = HDrholist[position]
    return rho

"""
Write large list of P and rho values to file
"""
newfile = open("SLYEOS_P_rho_list.txt","w")
for i in range(len(HDrholist)):
    P = expzeta(HDrholist[i])
    HDPlist.append(P)
    newfile.write("%e  %e\n" % (P, HDrholist[i]))
newfile.close()


"""
Generate solution to TOV equations with rho_core value
"""
rhoc = 1.75e14
M,R = tovsolve(rhoc)

""" Plots """
pylab.subplot(211)
pylab.plot(Rlist,Rholist,'-')
pylab.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
pylab.xlabel('Radius (cm)')
pylab.ylabel('Density g cm^-3')
pylab.grid()

pylab.subplot(212)
pylab.plot(Rlist,Plist,'-')
pylab.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
pylab.xlabel('Radius (cm)')
pylab.ylabel('Pressure (dyn cm^-2)')
pylab.grid()

#pylab.plot(HDrholist,Plist)
#pylab.xscale("log")
#pylab.yscale("log")

pylab.show()       
    
    
    
    
    
    
    
    
    
    
    