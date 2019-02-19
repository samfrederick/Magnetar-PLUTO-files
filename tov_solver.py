"""
tov.py
Date Created: 2-18-19

Use: This script solves the Tolman-Oppenheimer-Volkov (TOV) equations using 
a two component polytrope. 

Credit: 
This script is the work of Jolien Creighton and is detailed in a
presentation given on July 16th, 2012 entitled Lecture 2: Relativistic Stars.
The presentation can be accessed at 
http://www.ictp-saifr.org/schoolgr/Lecture2Creighton.pdf. 

Author's Note:
The code was modified from the presentation slides by Sam Frederick on 
February 18th, 2019 for potential use in research. Integration method was 
changed to scipy.integrate.odeint and various corrections were made to 
properly generate stuctural variable values in tovsolve. Ability to write these 
variables (radius, density, and pressure) to text files was added to assist in 
building strucutral models for neutron stars built off this code. 

One may enticipate adding a more complex EOS such as SLy EOS, which would 
replace the piecewise polytrope defined in the eos() function.  
"""
import pylab
import numpy as np
from scipy import integrate
from astropy import constants as const 
#from scipy.constants import pi, G, c, hbar, m_n

pi = np.pi
G = const.G.cgs.value
c = const.c.cgs.value
hbar = const.hbar.cgs.value
m_n = const.m_n.cgs.value
Msun = const.M_sun.cgs.value
#Msun = 1.98892e30
#rhoc = 2.0e15 

# piecewise polytrope equation of state
Gamma0 = 5.0/3.0 # low densities: soft non-relativistic degeneracy pressure
K0 = (3.0*pi**2)**(2.0/3.0)*hbar**2/(5.0*m_n**(8.0/3.0))
Gamma1 = 3.0 # high densities: stiffer equation of state
rho1 = 5.0e14
P1 = K0*rho1**Gamma0
K1 = P1/(rho1**Gamma1)

Upper  = 2.0e6 #cm 
Lower = 1e3 #cm 
dr = (Upper - Lower)/1999.0
r = pylab.arange(Lower,Upper,dr)

m = pylab.zeros_like(r)
P = pylab.zeros_like(r)
Pvals = list()
Rvals = list()
Rhovals = list()

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
    #print (Pval)
    rho = inveos(Pval)
    dPdr = -G*(rho + Pval/c**2)*(mval+4.0*pi*r**3*Pval/c**2)
    dPdr = dPdr/(r*(r-2.0*G*mval/c**2))
    dmdr = 4.0*pi*r**2*rho
    return pylab.array([dPdr,dmdr])

def tovsolve(rhoc):
    """ Solves TOV equations given a central density """
    m[0] = 4.0*pi*r[0]**3*rhoc
    P[0] = eos(rhoc)
    y = pylab.array([P[0],m[0]])
    i = 0 # integrate until density drops below zero 

    while P[i] > 0.0 and i < len(r) - 1:
        y= integrate.odeint(tov,y,[r[i],r[i+1]])
        y = y[-1]       
        P[i+1] = y[0]
        Pvals.append(y[0])
        Rvals.append(r[i])
        rhoi = inveos(y[0])
        Rhovals.append(rhoi)
        m[i+1] = y[1]
        i = i + 1
        if bool==1:
            newfile.write("%e  %e  %e\n" % (r[i], rhoi,y[0]))
    # return mass and radius of star (units of stellar mass, km respectively)  
    return m[i-1]/Msun,r[i-1]/1e6
    print "Rho_c: %.3e," % rhoc, "Stellar Mass: %.3f," % (m[i-1]/Msun),"Stellar Radius: %.3f" % (r[i-1]/1e6)

#M, R = tovsolve(rhoc)
#print (M,R)


''' Generate files for four intitial core density values and plot results '''
rhoc = pylab.logspace(14.9,15.5,num=4)
for i in range(len(rhoc)):
    bool = 1
    temprho = str(rhoc[i])
    newfile = open("EOS_rhoc_%.3e.txt" % rhoc[i],"w+")
    newfile.write("   Radius        Density      Pressure     \n")
    tovsolve(rhoc[i])
    newfile.close()
    pylab.plot(Rvals,Rhovals,'-')
#pylab.plot(Rhovals,Pvals)
pylab.plot(Rvals,Rhovals,'-')
#pylab.yscale("log")
#pylab.xscale("log")
pylab.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
pylab.xlabel('Radius (cm)')
pylab.ylabel('Density g cm^-3')
pylab.grid()
pylab.show()
#------------------------------------------------------------------------------#

''' plot mass-radius curve '''
rhoc = pylab.logspace(14.5,17) # logspace range of central densities
M = pylab.zeros_like(rhoc)
R = pylab.zeros_like(rhoc)
for i in range(len(rhoc)):
    bool = 0
    M[i], R[i] = tovsolve(rhoc[i])
pylab.plot(R,M)
pylab.xlabel('Radius (1e6 cm)')
pylab.ylabel('Mass (solar)')
pylab.grid()
pylab.show()   