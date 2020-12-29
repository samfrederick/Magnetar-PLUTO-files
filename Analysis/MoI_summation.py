# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 11:28:49 2019

@author: safrederick
"""
import numpy as np
import pylab as pt

M = (4.0*np.pi)/3.0
R = 1
I0 = (2.0/5.0)*M*R**2

rgrid = 50.0
tgrid = 63.0
pgrid = 80.0

rlist = np.linspace(0,R,rgrid)
thetalist = np.linspace(0,np.pi,tgrid)
philist = np.linspace(0,2*np.pi,pgrid)

test1 = list()
test2 = list ()
actual = (8*np.pi)/15.0


Izz = 0
Ixx = 0
rho = 1

for r in rlist:
    for theta in thetalist:
        for phi in philist:
            Izz += rho*((r*np.sin(theta)*np.sin(phi))**2+(r*np.cos(theta))**2)*(r**2)*np.sin(theta)*((1/rgrid)*((np.pi)/tgrid)*((2*np.pi)/pgrid))
            Ixx += rho*((r*np.sin(theta)*np.cos(phi))**2+(r*np.sin(theta)*np.sin(phi))**2)*(r**2)*np.sin(theta)*((1/rgrid)*((np.pi)/tgrid)*((2*np.pi)/pgrid))
            pxxerr = (np.abs(Ixx-actual))
            pyyerr = (np.abs(Ixx-actual))
            test1.append(pxxerr)
            test2.append(pyyerr)


print (1-Ixx/actual,1-Izz/actual)

ellip = (Izz-Ixx)/I0
print (Ixx/actual,Izz/actual,ellip)
pt.plot(test1,'o')
#pt.plot(test2,'o')
pt.yscale('log')
pt.show()            