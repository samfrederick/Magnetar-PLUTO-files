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
BMAX = 1e15
RMAX = 2.0

UNIT_DENSITY = 1.e8
UNIT_LENGTH = 1.e10
UNIT_VELOCITY = 1.e10

rin = np.linspace(0,1,99,endpoint=False)
rout = np.linspace(1,2,30)
r = np.concatenate((rin,rout),axis=0)
rghost = [2.01]
r = np.concatenate((r,rghost),axis=0)
print r
bx1list = list()
bx2list = list()
bx3list = list()

x2 = 1.7*CONST_PI # Theta set to pi/4 

Bx1 = 0.0
Bx2 = 0.0
Bx3 = 0.0

for x1 in r:
    if x1 < 1.0: # NEED TO ADD CONDITION FOR R == 0! 
    
        Bx1 = CONST_PI*CONST_PI*CONST_PI*x1*x1*x1 + 3*(CONST_PI*CONST_PI*x1*x1 -2)*np.sin(CONST_PI*x1)+6.0*CONST_PI*x1*np.cos(CONST_PI*x1)
        Bx1 = Bx1*(BMAX*np.cos(x2))/(CONST_PI*(CONST_PI*CONST_PI-6))
   
        Bx2 = -2*CONST_PI*CONST_PI*CONST_PI*x1*x1*x1+ 3*(CONST_PI*CONST_PI*x1*x1-2)*(np.sin(CONST_PI*x1)-CONST_PI*x1*np.cos(CONST_PI*x1))
        Bx2 = Bx2*(BMAX*np.sin(x2))/(2.0*CONST_PI*(CONST_PI*CONST_PI-6))
    
 
        Bx3 = (BMAX*np.sin(CONST_PI*x1)*np.sin(x2))/CONST_PI
    
    if x1 >= 1.0:    
        Bx1 = (BMAX*np.cos(x2))/(x1*x1*x1)
        Bx2 = (BMAX*np.sin(x2))/(2.0*x1*x1*x1)
        Bx3 = 0        
    if x1 > 2.0:
        Bx1 = (BMAX*np.cos(x2))/(RMAX*RMAX*RMAX)
        Bx2 = (BMAX*np.sin(x2))/(2.0*RMAX*RMAX*RMAX)
        Bx3 = 0.0       
    if x2 == 0:
        #need a condition for interior to star and exterior to star! B-field has different solution in these regions!
        if x1 < 1.0:
            Bx1 = CONST_PI*CONST_PI*CONST_PI*x1*x1*x1 + 3*(CONST_PI*CONST_PI*x1*x1 -2)*np.sin(CONST_PI*x1)+6.0*CONST_PI*x1*np.cos(CONST_PI*x1)
            Bx1 = Bx1*(BMAX*1)/(CONST_PI*(CONST_PI*CONST_PI-6))
        
            Bx2 = 0.0 
            Bx3 = 0.0 
        else:
            Bx1 = (BMAX*1)/(x1*x1*x1)
        
            Bx2 = 0.0 
            Bx3 = 0.0
    if x2 == 2*CONST_PI:
        #need a condition for interior to star and exterior to star! B-field has different solution in these regions!
        if x1 < 1.0:
            Bx1 = CONST_PI*CONST_PI*CONST_PI*x1*x1*x1 + 3*(CONST_PI*CONST_PI*x1*x1 -2)*np.sin(CONST_PI*x1)+6.0*CONST_PI*x1*np.cos(CONST_PI*x1)
            Bx1 = Bx1*(BMAX*-1)/(CONST_PI*(CONST_PI*CONST_PI-6))
        
            Bx2 = 0.0 
            Bx3 = 0.0 
        else:
            Bx1 = (BMAX*-1)/(x1*x1*x1)
        
            Bx2 = 0.0 
            Bx3 = 0.0        
    
    Bx1 = Bx1 / (np.sqrt(UNIT_DENSITY)*UNIT_VELOCITY)
    Bx2 = Bx2 / (np.sqrt(UNIT_DENSITY)*UNIT_VELOCITY)    
    Bx3 = Bx3 / (np.sqrt(UNIT_DENSITY)*UNIT_VELOCITY)

    
    bx1list.append(Bx1)
    bx2list.append(Bx2)
    bx3list.append(Bx3)
    

py.plot(r,bx1list,"bo-")
py.plot(r,bx2list,"ro-")
py.plot(r,bx3list,"go-")
#py.yscale("log")
py.xlim(0,2.1)
#py.ylim(-1,1)
#py.ylabel("Log(Density)")
#py.xlabel("Radius")
py.show()
