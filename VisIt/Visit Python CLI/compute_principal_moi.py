"""
Author: Sam Frederick
Date:  August 2019

This script automates the process of determining principal
MOI values (Ixx, Iyy, and Izz) for a region of the comput. 
domain defined by r <= R_star. 

A pseudocolor density plot is created, after which a clip 
operator sets the domain to be r <= R_star (1.0 in 
dimensionless units). The domain is the scaled up by a 
Transformation operator to the stellar radius in cm 
(X,Y,Z = 1e6 cm aka 10 km).

The DoublePrecisionQueryOverTime() function creates a text
file, to which MOI values are written in specified simulation
timesteps. Presently, the SetTimeSliderState is configured to 
iterate the timestep by 10 "cycles" or datafiles, which 
corresponds to a temporal timestep of 0.01 s as pluto.ini is 
presently configured. For iterating in steps of 10 cycles, the
range within the function should be set to the nearest rounded 
down integer corresponding to (Total Data File #) / 10. So if 
the number of datafiles is 1007, range should be 100. 


12-21-19 EDIT:
Added expressions for starting with 3d spherical hemisphere and 
converting to 2-d slice which is then revolved by Visit. Future 
itertion of this script may likely be intended for native 2-d 
simulation which is revolved by VisIT, making the runtime 
process for the simulation much quicker.      

"""



DeleteAllPlots()

print "Start\n---------------------"

# Create density plot with spherical clip operator 
AddPlot("Pseudocolor","rho")
p = PseudocolorAttributes()
p.colorTableName = "hot"
SetPlotOptions(p)

# (HEMISPHERE SPECIFIC CONDITION)
# Add Slice to convert 3-d to 2-d
#AddOperator("Slice")
#s = SliceAttributes()
#SetOperatorOptions(s)


# Revolve slice to 3-d sphere
#AddOperator("Revolve")
#r = RevolveAttributes()
#r.stopAngle = 180 # 180 for hemisphere 
#r.steps = 80
#SetOperatorOptions(r)

# Set spherical domain for r <= 1 (dimensionless radius)
AddOperator("Clip")
c = ClipAttributes()
c.quality = 1
c.funcType = 1
c.radius = 1
c.sphereInverse = 1
SetOperatorOptions(c)

# Transform cartesian coordinates to cm 
AddOperator("Transform")
# atts = GetOperatorOptions(1)
# print (atts)
t = TransformAttributes()
t.doScale = 1
t.scaleX = 1.0e6 # 10 km 
t.scaleY = 1.0e6
t.scaleZ = 1.0e6
SetOperatorOptions(t)

DrawPlots() 

# get the number of timesteps
nts = TimeSliderGetNStates()

# Write principal MOI values to file in specified simulation steps
# Function via http://visitusers.org/index.php?title=Query_over_time
def DoublePrecisionQueryOverTime():
   f = open("201206_MOI.txt", "w")
   f.write("              t                       Ixx                      Iyy                      Izz           \n")
   for ts in range(0, nts):
     TimeSliderSetState(ts) 
     Query("Time")
     t = GetQueryOutputValue()
     Query("Moment of Inertia")
     t2 = GetQueryOutputValue()
     str = "%25.15e %25.15e %25.15e %25.15e\n" %(t, t2[0],t2[4],t2[8]) # Print timestamp, Ixx, Iyy, and Izz
     f.write(str)
   f.close()


DoublePrecisionQueryOverTime()
