"""
Author: Sam Frederick
Date:  December 2020

VisIt CLI Python command
------------------------

Compute inertia tensor component values for stellar 
region defined by r <=R_star=1.0 on each timestep.

"""
# Clear any existing plots
DeleteAllPlots()

print "Computing Inertia Tensor\n"

# Create density plot with spherical clip operator 
AddPlot("Pseudocolor","rho")
p = PseudocolorAttributes()
p.colorTableName = "hot"
SetPlotOptions(p)

# Set spherical domain for r <= 1 (stellar radius)
AddOperator("Clip")
c = ClipAttributes()
c.quality = 1
c.funcType = 1
c.radius = 1
c.sphereInverse = 1
SetOperatorOptions(c)

# Transform unit coords to cm (R_star=10 km)
AddOperator("Transform")
t = TransformAttributes()
t.doScale = 1
t.scaleX = 1.0e6
t.scaleY = 1.0e6
t.scaleZ = 1.0e6
SetOperatorOptions(t)

DrawPlots() 

# Get the number of time steps
nstates = TimeSliderGetNStates()

# Write inertia tensor to file in specified simulation steps
# Function via http://visitusers.org/index.php?title=Query_over_time
def DoublePrecisionQueryOverTime():
    f = open("201206_InertiaTensor.csv", "w")
    f.write("t,Ixx,Ixy,Ixz,Iyx,Iyy,Iyz,Izx,Izy,Izz\n")
    for tstate in range(0, nstates):
        TimeSliderSetState(tstate)
        # Get data file time
        Query("Time")
        t = GetQueryOutputValue()
        # Compute inertia tensor for data file
        Query("Moment of Inertia")
        tensor = GetQueryOutputValue()
     
        # Write time and inertia tensor to file
        str = "%4.3f,%25.15e,%25.15e,%25.15e,%25.15e,%25.15e,%25.15e,%25.15e,%25.15e,%25.15e\n" \
	          %(t,tensor[0],tensor[1],tensor[2],tensor[3],
                tensor[4],tensor[5],tensor[6],tensor[7],tensor[8])
        f.write(str)

        # Print every tenth state to console
        if tstate % 10 == 0:
            nt "3.2f" % t

    f.close()


DoublePrecisionQueryOverTime()
