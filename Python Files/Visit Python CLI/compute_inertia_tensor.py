"""
Author: Sam Frederick
Date:  December 2020

VisIt CLI Python command
------------------------

Compute inertia tensor component values for stellar 
region defined by r <=R_star=1.0 on each timestep.

"""
AddWindow()

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


def InertiaTensorQuery(filename):
    """
    Write inertia tensor to file in specified simulation
    steps.

    Function modified via:
    http://visitusers.org/index.php?title=Query_over_time
    """
    f = open(filename, "w")
    f.write("t,Ixx,Ixy,Ixz,Iyx,Iyy,Iyz,Izx,Izy,Izz\n")

    for tstate in range(0, nstates):
     
        # Get data file time
        TimeSliderSetState(tstate)
        Query("Time")
        t = GetQueryOutputValue()

        # Compute inertia tensor for data file
        Query("Moment of Inertia")
        tsr = GetQueryOutputValue()  # tensor MOI values
     
        # Write time and inertia tensor to file
        str = "%4.3f,%25.15e,%25.15e,%25.15e,%25.15e,%25.15e,%25.15e,%25.15e,%25.15e,%25.15e\n" \
	          %(t,tsr[0],tsr[1],tsr[2],tsr[3],tsr[4],tsr[5],tsr[6],tsr[7],tsr[8])
        f.write(str)
     
        # Print every tenth state to console
        if tstate % 10 == 0:
            print "t=%3.2f"%(t)

    f.close()


InertiaTensorQuery(filename="201206_InertiaTensor.csv")
