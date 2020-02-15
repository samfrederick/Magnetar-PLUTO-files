

DeleteAllPlots()

print "Start\n---------------------"

# Create density plot with spherical clip operator
AddPlot("Pseudocolor","rho")
p = PseudocolorAttributes()
p.colorTableName = "hot"
SetPlotOptions(p)

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

# Write total mass to file in specified simulation steps
# Function via http://visitusers.org/index.php?title=Query_over_time
def DoublePrecisionQueryOverTime():
   f = open("Mass_Data.txt", "w")
   f.write("              t                       Mass                   \n")
   for i in range(1443): # For full analysis, set range to be TimeSliderGetNStates()
     SetTimeSliderState(i*10) # Select data files in multiples of 10
     Query("Time")
     t = GetQueryOutputValue()
     Query("Weighted Variable Sum")
     t2 = GetQueryOutputValue()
     str = "%25.15e %25.15e\n" %(t, t2) # Print timestamp, total mass
     f.write(str)
   f.close()


DoublePrecisionQueryOverTime()
