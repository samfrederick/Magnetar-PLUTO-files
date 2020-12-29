"""
Author: Sam Frederick
Date: 10-12-18

This script analyzes radial grid points surrounding a specified 
asymptote resulting from the condition number for a particular 
function with respect to radius in our computational domain. A 
tolerance is specified to see whether there exist any gridpoints 
significantly close to the asymptote which may result in poorly 
behaved computation around such points. The tolerance is determined 
by the threshold magnitude of the condition number below which we 
deem "well behaved" and vice versa. Traditionally, this threshold 
magnitude is of the order 10^8 or equivalent sufficiently large values. 

Here, we intentionally set the tolerance fairly low (1e-2) to show
even for such low tolerance, we do not expect poorly behaved points 
within our computational domain.

"""

import numpy as np

"""
This script is currently configured for an asymptote arising in the  
condition number for the theta component of the B-field.  
"""
rgridmin = 0
rgridmax = 2
divisions = 70
asymp = 0.765052
tolerance = 1e-2
count = 1


vals = np.linspace(rgridmin,rgridmax,divisions)
# print (vals)

for point in vals:
    if abs(point-asymp) < tolerance:
        if count == True:
            print ("The following values for r may blow up:") 
            count = False
        print ("r = %(point)f" % {"point": point }) 
        
if count == True:
    print ("No points are sufficiently close to the asymptote; " \
    "the computational domain is well behaved.")