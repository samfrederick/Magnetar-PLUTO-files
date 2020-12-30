#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 19:38:30 2020

@author: samfrederick
"""

def VolumeAvgTotalBfield():
    path = ('/Users/samfrederick/Documents/GitHub/Magnetar-PLUTO-files/'
            'Analysis/Bfield Analysis/201206_VolAvg_Total_Bfield.csv')
    DeleteAllPlots()
    f = open(path, 'r+')
    
    # Create scalar b-field plot 
    AddPlot("Pseudocolor","Scalar_B_Field_sqrd")
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

    DrawPlots()
    
    print 'Computing Volume Averaged Magnetic Field Strength'

    # Compute Volume Averaged Total B-field -------------
    #f.write("t,BTotal_v_avg\n")

    # Get the number of time steps
    nstates = TimeSliderGetNStates()


    i = 0
    for l in f:
        i += 1
    file_len = i - 1 # number of rows with data
        
    for tstate in range(file_len, nstates):
        
        TimeSliderSetState(tstate)
        Query("Time")
        t = GetQueryOutputValue()

        # Unit length used to avoid overflow
        Query("Volume")
        V = GetQueryOutputValue()

        Query("Weighted Variable Sum")
        Bhat = GetQueryOutputValue()
        B_avg = (Bhat/V)**(0.5)

        print t

        # Write time and volume averaged B magnitude to file
        str = ("%4.3f, %25.15e\n")\
              %(t, B_avg)
        f.write(str)
    f.close()
    print 'Computation completed'


def VolumeAvgPoloidalBfield():
    path = ('/Users/samfrederick/Documents/GitHub/Magnetar-PLUTO-files/'
            'Analysis/Bfield Analysis/201206_VolAvg_Poloidal_Bfield.csv')
    DeleteAllPlots()
    f = open(path, 'r+')
    
    # Create scalar b-field plot 
    AddPlot("Pseudocolor","Scalar_Poloidal_Field_sqrd")
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

    DrawPlots()
    
    print 'Computing Volume Averaged Magnetic Field Strength'

    # Compute Volume Averaged Poloidal B-field -------------
    #f.write("t,BPoloidal_v_avg\n")

    # Get the number of time steps
    nstates = TimeSliderGetNStates()


    i = 0
    for l in f:
        i += 1
    file_len = i - 1 # number of rows with data
        
    for tstate in range(file_len, nstates):
        

        TimeSliderSetState(tstate)
        Query("Time")
        t = GetQueryOutputValue()

        # Unit length used to avoid overflow
        Query("Volume")
        V = GetQueryOutputValue()

        Query("Weighted Variable Sum")
        Bhat = GetQueryOutputValue()
        B_avg = (Bhat/V)**(0.5)

        print t

        # Write time and volume averaged B magnitude to file
        str = ("%4.3f, %25.15e\n")\
              %(t, B_avg)
        f.write(str)
    f.close()
    print 'Computation completed'


def VolumeAvgToroidalBfield():
    path = ('/Users/samfrederick/Documents/GitHub/Magnetar-PLUTO-files/'
            'Analysis/Bfield Analysis/201206_VolAvg_Toroidal_Bfield.csv')
    DeleteAllPlots()
    f = open(path, 'r+')
    
    # Create scalar b-field plot 
    AddPlot("Pseudocolor","Scalar_Toroidal_Field_sqrd")
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

    DrawPlots()
    
    print 'Computing Volume Averaged Magnetic Field Strength'

    # Compute Volume Averaged Toroidal B-field -------------
    #f.write("t,BToroidal_v_avg\n")

    # Get the number of time steps
    nstates = TimeSliderGetNStates()


    i = 0
    for l in f:
        i += 1
    file_len = i - 1 # number of rows with data
        
    for tstate in range(file_len, nstates):
        

        TimeSliderSetState(tstate)
        Query("Time")
        t = GetQueryOutputValue()

        # Unit length used to avoid overflow
        Query("Volume")
        V = GetQueryOutputValue()

        Query("Weighted Variable Sum")
        Bhat = GetQueryOutputValue()
        B_avg = (Bhat/V)**(0.5)

        print t

        # Write time and volume averaged B magnitude to file
        str = ("%4.3f, %25.15e\n")\
              %(t, B_avg)
        f.write(str)
    f.close()
    print 'Computation completed'


VolumeAvgTotalBfield()
VolumeAvgPoloidalBfield()
VolumeAvgToroidalBfield()
