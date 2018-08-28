/*
Author: Sam Frederick
Date: 8-27-18

Simulates a basic disk with radius q1 centered around the origin in a 2D
Cartesian domain. Credit goes to the PLUTO User Manual, Chapter 5 Section 1,
Example #1 for listing these initial conditions for spatial coordinates.


 */
#include "pluto.h"

/* ******************************************************************* */
void Init (double *v, double x1, double x2, double x3)
{
    double r;

    r = sqrt(x1*x1 + x2*x2); /* magnitude of radius */
    v[VX1] = v[VX2] = 0.0; /* initilize vx and vy to zero */
    if (r < 1.0){ /* interior to disk set pressure and density uniform vals */
      v[RHO] = 10.0;
      v[PRS] = 30.0;
    }else{  /* exterior to disk set pressure and density to unit vals */
      v[RHO] = 1.0;
      v[PRS] = 1.0;
    }
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*!
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}
