/*
Author: Sam Frederick
Date: 8-27-18

Simulates a basic disk with radius q1 centered around the origin in a 2D
Cartesian domain. Credit goes to the PLUTO User Manual, Chapter 5 Section 1,
Example #1 for listing these initial conditions for spatial coordinates.


 */
#include "pluto.h"

#define RHO_C  10.0 /* Density at center of star */
#define K  5.0 /* Constant for polytropic pressure profile */

/* ******************************************************************* */
void Init (double *v, double x1, double x2, double x3)
{
    double r;

    /*r = sqrt((x1-0.5)*(x1-0.5)+(x2-0.5)*(x2-0.5)); /* magnitude of radius */
    v[VX1] = v[VX2] = 0.0; /* initilize vx and vy to zero */
    if ((x1 < 1.0) && (x1 != 0)){ /* interior to disk set pressure and density uniform vals */
      v[RHO] = (RHO_C*sin((CONST_PI*x1)/1.0)*1.0)/(x1*CONST_PI); /* Density for n = 1 polytrope*/
      v[PRS] = K*v[RHO]*v[RHO]; /* Pressure for n = 1 polytrope*/
    }else if(x1 = 0){
      v[RHO] = RHO_C; /* Density at star core */
      v[PRS] = K*v[RHO]*v[RHO];
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



/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*
 *
 *
 *********************************************************************** */
{
}



/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*
 *
 *********************************************************************** */
{
}
