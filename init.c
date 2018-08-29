/*
Author: Sam Frederick
Created: 8-27-18

Simulates a disk modeling a star with radius 10 centered around the origin in a
2D cylindrical domain.

 */
#include "pluto.h"

/* ******************************************************************** */
/*!
 * A Note on Constants:
 * These constants are defined relative to the 'UNIT_DENSITY','UNIT_LENGTH',
 * 'UNIT_VELOCITY' parameters in definitions.h. They are tentative, and will
 * very likely be updated to reflect accurate values. Changing their value
 * causes drastic variation in the stability of the resulting model. This
 * emphasizes the importance of determining accurate values for these parameters.
/* ******************************************************************** */
#define RHO_C  1 /* Density at center of star */
#define K  0.935 /* Constant for polytropic pressure profile */
#define G  6e-8 /* Gravitational constant */
#define M_STAR  1.e5 /* Mass of Star */
#define R 10.0 /*Radius of Star*/

/* ******************************************************************* */
void Init (double *v, double x1, double x2, double x3)
{
    double r;

    v[VX1] = v[VX2] = 0.0; /* initilize vx and vy to zero */
    if ((x1 < 1.0) && (x1 != 0)){ /* set pressure and density inside disk */
      v[RHO] = (RHO_C*sin((CONST_PI*x1)/1.0)*1.0)/(x1*CONST_PI); /* Density for n = 1 polytrope*/
      v[PRS] = K*v[RHO]*v[RHO]; /* Pressure for n = 1 polytrope*/
    }else if(x1 = 0){ /* Density at center, this may be causing errors in simulation */
      v[RHO] = RHO_C; /* Density at star core */
      v[PRS] = K*RHO_C*RHO_C;
    }else{  /* exterior to disk set pressure and density to unit vals */
      v[RHO] = 1e-5;
      v[PRS] = 1e-5;
    }
}



/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*!
 * Assign initial condition by looping over the computatiohttps://photos.google.com/nal domain.
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
/*#if (BODY_FORCE & POTENTIAL)*/
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 * A graviational poential is assigned to the problem. This potential comes
 * from the work of Kuhn's Thesis, with a note to a minor typo to Kuhn's
 * potential at r = 0 (referred to as phi nought).
 *************************************************************************** */
{
  double phi;

  if ((x1 < 1) && (x1 != 0)) { /* Potential interior to star except r = 0 */
    phi = 4*G*RHO_C*(-(R*R*R*sin((CONST_PI*x1)/R))/(CONST_PI*CONST_PI*x1)-(M_STAR)/(4*R*RHO_C));
  }
  if (x1 >= 1){ /* Potential exterior to star */
    phi = -(G*M_STAR) / x1;
  }
  if (x1 = 0){ /* Potential at r = 0 */
    phi = 4*G*RHO_C*(-(R*R)/(CONST_PI)-(M_STAR)/(4*R*RHO_C));
  }
  return phi;
}
