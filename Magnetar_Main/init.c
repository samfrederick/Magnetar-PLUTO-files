/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful
  for problem configuration.
  It is automatically searched for by the makefile.

  \ author: Sam Frederick, Davidson College '19
  \ date created: 8-27-18, work in progress.
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*

FOREWORD (9-7-18):

Units defined in cgs. Parameter values are first computed using cgs constants
and are subsequently converted into 'dimensionless units' with considerations
for minimizing computational over/underflow and in supporting the robustness of
chosen numerical methods for computation.

The radius r, which will be referred to in the code as 'x1' is a value that
ranges from {0 < r < C} where C is the outermost boundary of the computational
domain specified in line 3 of the pluto.ini file. The domain of the star is
defined by {0 < r < 1}, where R = 1 is the radius of the star. There are numerous
points in the code where we make reference in both equations of state and
potential to r / R. We must be quite careful in our understanding of this value,
as we can not simply use the r specified prior as 'x1' in this ratio. Say we're
computing the value of this ratio at r = 8e5 cm  = 8000 km, where the radius of
the star is 1e6 cm =10000 km. We want this ratio to be 8000 km / 10000 km
= 0.8. NOT 0.8 / 100000 km as I mistakenly have previously specified. This means
that if we seek to use the domain {0 < x1 < 2}, we need to multiply this value
by 10000 to 'normalize' the ratio r / R such that (10000*x1)/ R = r / R.


/* ********************************************************************* */

/* PHYSICAL CONSTANTS */
#define G_CONST  6.67e-8 /*cm^3 g^-1 s^-2 */        /* Gravitational constant */
#define M_STAR  2.785e33 /* g             */        /* Mass of Star           */
#define RHO_C 2.2e15     /* g cm^-3       */        /* Core density of star   */
#define R 1.e6           /* cm            */        /* Radius of Star         */
#define K 4.25e4         /* cm^5 g^-1 s^-2*/        /* Coeff. for density eqn.*/
#define BMAX 1e17        /* gauss  := g^1/2 * cm^-1/2 * s                     */
#define Lambda 2.3619    /* First Eigenvalue for mixed Haskell Equations      */

/* HARD CODED VALUES */
#define RPOT 1.857595e20 /* Magnitude of the gravitational potential at r = R */
#define GPRSQ 1.4674e20  /* G x rho_c x R^2                                   */

/* COMPUTATIONAL CONSTANTS*/
#define h 1e-4           /* Finite difference step-size                       */
#define VACUUM 1e6       /* g cm^-3, Vacuum 'density', for comput. stability  */


/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*!
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical}
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
  /* */
  g_gamma = 2.0; /* Ratio of specific heats for an ideal gas, specifically
  defined for N = 1 polytrope [USING IDEAL EOS] */

  /* Global Initialization of state variables */
  v[RHO] = VACUUM / UNIT_DENSITY; /* Vacuum density surrounding star */
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;


  /* Comment out velocity perturbations, 2-10-19
  if ((x1< 1.0) && (x2 < (65*CONST_PI)/100) && (x2 > (35*CONST_PI)/100)){
    v[VX2] = ((3*cos(x2)*cos(x2)-1)*(2*cos(2*x3)))/(sin(x2));
    v[VX3] = cos(x2)*sin(x2)*cos(2*x3);
  }else{
    v[VX2] = 0.0;
    v[VX3] = 0.0;
  }
  */


  /* Velocity component normalization */
  v[VX1] = v[VX1] / UNIT_VELOCITY;
  v[VX2] = v[VX2] / UNIT_VELOCITY;
  v[VX3] = v[VX3] / UNIT_VELOCITY;

  /* Assign "vacuum pressure" exterior to star */
  v[PRS] = (K*VACUUM*VACUUM)/(UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);
  v[TRC] = 0.0; /* Tracer (passive scalar, Q) */


  /* Assign physical attributes for stellar interior */
  if ((x1 < 1.0) && (x1!= 0))
  {
    /* Assign B-field Component Values (Haskell et al. 2008) */
    v[BX1] = (2*A(x1)*cos(x2))/((x1*R)*(x1*R));
    v[BX2] = (-dA(x1)*sin(x2))/(x1*R*R);
    v[BX3] = (Lambda*CONST_PI*A(x1)*sin(x2))/(x1*R*R);

    /* Normalization */
    v[BX1] = v[BX1] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
    v[BX2] = v[BX2] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
    v[BX3] = v[BX3] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);

    /* Calcuate values for pressure and density using N = 1 polytrope EOS */
    v[RHO] = (RHO_C*sin(CONST_PI*x1))/(x1*CONST_PI) + VACUUM;
    v[PRS] = K*v[RHO]*v[RHO];

    /* Normalize values for density and pressure */
    v[RHO] = v[RHO] / UNIT_DENSITY; /* Converting to UNITLESS computational values */
    v[PRS] = v[PRS] / (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);


  }
  /* Assign physical attributes for stellar core */
  else if(x1 == 0){
     v[RHO] = RHO_C + VACUUM; /* Density at stellar core  */
     v[PRS] = K*RHO_C*RHO_C;  /* Pressure at stellar core */

     /* Normalize values for density and pressure */
     v[RHO] = v[RHO] / UNIT_DENSITY;
     v[PRS] = v[PRS] / (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);

     /* Assign B-field components at core */
     v[BX1] = -2*CONST_PI*CONST_PI*BMAX*cos(x2)*(2*CONST_PI*pow(Lambda,3)+
     (1-3*Lambda*Lambda)*sin(CONST_PI*Lambda)+CONST_PI*(3*Lambda*Lambda-1)*Lambda*cos(CONST_PI*Lambda));
     v[BX1] = v[BX1] / (3*pow((Lambda*Lambda-1),2)*(CONST_PI*Lambda*cos(CONST_PI*Lambda)-sin(CONST_PI*Lambda)));

     v[BX2] = 2*CONST_PI*CONST_PI*BMAX*sin(x2)*(2*CONST_PI*pow(Lambda,3)+
     (1-3*Lambda*Lambda)*sin(CONST_PI*Lambda)+CONST_PI*(3*Lambda*Lambda-1)*Lambda*cos(CONST_PI*Lambda));
     v[BX2]  = v[BX2]  / (3*pow((Lambda*Lambda-1),2)*(CONST_PI*Lambda*cos(CONST_PI*Lambda)-sin(CONST_PI*Lambda)));
     v[BX3] = 0;

     /* B-field Normalization */
     v[BX1] = v[BX1] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
     v[BX2] = v[BX2] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
     v[BX3] = v[BX3] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
  }
}
/* ********************************************************************* */
double A(double x1)
/*
Stream Function via Haskell et al. 2008 p. 540  Mixed Field Equations

The stream function is a mathematical tool appearing in the definition of the
field component equations.

To be specific, the stream function is defined as S(r, theta) = A(r)sin^2(theta).
Therefore, what we are computing is actually the stream function divided by
sin^2(theta). We make this choice because A(r), not S(r, theta) appers explicitly
in the expressions for the field component functions on pg. 540 of Haskell et al.

/* ********************************************************************* */
{
    double Aval;
    if (x1 != 0)
    {
      Aval =
      ((BMAX*R*R)/((Lambda*Lambda-1)*(Lambda*Lambda-1)*CONST_PI*x1))*
      (2*CONST_PI*((Lambda*CONST_PI*x1*cos(Lambda*CONST_PI*x1)-sin(Lambda*CONST_PI*x1))/
      (CONST_PI*Lambda*cos(CONST_PI*Lambda)-sin(CONST_PI*Lambda))) +
      ((1-Lambda*Lambda)*(CONST_PI*x1)*(CONST_PI*x1)-2)*sin(CONST_PI*x1) +
      2*CONST_PI*x1*cos(CONST_PI*x1));
    }
    else
    {
        Aval = 0;
    }
    return Aval;

}
/* ********************************************************************* */
double dA(double x1)
{
/*
Radial Derivative of Stream Function as defined above

/* ********************************************************************* */
double D_Aval;
if (x1 != 0)
  {
  D_Aval = (A(x1+h) - A(x1-h))/(2*h);  /* Central Difference Approx */
  }
else
  {
  D_Aval = 0;      /* Forward Differnce Approx */
  }


  return D_Aval;
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
// double simpsonzz (int no, double min, double max)
// /*!
// *  Perform simpsons approximation for integration of Izzint
// * Credit given to the following for publicly posting code:
// *								                                                       *
// *  taken from: "Projects in Computational Physics" by Landau and Paez  *
// *	       copyrighted by John Wiley and Sons, New York                  *
// *                                                                      *
// *  written by: students in PH465/565, Computational Physics,           *
// *	       at Oregon State University                                    *
// *              code copyrighted by RH Landau                           *
// *  supported by: US National Science Foundation, Northwest Alliance    *
// *                for Computational Science and Engineering (NACSE),    *
// *                US Department of Energy 	                             *
// *								                                                       *
// *  UNIX (DEC OSF, IBM AIX): cc integ.c -lm                             *
// *    			                                                             *
// */
// {
// int n;
// double interval, sum=0.0, x;
//
// interval = ((max-min))/(no-1);
//
// for (n=2; n<no; n+=2)           /* loop for odd points      */
//    {
//        x = interval * (n-1);
//        sum += 4 * Izzint(x);
//    }
//    for (n=3; n<no; n+=2)        /* loop for even points     */
//    {
//       x = interval * (n-1);
//       sum += 2 * Izzint(x);
//    }
//    sum +=  Izzint(min) + Izzint(max);	 	  /* add first and last value */
//    sum *= interval/3.;        	/* then multilpy by interval*/
//
//    return (sum);
//
// }
// /* ********************************************************************* */
// double simpsonxx (int no, double min, double max)
// /*!
// *  Perform simpsons approximation for integration of Ixxint
// * Credit given to the following for publicly posting code:
// *								                                                       *
// *  taken from: "Projects in Computational Physics" by Landau and Paez  *
// *	       copyrighted by John Wiley and Sons, New York                  *
// *                                                                      *
// *  written by: students in PH465/565, Computational Physics,           *
// *	       at Oregon State University                                    *
// *              code copyrighted by RH Landau                           *
// *  supported by: US National Science Foundation, Northwest Alliance    *
// *                for Computational Science and Engineering (NACSE),    *
// *                US Department of Energy 	                             *
// *								                                                       *
// *  UNIX (DEC OSF, IBM AIX): cc integ.c -lm                             *
// *    			                                                             *
// **************************************************************************/
// {
// int n;
// double interval, sum=0.0, x;
//
// interval = ((max-min))/(no-1);
//
// for (n=2; n<no; n+=2)           /* loop for odd points      */
//    {
//        x = interval * (n-1);
//        sum += 4 * Ixxint(x);
//    }
//    for (n=3; n<no; n+=2)        /* loop for even points     */
//    {
//       x = interval * (n-1);
//       sum += 2 * Ixxint(x);
//    }
//    sum +=  Ixxint(min) + Ixxint(max);	 	  /* add first and last value */
//    sum *= interval/3.;        	/* then multilpy by interval*/
//
//    return (sum);
//
// }
// /* ********************************************************************* */
// double Izzint(const Data *d, Grid *grid)
// /*!
//  *  Integrand for Izz moment of inertia
//  *
//  *
//  *********************************************************************** */
// {
//   int i,j,k;
//   double  *x1, *x2, *x3,rho,prs,dv;
//
//   x1 = grid->x[IDIR];
//   x2 = grid->x[JDIR];
//   x3 = grid->x[KDIR];
//
//   rho = d->Vc[RHO][k][j][i];
//
//   dv = (1/100.0)(CONST_PI/20.0)(CONST_PI/(40.0));
//
//   return ((rho*((x1[i]*sin(x2[j])*cos(x3[k])*(x1[i]*sin(x2[j])*sin(x3[k])))+
//   ((x1[i]*sin(x2[j])*sin(x3[k]))*(x1[i]*sin(x2[j])*sin(x3[k]))))*x1[i]*x1[i]*sin(x2[j]))*dv);
//
// }
// /* ********************************************************************* */
// double Ixxint(const Data *d, Grid *grid)
// /*!
//  *  Integrand for Ixx moment of inertia
//  *
//  *
//  *********************************************************************** */
// {
//   int i,j,k;
//   double  *x1, *x2, *x3,rho,prs,dv;
//
//   x1 = grid->x[IDIR];
//   x2 = grid->x[JDIR];
//   x3 = grid->x[KDIR];
//
//   rho = d->Vc[RHO][k][j][i];
//
//   dv = (1/100.0)(CONST_PI/20.0)(CONST_PI/(40.0));
//
//   return ((rho*((x1[i]*sin(x2[j])*sin(x3[k])*(x1[i]*sin(x2[j])*sin(x3[k])))+
//   ((x1[i]*cos(x2[j]))*(x1[i]*cos(x2[j]))))*x1[i]*x1[i]*sin(x2[j]))*dv);
//
// }
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*!
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures
 *
 *********************************************************************** */
{

int i,j,k;
double  *x1, *x2, *x3,rho,prs,Ixx,Izz,diff,Bt,Btsquare;
double *dr, *dt, *dp, dV;
double rgrid, tgrid, pgrid,Bphi;

x1 = grid->x[IDIR];
x2 = grid->x[JDIR];
x3 = grid->x[KDIR];

dr = grid->dx[IDIR];
dt = grid->dx[JDIR];
dp = grid->dx[KDIR];

rho = d->Vc[RHO][k][j][i];
prs = d->Vc[PRS][k][j][i];
Bphi = d->Vc[BX3][k][j][i];



rgrid = 100.0;
tgrid = 40.0;
pgrid = 20.0;

//Ixx = simpsonxx(100,0,1)
Ixx = 0.0;
//Izz = simpsonzz(100,0,1)
Izz = 0.0;
Bt = 0.0;

 DOM_LOOP(k,j,i){

     dV = x1[i]*x1[i]*sin(x2[j])*dr[i]*dt[j]*dp[k];

     if (rho > (VACUUM / UNIT_DENSITY) ){
//
//       //printf("Stuff is happening! %.5e\n",rho);
//


       Bt = Bt + (d->Vc[BX3][k][j][i]*d->Vc[BX3][k][j][i])*dV;
       Izz = Izz +  d->Vc[PRS][k][j][i]*(pow(x1[i]*sin(x2[j])*sin(x3[k]),2)+pow(x1[i]*cos(x2[j]),2))*dV;
       Ixx = Ixx +  d->Vc[PRS][k][j][i]*(pow(x1[i]*sin(x2[j])*cos(x3[k]),2)+pow(x1[i]*sin(x2[j])*sin(x3[k]),2))*dV;
//       fprintf(f,"Density at r =%f, theta = %f,phi = %f\n",x1[i],x2[j],x3[k]);
//
//       //fprintf(f,"%.5e\n",rho);
//       // if (f == NULL ){
//       //     printf("Error opening file!\n");
//       //     exit(1);
//       // }
//       // fprintf(f,"Radius: %12.6f Density: %12.6e  Pressure: %12.6e \n",x1[i],rho,prs);
//
//
 }
}

Btsquare = Bt;
diff = Izz-Ixx;


if (prank == 0){
  char fname[512];
  static double tpos = -1.0;
  FILE *f;

if (g_stepNumber==0){
   f = fopen("ANALYSIS.txt","w");
   fprintf(f,"%7s  %12s  %12s  %12s     %12s\n","t","dt","Izz","Ixx","Izz-Ixx");
}else{
  if (tpos < 0.0){
    char sline[512];
    f = fopen("ANALYSIS.txt","r");
    while (fgets(sline,512,f)) {}
    sscanf(sline, "%lf\n",&tpos);
    fclose(f);
  }
  f = fopen("ANALYSIS.txt","a");
}
if (g_time > tpos){
  fprintf(f,"%12.6e  %12.6e  %12.6e  %12.6e  %12.6e \n",g_time,g_dt,Izz,Ixx,diff);
}
fclose(f);
}

}
/* ********************************************************************* */
 // void BackgroundField (double x1, double x2, double x3, double *B0)
 // /*!
 //  * Define the component of a static, curl-free background
 //  * magnetic field.
 //  *
 //  * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 //  * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 //  * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 //  * \param [out] B0 array containing the vector components of the background
 //  *                 magnetic field
 //  *********************************************************************** */
 // {

 // }


/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*!
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and
 *                    staggered magnetic fields (d->Vs, when used) to
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END,
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

    if (side == 0) {
      TOT_LOOP(k,j,i){
        if ((x1[i] >= 1.0) /*&& (d->Vc[BX3][k][j][i]) != 0.0*/){
          d->Vc[BX3][k][j][i] = 0.0;
          d->Vc[BX2][k][j][i] = 0.0;
          d->Vc[BX1][k][j][i] = 0.0;
          d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;
        }
        if ((x1[i] > .98) && (x1[i] < 1.01)){
          d->Vc[BX2][k][j][i] =  (d->Vc[BX2][k][j][i]+ d->Vc[BX2][k][j][i-1])/2;
          d->Vc[BX1][k][j][i] =  (d->Vc[BX1][k][j][i]+ d->Vc[BX1][k][j][i-1])/2;
        }
        if (d->Vc[RHO][k][j][i] < (VACUUM) / UNIT_DENSITY){
            /* Replace negative values with vacuum density */
              d->Vc[RHO][k][j][i] = (VACUUM) / UNIT_DENSITY;
        }
        if (d->Vc[PRS][k][j][i] < (K*VACUUM*VACUUM)/(UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY)){
            /* Replace negative values with vacuum pressure */
            //printf("Low Pressure %f\n",d->Vc[PRS][k][j][i]);
            d->Vc[PRS][k][j][i] = (K*VACUUM*VACUUM)/(UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);
            //printf("Pressure corrected at %f\n",x1[i] );
        }
        if (x1[i] <= 1*(RMAX-RMIN)/(RGRID)){
          /* Constraint on core density and pressure; time independent values */
          d->Vc[RHO][k][j][i] = (RHO_C + VACUUM) / UNIT_DENSITY;
          d->Vc[PRS][k][j][i] = (K*RHO_C*RHO_C)/ (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);
          d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY; /* These values are TIME INDEPENDENT */
         }
        //if (d->Vc[RHO][k][j][i] > 2.5e15/UNIT_DENSITY){
        //   d->Vc[RHO][k][j][i] = 2.5e15/UNIT_DENSITY;
        // }
        if (x1[i] > 1.98 ){
          /* Set density/pressure fixed at outermost boundary of simulation */
            d->Vc[RHO][k][j][i] = (VACUUM) / UNIT_DENSITY;
            d->Vc[PRS][k][j][i] = (K*VACUUM*VACUUM)/ (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);
            d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY; /* These values are TIME INDEPENDENT */
        }
        // if (x2[j] > 1.56) && (x2[j] < 1.58){
        //   d->Vc[BX1][k][j][i] = 0.0;
        // }
      }
    }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){
        d->Vc[RHO][k][j][i] = RHO_C + VACUUM;/* Density at star core */
        d->Vc[PRS][k][j][i] = K*RHO_C*RHO_C;

        /* Normalize values for density and pressure */
        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i] / UNIT_DENSITY;
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][i] / (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);

        /* B-field values */
        d->Vc[BX1][k][j][i] = -2*CONST_PI*CONST_PI*BMAX*cos(x2[j])*(2*CONST_PI*pow(Lambda,3)+
        (1-3*Lambda*Lambda)*sin(CONST_PI*Lambda)+CONST_PI*(3*Lambda*Lambda-1)*Lambda*cos(CONST_PI*Lambda));
        d->Vc[BX1][k][j][i] = d->Vc[BX1][k][j][i] / (3*pow((Lambda*Lambda-1),2)*(CONST_PI*Lambda*cos(CONST_PI*Lambda)-sin(CONST_PI*Lambda)));

        d->Vc[BX2][k][j][i] = 2*CONST_PI*CONST_PI*BMAX*sin(x2[j])*(2*CONST_PI*pow(Lambda,3)+
        (1-3*Lambda*Lambda)*sin(CONST_PI*Lambda)+CONST_PI*(3*Lambda*Lambda-1)*Lambda*cos(CONST_PI*Lambda));
        d->Vc[BX2][k][j][i] = d->Vc[BX2][k][j][i] / (3*pow((Lambda*Lambda-1),2)*(CONST_PI*Lambda*cos(CONST_PI*Lambda)-sin(CONST_PI*Lambda)));

        d->Vc[BX3][k][j][i] = 0.0;

        /* Normalize values for Bfield */
        d->Vc[BX1][k][j][i] = d->Vc[BX1][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
        d->Vc[BX2][k][j][i] = d->Vc[BX2][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
        d->Vc[BX3][k][j][i] = d->Vc[BX3][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);


        }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){

        /* B-field values */
        d->Vc[BX1][k][j][i] = 0.0;
        d->Vc[BX2][k][j][i] = 0.0;
        d->Vc[BX3][k][j][i] = 0.0;

        /* Normalize values for Bfield */
        d->Vc[BX1][k][j][i] = d->Vc[BX1][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
        d->Vc[BX2][k][j][i] = d->Vc[BX2][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
        d->Vc[BX3][k][j][i] = d->Vc[BX3][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);

        /* Set and normalize density and pressure and vacuum values */
        d->Vc[RHO][k][j][i] = (VACUUM) / UNIT_DENSITY;
        d->Vc[PRS][k][j][i] = (K*VACUUM*VACUUM)/ (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);

      }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){
        if (x1[i] < 1.0){

          /* B-field values interior to star along theta axis*/
          d->Vc[BX1][k][j][i] = (-CONST_PI*CONST_PI*(Lambda*Lambda-1)*x1[i]*x1[i]-2)*sin(CONST_PI*x1[i])+
          ((2*CONST_PI)/(CONST_PI*Lambda*cos(CONST_PI*Lambda)-sin(CONST_PI*Lambda)))*
          (CONST_PI*Lambda*x1[i]*cos(CONST_PI*Lambda*x1[i])-sin(CONST_PI*Lambda*x1[i]))+
          2*CONST_PI*x1[i]*cos(CONST_PI*x1[i]);
          d->Vc[BX1][k][j][i] = d->Vc[BX1][k][j][i]*((2*BMAX)/(CONST_PI*pow((Lambda*Lambda-1),2)*pow(x1[i],3)));
          d->Vc[BX2][k][j][i] = 0.0; /* No theta bfield component at central (theta) axis*/
          d->Vc[BX3][k][j][i] = 0.0;

          /* Normalize values for Bfield */
          d->Vc[BX1][k][j][i] = d->Vc[BX1][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
          d->Vc[BX2][k][j][i] = d->Vc[BX2][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
          d->Vc[BX3][k][j][i] = d->Vc[BX3][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
        }
      else{

        /* B-field values exterior to star along theta axis*/
        d->Vc[BX1][k][j][i] = 0.0;
        d->Vc[BX2][k][j][i] = 0.0; /* No theta bfield component at central (theta) axis*/
        d->Vc[BX3][k][j][i] = 0.0;

        /* Normalize values for Bfield */
        d->Vc[BX1][k][j][i] = d->Vc[BX1][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
        d->Vc[BX2][k][j][i] = d->Vc[BX2][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
        d->Vc[BX3][k][j][i] = d->Vc[BX3][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
      }

      }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){
        if (x1[i] < 1.0){

          /* B-field values interior to star along theta axis*/
          d->Vc[BX1][k][j][i] = (-CONST_PI*CONST_PI*(Lambda*Lambda-1)*x1[i]*x1[i]-2)*sin(CONST_PI*x1[i])+
          ((2*CONST_PI)/(CONST_PI*Lambda*cos(CONST_PI*Lambda)-sin(CONST_PI*Lambda)))*
          (CONST_PI*Lambda*x1[i]*cos(CONST_PI*Lambda*x1[i])-sin(CONST_PI*Lambda*x1[i]))+
          2*CONST_PI*x1[i]*cos(CONST_PI*x1[i]);
          d->Vc[BX1][k][j][i] = d->Vc[BX1][k][j][i]*((-2*BMAX)/(CONST_PI*pow((Lambda*Lambda-1),2)*pow(x1[i],3)));
          d->Vc[BX2][k][j][i] = 0.0; /* No theta bfield component at central (theta) axis*/
          d->Vc[BX3][k][j][i] = 0.0;

          /* Normalize values for Bfield */
          d->Vc[BX1][k][j][i] = d->Vc[BX1][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
          d->Vc[BX2][k][j][i] = d->Vc[BX2][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
          d->Vc[BX3][k][j][i] = d->Vc[BX3][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
        }
      else{

        /* B-field values exterior to star along theta axis*/
        d->Vc[BX1][k][j][i] = 0.0;
        d->Vc[BX2][k][j][i] = 0.0; /* No theta bfield component at central (theta) axis*/
        d->Vc[BX3][k][j][i] = 0.0;

        /* Normalize values for Bfield */
        d->Vc[BX1][k][j][i] = d->Vc[BX1][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
        d->Vc[BX2][k][j][i] = d->Vc[BX2][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
        d->Vc[BX3][k][j][i] = d->Vc[BX3][k][j][i] / (sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
      }
      }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  // if (side == X3_BEG){  /* -- X3_BEG boundary -- */
  //   if (box->vpos == CENTER) {
  //     BOX_LOOP(box,k,j,i){  }
  //   }else if (box->vpos == X1FACE){
  //     BOX_LOOP(box,k,j,i){  }
  //   }else if (box->vpos == X2FACE){
  //     BOX_LOOP(box,k,j,i){  }
  //   }else if (box->vpos == X3FACE){
  //     BOX_LOOP(box,k,j,i){  }
  //   }
  // }
  //
  // if (side == X3_END){  /* -- X3_END boundary -- */
  //   if (box->vpos == CENTER) {
  //     BOX_LOOP(box,k,j,i){  }
  //   }else if (box->vpos == X1FACE){
  //     BOX_LOOP(box,k,j,i){  }
  //   }else if (box->vpos == X2FACE){
  //     BOX_LOOP(box,k,j,i){  }
  //   }else if (box->vpos == X3FACE){
  //     BOX_LOOP(box,k,j,i){  }
  //   }
  // }
}
//
// #if BODY_FORCE != NO
// /* ********************************************************************* */
// void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
// /*!
//  * Prescribe the acceleration vector as a function of the coordinates
//  * and the vector of primitive variables *v.
//  *
//  * \param [in] v  pointer to a cell-centered vector of primitive
//  *                variables
//  * \param [out] g acceleration vector
//  * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
//  * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
//  * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
//  *
//  *********************************************************************** */
// {
//   g[IDIR] = 0.0;
//   g[JDIR] = 0.0;
//   g[KDIR] = 0.0;
// }
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  /* similar normalization of computaitonal radius 'r' to be in proportion with R*/
   double phi;

   if ((x1 <= 1) && (x1 != 0)) { /* Potential interior to star except r = 0 */
     phi = (-4*GPRSQ*sin(CONST_PI*x1))/(CONST_PI*CONST_PI*x1) - RPOT ; /* Factor of R has been taken out for first term
     denominator because of normalization*/
     phi = phi/(UNIT_VELOCITY*UNIT_VELOCITY);
   }
   if (x1 > 1){ /* Potential exterior to star */
     phi = -(G_CONST*M_STAR) / (R*x1);
     phi = phi/(UNIT_VELOCITY*UNIT_VELOCITY);
   }
   if (x1 == 0){ /* Potential at r = 0 */
     phi = 4*G_CONST*RHO_C*(-(R*R)/(CONST_PI)-(M_STAR)/(4*R*RHO_C));
     phi = phi/(UNIT_VELOCITY*UNIT_VELOCITY);
   }

  return phi;
}
