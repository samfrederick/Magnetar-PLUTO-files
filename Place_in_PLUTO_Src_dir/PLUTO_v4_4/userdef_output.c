#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k;
  double ***Mr, ***Mp, ***Mt;
  double ***p,***rho;
  double ***bx1,***bx2,***bx3;
  //double ***Temp;

  //Temp = GetUserVar("T");
  
  #if PHYSICS == MHD
  /*
  Mr = GetUserVar("Mr");
  Mp = GetUserVar("Mp");
  Mt = GetUserVar("Mt");

  rho = d->Vc[RHO];
  p   = d->Vc[PRS];
  bx1 = d->Vc[BX1];
  bx2 = d->Vc[BX2];
  bx3 = d->Vc[BX3];

  
  DOM_LOOP(k,j,i){
  }
  //Temp[k][j][i] = p[k][j][i] / rho[k][j][i];
  Mr[k][j][i] = bx1[k][j][i]*(sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
  Mt[k][j][i] = bx2[k][j][i]*(sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
  Mp[k][j][i] = bx3[k][j][i]*(sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
  */
  #endif

  double ***gravpot;
  double r; 
  double *x1, *x2, *x3; 

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];
  
  //i = g_i; j = g_j; k = g_k; 

 #if BODY_FORCE != NO
 gravpot = GetUserVar("gpot");
  
 DOM_LOOP(k,j,i){
    //r = (i - 1.0) / RGRID; // radial grid integer divided by number of grid points
    gravpot[k][j][i] = gpot(x1[i], x2[j], x3[k])*(UNIT_VELOCITY*UNIT_VELOCITY); 
  }
 
 #endif 

  
  double rho_r, prs_r, h;
  double del_prs, del_phi, ***hstatic_cond;
  double rho_pls_dr, prs_pls_dr;
  double rho_mns_dr, prs_mns_dr;

  hstatic_cond = GetUserVar("hstatic_condition");
  h = 1e-2;

  DOM_LOOP(k,j,i){
		
	  //r = (i - 1.0) / 100.0; // radial grid integer divided by number of grid points
	  rho_r = d->Vc[RHO][k][j][i];
	  prs_r = d->Vc[PRS][k][j][i];	


	  prs_pls_dr = d->Vc[PRS][k][j][i+1];
	  prs_mns_dr = d->Vc[PRS][k][j][i-1];	

	  if (i != 0)
	  {
	     
	     del_prs = (prs_pls_dr - prs_mns_dr) / (2*h);
	     del_phi = (gpot(x1[i], x2[j], x3[k]) - gpot(x1[i-1], x2[j-1], x3[k-1])) / (2*h);
	     hstatic_cond[k][j][i] = (del_prs / rho_r) + del_phi; 
	  }
	else
	  {
	     hstatic_cond[k][j][i]  = 0;  
	  }
	  
  }

  
}
/* ************************************************************* */
void ChangeOutputVar ()
/*
 *
 *
 *************************************************************** */
{
  Image *image;

#ifdef PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}
