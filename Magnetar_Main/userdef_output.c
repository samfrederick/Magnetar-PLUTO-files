#include "pluto.h"

/*
Author: Sam Frederick
Date: 9-13-19

Created to assign b-field component values to the user-defined 
variables BR, BT, and BP. These can then be visualized along with
state variables in VisIt. 

This file should replace the default userdef_output.c under the 
primary PLUTO folder (PLUTO/Src). 

*/

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
  double ***BR, ***BT, ***BP;

  BR = GetUserVar("BR");
  BT = GetUserVar("BT");
  BP = GetUserVar("BP");
	
  

  DOM_LOOP(k,j,i){
    BR[k][j][i] = d->Vc[BX1][k][j][i]*(sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
    BT[k][j][i] = d->Vc[BX2][k][j][i]*(sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
    BP[k][j][i] = d->Vc[BX3][k][j][i]*(sqrt(4*CONST_PI*UNIT_DENSITY)*UNIT_VELOCITY);
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





