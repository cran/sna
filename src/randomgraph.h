/*
######################################################################
#
# randomgraph.h
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 12/27/04
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains headers for randomgraph.c.
#
######################################################################
*/
#ifndef RANDOMGRAPH_H
#define RANDOMGRAPH_H

/*DECLARATIONS/INCLUSIONS---------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>


/*INTERNAL ROUTINES---------------------------------------------------------*/


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

void udrewire_R(double *g, double *pn, double *pnv, double *pp);
void wsrewire_R(double *gi, double *go, double *pn, double *pnv, double *pp);


#endif
