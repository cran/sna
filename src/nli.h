/*
######################################################################
#
# nli.h
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 11/21/04
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains headers for nli.c.
#
######################################################################
*/
#ifndef NLI_H
#define NLI_H

/*DECLARATIONS/INCLUSIONS---------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>


/*INTERNAL ROUTINES---------------------------------------------------------*/


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

void betweenness_R(double *g, double *pn, double *bet, double *gd, 
double *sigma);

void stresscent_R(double *g, double *pn, double *stress, double *gd, 
double *sigma);


#endif
