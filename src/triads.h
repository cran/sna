/*
######################################################################
#
# triads.h
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 11/13/04
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains headers for triads.c.
#
######################################################################
*/
#ifndef TRIADS_H
#define TRIADS_H

/*DECLARATIONS/INCLUSIONS---------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>


/*INTERNAL ROUTINES---------------------------------------------------------*/
int triad_classify(int *g, int gn, int i, int j, int k);


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

void triad_census_R(int *g, int *n, double *t);
void triad_classify_R(int *g, int *tt);


#endif
