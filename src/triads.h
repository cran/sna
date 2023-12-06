/*
######################################################################
#
# triads.h
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 6/26/11
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
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
#include "utils.h"


/*INTERNAL ROUTINES---------------------------------------------------------*/

int triad_classify(int *g, int gn, int i, int j, int k, int gm);

int triad_classify_el(snaNet *g, int i, int j, int k, int gm, int checkmissing);


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

void transitivity_R(double *mat, int *n, int *m, double *t, int *meas, int *checkna);

void triad_census_R(double *g, int *n, int *m, double *t, int *gm, int *checkna);

void triad_classify_R(int *g, int *tt, int *gm);


#endif
