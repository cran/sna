/*
######################################################################
#
# geodist.h
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 11/21/04
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains headers for geodist.c.
#
######################################################################
*/
#ifndef GEODIST_H
#define GEODIST_H

/*DECLARATIONS/INCLUSIONS---------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>


/*INTERNAL ROUTINES---------------------------------------------------------*/


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

void geodist_R(double *g, double *pn, double *gd, double *sigma);


#endif
