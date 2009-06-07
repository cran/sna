/*
######################################################################
#
# gli.h
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 1/18/09
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains headers for gli.c.
#
######################################################################
*/

#ifndef GLI_H
#define GLI_H

/*DECLARATIONS/INCLUSIONS---------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "utils.h"


/*INTERNAL ROUTINES---------------------------------------------------------*/


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

void brokerage_R(double *g, int *pn, int *pm, int *cl, double *brok);

void lubness_con_R(double *g, double *pn, int *r, double *viol);


#endif
