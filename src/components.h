/*
######################################################################
#
# components.h
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 11/26/04
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains headers for components.c.
#
######################################################################
*/
#ifndef COMPONENTS_H
#define COMPONENTS_H

/*DECLARATIONS/INCLUSIONS---------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>


/*INTERNAL ROUTINES---------------------------------------------------------*/


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

void component_dist_R(double *g, double *pn, double *memb);


#endif
