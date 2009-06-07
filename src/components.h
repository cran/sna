/*
######################################################################
#
# components.h
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 3/26/09
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
#include "utils.h"


/*INTERNAL ROUTINES---------------------------------------------------------*/


slelement *BFS(snaNet *g, int *n, int v, int transpose);

int numStrongComponents(snaNet *g, int *n);

slelement *strongComponentByVertex(snaNet *g, int *n, int v);

int *strongComponents(snaNet *g, int *n);

void strongComponentsRecurse(snaNet *g, int *n, int v, int *rindex, int *index, int *ccount, element *dfs);

int *undirComponents(snaNet *g);

void undirComponentsRecurse(snaNet *g,int v,int *memb);


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

void component_dist_R(double *g, double *pn, double *memb);


#endif
