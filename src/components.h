/*
######################################################################
#
# components.h
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 7/19/16
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
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
#include <Rdefines.h>
#include <Rinternals.h>
#include "utils.h"


/*INTERNAL ROUTINES---------------------------------------------------------*/


slelement *BFS(snaNet *g, int *n, int v, int transpose);

element *BFS_unord(snaNet *g, int *n, int v, int transpose);

int numStrongComponents(snaNet *g, int *n);

slelement *strongComponentByVertex(snaNet *g, int *n, int v);

int *strongComponents(snaNet *g, int *n);

void strongComponentsRecurse(snaNet *g, int *n, int v, int *rindex, int *index, int *ccount, element *dfs);

int *undirComponents(snaNet *g);

void undirComponentsRecurse(snaNet *g,int v,int *memb);

void undirComponentsNoRecurse(snaNet *g, int *memb);


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

void component_dist_R(double *g, double *pn, double *memb);

void compsizes_R(double *mat, int *n, int *m, int *csizes);

SEXP reachability_R(SEXP mat, SEXP sn, SEXP sm);

void undirComponents_R(double *mat, int *n, int *m, int *memb);

#endif
