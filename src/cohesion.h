/*
######################################################################
#
# cohesion.h
#
# copyright (c) 2007, Carter T. Butts <buttsc@uci.edu>
# Last Modified 5/1/09
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
#
# Part of the R/sna package
#
# This file contains headers for cohesion.c.
#
######################################################################
*/
#ifndef COHESION_H
#define COHESION_H

/*DECLARATIONS/INCLUSIONS---------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "utils.h"
#include "components.h"


/*INTERNAL ROUTINES---------------------------------------------------------*/

void bicomponentRecurse(snaNet *g, element *complist, element *estack, int *parent, int *num, int *back, int *dfn, int v);

slelement *cliqueFirstChild(snaNet *g, slelement *cl);

void cliqueRecurse(snaNet *g, slelement *k, int parind, element **clist, double *ccount, int *compmemb);

void cutpointUndirRecurse(snaNet *g, int *cpstatus, int *minvis, int *visdep, int depth, int v, int src);


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

SEXP bicomponents_R(SEXP net, SEXP sn, SEXP sm);

SEXP cliques_R(SEXP net, SEXP sn, SEXP sm, SEXP stabulatebyvert, SEXP scomembership, SEXP senumerate);

void cutpointsDir_R(double *mat, int *n, int *m, int *cpstatus);

void cutpointsUndir_R(double *mat, int *n, int *m, int *cpstatus);

void kcores_R(double *mat, int *n, int *m, double *corevec, int *dtype, int *pdiag, int *pigeval);

#endif
