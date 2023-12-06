/*
######################################################################
#
# geodist.h
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 4/26/09
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
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
#include <Rdefines.h>
#include <Rinternals.h>
#include "utils.h"


/*INTERNAL ROUTINES---------------------------------------------------------*/

void spsp(int ego, snaNet *g, double *gd, double *sigma, element **pred, int *npred, int checkna);

void spsp_val(int ego, snaNet *g, double *gd, double *sigma, element **pred, int *npred, int checkna);


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

void geodist_adj_R(double *g, double *pn, double *gd, double *sigma);

SEXP geodist_R(SEXP mat, SEXP sn, SEXP m, SEXP scheckna, SEXP scalcsig, SEXP scalcpred);

SEXP geodist_val_R(SEXP mat, SEXP sn, SEXP sm, SEXP scheckna, SEXP scalcsig, SEXP scalcpred);

void maxflow_EK_R(double *g,int *pn,int *psource,int *psink,double *flow);

#endif
