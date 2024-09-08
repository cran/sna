/*
######################################################################
#
# randomgraph.h
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 8/01/24
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
#
# Part of the R/sna package
#
# This file contains headers for randomgraph.c.
#
######################################################################
*/
#ifndef RANDOMGRAPH_H
#define RANDOMGRAPH_H

/*DECLARATIONS/INCLUSIONS---------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include "utils.h"

/*Homogeneity classes for rgbern*/
#define BERNHOM 0               /*Homogeneous throughout the adjacency matrix*/
#define BERNROW 1               /*Homogeneous within rows*/
#define BERNCOL 2               /*Homogeneous within columns*/
#define BERNHET 3               /*Heteogeneous*/


/*INTERNAL ROUTINES---------------------------------------------------------*/


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

void bn_cftp_R(int *g, int *pn, double *pi, double *sigma, double *rho, double *d, int *pmaxiter, int *sibdichot);

void bn_mcmc_R(int *g, double *pn, double *pdraws, double *pburn, int *pthin, double *pi, double *sigma, double *rho, double *d, double *delta, double *epsilon, int *sibdichot, double *maxedge);

SEXP rgbern_R(SEXP sn, SEXP stp, SEXP sdirected, SEXP sloops, SEXP spmode);

void udrewire_R(double *g, double *pn, double *pnv, double *pp);

void wsrewire_R(double *gi, double *go, double *pn, double *pnv, double *pp);


#endif
