/*
######################################################################
#
# nli.h
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 3/29/09
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
#
# Part of the R/sna package
#
# This file contains headers for nli.c.
#
######################################################################
*/
#ifndef NLI_H
#define NLI_H

/*DECLARATIONS/INCLUSIONS---------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "utils.h"
#include "geodist.h"

/*Definitions for various measures to be computed by betweenness_R (mostly
based on Brandes (2008)); note that some are not forms of betweenness, but
can be calculated using that routine.*/
#define BETSTANDARD     0         /*"Standard" form betweenness (a la Freeman)*/
#define BETWENDPTS      1                    /*Betweenness including endpoints*/
#define BETPROXIMALSRC  2                        /*Proximal source betweenness*/
#define BETPROXIMALTAR  3                        /*Proximal target betweenness*/
#define BETPROXIMALSUM  4                         /*Total proximal betweenness*/
#define BETLENSCALED    5                          /*Length-scaled betweenness*/
#define BETLINSCALED    6                        /*Linearly-scaled betweenness*/
#define BETSTRESS       7     /*Shimbel's stress centrality (not betweenness!)*/
#define BETLOAD       8/*Goh's load centrality (must be given transpose graph)*/
      
      
/*INTERNAL ROUTINES---------------------------------------------------------*/


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

SEXP betweenness_R(SEXP mat, SEXP sn, SEXP sm, SEXP smeasure, SEXP sprecomp, SEXP signoreevals, SEXP sgd, SEXP ssigma, SEXP spred);

void degree_R(double *g, int *pm, int *cmode, int *diag, int *igeval, double *d);

void evcent_R(double *mat, int *n, int *m, double *ev, double *tol, int *maxiter, int *checkna, int *ignoreeval);

void stresscent_R(double *g, double *pn, double *stress, double *gd, 
double *sigma);

void gilschmidt_R(double *mat, int *n, int *m, double *scores, int *normalize);

#endif
