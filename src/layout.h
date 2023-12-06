/*
######################################################################
#
# layout.h
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 11/21/10
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
#
# Part of the R/sna package
#
# This file contains headers for layout.c.
#
######################################################################
*/
#ifndef LAYOUT_H
#define LAYOUT_H

/*DECLARATIONS/INCLUSIONS---------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include "utils.h"

/*Macro: is a within the [b,c] interval?*/
#define ININT(a,b,c) ((a-b)*(a-c)<=0.0)

/*INTERNAL ROUTINES---------------------------------------------------------*/

double angdist(double a, double b, double ilen);

double poldist(double ra,double ta,double rb,double tb);

double pollinedist(double ra,double ta,double rb,double tb, double rc, double tc);

int poledgecross(double ra, double ta, double rb, double tb, double rc, double tc, double rd, double td);


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

void gplot_layout_target_R(int *d, double *pn, int *pniter, double *elen, double *radii, int *core, double *pdisconst, double *pcrossconst, double *prepconst, double *pminpdis, double *pinitemp, double *pcoolexp, double *pmaxdelta, double *theta);

void gplot_layout_fruchtermanreingold_R(double *d, double *pn, double *pm, 
int *pniter, double *pmaxdelta, double *pvolume, double *pcoolexp, double 
*prepulserad, int *pncell, double *pcjit, double *pcppr, double *pcpcr, double
*pcccr, double *x, double *y);

void gplot_layout_fruchtermanreingold_old_R(double *d, int *pn, int *pm, int *pniter, double *pmaxdelta, double *pvolume, double *pcoolexp, double *prepulserad, double *x, double *y);  /*Deprecated code, to be removed*/

void gplot_layout_kamadakawai_R(int *pn, int *pniter, double *elen, double *pinitemp, double *pcoolexp, double *pkkconst, double *psigma, double *x, double *y);

void gplot3d_layout_fruchtermanreingold_R(double *d, int *pn, int *pm, int *pniter, double *pmaxdelta, double *pvolume, double *pcoolexp, double *prepulserad, double *x, double *y, double *z);

void gplot3d_layout_kamadakawai_R(double *pn, int *pniter, double *elen, double *pinitemp, double *pcoolexp, double *pkkconst, double *psigma, double *x, double *y, double *z);

#endif
