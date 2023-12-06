/*
######################################################################
#
# likelihood.h
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 3/10/05
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
#
# Part of the R/sna package
#
# This file contains headers for likelihood.c.
#
######################################################################
*/
#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

/*DECLARATIONS/INCLUSIONS---------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>


/*INTERNAL ROUTINES---------------------------------------------------------*/

double bn_lpka(long int k,double pi, double sigma, double rho, double d);
double bn_lpkm(long int k,double pi, double sigma, double rho, double d);
double bn_lpkn(long int k,double pi, double sigma, double rho, double d);
double bn_lpt(int xy, int yx, int yz, int zy, int xz, int zx, long int kxy, long int kyz, long int kxz, double pi, double sigma, double rho, double d);
double bn_lpt_M(long int m,double pi,double sigma,double rho,double d);
double bn_lpt_a(long int m,double pi,double sigma,double rho,double d);
double bn_lpt_N(long int m,double pi,double sigma,double rho,double d);
double bn_lpt_Mp1(long int m,double pi,double sigma,double rho,double d);
double bn_lpt_ap1(long int m,double pi,double sigma,double rho,double d);
double bn_lpt_Np1(long int m,double pi,double sigma,double rho,double d);
double bn_lpt_M1(double pi,double sigma,double rho,double d);
double bn_lpt_a1(double pi,double sigma,double rho,double d);
double bn_lpt_N1(double pi,double sigma,double rho,double d);
double bn_lpt_Sr(double pi,double sigma,double rho,double d);
double bn_lpt_1mSr(double pi,double sigma,double rho,double d);


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

void bn_dyadstats_R(int *g, double *pn, double *stats);
void bn_triadstats_R(int *g, double *pn, double *stats);
void bn_lpl_dyad_R(double *stats, double *psr, double *pi, double *sigma, double *rho, double *d, double *lpl);
void bn_lpl_triad_R(int *g, double *stats, double *pn, double *pi, double *sigma, double *rho, double *d, double *lpl);
void bn_ptriad_R(double *pi, double *sigma, double *rho, double *d, double *pt);

#endif
