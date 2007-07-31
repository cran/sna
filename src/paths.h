/*
######################################################################
#
# paths.h
#
# copyright (c) 2007, Carter T. Butts <buttsc@uci.edu>
# Last Modified 4/27/07
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains headers for paths.c.
#
######################################################################
*/
#ifndef PATHS_H
#define PATHS_H

/*DECLARATIONS/INCLUSIONS---------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>


/*INTERNAL ROUTINES---------------------------------------------------------*/

void edgewisePathRecurse(snaNet *g, int src, int dest, int curnode, int  *availnodes, int availcount, int *usednodes, int curlen, double *count, double *cpcount, double *dpcount, int maxlen, int directed, int byvertex, int copaths, int dyadpaths);

void edgewiseCycleCensus(snaNet *g, int src, int dest, double *count, double *cccount, int maxlen, int directed, int byvertex, int cocycles);

void dyadPathCensus(snaNet *g, int src, int dest, double *count, double *cpcount, double *dpcount, int maxlen, int directed, int byvertex, int copaths, int dyadpaths);


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

void cycleCensus_R(int *g, int *pn, double *count, double *cccount, int *pmaxlen, int *pdirected, int *pbyvertex, int *pcocycles);

void pathCensus_R(double *g, int *pn, double *count, double *cpcount, double *dpcount, int *pmaxlen, int *pdirected, int *pbyvertex, int *pcopaths, int *pdyadpaths);

#endif
