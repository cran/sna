/*
######################################################################
#
# nli.c
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 11/30/04
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines related to the computation of node level
# indices (NLIs).
#
######################################################################
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "nli.h"


void betweenness_R(double *g, double *pn, double *bet, double *gd, 
double *sigma)
/*
Compute betweenness stats for the graph in g.  Geodesic distances are assumed to
have been stored in gd, and the path counts in sigma (both being nxn matrices).
It is also assumed that bet has been initialized to 0.
*/
{
  long int n,i,j,k;

  /*Set up stuff*/
  n=*pn;
  /*Cycle through each triad, accumulating paths*/
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      for(k=0;k<n;k++){
        if((j!=i)&&(k!=i)&&(j!=k)&&(gd[j+k*n]>=gd[j+i*n]+gd[i+k*n])&& (sigma[j+k*n]>0.0))
          bet[i]+=sigma[j+i*n]*sigma[i+k*n]/sigma[j+k*n];
      }
    }
  }
}


void stresscent_R(double *g, double *pn, double *stress, double *gd, 
double *sigma)
/*
Compute stress centrality for the graph in g.  Geodesic distances are assumed to
have been stored in gd, and the path counts in sigma (both being nxn matrices).
It is also assumed that stress has been initialized to 0.
*/
{
  long int n,i,j,k;

  /*Set up stuff*/
  n=*pn;
  /*Cycle through each triad, accumulating paths*/
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      for(k=0;k<n;k++){
        if((j!=i)&&(k!=i)&&(j!=k)&&(gd[j+k*n]>=gd[j+i*n]+gd[i+k*n]))
          stress[i]+=sigma[j+i*n]*sigma[i+k*n];
      }
    }
  }
}

