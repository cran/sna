/*
######################################################################
#
# gli.c
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 11/26/04
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines related to the computation of graph 
# level indices (GLIs).
#
######################################################################
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "gli.h"

void lubness_con_R(double *g, double *pn, int *r, double *viol)
/*Calculate lubness violations for a weakly connected graph, g, having
strong reachability matrix r (both must be nxn).  This is obviously intended
to be called on a component by component basis.*/
{
  long int i,j,k,l,n,nub,*ub;
  char lub,notlub;

  /*Set things up*/
  *viol=0.0;
  n=*pn;
  ub=(long int *)R_alloc(n,sizeof(long int));
  /*Accumulate LUBness violations*/
  if(n>2){  /*No violations unless n>2, given weak connectivity*/
    for(i=0;i<n;i++)      /*Walk vertex pairs*/
      for(j=i+1;j<n;j++){
        /*Accumulate upper bounds*/
        nub=0;
        for(k=0;k<n;k++)
          if(r[k+i*n]&&r[k+j*n])
            ub[nub++]=k;
        /*Seek a least upper bound*/
        lub=0;
        for(k=0;(k<nub)&&(!lub);k++){
          notlub=0;
          for(l=0;(l<nub)&&(!notlub);l++)
            if(!r[ub[k]+ub[l]*n])
              notlub++;
          if(!notlub)
            lub++;
        }
        /*Aggregate the violation count*/
        if(!lub)
          (*viol)++;
      }
  }  
}
