/*
######################################################################
#
# triads.c
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 4/24/05
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines related to the classification and
# counting of triads.
#
######################################################################
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "triads.h"

int triad_classify(int *g, int gn, int i, int j, int k, int gm)
/*
If gm=1, compute a Holland and Leinhardt classification for the {i,j,k} triad
in graph g; otherwise, compute the four-state undirected classification.  Note
that this routine assumes dichotomized data.

This routine is not intended to be called from R.
*/
{
  int m,a,n,di,dj,dk;
  
  /*If we are in the undirected case, take the easy way out*/
  if(!gm)
    return g[i+j*gn]+g[j+k*gn]+g[i+k*gn];
      
  /*Get MAN information*/
  m=g[i+j*gn]*g[j+i*gn]+g[i+k*gn]*g[k+i*gn]+g[j+k*gn]*g[k+j*gn];
  n=(1-g[i+j*gn])*(1-g[j+i*gn])+(1-g[i+k*gn])*(1-g[k+i*gn])+(1-g[j+k*gn])*(1-g[k+j*gn]);
  a=3-m-n;

  /*Now classify, using dyad census as a first cut*/
  if(n==3)                    /*003*/
    return 0;
  else if((a==1)&&(n==2))     /*012*/
    return 1; 
  else if((m==1)&&(n==2))     /*102*/
    return 2; 
  else if((a==2)&&(n==1)){    /*021*/
    di=g[i+j*gn]+g[i+k*gn];
    if(di==2)
      return 3;                 /*021D*/
    dj=g[j+i*gn]+g[j+k*gn];
    if(dj==2)
      return 3;                 /*021D*/
    dk=g[k+i*gn]+g[k+j*gn];
    if(dk==2)
      return 3;                 /*021D*/
    di=g[j+i*gn]+g[k+i*gn];
    if(di==2)
      return 4;                 /*021U*/
    dj=g[i+j*gn]+g[k+j*gn];
    if(dj==2)
      return 4;                 /*021U*/
    dk=g[i+k*gn]+g[j+k*gn];
    if(dk==2)
      return 4;                 /*021U*/
    return 5;                   /*021C*/
  }else if((m==1)&&(n==1)){   /*111*/
    di=g[j+i*gn]+g[k+i*gn];
    if((di==0)||(di==2))
      return 6;                 /*111D*/
    dj=g[i+j*gn]+g[k+j*gn];
    if((dj==0)||(dj==2))
      return 6;                 /*111D*/
    return 7;                   /*111U*/
  }else if(a==3){            /*030*/
    di=g[j+i*gn]+g[k+i*gn];
    if((di==2)||(di==0))
      return 8;                 /*030T*/ 
    dj=g[i+j*gn]+g[k+j*gn];
    if((dj==2)||(dj==0))
      return 8;                 /*030T*/ 
    return 9;                   /*030C*/ 
  }else if((m==2)&&(n==1))    /*201*/
    return 10; 
  else if((m==1)&&(a==2)){    /*120*/
    di=g[j+i*gn]+g[k+i*gn];
    if(di==0)
      return 11;                /*120D*/
    dj=g[i+j*gn]+g[k+j*gn];
    if(dj==0)
      return 11;                /*120D*/
    dk=g[i+k*gn]+g[j+k*gn];
    if(dk==0)
      return 11;                /*120D*/
    di=g[i+j*gn]+g[i+k*gn];
    if(di==0)
      return 12;                /*120U*/
    dj=g[j+i*gn]+g[j+k*gn];
    if(dj==0)
      return 12;                /*120U*/
    dk=g[k+i*gn]+g[k+j*gn];
    if(dk==0)
      return 12;                /*120U*/
    return 13;                  /*120C*/
  }else if((m==2)&&(a==1))    /*210*/
    return 14; 
  else                        /*300*/
    return 15; 
}


void triad_classify_R(int *g, int *tt, int *gm)
/*
Given a triadic adjacency matrix, classify the triad in question.  Note that this routine assumes dichotomized data.  (This is a wrapper for triad_classify.)

This routine may be called from R using .C.
*/
{

/*Perform the classification*/
*tt=triad_classify(g,3,0,1,2,*gm);
}


void triad_census_R(int *g, int *n, double *t, int *gm)
/*
Compute a Holland and Leinhardt triad census for the graph with adjacency
matrix g.  Note that this routine assumes dichotomized data.

This routine may be called from R using .C.
*/
{
  int i,j,k;

  /*Clear out triad structure*/
  for(i=0;i<4+(*gm)*12;i++)     /*Vector length depends on gm*/
    t[i]=0.0;
  /*Get the triad counts*/
  for(i=0;i<*n;i++)
    for(j=i+1;j<*n;j++)
      for(k=j+1;k<*n;k++)
        t[triad_classify(g,*n,i,j,k,*gm)]++;
}
