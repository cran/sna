/*
######################################################################
#
# randomgraph.c
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 12/27/04
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines related to the generation of random
# graphs.
#
######################################################################
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "randomgraph.h"

void udrewire_R(double *g, double *pn, double *pnv, double *pp)
/*Perform a uniform rewiring process on the adjacency array pointed
to by *g.  It is assumed that g contains a *pn x *pnv *pnv array, whose dyads
are rewired (symmetrically) with uniform probability *pp.*/
{
  long int n,nv,i,j,k,h,t;
  double p,r,tempht,tempth;
  
  /*Take care of preliminaries*/
  n=(long int)*pn;
  nv=(long int)*pnv;
  p=*pp;
  GetRNGstate();

  /*Rewire the array*/
  for(i=0;i<n;i++){
    for(j=0;j<nv;j++){
      for(k=j+1;k<nv;k++){
        /*Rewire w/prob p*/
        if(runif(0.0,1.0)<p){
          t=j;  /*Save the head, tail*/
          h=k;
          if(runif(0.0,1.0)<0.5){   /*Switch head or tail w/50% prob*/
            while((h==j)||(h==k)) /*Draw h until legal*/
              h=(long int)floor(runif(0.0,1.0)*nv);
          }else{
            while((t==j)||(t==k)) /*Draw t until legal*/
              t=(long int)floor(runif(0.0,1.0)*nv);
          }
          /*Swap the dyad states*/
          tempth=g[i+t*n+h*n*nv];
          tempht=g[i+h*n+t*n*nv];
          g[i+t*n+h*n*nv]=g[i+j*n+k*n*nv];
          g[i+h*n+t*n*nv]=g[i+k*n+j*n*nv];
          g[i+j*n+k*n*nv]=tempth;
          g[i+k*n+j*n*nv]=tempht;
       }
      }
    }    
  }
  /*Reset the random number generator*/
  PutRNGstate();
}

void wsrewire_R(double *gi, double *go, double *pn, double *pnv, double *pp)
/*Perform a Watts-Strogatz rewiring process on the adjacency array pointed
to by *gi, storing the results in *go.  It is assumed that gi contains a 
*pn x *pnv *pnv array, whose non-null dyads are rewired (symmetrically) with
uniform probability *pp.  *go should be a copy of *gi.*/
{
  long int n,nv,i,j,k,h,t;
  double p,r,tempht,tempth;
  char flag;
  
  /*Take care of preliminaries*/
  n=(long int)*pn;
  nv=(long int)*pnv;
  p=*pp;
  GetRNGstate();

  /*Rewire the array*/
  for(i=0;i<n;i++){
    for(j=0;j<nv;j++){
      for(k=j+1;k<nv;k++){
        /*If the original dyad is non-null, rewire it w/prob p*/
        if(((gi[i+j*n+k*n*nv]!=0.0)||(gi[i+j*n+k*n*nv]!=0.0)) &&(runif(0.0,1.0)<p)){
          flag=0;
          while(!flag){
            t=j;  /*Save the head, tail*/
            h=k;
            if(runif(0.0,1.0)<0.5){   /*Switch head or tail w/50% prob*/
              h=(long int)floor(runif(0.0,1.0)*nv);
              if((h!=j)&&(h!=k)&&(go[i+t*n+h*n*nv]==0.0)&& (go[i+h*n+t*n*nv]==0.0)) /*Is h legal?*/
                flag++;
            }else{
              t=(long int)floor(runif(0.0,1.0)*nv);
              if((t!=j)&&(t!=k)&&(go[i+t*n+h*n*nv]==0.0)&& (go[i+h*n+t*n*nv]==0.0)) /*Is t legal?*/
                flag++;
            }
          }
          /*Swap the dyad states*/
          tempth=go[i+t*n+h*n*nv];
          tempht=go[i+h*n+t*n*nv];
          go[i+t*n+h*n*nv]=go[i+j*n+k*n*nv];
          go[i+h*n+t*n*nv]=go[i+k*n+j*n*nv];
          go[i+j*n+k*n*nv]=tempth;
          go[i+k*n+j*n*nv]=tempht;
        }
      }
    }    
  }
  /*Reset the random number generator*/
  PutRNGstate();
}
