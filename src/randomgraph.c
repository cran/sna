/*
######################################################################
#
# randomgraph.c
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 4/15/06
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

void bn_mcmc_R(int *g, double *pn, double *pdraws, double *pburn, int *pthin, double *pi, double *sigma, double *rho, double *d)
{
  long int n,i,j,k,x,draws,burn,bc,*parents;
  double lne,lnpar,lnsib,lndblr,ep;
  int thin,tc;
  
  /*Initialize various things*/
  n=(long int)*pn;
  draws=(long int)*pdraws;
  burn=(long int)*pburn;
  thin=(int)*pthin;
  GetRNGstate();
  parents=(long int *)R_alloc(n*n,sizeof(long int));
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      g[i*draws+j*n*draws]=0;
      parents[i+j*n]=0;  /*Eventually fix symmetry condition*/
    }    
  }  
  lne=log(1.0-*d);
  lnpar=log(1.0-*pi);
  lnsib=log(1.0-*sigma);
  lndblr=log(1.0-*rho);

  /*Run the MCMC loop*/
  bc=0;
  tc=0;
  for(i=1;i<draws;i++){
    /*
    if(bc<burn)
      Rprintf("Burn-in Iteration %ld\n",bc);
    else
      Rprintf("Iteration %ld\n",i);
    Rprintf("\tCopying old state\n");
    */
    /*Copy our old state into the next iteration slot*/
    for(j=0;j<n;j++)
      for(k=0;k<n;k++)
        g[i+j*draws+k*n*draws]=g[i-1+j*draws+k*n*draws];
    /*Rprintf("\tEdge selection\n");*/
    /*Select an edge to redraw*/
    j=(long int)floor(runif(0.0,1.0)*n);
    do{
      k=(long int)floor(runif(0.0,1.0)*n);
    }while(j==k);
    /*Rprintf("\tDrawing and updating\n");*/
    /*Redraw the edge from the full conditional*/
    ep=1.0-exp(lne+g[i+k*draws+j*n*draws]*lnpar+parents[j+k*n]*lnsib+ g[i+k*draws+j*n*draws]*parents[j+k*n]*lndblr);
    if(runif(0.0,1.0)<=ep){
      g[i+j*draws+k*n*draws]=1;    /*Set the edge*/
      /*If something has changed update the parent count*/
      if(g[i-1+j*draws+k*n*draws]==0){
        /*Rprintf("\tAdded edge update\n");*/
        for(x=0;x<n;x++)
          if((g[i+j*draws+x*n*draws])&&(j!=x)&&(x!=k)){
            parents[k+x*n]++;
            parents[x+k*n]++;
            /*
            if(parents[k+x*n]>n-2)
              Rprintf("Parent overflow: iter=%ld j=%ld k=%ld x=%ld",i,j,k,x);
            */
          }
      }
    }else{
      g[i+j*draws+k*n*draws]=0;   /*Unset the edge*/
      /*If something has changed update the parent count*/
      if(g[i-1+j*draws+k*n*draws]==1){
        /*Rprintf("\tDeleted edge update\n");*/
        for(x=0;x<n;x++)
          if((g[i+j*draws+x*n*draws])&&(j!=x)&&(x!=k)){
            parents[k+x*n]--;
            parents[x+k*n]--;
            /*
            if(parents[k+x*n]<0)
              Rprintf("Parent underflow: iter=%ld j=%ld k=%ld x=%ld",i,j,k,x);
            */
          }
      }
    }
    /*Burn-in check*/
    if(bc<burn){
      /*Save the state*/
      for(j=0;j<n;j++)
        for(k=0;k<n;k++)
          g[i-1+j*draws+k*n*draws]=g[i+j*draws+k*n*draws];
      i--;   
      bc++;
    }else{ 
      /*Rprintf("\tThinning?\n");*/
      if(tc%thin!=0){
        /*Rprintf("\tThinning\n");*/
        /*Save the state*/
        for(j=0;j<n;j++)
          for(k=0;k<n;k++)
            g[i-1+j*draws+k*n*draws]=g[i+j*draws+k*n*draws];
        i--;   
      }
      tc++;
    }
  }
  /*Save the RNG state*/
  PutRNGstate();
}

void udrewire_R(double *g, double *pn, double *pnv, double *pp)
/*Perform a uniform rewiring process on the adjacency array pointed
to by *g.  It is assumed that g contains a *pn x *pnv *pnv array, whose dyads
are rewired (symmetrically) with uniform probability *pp.*/
{
  long int n,nv,i,j,k,h,t;
  double p,tempht,tempth;
  
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
  double p,tempht,tempth;
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
