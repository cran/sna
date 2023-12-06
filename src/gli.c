/*
######################################################################
#
# gli.c
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 6/25/09
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
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


void brokerage_R(double *g, int *pn, int *pm, int *cl, double *brok)
/*Calculate Gould-Fernandez brokerage scores for the vertices of g, based on
the vertex class vector cl.  Scores are recorded in a *pn x 5 matrix, brok,
whose columns are (in order) counts of coordinator, representative, 
gatekeeper, itinerant, and liason broker roles for each vertex.*/
{
  int n,i,j,k;
  slelement *ep,*ep2;
  snaNet *net;

  /*Set things up*/
  n=*pn;
  for(i=0;i<n;i++)
    for(j=0;j<5;j++)
      brok[i+n*j]=0.0;
  GetRNGstate();
  net=elMatTosnaNet(g,pn,pm);
  PutRNGstate();

  /*Calculate those scores!*/
  for(i=0;i<n;i++){                /*Walk the egos*/
    for(ep=snaFirstEdge(net,i,1);ep!=NULL;ep=ep->next[0])    /*ego->alter*/
      if(ep->val!=(double)i){
        for(ep2=snaFirstEdge(net,(int)(ep->val),1);ep2!=NULL;ep2=ep2->next[0])  /*alt->alt*/
          if((ep2->val!=(double)i)&&(ep2->val!=ep->val)){  /*Found 2-path?*/
            if(!snaIsAdjacent(i,(int)(ep2->val),net,0)){   /*Found broker?*/
              j=(int)(ep->val);
              k=(int)(ep2->val);
              /*Classify by type*/
              if(cl[j]==cl[i]){
                if(cl[j]==cl[k])
                  brok[j]++;         /*Type 0: Within-group (wI) [i j k]*/
                else
                  brok[j+2*n]++;     /*Type 2: Representative (bIO) [i j] [k]*/
              }else if(cl[j]==cl[k]){
                brok[j+3*n]++;       /*Type 3: Gatekeeping (bOI) [i] [j k]*/
              }else if(cl[i]==cl[k]){
                brok[j+n]++;         /*Type 2: Itinerant (WO) [j] [i k]*/
              }else
                brok[j+4*n]++;       /*Type 4: Liason (bO) [i] [j] [k]*/
            }
          }
      }
  }
}


void connectedness_R(double *mat, int *n, int *m, double *con)
/*Compute Krackhardt's connectedness for the graph in mat (which must be pre-
symmetrized using mode=="weak", since the measure is semipath based).*/
{
  snaNet *g;
  int i,*memb,*csize;

  /*Calculate the weak components of g*/  
  GetRNGstate();
  g=elMatTosnaNet(mat,n,m);
  PutRNGstate();
  memb=undirComponents(g);

  /*Tabulate the component sizes*/
  csize=(int *)R_alloc(memb[0],sizeof(int));
  for(i=0;i<memb[0];i++)
    csize[i]=0;
  for(i=0;i<*n;i++)
    csize[memb[i+1]-1]++;
  
  /*Compute the connectedness score (density of the reachability matrix)*/
  *con=0.0;
  for(i=0;i<memb[0];i++)
    *con+=csize[i]*(csize[i]-1.0)/2.0;
  *con/=(*n)*(*n-1.0)/2.0;
}


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
