/*
######################################################################
#
# geodist.c
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 11/21/04
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines related to the computation of geodesics.
#
######################################################################
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "geodist.h"

void geodist_R(double *g, double *pn, double *gd, double *sigma)
/*
Compute geodesics for the graph in g.  The geodesic distances are stored in
gd, and the path counts in sigma (both being nxn matrices).  Note that these
should be initialized to all infs and all 0s, respectively.
*/
{
  char *visited;
  long int n,v,i,nod,s1count;

  /*Set up stuff*/
  n=*pn;
  /*Allocate memory for visited list*/
  visited=(char *)R_alloc(n,sizeof(char));
  /*Cycle through each node, performing a BFS*/
  for(v=0;v<n;v++){
    /*Clear the visit list*/
    for(i=0;i<n;i++)
      visited[i]=0; 
    s1count=0;
    /*Start with the source node*/
    nod=v;
    visited[nod]=1;
    s1count++;
    sigma[v+v*n]=1.0;
    gd[v+v*n]=0.0;
    /*Now, conduct the trace*/
    while(s1count>0){
      while(s1count>0){
        /*Find the next visitable node, and change its state*/
        for(nod=0;visited[nod]!=1;nod++); /*Only OK b/c s1count>0*/
        visited[nod]=3;
        s1count--;
        for(i=0;i<n;i++)   /*Walk the unvisited neighborhood*/
          if((g[nod+i*n]!=0.0)&&((visited[i]==0)||(visited[i]==2))){
            if(visited[i]==0)  /*If j is unvisited, visit it next time*/
              visited[i]=2;
            if(gd[v+i*n]-gd[v+nod*n]>=g[nod+i*n]){
              gd[v+i*n]=gd[v+nod*n]+g[nod+i*n];  /*Geodist is nod's+g*/
              sigma[v+i*n]+=sigma[v+nod*n];      /*Add to path count*/
            }
          }
      } /*Continue until we run out of nodes for this iteration*/
      /*Mark all "to-be-visited" nodes as visitable*/
      for(i=0;i<n;i++)
        if(visited[i]==2){
          visited[i]=1;
          s1count++;
        }
    } /*Keep going until all nodes are accounted for*/
  } 
}

