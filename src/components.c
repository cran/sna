/*
######################################################################
#
# components.c
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 11/26/04
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines related to the identification of 
# components.
#
######################################################################
*/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "components.h"

void component_dist_R(double *g, double *pn, double *memb)
/*
Determine component memberships in g.  The memberships are stored in memb,
which must be a zero-initialized vector of length *pn.
*/
{
  char *visited;
  long int n,v,nod,i,s1count;
  double comp=0.0;

  /*Set up stuff*/
  n=*pn;
  /*Allocate memory for visited list*/
  visited=(char *)R_alloc(n,sizeof(char));
  /*Cycle through each node, performing a BFS*/
  for(v=0;v<n;v++){
    if(memb[v]==0.0){   /*Ignore nodes w/known membership*/
      comp++;           /*Increment component counter*/
      for(i=0;i<n;i++)  /*Mark all nodes unvisited*/
        visited[i]=0;
      s1count=0;
      visited[v]++;     /*Mark v as "to be visited"*/
      s1count++;
      memb[v]=comp;     /*v belongs to new component*/
      while(s1count){
        while(s1count){
          /*Find next node to be visited, change state*/
          for(nod=v;visited[nod]!=1;nod++); /*Only OK b/c s1count>0*/
          visited[nod]=3;               /*Mark as visited*/
          s1count--;
          memb[nod]=comp;               /*Set membership to comp*/
          for(i=v+1;i<n;i++)            /*Walk the unvisited neighborhood*/
            if((g[nod+i*n]!=0.0)&&(visited[i]==0)){
              visited[i]=2;               /*Visit this next time*/
            }
        } /*Continue until we run out of nodes at this level*/
        /*Mark all "to-be-visited" nodes as visitable*/
        for(i=v+1;i<n;i++)
          if(visited[i]==2){
            visited[i]=1;
            s1count++;
          }
      } /*Keep going until all nodes are accounted for*/
    }
  }
}
