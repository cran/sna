/*
######################################################################
#
# components.c
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 03/26/09
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


slelement *BFS(snaNet *g, int *n, int v, int transpose)
/*This is a boring old BFS, included as a utility.  It starts with vertex v, and proceeds forward (transpose=0) or backward (transpose=1), building a skip list of reached vertices.  A pointer to the resulting list is returned.*/
{
  int i;
  char *vis;
  element *tovis=NULL,cur;
  slelement *reach=NULL,*ep;
  
  /*Initialize visit list*/  
  vis=(char *)R_alloc(*n,sizeof(char));
  for(i=0;i<*n;i++)
    vis[i]=0;

  /*Run the BFS*/
  tovis=push(tovis,(double)v,NULL);
  vis[v]=1;
  while(tovis!=NULL){
    /*Pop the node to visit, and add to list*/
    cur=pop(tovis);
    tovis=cur.next;
    reach=slistInsert(reach,cur.val,NULL);
    /*Add the neighbors of this node, if not already marked*/
    if(!transpose){
      for(ep=snaFirstEdge(g,(int)(cur.val),1);ep!=NULL;ep=ep->next[0])
        if(!vis[(int)(ep->val)]){
          tovis=push(tovis,ep->val,NULL);
          vis[(int)(ep->val)]++;
        }
    }else{
      for(ep=snaFirstEdge(g,(int)(cur.val),0);ep!=NULL;ep=ep->next[0])
        if(!vis[(int)(ep->val)]){
          tovis=push(tovis,ep->val,NULL);
          vis[(int)(ep->val)]++;
        }
    }
  }
  
  /*Return the reach list.  R will deallocate other memory (eventually).*/
  return reach;
}


int numStrongComponents(snaNet *g, int *n)
/*Return the number of strong components in g.*/
{
  int *comp,i,ccount=*n;

  /*Get the components*/
  comp=strongComponents(g,n);

  /*Count 'em (should be equal to n-min component number)*/
  for(i=0;i<*n;i++)
    ccount=MIN(ccount,comp[i]);

  return *n-ccount;
}


slelement *strongComponentByVertex(snaNet *g, int *n, int v)
/*Use a lame BFS algorithm to return the maximal strong component to which the specified vertex belongs (as an slist).  This is placed here as a utility, rather than as an R-callable routine.*/
{
  slelement *olist=NULL,*ilist=NULL,*comp=NULL;
  
  /*Get out/in lists*/
  olist=BFS(g,n,v,0);
  ilist=BFS(g,n,v,1);
  
  /*Find the intersection*/
  olist=olist->next[0];
  ilist=ilist->next[0];
  while((olist!=NULL)&&(ilist!=NULL)){
    if(olist->val==ilist->val){
      comp=slistInsert(comp,olist->val,NULL);
      ilist=ilist->next[0];
      olist=olist->next[0];
    }else{
      if(olist->val<ilist->val){
        olist=olist->next[0];
      }else{
        ilist=ilist->next[0];
      }
    }
  }
  
  /*Return the result*/
  return comp;
}


int *strongComponents(snaNet *g, int *n)
/*This function uses a variant of Tarjan's DFS algorithm published in a technical report by David J. Pearce to find the strongly connected components of g (which are returned via an index vector).  Pearce's algorithm has the same running time as Tarjan's, but is (very) slightly more space efficient.*/
{
  int i,*index,*ccount,*rindex;
  element *dfs;
  
  /*Initialize everything*/
  dfs=(element *)R_alloc(1,sizeof(element));
  rindex=(int *)R_alloc(*n,sizeof(int));
  index=(int *)R_alloc(1,sizeof(int));
  ccount=(int *)R_alloc(1,sizeof(int));
  for(i=0;i<*n;i++)
    rindex[i]=0;
  ccount[0]=*n-1;
  index[0]=1;
  dfs->next=NULL;

  /*Find the components*/
  for(i=0;i<*n;i++)
    if(rindex[i]==0)
      strongComponentsRecurse(g,n,i,rindex,index,ccount,dfs);

//  for(i=0;i<*n;i++)
//    Rprintf("%d ",rindex[i]);
//  Rprintf("\n");

  /*Return the result*/
  return rindex;
}


void strongComponentsRecurse(snaNet *g, int *n, int v, int *rindex, int *index, int *ccount, element *dfs)
{
  char root=1;
  element w;
  slelement *ep;
  
  /*Set index for v, and increment*/
  rindex[v]=*index;
  (*index)++;

  /*Visit v's unvisited out-neighbors*/
  for(ep=snaFirstEdge(g,v,1);ep!=NULL;ep=ep->next[0]){
    if(rindex[(int)(ep->val)]==0)
      strongComponentsRecurse(g,n,(int)(ep->val),rindex,index,ccount,dfs);
    if(rindex[(int)(ep->val)]<rindex[v]){
      rindex[v]=rindex[(int)ep->val];
      root=0;
    }
  }
  
  /*If v is the root of its tree, pop its subtree and mark as one component*/
  if(root){
    (*index)--;
    while((dfs->next!=NULL)&&(rindex[v]<=rindex[(int)(dfs->next->val)])){
      w=pop(dfs->next);
      dfs->next=w.next;
      rindex[(int)(w.val)]=*ccount;
      (*index)--;
    }
    rindex[v]=*ccount;
    (*ccount)--;
  }else                      /*Otherwise, add v to the stack*/
    dfs->next=push(dfs->next,(double)v,NULL);
}


int *undirComponents(snaNet *g)
/*Returns the components of an undirected graph g as a vector of integers.  The vector length is g->n+1, with the first element being the number of components.  (Yes, this is slightly different from strongComponents, above, and the two should be harmonized eventually.)*/
{
  int i,*memb;
  
  /*Initialize*/
  memb=(int *)R_alloc(g->n+1,sizeof(int));
  for(i=0;i<g->n+1;i++)
    memb[i]=0;
  
  /*Perform a DFS to find the components*/
  for(i=0;i<g->n;i++)
    if(memb[i+1]==0){
      memb[0]++;
      undirComponentsRecurse(g,i,memb);
    }
  
  return memb;
}


void undirComponentsRecurse(snaNet *g,int v,int *memb)
/*Recursive computation for components of an undirected graph (via DFS).  v is the current vertex, and memb is a membership vector with the i+1th element corresponding to the ith vertex (the first element contains the current number of components).  We assume that we are being called on the "current" component, and thus memb[0] is taken as the assignment for all reachable vertices; memb[i]==0 is used to indicate an unreached vertex.*/
{
  slelement *sp;
  
  memb[v+1]=memb[0];   /*Label the current vertex*/
  /*If any neighbors, recurse accordingly*/
  if(g->outdeg[v]>0)
    for(sp=g->oel[v]->next[0];sp!=NULL;sp=sp->next[0])
      if(memb[(int)(sp->val)+1]==0)
        undirComponentsRecurse(g,(int)(sp->val),memb);
}


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
