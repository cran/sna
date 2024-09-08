/*
######################################################################
#
# paths.c
#
# copyright (c) 2007, Carter T. Butts <buttsc@uci.edu>
# Last Modified 9/7/24
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
#
# Part of the R/sna package
#
# This file contains routines related to cycle and path counting.
#
######################################################################
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "utils.h"
#include "paths.h"


/*INTERNAL ROUTINES---------------------------------------------------------*/

void edgewisePathRecurse(snaNet *g, int src, int dest, int curnode, int *availnodes, int availcount, int *usednodes, int curlen, double *count, double *cpcount, double *dpcount, int maxlen, int directed, int byvertex, int copaths, int dyadpaths)
/*Recursively count the paths from src to dest.  (This is an adaptation of the routine I wrote for statnet.)  count should be vector of path counts (starting with length 1) if !byvertex, or else a matrix of dimension (maxlen-1)x(n+1) whose first column contains aggregate counts and i+1th column contains counts for vertex i.  If copaths==1, cpcount should contain an nxn matrix containing path co-membership counts.  If copaths=2, cpcount should contain a (maxlen-1)xnxn array containing path co-membership counts at each length.  This is ignored if copaths==0.  (Odd note: if this is being used to construct a path census, rather than a cycle census, maxlen should be one greater than the true maximum length.  It looks wrong, but it isn't.  All of the maxlen-1 stuff in here is also correct, even though (e.g., for dyad paths) it may appear otherwise.)*/
{
  int *newavail,i,j,k,newavailcount,*newused,n;

  /*Rprintf("\t\t\tRecursion: src=%d, dest=%d, curnode=%d, curlen=%d, availcount=%d\n",src,dest, curnode,curlen,availcount);*/

  n=g->n;
  /*Rprintf("N=%d\n",n);*/
  /*If we've found a path to the destination, increment the census vector*/ 
  if(directed||(curnode<dest)){
    if(snaIsAdjacent(curnode,dest,g,2)){
      /*Rprintf("\t\t\t\t%d is adjacent to target (%d)\n",curnode,dest);*/
      count[curlen]++;                       /*Basic update*/
      if(byvertex){                          /*Update path incidence counts*/
        for(j=0;j<curlen;j++)
          count[curlen+(1+usednodes[j])*(maxlen-1)]++;
        count[curlen+(1+curnode)*(maxlen-1)]++;
        count[curlen+(1+dest)*(maxlen-1)]++;
      }
      if(copaths==1){                        /*Update copath incidence counts*/
        for(j=0;j<curlen;j++){
          for(k=j;k<curlen;k++){
            cpcount[usednodes[j]+usednodes[k]*n]++;
            if(k!=j)
              cpcount[usednodes[k]+usednodes[j]*n]++;
          }
          cpcount[usednodes[j]+dest*n]++;
          cpcount[dest+usednodes[j]*n]++;
          cpcount[usednodes[j]+curnode*n]++;
          cpcount[curnode+usednodes[j]*n]++;
        }
        cpcount[curnode+dest*n]++;
        cpcount[dest+curnode*n]++;
        cpcount[curnode+curnode*n]++;
        cpcount[dest+dest*n]++;
      }
      if(copaths==2){                        /*Update copath counts using len*/
        for(j=0;j<curlen;j++){
          for(k=j;k<curlen;k++){
            cpcount[curlen+usednodes[j]*(maxlen-1)+ usednodes[k]*(maxlen-1)*n]++;
            if(k!=j)
              cpcount[curlen+usednodes[k]*(maxlen-1)+ usednodes[j]*(maxlen-1)*n]++;
          }
          cpcount[curlen+usednodes[j]*(maxlen-1)+dest*(maxlen-1)*n]++;
          cpcount[curlen+dest*(maxlen-1)+usednodes[j]*(maxlen-1)*n]++;
          cpcount[curlen+usednodes[j]*(maxlen-1)+curnode*(maxlen-1)*n]++;
          cpcount[curlen+curnode*(maxlen-1)+usednodes[j]*(maxlen-1)*n]++;
        }
        cpcount[curlen+curnode*(maxlen-1)+dest*(maxlen-1)*n]++;
        cpcount[curlen+dest*(maxlen-1)+curnode*(maxlen-1)*n]++;
        cpcount[curlen+dest*(maxlen-1)+dest*(maxlen-1)*n]++;
        cpcount[curlen+curnode*(maxlen-1)+curnode*(maxlen-1)*n]++;
      }
      if(dyadpaths==1){                      /*Update dyadic path counts*/
        dpcount[src+dest*n]++;
        if(!directed)
          dpcount[dest+src*n]++;
      }
      if(dyadpaths==2){                  /*Update dyadic path counts using len*/
        dpcount[curlen+src*(maxlen-1)+dest*(maxlen-1)*n]++;
        if(!directed)
          dpcount[curlen+dest*(maxlen-1)+src*(maxlen-1)*n]++;
      }
    }
  }else{
    if(snaIsAdjacent(dest,curnode,g,2)){
      count[curlen]++;                       /*Basic update*/
      if(byvertex){                          /*Update path incidence counts*/
        for(j=0;j<curlen;j++)
          count[curlen+(1+usednodes[j])*(maxlen-1)]++;
        count[curlen+(1+curnode)*(maxlen-1)]++;
        count[curlen+(1+dest)*(maxlen-1)]++;
      }
      if(copaths==1){                       /*Update copath incidence counts*/
        for(j=0;j<curlen;j++){
          for(k=j;k<curlen;k++){
            cpcount[usednodes[j]+usednodes[k]*n]++;
            if(k!=j)
              cpcount[usednodes[k]+usednodes[j]*n]++;
          }
          cpcount[usednodes[j]+dest*n]++;
          cpcount[dest+usednodes[j]*n]++;
          cpcount[usednodes[j]+curnode*n]++;
          cpcount[curnode+usednodes[j]*n]++;
        }
        cpcount[curnode+dest*n]++;
        cpcount[dest+curnode*n]++;
        cpcount[curnode+curnode*n]++;
        cpcount[dest+dest*n]++;
      }
      if(copaths==2){                      /*Update copath counts using len*/
        for(j=0;j<curlen;j++){
          for(k=j;k<curlen;k++){
            cpcount[curlen+usednodes[j]*(maxlen-1)+ usednodes[k]*(maxlen-1)*n]++;
            if(k!=j)
              cpcount[curlen+usednodes[k]*(maxlen-1)+ usednodes[j]*(maxlen-1)*n]++;
          }
          cpcount[curlen+usednodes[j]*(maxlen-1)+dest*(maxlen-1)*n]++;
          cpcount[curlen+dest*(maxlen-1)+usednodes[j]*(maxlen-1)*n]++;
          cpcount[curlen+usednodes[j]*(maxlen-1)+curnode*(maxlen-1)*n]++;
          cpcount[curlen+curnode*(maxlen-1)+usednodes[j]*(maxlen-1)*n]++;
        }
        cpcount[curlen+curnode*(maxlen-1)+dest*(maxlen-1)*n]++;
        cpcount[curlen+dest*(maxlen-1)+curnode*(maxlen-1)*n]++;
        cpcount[curlen+curnode*(maxlen-1)+curnode*(maxlen-1)*n]++;
        cpcount[curlen+dest*(maxlen-1)+dest*(maxlen-1)*n]++;
      }
      if(dyadpaths==1){                      /*Update dyadic path counts*/
        dpcount[src+dest*n]++;
        if(!directed)
          dpcount[dest+src*n]++;
      }
      if(dyadpaths==2){                  /*Update dyadic path counts using len*/
        dpcount[curlen+src*(maxlen-1)+dest*(maxlen-1)*n]++;
        if(!directed)
          dpcount[curlen+dest*(maxlen-1)+src*(maxlen-1)*n]++;
      }
    }
  }
  
  /*If possible, keep searching for novel paths*/
  if((availcount>0)&&(curlen<maxlen-2)){
    if(availcount>1){    /*Remove the current node from the available list*/
      /*Rprintf("\t\t\tRemoving %d from available node list (availcount=%d)\n", curnode,availcount);*/
      if((newavail=(int *)R_Calloc((size_t)(availcount-1),int))==NULL){
        Rprintf("Unable to allocate %ld bytes for available node list in edgewisePathRecurse.  Trying to terminate recursion gracefully, but your path count is probably wrong.\n",(long int)(sizeof(int)*(availcount-1)));
        return;
      }
      j=0;
      for(i=0;i<availcount;i++)      /*Create the reduced list, fur passin'*/
        if(availnodes[i]!=curnode)
          newavail[j++]=availnodes[i];
      /*Rprintf("\t\t\tBuilt newavail without apparent issue\n");*/
    }else
      newavail=NULL;                 /*Set to NULL if we're out of nodes*/
    newavailcount=availcount-1;      /*Decrement the available count*/
    if(byvertex||copaths||dyadpaths){  /*Add the current node to the used list*/
      if((newused=(int *)R_Calloc((size_t)(curlen+1), int))==NULL){
        Rprintf("Unable to allocate %ld bytes for used node list in edgewisePathRecurse.  Trying to terminate recursion gracefully, but your path count is probably wrong.\n",(long int)(sizeof(int)*(curlen+1)));
        return;
      }
      for(i=0;i<curlen;i++)
        newused[i]=usednodes[i];
      newused[curlen]=curnode;
    }else
      newused=NULL;
    /*Recurse on all available nodes*/
    /*Rprintf("\t\t\tAbout to recurse on available nodes (newavail=%d)\n", newavailcount);*/
    for(i=0;i<newavailcount;i++)
      if(directed||(curnode<newavail[i])){
        if(snaIsAdjacent(curnode,newavail[i],g,2))
          edgewisePathRecurse(g,src,dest,newavail[i],newavail,newavailcount,
            newused,curlen+1,count,cpcount,dpcount,maxlen,directed,byvertex,
            copaths,dyadpaths);
      }else{
        if(snaIsAdjacent(newavail[i],curnode,g,2))
          edgewisePathRecurse(g,src,dest,newavail[i],newavail,newavailcount,
            newused,curlen+1,count,cpcount,dpcount,maxlen,directed,byvertex,
            copaths,dyadpaths);
      }
    /*Free the available node and used node lists*/
    /*Rprintf("\t\t\tDone with available node recursion, freeing\n");*/
    if(newavail!=NULL){
      /*Rprintf("\t\t\t\tFreeing newavail; count=%d\n",newavailcount);*/
      R_Free(newavail);
    }
    if(newused!=NULL){
      /*Rprintf("\t\t\t\tFreeing newused; count=%d\n",curlen);*/
      R_Free(newused);
    }
  }

  /*Free the used node list*/
  /*Rprintf("\t\t\tBacking out for src=%d, dest=%d, curnode=%d\n",src,dest, curnode);*/
  /*if(usednodes!=NULL)
    free((void *)usednodes);*/
  /*Check for interrupts (if recursion is taking way too long...)*/
  R_CheckUserInterrupt();
}


void edgewiseCycleCensus(snaNet *g, int src, int dest, double *count, double *cccount, int maxlen, int directed, int byvertex, int cocycles)
/*Count the number of cycles associated with the (src,dest) edge in g, assuming that this edge exists.  The byvertex and cocycles flags indicate whether cycle counts should be broken down by participating vertex, and whether a cycle co-membership matrix should be returned (respectively).  In either case, count and cccount must be structured per count and pccount in edgewisePathRecurse.*/
{
  int n,i,j,*availnodes,*usednodes;

  /*Set things up*/
  n=g->n;
  usednodes=NULL;

  /*First, check for a 2-cycle (but only if directed)*/
  /*Rprintf("\t\tChecking for (%d,%d) edge\n",dest,src);*/
  if(directed&&snaIsAdjacent(dest,src,g,2)){
    count[0]++;
    if(byvertex){
      count[(1+src)*(maxlen-1)]++;
      count[(1+dest)*(maxlen-1)]++;
    }
    if(cocycles==1){
      cccount[src+dest*n]++;
      cccount[dest+src*n]++;
      cccount[src+src*n]++;
      cccount[dest+dest*n]++;
    }
    if(cocycles==2){
      cccount[src*(maxlen-1)+dest*(maxlen-1)*n]++;
      cccount[dest*(maxlen-1)+src*(maxlen-1)*n]++;
      cccount[src*(maxlen-1)+src*(maxlen-1)*n]++;
      cccount[dest*(maxlen-1)+dest*(maxlen-1)*n]++;
    }
  }
  if(n==2)
    return;                 /*Failsafe for graphs of order 2*/
  
  /*Perform the recursive path count*/
  if((availnodes=(int *)R_Calloc((size_t)(n-2),int))==NULL){
    Rprintf("Unable to allocate %ld bytes for available node list in edgewiseCycleCensus.  Exiting.\n",(long int)(sizeof(int)*(n-2)));
    return;
  }
  j=0;                             /*Initialize the list of available nodes*/
  for(i=0;i<n;i++)
    if((i!=src)&&(i!=dest))
      availnodes[j++]=i;
  if(byvertex||cocycles){          /*Initialize the list of already used nodes*/
    if((usednodes=(int *)R_Calloc(1,int))==NULL){
      Rprintf("Unable to allocate %ld bytes for used node list in edgewiseCycleCensus.  Exiting.\n",(long int)(sizeof(int)));
      return;
    }
    usednodes[0]=dest;
  }
  /*Rprintf("\t\tBeginning recursion\n");*/
  for(i=0;i<n-2;i++)               /*Recurse on each available vertex*/
    if(directed||(dest<availnodes[i])){
      if(snaIsAdjacent(dest,availnodes[i],g,2))
        edgewisePathRecurse(g,dest,src,availnodes[i],availnodes,n-2,usednodes,1,
          count,cccount,NULL,maxlen,directed,byvertex,cocycles,0);
    }else{
      if(snaIsAdjacent(availnodes[i],dest,g,2))
        edgewisePathRecurse(g,dest,src,availnodes[i],availnodes,n-2,usednodes,1,
          count,cccount,NULL,maxlen,directed,byvertex,cocycles,0);
    }
  /*Rprintf("\t\tReturned from recursion; freeing memory\n");*/
  if(availnodes!=NULL)
    R_Free(availnodes);  /*Free the available node list*/
  if(usednodes!=NULL)
    R_Free(usednodes); /*Free the used node list, if needed*/
}


void dyadPathCensus(snaNet *g, int src, int dest, double *count, double *cpcount, double *dpcount, int maxlen, int directed, int byvertex, int copaths, int dyadpaths)
/*Calculate a path census.  Note that edgewisePathRecurse is called with maxlen of 1 greater than the true value.  This is strange, but correct!*/
{
  int n,i,j,*availnodes,*usednodes;

  /*Set things up*/
  n=g->n;
  usednodes=NULL;

  if(n<2)
    return;                 /*Failsafe for graphs of order 2*/

  /*Check for a 1-path (i.e., edge)*/
  /*Rprintf("\t\tChecking for (%d,%d) edge\n",src,dest);*/
  if(snaIsAdjacent(src,dest,g,2)||((!directed)&&snaIsAdjacent(dest,src,g,2))){
    count[0]++;
    if(byvertex){
      count[(1+src)*maxlen]++;
      count[(1+dest)*maxlen]++;
    }
    if(copaths==1){
      cpcount[src+dest*n]++;
      cpcount[dest+src*n]++;
      cpcount[src+src*n]++;
      cpcount[dest+dest*n]++;
    }
    if(copaths==2){
      cpcount[src*maxlen+dest*maxlen*n]++;
      cpcount[dest*maxlen+src*maxlen*n]++;
      cpcount[src*maxlen+src*maxlen*n]++;
      cpcount[dest*maxlen+dest*maxlen*n]++;
    }
    if(dyadpaths==1){                      /*Update dyadic path counts*/
      dpcount[src+dest*n]++;
      if(!directed)
        dpcount[dest+src*n]++;
    }
    if(dyadpaths==2){                  /*Update dyadic path counts using len*/
      dpcount[src*maxlen+dest*maxlen*n]++;
      if(!directed)
        dpcount[dest*maxlen+src*maxlen*n]++;
    }
  }
  
  /*Perform the recursive path count*/
  if((availnodes=(int *)R_Calloc((size_t)(n-2),int))==NULL){
    Rprintf("Unable to allocate %ld bytes for available node list in dyadPathCensus.  Exiting.\n",(long int)sizeof(int)*(n-2));
    return;
  }
  j=0;                             /*Initialize the list of available nodes*/
  for(i=0;i<n;i++)
    if((i!=src)&&(i!=dest))
      availnodes[j++]=i;
  if(byvertex||copaths){          /*Initialize the list of already used nodes*/
    if((usednodes=(int *)R_Calloc(1,int))==NULL){
      Rprintf("Unable to allocate %ld bytes for used node list in edgewiseCycleCensus.  Exiting.\n",(long int)sizeof(int));
      return;
    }
    usednodes[0]=src;
  }
  for(i=0;i<n-2;i++)               /*Recurse on each available vertex*/
    if(directed||(dest<availnodes[i])){  
      if(snaIsAdjacent(src,availnodes[i],g,2))
        edgewisePathRecurse(g,src,dest,availnodes[i],availnodes,n-2,usednodes,1,
          count,cpcount,dpcount,maxlen+1,directed,byvertex,copaths,dyadpaths);
    }else{
      if(snaIsAdjacent(availnodes[i],src,g,2))
        edgewisePathRecurse(g,src,dest,availnodes[i],availnodes,n-2,usednodes,1,
          count,cpcount,dpcount,maxlen+1,directed,byvertex,copaths,dyadpaths);
    }
  R_Free(availnodes);  /*Free the available node list*/
  if(usednodes!=NULL)
    R_Free(usednodes); /*Free the used node list, if needed*/
}


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

void cycleCensus_R(int *g, int *pn, int *pm, double *count, double *cccount, int *pmaxlen, int *pdirected, int *pbyvertex, int *pcocycles)
/*Count the number of cycles associated with the (src,dest) edge in g, assuming that this edge exists.  The byvertex and cocycles flags indicate whether cycle counts should be broken down by participating vertex, and whether cycle co-membership counts should be returned (respectively).  In either case, count and cccount must be structured per count and pccount in edgewisePathRecurse.*/
{
  int i,r,c,n,m;
  double *dval;
  snaNet *ng;

  GetRNGstate();
  /*Allocate memory for the new graph object*/
  /*Rprintf("Initializing ng\n");*/
  n=(*pn);
  m=(*pm);
  ng=(snaNet *)R_alloc(1,sizeof(struct snaNettype));
  ng->n=(*pn);
  ng->indeg=(int *)R_alloc(n,sizeof(int));
  ng->outdeg=(int *)R_alloc(n,sizeof(int));
  ng->iel=(slelement **)R_alloc(n,sizeof(slelement *));
  ng->oel=(slelement **)R_alloc(n,sizeof(slelement *));

  /*Initialize the graph*/
  for(i=0;i<n;i++){
    ng->indeg[i]=0;
    ng->outdeg[i]=0;
    ng->iel[i]=NULL;
    ng->oel[i]=NULL;
  }

  /*Walk the graph, adding edges and accumulating cycles*/
  /*Rprintf("Building graph/accumulating cycles\n\tn=%d,%d\n",n,ng->n);*/
  for(i=0;i<m;i++)
    if((!IISNA(g[i+2*m]))&&((*pdirected)||(g[i]<g[i+m]))){
      r=g[i]-1;
      c=g[i+m]-1;
      /*First, accumulate the cycles to be formed by the (r,c) edge*/
      /*Rprintf("\tEdge at (%d,%d); counting cycles\n",r+1,c+1);*/
      edgewiseCycleCensus(ng,r,c,count,cccount,*pmaxlen,*pdirected, 
        *pbyvertex,*pcocycles);
      /*for(k=0;k<*pmaxlen-1;k++)
      Rprintf("%d:%d ",k+2,(int)(count[k]));
      Rprintf("\n");*/
      /*Next, add the (r,c) edge to the graph*/
      /*Rprintf("\tGot cycles, now adding edge to ng\n");*/
      dval=(double *)R_alloc(1,sizeof(double));   /*Create iel element*/
      dval[0]=(double)g[i+2*m];
      ng->iel[c]=slistInsert(ng->iel[c],(double)r,(void *)dval);
      ng->indeg[c]++;
      dval=(double *)R_alloc(1,sizeof(double));   /*Create oel element*/
      dval[0]=(double)g[i+2*m];
      ng->oel[r]=slistInsert(ng->oel[r],(double)c,(void *)dval);
      ng->outdeg[r]++;
      if(!(*pdirected)){
        dval=(double *)R_alloc(1,sizeof(double));   /*Create iel element*/
        dval[0]=(double)g[i+2*m];
        ng->iel[r]=slistInsert(ng->iel[r],(double)c,(void *)dval);
        ng->indeg[r]++;
        dval=(double *)R_alloc(1,sizeof(double));   /*Create oel element*/
        dval[0]=(double)g[i+2*m];
        ng->oel[c]=slistInsert(ng->oel[c],(double)r,(void *)dval);
        ng->outdeg[c]++;
      }
    }
    
  PutRNGstate();
}


void pathCensus_R(double *g, int *pn, int *pm, double *count, double *cpcount, double *dpcount, int *pmaxlen, int *pdirected, int *pbyvertex, int *pcopaths, int *pdyadpaths)
/*Conduct a census of paths in g, out to length maxlen.  The byvertex and copaths flags indicate whether path counts should be broken down by participating vertex, and whether path co-membership counts should be returned (respectively).  In either case, count and pccount must be structured per count and pccount in edgewisePathRecurse.*/
{
  int i,j,n;
  snaNet *ng;

  /*Create the new graph object*/
  n=(*pn);
  GetRNGstate();
  ng=elMatTosnaNet(g,pn,pm);

  /*Walk the graph, counting paths associated with each pair*/
  for(i=0;i<n;i++)
    for(j=(!(*pdirected))*(i+1);j<n;j++)
      if(i!=j){
        dyadPathCensus(ng,i,j,count,cpcount,dpcount,*pmaxlen,*pdirected, 
          *pbyvertex,*pcopaths,*pdyadpaths);
      }
      
  PutRNGstate();
}

