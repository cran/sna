/*
######################################################################
#
# cohesion.c
#
# copyright (c) 2007, Carter T. Butts <buttsc@uci.edu>
# Last Modified 8/28/09
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
#
# Part of the R/sna package
#
# This file contains routines related to the identification of 
# cohesive subgroups.
#
######################################################################
*/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "cohesion.h"

/*INTERNAL ROUTINES---------------------------------------------------------*/

void bicomponentRecurse(snaNet *g, element *complist, element *estack, int *parent, int *num, int *back, int *dfn, int v)
/*Perform a depth-first recusion to identify bicomponents.  This call implements
a visit to v.  Components are stored as vertex lists, in an element list with 
the following structure:

    dp-> (last component)
    |
  head -next-> component1 -next> ...
   |           |      |
  val         val    dp-> elem -next-> ...
(#comps)    (#vert)        |
                          val (vert ID)

Bicomponents are accumulated through their edges; these are stored in
*estack until needed.  val is used to store both endpoints, via the 
following encoding scheme:

(i,j) -> k=i+n*j; k -> (i=k%n,j=floor(k/n))

estack itself is of the form

head -next-> edge1 -next-> ...
               |
             val (i+n*j)

The use of the fixed head element allows the stack to be passed around during
recursion without ill effects.  Just remember that the stack's lead pointer
is estack->next.

(Note: ultimately, it would be nice to handle biblocks (the bicomponents in
digraphs) and bijoin points (articulation points in digraphs bijoin point).)
*/
{
  element *es,*es2,*cp;
  slelement *ep;
  int vj,flag,n,vert;

  //Rprintf("\tEntered bicomp recursion for %d.\n",v);
  n=g->n;
  es=estack;
  /*Update counters*/
  num[v]=back[v]=(++(*dfn));
  /*Recursively seek cutpoints*/
  for(ep=snaFirstEdge(g,v,1);ep!=NULL;ep=ep->next[0]){     /*Walk v's neighbors*/
    if(((int)ep->val!=v)&&((int)ep->val!=parent[v])){
      vj=(int)ep->val;
      if(num[vj]==0){                     /*We're seeing vj for the first time*/
        es->next=push(es->next,v+(double)vj*n,NULL);   /*This is a tree edge*/
        parent[vj]=v;
        bicomponentRecurse(g,complist,es,parent,num,back,dfn,vj);
        //Rprintf("Returned from bicomp recursion (%d to %d).\n",v,vj);
        //Rprintf("\tback[%d]=%d, num[%d]=%d\n",vj,back[vj],v,num[v]);
        if(back[vj]>=num[v]){
          /*Create/insert new component*/
          cp=(element *)R_alloc(1,sizeof(element));
          cp->dp=NULL;
          cp->next=NULL;
          cp->val=0.0;
          if((int)complist->val==0){
            complist->next=cp;
          }else{
            ((element *)(complist->dp))->next=cp;
          }
          complist->dp=(void *)cp;
          complist->val++;
          /*Add vertices to component*/
          flag=0;
          //Rprintf("\tFound backedge with (%d,%d); popping stack.\n",v,vj);
          for(es2=es->next;(es2!=NULL)&&(!flag);es2=es2->next){
            //Rprintf("\t\t(%d,%d)\n",(int)(es2->val) % n, (int)floor((es2->val)/((double)n)));
            if(es2->val==(double)v+(double)vj*n) /*Check to see if we've pulled last edge*/
              flag++;
            if(!flag){
              vert=(int)fmod(es2->val,(double)n);
              if(!isinstack((element *)(cp->dp),vert)){
                cp->dp=(void *)listInsert((element *)(cp->dp),vert,NULL);
                cp->val++;
              }
              vert=(int)floor((es2->val)/((double)n));
              if(!isinstack((element *)(cp->dp),vert)){
                cp->dp=(void *)listInsert((element *)(cp->dp),vert,NULL);
                cp->val++;
              }
            }
          }
          es->next=es2; /*Remove popped edges from stack (R will deallocate)*/
        }else{
          //Rprintf("\tSetting back[%d]=%d\n",v,MIN(back[v],back[vj]));
          back[v]=MIN(back[v],back[vj]);
        }
      }else if((num[vj]<num[v])&&(vj!=parent[v])){       /*This is a back edge*/
        //Rprintf("(%d,%d) is a back edge.  Pushing and setting back[%d]=%d\n",v,vj,v,MIN(num[vj],back[v]));
        es->next=push(es->next,(double)v+(double)vj*n,NULL);
        back[v]=MIN(num[vj],back[v]);
      }
    }
  }  
}


slelement *cliqueFirstChild(snaNet *g, slelement *cl)
/*Given the "seed" clique in cl, find the (lexically) first maximal clique
containing cl.  Note that cl is "expanded" into the new clique, and is thus
changed in place.  A pointer to cl is superfluously returned.*/
{
  slelement *i,*ep,*ep2;
     
  /*Check to see if we were called with an empty list*/
  if((cl==NULL)||(cl->val==0.0))
    return cl;
  /*Index the first existing clique entry*/
  ep=cl->next[0];
  /*If the first element is an isolate, return; else, search within the*/
  /*neighbors of the first element (since this is much smaller than n!).*/
  if(g->outdeg[(int)(ep->val)]==0)
    return cl;
  i=g->oel[(int)(ep->val)]->next[0];
  while(i!=NULL){
    /*If i is equal to the currently indexed clique element, advance both*/
    while((ep!=NULL)&&(i->val==ep->val)){
      i=i->next[0];
      ep=ep->next[0];
    }
    /*Otherwise, test for adjacency*/
    if(i!=NULL){
      for(ep2=cl->next[0];(ep2!=NULL)&& snaIsAdjacent((int)(i->val),(int)(ep2->val),g,2); ep2=ep2->next[0]);
      /*If we got to the end, add i->val to the clique*/
      if(ep2==NULL)
        cl=slistInsert(cl,i->val,NULL);
    }
    /*Increment i*/
    i=i->next[0];
  }
  
  /*Return a pointer to the updated clique (should not have changed)*/
  return cl;
}


void cliqueRecurse(snaNet *g, slelement *k, int parind, element **clist, double *ccount, int *compmemb)
/*Recursive clique enumeration, using the algorithm of Makino and Uno (2004) (with some implementation differences).  k should contain (as an slist) the clique to evaluate -- it is stored, and all lexical descendents are then computed.  Cliques are stored in clist, which must be a vector of stacks of length n.  clist[i] is then the stack of cliques of length i+1; the cliques themselves are stored as slists within said stacks.  The vector ccount should likewise be an n-vector, and is used to store the number of cliques of each size so far enumerated.  It should be initialized to zero, obviously.  Finally, parind is the "parent index" of k, which is used in the recursion process.  To enumerate all cliques within a component, call this function with the lexically first clique in the component, and a parind of the lead vertex number.  (Because this algorithm doesn't work properly across components, one must also pass a vector of component memberships.  compmemb is assumed to be such an indicator vector, shifted forward by 1 (so that the ith vertex in *g reckoning is in the i+1th position), as returned by undirComponents.)*/
{
  int i,j,flag,cm;
  slelement *kc,*ep,*ep2;

  /*0th step: see if we need to exit*/
  R_CheckUserInterrupt();
  
  /*First, store k and update ccount*/
  clist[(int)(k->val)-1]=push(clist[(int)(k->val)-1],k->val,(void *)k);
  ccount[(int)(k->val)-1]++;
  
  /*Try finding "child" cliques using all non-k members > parind*/
  cm=compmemb[(int)(k->next[0]->val)+1];
  for(i=parind+1;i<g->n;i++)
    if((cm==compmemb[i+1])&&(!isInSList(k,(double)i))){
      /*Create candidate seed*/
      kc=slistInsert(NULL,(double)i,NULL);  /*i U {k[k<=i] int N(i)}*/
      for(ep=k->next[0];(ep!=NULL)&&(ep->val<=(double)i);ep=ep->next[0])
        if(snaIsAdjacent(i,(int)(ep->val),g,2))
          kc=slistInsert(kc,ep->val,NULL);
      /*Now, test to see if kc produces a child of k*/
      j=flag=0;                             /*"Lemma 3" condition*/
      ep=kc->next[0];
      while((j<i)&&(!flag)){
        /*If j is not in the same component as k, increment*/
        while((cm!=compmemb[j+1])&&(j<i))
          j++;   
        /*If j is equal to a currently indexed seed element, advance both*/
        while((ep!=NULL)&&(j==(int)ep->val)){
          j++;
          ep=ep->next[0];
        }
        /*Test to see if j adjacent to all seed elements*/
        if(j<i){
          flag=1;
          for(ep2=kc->next[0];flag&&(ep2!=NULL);ep2=ep2->next[0])
            if(!snaIsAdjacent(j,(int)(ep2->val),g,2))
              flag--;
        }
        /*Increment j*/
        j++;
      }
      j=0;                                  /*"Lemma 4" condition*/
      ep=k->next[0];
      while((j<i)&&(!flag)){
        /*If j is not in the same component as k, increment*/
        while((cm!=compmemb[j+1])&&(j<i))
          j++;   
        /*If j is equal to a currently indexed k element, advance both*/
        while((ep!=NULL)&&(j==(int)ep->val)){
          j++;
          ep=ep->next[0];
        }
        /*Test to see if j adjacent to all k<=j*/
        if(j<i){
          flag=1;
          for(ep2=k->next[0];flag&&(ep2!=NULL)&&(ep2->val<=(double)j); ep2=ep2->next[0])
            if(!snaIsAdjacent(j,(int)(ep2->val),g,2))
              flag--;
        }
        /*If we got this far, test to see if j adjacent to non-i seeds*/
        if((j<i)&&flag){
          for(ep2=kc->next[0];flag&&(ep2!=NULL);ep2=ep2->next[0])
            if(((int)(ep2->val)!=i)&&(!snaIsAdjacent(j,(int)(ep2->val),g,2)))
              flag--;
        }
        /*Increment j*/
        j++;
      }
      /*If both tests were passed, then create the child clique and recurse*/
      if(!flag){
        kc=cliqueFirstChild(g,kc);  /*Get first clique containing seed nodes*/
        cliqueRecurse(g,kc,i,clist,ccount,compmemb);/*Recurse on kc, parind=i*/
      }
    }
}


void cutpointUndirRecurse(snaNet *g, int *cpstatus, int *minvis, int *visdep, int depth, int v, int src)
/*Perform a depth-first recusion to identify cutpoints in undirected graphs.
This call implements a visit from src to v; for the search root, use src=-1.  On
the initial pass, visdep and mindis should be 0, as should depth and cpstatus.*/
{
  slelement *ep;
  int vj,ccount=0;

  /*Update counters*/
  minvis[v]=visdep[v]=(++depth);
  /*Recursively seek cutpoints*/
  for(ep=snaFirstEdge(g,v,1);ep!=NULL;ep=ep->next[0])  /*Walk v's neighbors*/
    if((int)ep->val!=src){
      vj=(int)ep->val;
      if(visdep[vj]==0){          /*Is vj unvisited?*/
        if((visdep[v]==1)&&(++ccount>1)) /*Root node w/>1 child*/
          cpstatus[v]=1;
        cutpointUndirRecurse(g,cpstatus,minvis,visdep,depth,vj,v);
        minvis[v]=MIN(minvis[v],minvis[vj]);
        if((visdep[v]!=1)&&(minvis[vj]>=visdep[v]))
          cpstatus[v]=1;
      }else
        minvis[v]=MIN(minvis[v],visdep[vj]);  /*If shorter path, note it*/
    }
}


/*R-CALLABLE ROUTINES-------------------------------------------------------*/

SEXP bicomponents_R(SEXP net, SEXP sn, SEXP sm)
{
  snaNet *g;
  int *parent,*num,*back,*dfn,i,j,n,count,pc=0;
  element *complist,*ep,*ep2,*es;
  SEXP bicomps,bcl,memb,outlist;

  /*Coerce what needs coercin'*/
  //Rprintf("Initial coercion\n");
  PROTECT(sn=coerceVector(sn,INTSXP)); pc++;
  PROTECT(sm=coerceVector(sm,INTSXP)); pc++;
  PROTECT(net=coerceVector(net,REALSXP)); pc++;
  n=INTEGER(sn)[0];

  /*Initialize sna internal network*/
  GetRNGstate();
  g=elMatTosnaNet(REAL(net),INTEGER(sn),INTEGER(sm));

  /*Calculate the sorting stat*/
  parent=(int *)R_alloc(n,sizeof(int));
  num=(int *)R_alloc(n,sizeof(int));
  back=(int *)R_alloc(n,sizeof(int));
  dfn=(int *)R_alloc(1,sizeof(int));
  for(i=0;i<n;i++){
    parent[i]=-1;
    num[i]=0;
    back[i]=n+1;
  }
  *dfn=0;

  /*Initialize the component list*/
  complist=(element *)R_alloc(1,sizeof(element));
  complist->val=0.0;
  complist->next=NULL;
  complist->dp=NULL;

  /*Walk the graph components, finding bicomponents*/
  es=(element *)R_alloc(1,sizeof(element));
  for(i=0;i<n;i++)
    if(num[i]==0){
      es->next=NULL;
      bicomponentRecurse(g,complist,es,parent,num,back,dfn,i);
    }

  /*Transfer information from complist to output vector*/
  //Rprintf("Gathering components...\n");
  count=(int)complist->val;
  PROTECT(outlist=allocVector(VECSXP,2)); pc++;
  PROTECT(bicomps=allocVector(VECSXP,count)); pc++;
  PROTECT(memb=allocVector(INTSXP,n)); pc++;
  for(i=0;i<n;i++)  /*Initialize memberships, since some have none*/
    INTEGER(memb)[i]=-1;
  ep=complist->next;
  for(i=0;i<count;i++){
    PROTECT(bcl=allocVector(INTSXP,(int)ep->val));
    j=0;
    for(ep2=(element *)ep->dp;ep2!=NULL;ep2=ep2->next){
      INTEGER(bcl)[j++]=(int)ep2->val+1;
      INTEGER(memb)[(int)ep2->val]=i+1;
    }
    SET_VECTOR_ELT(bicomps,i,bcl);
    UNPROTECT(1);
    ep=ep->next;
  }
  SET_VECTOR_ELT(outlist,0,bicomps); 
  SET_VECTOR_ELT(outlist,1,memb); 

  /*Unprotect and return*/
  PutRNGstate();
  UNPROTECT(pc);
  return outlist;
}


SEXP cliques_R(SEXP net, SEXP sn, SEXP sm, SEXP stabulatebyvert, SEXP scomembership, SEXP senumerate)
/*Maximal clique enumeration as an R-callable (.Call) function.  net should be an sna edgelist (w/n vertices and m/2 edges), and must be pre-symmetrized.  stabulatebyvert should be 0 if no tabulation is to be performed, or 1 for vertex-level tabulation of clique membership.  scomembership should be 0 for no co-membership tabulation, 1 for aggregate vertex-by-vertex tabulation, and 2 for size-by-vertex-by-vertex tabulation.  Finally, senumerate should be 1 iff the enumerated clique list should be returned.  (The current algorithm enumerates them internally, regardless.  This is b/c I am lazy, and didn't fold all of the tabulation tasks into the recursion process.  Life is hard.)*/
{
  int n,tabulate,comemb,enumerate,*gotcomp,*compmemb,i,j,k,maxcsize,pc=0;
  double *ccount,*pccountvec,*pcocliquevec=NULL;
  snaNet *g;
  slelement *sep,*sep2,*k0;
  element **clist,*ep;
  SEXP smaxcsize,ccountvec,outlist,cliquevec=R_NilValue;
  SEXP temp=R_NilValue,sp=R_NilValue,cocliquevec=R_NilValue;

  /*Coerce what needs coercin'*/
  PROTECT(sn=coerceVector(sn,INTSXP)); pc++;
  PROTECT(net=coerceVector(net,REALSXP)); pc++;
  PROTECT(stabulatebyvert=coerceVector(stabulatebyvert,INTSXP)); pc++;
  PROTECT(scomembership=coerceVector(scomembership,INTSXP)); pc++;
  PROTECT(senumerate=coerceVector(senumerate,INTSXP)); pc++;
  n=INTEGER(sn)[0];
  tabulate=INTEGER(stabulatebyvert)[0];
  comemb=INTEGER(scomembership)[0];
  enumerate=INTEGER(senumerate)[0];

  /*Pre-allocate what needs pre-allocatin'*/
  ccount=(double *)R_alloc(n,sizeof(double));
  PROTECT(smaxcsize=allocVector(INTSXP,1)); pc++;
  clist=(element **)R_alloc(n,sizeof(element *));
  for(i=0;i<n;i++){
    ccount[i]=0.0;
    clist[i]=NULL;
  }
    
  /*Initialize sna internal network*/
  GetRNGstate();
  g=elMatTosnaNet(REAL(net),INTEGER(sn),INTEGER(sm));

  /*Calculate the components of g*/
  compmemb=undirComponents(g);

  /*Accumulate cliques across components*/
  gotcomp=(int *)R_alloc(compmemb[0],sizeof(int));
  for(i=0;i<compmemb[0];i++)
    gotcomp[i]=0;
  for(i=0;i<n;i++)                   /*Move through vertices in order*/
    if(!gotcomp[compmemb[i+1]-1]){   /*Take first vertex of each component*/
      gotcomp[compmemb[i+1]-1]++;              /*Mark component as visited*/
      /*Get the first maximal clique in this component*/
      k0=slistInsert(NULL,(double)i,NULL);
      k0=cliqueFirstChild(g,k0);
      /*Recursively enumerate all cliques within the component*/
      cliqueRecurse(g,k0,i,clist,ccount,compmemb);
    }
  PutRNGstate();
  
  /*Find the maximum clique size (to cut down on subsequent memory usage)*/
  INTEGER(smaxcsize)[0]=n+1;
  for(i=n-1;(i>=0)&(INTEGER(smaxcsize)[0]==n+1);i--)
    if(ccount[i]>0.0)
      INTEGER(smaxcsize)[0]=i+1;
  maxcsize=INTEGER(smaxcsize)[0];

  /*Allocate memory for R return value objects*/
  if(tabulate){
    PROTECT(ccountvec=allocVector(REALSXP,maxcsize*(1+n))); pc++;
    for(i=0;i<maxcsize*(1+n);i++)
      REAL(ccountvec)[i]=0.0;
  }else{
    PROTECT(ccountvec=allocVector(REALSXP,maxcsize)); pc++;
    for(i=0;i<maxcsize;i++)
      REAL(ccountvec)[i]=0.0;
  }
  pccountvec=REAL(ccountvec);
  switch(comemb){
    case 0:
      cocliquevec=R_NilValue;
      pcocliquevec=NULL;
      break;
    case 1:
      PROTECT(cocliquevec=allocVector(REALSXP,n*n)); pc++;
      for(i=0;i<n*n;i++)
        REAL(cocliquevec)[i]=0.0;
      pcocliquevec=REAL(cocliquevec);
      break;
    case 2:
      PROTECT(cocliquevec=allocVector(REALSXP,maxcsize*n*n)); pc++;
      for(i=0;i<maxcsize*n*n;i++)
        REAL(cocliquevec)[i]=0.0;
      pcocliquevec=REAL(cocliquevec);
      break;
  }
  if(enumerate){
    PROTECT(cliquevec=allocVector(VECSXP,maxcsize)); pc++;
    for(i=0;i<maxcsize;i++){
      if(ccount[i]==0.0)
        SET_VECTOR_ELT(cliquevec,i,R_NilValue);
      else{
        PROTECT(temp=allocVector(VECSXP,(int)(ccount[i])));
        SET_VECTOR_ELT(cliquevec,i,temp);
        UNPROTECT(1);
      }
    }
  }

  /*Tabulate, enumerate, and other good things*/
  for(i=0;i<maxcsize;i++){
    pccountvec[i+tabulate*maxcsize*n]=ccount[i];
    if(ccount[i]>0.0){
      if(enumerate)
        sp=VECTOR_ELT(cliquevec,i);
      /*Walk through every clique of size i+1*/
      for(j=0,ep=clist[i];ep!=NULL;ep=ep->next){
        if(enumerate)
          PROTECT(temp=allocVector(INTSXP,i+1));
        /*Walk through every clique member*/
        for(k=0,sep=((slelement *)(ep->dp))->next[0];sep!=NULL; sep=sep->next[0]){
          if(enumerate)        /*Add to enumeration list*/
            INTEGER(temp)[k++]=(int)(sep->val)+1;
          if(tabulate)         /*Add to vertex-by-size tabulation*/
            pccountvec[i+maxcsize*((int)(sep->val))]++;
          switch(comemb){      /*Add co-membership information*/
            case 0:                /*Case 0 - do nothing*/
              break;
            case 1:                /*Case 1 - just co-membership*/
              for(sep2=((slelement *)(ep->dp))->next[0];sep2!=sep; sep2=sep2->next[0]){
                pcocliquevec[((int)(sep->val))+n*((int)(sep2->val))]++;
                pcocliquevec[((int)(sep2->val))+n*((int)(sep->val))]++;
              }
              pcocliquevec[((int)(sep->val))+n*((int)(sep->val))]++;
              break;
            case 2:                /*Case 2 - co-membership by size*/
              for(sep2=((slelement *)(ep->dp))->next[0];sep2!=sep; sep2=sep2->next[0]){
                pcocliquevec[i+maxcsize*((int)(sep->val))+ maxcsize*n*((int)(sep2->val))]++;
                pcocliquevec[i+maxcsize*((int)(sep2->val))+ maxcsize*n*((int)(sep->val))]++;
              }
              pcocliquevec[i+maxcsize*((int)(sep->val))+ maxcsize*n*((int)(sep->val))]++;
              break;
          }
        }
        if(enumerate){
          SET_VECTOR_ELT(sp,j++,temp);
          UNPROTECT(1);
        }
      }
    }
  }
  
  /*Prepare and return the results*/
  PROTECT(outlist=allocVector(VECSXP,4)); pc++;
  SET_VECTOR_ELT(outlist,0,smaxcsize);
  SET_VECTOR_ELT(outlist,1,ccountvec);
  SET_VECTOR_ELT(outlist,2,cocliquevec);
  SET_VECTOR_ELT(outlist,3,cliquevec);
  UNPROTECT(pc);
  return outlist;
}


void cutpointsDir_R(double *mat, int *n, int *m, int *cpstatus)
/*Compute (strong) cutpoints in a directed graph.  mat should be the edgelist matrix (of order n), and cpstatus should be a zero-initialized vectors to contain the cutpoint status (0=not a cutpoint, 1=cutpoint).  Lacking a good algorithm, I've used something horribly slow and ugly -- nevertheless, it will get the job done for graphs of typical size.  Although this should work fine with undirected graphs, it will be hideously slow...use the undirected variant wherever possible.*/
{
  snaNet *g;
  int i,j,ccount,ccountwoi,tempideg,tempodeg;
  slelement *sep,*tempiel,*tempoel,**tempentries;

  //Rprintf("Now in cutpointsDir_R.  Setting up snaNet\n");
  /*Initialize sna internal network*/
  GetRNGstate();
  g=elMatTosnaNet(mat,n,m);
  for(i=0;i<*n;i++){
    cpstatus[i]=0;
  }
  
  /*Walk the vertices, finding cutpoints by brute force*/
  ccount=numStrongComponents(g,n);
  //Rprintf("Original number of components: %d\n",ccount);
  for(i=0;i<*n;i++)
    if((g->indeg[i]>0)&&(g->outdeg[i]>0)){  /*Must be internal to a path*/
      //Rprintf("\tEntering with %d\n",i);
      /*Temporarily make i an isolate*/
      //Rprintf("\tMoving out %d's edges\n",i);
      tempideg=g->indeg[i];
      tempodeg=g->outdeg[i];
      tempiel=g->iel[i];
      tempoel=g->oel[i];
      g->indeg[i]=0;
      g->outdeg[i]=0;
      g->iel[i]=NULL;
      g->oel[i]=NULL;
      tempentries=(slelement **)R_alloc(tempideg,sizeof(slelement *));
      //Rprintf("\tMoving out %d edges pointing to %d\n",tempideg,i);
      if(tempiel==NULL)
        sep=NULL;
      else
        sep=tempiel->next[0];
      for(j=0;sep!=NULL;sep=sep->next[0]){  /*Remove edges pointing to i*/
        //Rprintf("\t\t%d, about to do slistDelete\n",j);
        tempentries[j++]=slistDelete(g->oel[(int)(sep->val)],(double)i);
        //Rprintf("\t\tSending vertex is %d\n",(int)(sep->val));
        //Rprintf("\t\t%d, about to do decrement outdegrees\n",j);
        /*Decrement outdegree*/
         //Rprintf("\t\toutdegree is %d\n", g->outdeg[(int)(sep->val)]);
        g->outdeg[(int)(sep->val)]--;
        //Rprintf("\t\tnew outdegree is %d\n", g->outdeg[(int)(sep->val)]);
        //Rprintf("\t%d -> %d [%.1f]\n",(int)(sep->val), (int)(tempentries[j-1]->val), *((double *)(tempentries[j-1]->dp)));
        //Rprintf("\t\tfinished tracer\n");
      }
      /*Recalculate components (told you this was ugly!)*/
      ccountwoi=numStrongComponents(g,n)-1;  /*Remove 1 for i*/
      //Rprintf("\tNumber of components w/out %d: %d\n",i,ccountwoi);
      if(ccountwoi>ccount)
        cpstatus[i]++;
      /*Restore i to its former glory*/
      g->indeg[i]=tempideg;
      g->outdeg[i]=tempodeg;
      g->iel[i]=tempiel;
      g->oel[i]=tempoel;
      //Rprintf("\tRestoring edges to %d\n",i);
      if(tempiel==NULL)
        sep=NULL;
      else
        sep=tempiel->next[0];
      for(j=0;sep!=NULL;sep=sep->next[0]){  /*Restore edges->i*/
        g->oel[(int)(sep->val)]=slistInsert(g->oel[(int)(sep->val)],(double)i, tempentries[j++]->dp);
        /*Increment outdegree*/
        g->outdeg[(int)(sep->val)]++;
        //Rprintf("\t\tnew outdegree is %d\n", g->outdeg[(int)(sep->val)]);
        //Rprintf("\t%d -> %d [%.1f]\n",(int)(sep->val), (int)(tempentries[j-1]->val), *(double*)(tempentries[j-1]->dp));
      }
    }    
  PutRNGstate();
}


void cutpointsUndir_R(double *mat, int *n, int *m, int *cpstatus)
/*Compute cutpoints in an undirected graph.  mat should be edgelist matrix (of order n, w/m edges), and cpstatus should be a zero-initialized vectors to contain the cutpoint status (0=not a cutpoint, 1=cutpoint).  The standard DFS method is used here -- this will _not_ work for non-mutual digraphs, but is much faster than the currently implemented digraph method.  Use this whenever appropriate!*/
{
  snaNet *g;
  int *minvis,*visdep,i;

//  Rprintf("Now in cutpointsUndir_R.  Setting up snaNet\n");
  /*Initialize sna internal network*/
  GetRNGstate();
  g=elMatTosnaNet(mat,n,m);

  /*Initialize cutpoint/visit structures*/
//  Rprintf("Initializing\n");
  minvis=(int *)R_alloc(*n,sizeof(int));
  visdep=(int *)R_alloc(*n,sizeof(int));
  for(i=0;i<*n;i++){
    cpstatus[i]=minvis[i]=visdep[i]=0;
  }

  /*Walk the graph components, finding cutpoints*/
//  Rprintf("Finding cutpoints\n");
  for(i=0;i<*n;i++)
    if(visdep[i]==0)
      cutpointUndirRecurse(g,cpstatus,minvis,visdep,0,i,-1);
//  Rprintf("Returning\n");
  PutRNGstate();
}


void kcores_R(double *mat, int *n, int *m, double *corevec, int *dtype, int *pdiag, int *pigeval)
/*Compute k-cores for an input graph.  Cores to be computed can be based on
indegree (dtype=0), outdegree (dtype=1), or total degree (dtype=2).  Algorithm
used is based on Batagelj and Zaversnik (2002), with some pieces filled in.
It's quite fast -- for large graphs, far more time is spent processing the
input than computing the k-cores!  When processing edge values, igeval
determines whether edge values should be ignored (0) or used (1); missing edges
are not counted in either case.  When diag=1, diagonals are used; else they are
also omitted.*/
{
  int i,j,k,temp,*ord,*nod,diag,igev;
  double *stat;
  snaNet *g;
  slelement *ep;

  diag=*pdiag;
  igev=*pigeval;

  /*Initialize sna internal network*/
  GetRNGstate();
  g=elMatTosnaNet(mat,n,m);
  PutRNGstate();

  /*Calculate the sorting stat*/
  stat=(double *)R_alloc(*n,sizeof(double));
  switch(*dtype){
    case 0:  /*Indegree*/
      for(i=0;i<*n;i++){
        stat[i]=0.0;
        for(ep=snaFirstEdge(g,i,0);ep!=NULL;ep=ep->next[0])
          if(((diag)||(i!=(int)(ep->val)))&&((ep->dp!=NULL)&&(!ISNAN(*(double *)(ep->dp)))))
          stat[i]+= igev ? *((double *)(ep->dp)) : 1.0;
      }
    break;
    case 1:  /*Outdegree*/
      for(i=0;i<*n;i++){
        stat[i]=0.0;
        for(ep=snaFirstEdge(g,i,1);ep!=NULL;ep=ep->next[0])
          if(((diag)||(i!=(int)(ep->val)))&&((ep->dp!=NULL)&&(!ISNAN(*(double *)(ep->dp)))))
          stat[i]+= igev ? *((double *)(ep->dp)) : 1.0;
      }
    break;
    case 2:  /*Total degree*/
      for(i=0;i<*n;i++){
        stat[i]=0.0;
        for(ep=snaFirstEdge(g,i,0);ep!=NULL;ep=ep->next[0])
          if(((diag)||(i!=(int)(ep->val)))&&((ep->dp!=NULL)&&(!ISNAN(*(double *)(ep->dp)))))
          stat[i]+= igev ? *((double *)(ep->dp)) : 1.0;
        for(ep=snaFirstEdge(g,i,1);ep!=NULL;ep=ep->next[0])
          if(((diag)||(i!=(int)(ep->val)))&&((ep->dp!=NULL)&&(!ISNAN(*(double *)(ep->dp)))))
          stat[i]+= igev ? *((double *)(ep->dp)) : 1.0;
      }
    break;
  }
    
  /*Set initial core/order values*/
  ord=(int *)R_alloc(*n,sizeof(int));
  nod=(int *)R_alloc(*n,sizeof(int));
  for(i=0;i<*n;i++){
    corevec[i]=stat[i];
    ord[i]=nod[i]=i;
  }

  /*Heap reminder: i->(2i+1, 2i+2); parent at floor((i-1)/2); root at 0*/
  /*Build a heap, based on the stat vector*/
  for(i=1;i<*n;i++){
    j=i;
    while(j>0){
      k=(int)floor((j-1)/2);    /*Parent node*/
      if(stat[nod[k]]>stat[nod[j]]){ /*Out of order -- swap*/
        temp=nod[k];
        nod[k]=nod[j];
        nod[j]=temp;
        ord[nod[j]]=j;
        ord[nod[k]]=k;
      }
      j=k;                 /*Move to parent*/
    }
  }

  /*Heap test
  for(i=0;i<*n;i++){
    Rprintf("Pos %d (n=%d, s=%.0f, check=%d): ",i,nod[i],stat[nod[i]],ord[nod[i]]==i);
    j=(int)floor((i-1)/2.0);
    if(j>=0)
      Rprintf("Parent %d (n=%d, s=%.0f), ",j,nod[j],stat[nod[j]]);
    else
      Rprintf("No Parent (root), ");
    j=2*i+1;
    if(j<*n)
      Rprintf("Lchild %d (n=%d, s=%.0f), ",j,nod[j],stat[nod[j]]);
    else
      Rprintf("No Lchild, ");
    j=2*i+2;
    if(j<*n)
      Rprintf("Rchild %d (n=%d, s=%.0f)\n",j,nod[j],stat[nod[j]]);
    else
      Rprintf("No Rchild\n");
  }
  */

  /*Now, find the k-cores*/
  for(i=*n-1;i>=0;i--){
    /*Rprintf("Stack currently spans positions 0 to %d.\n",i);*/
    corevec[nod[0]]=stat[nod[0]];  /*Coreness of top element is fixed*/
    /*Rprintf("Pulled min vertex (%d): coreness was %.0f\n",nod[0],corevec[nod[0]]);*/
    /*Swap root w/last element (and re-heap) to remove it*/
    temp=nod[0];
    nod[0]=nod[i];
    nod[i]=temp;
    ord[nod[0]]=0;
    ord[nod[i]]=i;
    j=0;
    while(2*j+1<i){
      k=2*j+1;                                   /*Get first child*/
      if((k<i-1)&&(stat[nod[k+1]]<stat[nod[k]]))   /*Use smaller child node*/
        k++;
      if(stat[nod[k]]<stat[nod[j]]){               /*If child smaller, swap*/
        temp=nod[j];
        nod[j]=nod[k];
        nod[k]=temp;
        ord[nod[j]]=j;
        ord[nod[k]]=k;
      }else
        break;
      j=k;                                   /*Move to child, repeat*/
    }
    /*Having removed top element, adjust its neighbors downward*/
    switch(*dtype){
      case 0:  /*Indegree -> update out-neighbors*/
         /*Rprintf("Reducing indegree of %d outneighbors...\n",g->outdeg[nod[i]]);*/
         for(ep=snaFirstEdge(g,nod[i],1);ep!=NULL;ep=ep->next[0]){
           j=(int)ep->val;
           if(ord[j]<i){                 /*Don't mess with removed nodes!*/
             /*Adjust stat*/
             /*Rprintf("\t%d: %.0f ->",j,stat[j]);*/
             stat[j]=MAX(stat[j]-*((double *)(ep->dp)),corevec[nod[i]]);
             /*Rprintf(" %.0f\n",stat[j]);*/
             /*Percolate heap upward (stat can only go down!)*/
             j=ord[j];
             while(floor((j-1)/2)>=0){
               k=floor((j-1)/2);               /*Parent node*/
               if(stat[nod[k]]>stat[nod[j]]){   /*If parent greater, swap*/
                 temp=nod[j];
                 nod[j]=nod[k];
                 nod[k]=temp;
                 ord[nod[j]]=j;
                 ord[nod[k]]=k;
               }else
                 break;
               j=k;                             /*Repeat w/new parent*/
             }
           }
         }
         break;
      case 1:  /*Outdegree -> update in-neighbors*/
         for(ep=snaFirstEdge(g,nod[i],0);ep!=NULL;ep=ep->next[0]){
           j=(int)ep->val;
           if(ord[j]<i){                 /*Don't mess with removed nodes!*/
             /*Adjust stat*/
             /*Rprintf("\t%d: %.0f ->",j,stat[j]);*/
             stat[j]=MAX(stat[j]-*((double *)(ep->dp)),corevec[nod[i]]);
             /*Rprintf(" %.0f\n",stat[j]);*/
             /*Percolate heap upward (stat can only go down!)*/
             j=ord[j];
             while(floor((j-1)/2)>=0){
               k=floor((j-1)/2);               /*Parent node*/
               if(stat[nod[k]]>stat[nod[j]]){   /*If parent greater, swap*/
                 temp=nod[j];
                 nod[j]=nod[k];
                 nod[k]=temp;
                 ord[nod[j]]=j;
                 ord[nod[k]]=k;
               }else
                 break;
               j=k;                             /*Repeat w/new parent*/
             }
           }
         }
         break;
      case 2:  /*Total degree -> update all neighbors*/
         for(ep=snaFirstEdge(g,nod[i],1);ep!=NULL;ep=ep->next[0]){
           j=(int)ep->val;
           if(ord[j]<i){                 /*Don't mess with removed nodes!*/
             /*Adjust stat*/
             /*Rprintf("\t%d: %.0f ->",j,stat[j]);*/
             stat[j]=MAX(stat[j]-*((double *)(ep->dp)),corevec[nod[i]]);
             /*Rprintf(" %.0f\n",stat[j]);*/
             /*Percolate heap upward (stat can only go down!)*/
             j=ord[j];
             while(floor((j-1)/2)>=0){
               k=floor((j-1)/2);               /*Parent node*/
               if(stat[nod[k]]>stat[nod[j]]){   /*If parent greater, swap*/
                 temp=nod[j];
                 nod[j]=nod[k];
                 nod[k]=temp;
                 ord[nod[j]]=j;
                 ord[nod[k]]=k;
               }else
                 break;
               j=k;                             /*Repeat w/new parent*/
             }
           }
         }
         for(ep=snaFirstEdge(g,nod[i],0);ep!=NULL;ep=ep->next[0]){
           j=(int)ep->val;
           if(ord[j]<i){                 /*Don't mess with removed nodes!*/
             /*Adjust stat*/
             /*Rprintf("\t%d: %.0f ->",j,stat[j]);*/
             stat[j]=MAX(stat[j]-*((double *)(ep->dp)),corevec[nod[i]]);
             /*Rprintf(" %.0f\n",stat[j]);*/
             /*Percolate heap upward (stat can only go down!)*/
             j=ord[j];
             while(floor((j-1)/2)>=0){
               k=floor((j-1)/2);               /*Parent node*/
               if(stat[nod[k]]>stat[nod[j]]){   /*If parent greater, swap*/
                 temp=nod[j];
                 nod[j]=nod[k];
                 nod[k]=temp;
                 ord[nod[j]]=j;
                 ord[nod[k]]=k;
               }else
                 break;
               j=k;                             /*Repeat w/new parent*/
             }
           }
         }
         break;
    }
  }
}
