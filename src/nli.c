/*
######################################################################
#
# nli.c
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 7/18/16
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
#
# Part of the R/sna package
#
# This file contains routines related to the computation of node level
# indices (NLIs).
#
######################################################################
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "nli.h"


SEXP betweenness_R(SEXP mat, SEXP sn, SEXP sm, SEXP smeasure, SEXP sprecomp, SEXP signoreeval, SEXP sgd, SEXP ssigma, SEXP spred)
/*
Compute betweenness (and some related measures) for the network in mat.  If sprecomp==TRUE, then sgd, ssigma, and spred are taken to hold geodesic distances, path counts, and predecessor lists (as returned by geodist_R or geodist_val_r); else, these are computed on the fly (which, BTW, saves memory, though it prohibits reuse).  If signoreevals==TRUE, then edge values are not used when computing paths (irrelevant if called with precomputed geodesics).
*/
{
  int n,i,j,k,wv,precomp,*npred,ignoreeval,measure,pc=0;
  double *gd, *sigma,*bet,*delta;
  element **pred,*w,*v;
  snaNet *g;
  SEXP sbet,lp,vp;

  /*Coerce inputs*/
  PROTECT(mat=coerceVector(mat,REALSXP)); pc++;
  PROTECT(sn=coerceVector(sn,INTSXP)); pc++;
  PROTECT(sm=coerceVector(sm,INTSXP)); pc++;
  PROTECT(sprecomp=coerceVector(sprecomp,INTSXP)); pc++;
  PROTECT(smeasure=coerceVector(smeasure,INTSXP)); pc++;
  PROTECT(signoreeval=coerceVector(signoreeval,INTSXP)); pc++;
  n=INTEGER(sn)[0];
  precomp=INTEGER(sprecomp)[0];
  measure=INTEGER(smeasure)[0];
  ignoreeval=INTEGER(signoreeval)[0];
  if(precomp){
    PROTECT(sgd=coerceVector(sgd,REALSXP)); pc++;
    PROTECT(ssigma=coerceVector(ssigma,REALSXP)); pc++;
  }

  /*Allocate memory*/
  PROTECT(sbet=allocVector(REALSXP,n)); pc++;
  npred=(int *)R_alloc(n,sizeof(int));
  pred=(element **)R_alloc(n,sizeof(element *));
  gd=(double *)R_alloc(n,sizeof(double));
  sigma=(double *)R_alloc(n,sizeof(double));
  delta=(double *)R_alloc(n,sizeof(double));
  bet=REAL(sbet);
  
  /*Set up stuff*/
  GetRNGstate();
  g=elMatTosnaNet(REAL(mat),INTEGER(sn),INTEGER(sm));
  PutRNGstate();
  for(i=0;i<n;i++)
    bet[i]=0.0;

  /*Calculate betweenness*/
  for(i=0;i<n;i++){
    R_CheckUserInterrupt();
    /*Get geodesic information*/
    if(precomp){                   /*Use pre-computed values*/
      lp=VECTOR_ELT(spred,i);
      for(j=0;j<n;j++){
        gd[j]=REAL(sgd)[i+j*n];
        sigma[j]=REAL(ssigma)[i+j*n];
        pred[j]=NULL;
        vp=PROTECT(coerceVector(VECTOR_ELT(lp,j),REALSXP));
        npred[j]=length(vp);
        for(k=npred[j]-1;k>=0;k--)
          pred[j]=push(pred[j],REAL(vp)[k]-1.0,NULL);
        UNPROTECT(1);
      }
    }else{                         /*Compute on the fly*/
      if(ignoreeval)
        spsp(i,g,gd,sigma,pred,npred,1);
      else
        spsp_val(i,g,gd,sigma,pred,npred,1);
    }
    /*Accumulate betweenness incremements*/
    switch(measure){
      case BETSTANDARD:    /*"Standard" form betweenness (a la Freeman)*/
        for(j=0;j<n;j++)
          delta[j]=0.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next)
            delta[(int)(v->val)]+=sigma[(int)(v->val)]/sigma[wv] * (1.0+delta[wv]);
          if(i!=wv)
            bet[wv]+=delta[wv];
        }
        break;
      case BETWENDPTS:   /*Betweenness including endpoints*/
        bet[i]+=npred[i]-1.0;
        for(j=0;j<n;j++)
          delta[j]=0.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next)
            delta[(int)(v->val)]+=sigma[(int)(v->val)]/sigma[wv] * (1.0+delta[wv]);
          if(i!=wv)
            bet[wv]+=delta[wv]+1.0;
        }
        break;
      case BETPROXIMALSRC:   /*Proximal source betweenness*/
        for(j=0;j<n;j++)
          delta[j]=0.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next){
            if((int)(v->val)!=i)
              bet[(int)(v->val)]+=sigma[(int)(v->val)]/sigma[wv];
          }
        }
        break;
      case BETPROXIMALTAR:   /*Proximal target betweenness*/
        for(j=0;j<n;j++)
          delta[j]=0.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next){
            delta[(int)(v->val)]+=sigma[(int)(v->val)]/sigma[wv] * (1.0+delta[wv]);
            if((int)(v->val)==i)
              bet[wv]+=delta[wv];
          }
        }
        break;
      case BETPROXIMALSUM:  /*Total proximal betweenness*/
        for(j=0;j<n;j++)
          delta[j]=0.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next){
            delta[(int)(v->val)]+=sigma[(int)(v->val)]/sigma[wv] * (1.0+delta[wv]);
            if((int)(v->val)!=i)
              bet[(int)(v->val)]+=sigma[(int)(v->val)]/sigma[wv];
            else
              bet[wv]+=delta[wv];
          }
        }
        break;
      case BETLENSCALED:   /*Length-scaled betweenness*/
        for(j=0;j<n;j++)
          delta[j]=0.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next)
            delta[(int)(v->val)]+=sigma[(int)(v->val)]/sigma[wv] * (1.0/gd[wv]+delta[wv]);
          if(i!=wv)
            bet[wv]+=delta[wv];
        }
        break;
      case BETLINSCALED:   /*Linearly-scaled betweenness*/
        for(j=0;j<n;j++)
          delta[j]=0.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next)
            delta[(int)(v->val)]+=sigma[(int)(v->val)]/sigma[wv] * (1.0/gd[wv]+delta[wv]);
          if(i!=wv)
            bet[wv]+=gd[wv]*delta[wv];
        }
        break;
      case BETSTRESS:   /*Shimbel's stress centrality (not betweenness!)*/
        for(j=0;j<n;j++)
          delta[j]=0.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next)
            delta[(int)(v->val)]+=1.0+delta[wv];
          if(i!=wv)
            bet[wv]+=sigma[wv]*delta[wv];
        }
        break;
      case BETLOAD:    /*Goh's load centrality (must be given transpose graph)*/
        for(j=0;j<n;j++)
          delta[j]=1.0;
        while(pred[i]!=NULL){
          w=pred[i];
          wv=(int)(w->val);
          pred[i]=pred[i]->next;
          for(v=pred[wv];v!=NULL;v=v->next)
            delta[(int)(v->val)]+=delta[wv]/(double)npred[wv];
          bet[wv]+=delta[wv];
        }
        break;
    }
  }

  /*Unprotect and return*/
  UNPROTECT(pc);
  return sbet;
}


void degree_R(double *g, int *pm, int *cmode, int *diag, int *igeval, double *d)
/*
Compute degree centralities for the graph in g (assumed to be in sna edgelist form).  cmode should be 0 for indegree, 1 for outdegree, and 2 for total degree, with diag=1 if the diagonal should be considered, and igeval=1 if edge values should be ignored.  Missing edges are omitted in the sum (and corresponding edge counts, where edge values are ignored).  The centrality scores themselves are stored in d, which is assumed to be initialized to 0.
*/
{
  int m=*pm,i;
  
  for(i=0;i<m;i++)
    if(!ISNAN(g[i+2*m])){
      if(g[i]==g[i+m]){
        if(*diag){
          d[(int)g[i]-1]+= (*igeval) ? 1.0 : g[i+2*m];
        }
      }else{
        switch(*cmode){
          case 0:             /*Indegree*/
            d[(int)g[i+m]-1]+= (*igeval) ? 1.0 : g[i+2*m];
            break;
          case 1:             /*Outdegree*/
            d[(int)g[i]-1]+= (*igeval) ? 1.0 : g[i+2*m];
            break;
          case 2:             /*Total Degree*/
            d[(int)g[i]-1]+= (*igeval) ? 1.0 : g[i+2*m];
            d[(int)g[i+m]-1]+= (*igeval) ? 1.0 : g[i+2*m];
            break;
        }
      }
    }
}


void evcent_R(double *mat, int *n, int *m, double *ev, double *tol, int *maxiter, int *checkna, int *ignoreeval)
/*Compute the eigenvector centrality for mat using the power method.  Should
be faster than eigen() for large, sparse graphs (although not perhaps as 
stable).  The convergence tolerance (L2 norm for successive differences) is
given in tol.  Missing data behavior is controlled by checkna; if checkna==0,
missing data is ignored (i.e., edges are counted regardless), while checkna==1
or 2 results in missing edges being removed.  Likewise, edge values are ignored 
if ignoreeval==1 (i.e., all edges are treated as 1.0).
*/
{
  snaNet *g;
  int i,j;
  double diff,norm,*ev2;
  slelement *ep;
    
  /*Set things up*/  
  GetRNGstate();
  g=elMatTosnaNet(mat,n,m);
  PutRNGstate();
  ev2=(double *)R_alloc(g->n,sizeof(double));
  for(i=0;i<*n;i++)
    ev[i]=1.0/sqrt((double)(g->n));
    
  /*Iterate until convergence*/
  diff=1.0;
  j=0;
  while((sqrt(diff)>(*tol))&&(j<(*maxiter))){
    j++;
    R_CheckUserInterrupt();
    /*Calculate unnormalized values*/
    for(i=0;i<*n;i++){
      ev2[i]=0.0;
      for(ep=snaFirstEdge(g,i,1);ep!=NULL;ep=ep->next[0])
        if((!(*checkna))||((ep->dp!=NULL)&&(!ISNAN(*(double *)(ep->dp))))){
          if(*ignoreeval)
            ev2[i]+=ev[(int)(ep->val)];
          else
            ev2[i]+=(*(double *)(ep->dp))*ev[(int)(ep->val)];
        }
    }
    /*Obtain norm*/
    norm=0.0;
    for(i=0;i<*n;i++)
      norm+=ev2[i]*ev2[i];
    norm=sqrt(norm);
    /*Store normalized values*/
    diff=0.0;
    for(i=0;i<*n;i++){
      ev2[i]/=norm;
      diff+=(ev[i]-ev2[i])*(ev[i]-ev2[i]);
      ev[i]=ev2[i];
    }
  }
  if(j==(*maxiter))
    warning("Maximum iterations exceeded in evcent_R without convergence.  This matrix may be pathological - increase maxiter or try eigen().\n");
}


void gilschmidt_R(double *mat, int *n, int *m, double *scores, int *normalize)
{
  snaNet *g;
  double *gd,*sigma;
  element **pred,*ptr;
  int i,*npred;
    
  /*Initialize sna internal network*/
  GetRNGstate();
  g=elMatTosnaNet(mat,n,m);
  PutRNGstate();

  /*Initialize various other things*/
  gd=(double *)R_alloc(*n,sizeof(double));
  sigma=(double *)R_alloc(*n,sizeof(double));
  pred=(element **)R_alloc(*n,sizeof(element *));
  npred=(int *)R_alloc(*n,sizeof(int));

  /*Now, find the unnormalized GS score for each vertex*/
  for(i=0;i<*n;i++){
    scores[i]=0.0;
    spsp(i,g,gd,sigma,pred,npred,0);
    for(ptr=pred[i];ptr!=NULL;ptr=ptr->next)  /*Walk those ego can reach*/
      if((int)ptr->val!=i)
        scores[i]+=1.0/gd[(int)ptr->val];     /*Increment as 1/dist*/
    if(*normalize)
      scores[i]/=npred[i]-1.0;
  }
}


void stresscent_R(double *g, double *pn, double *stress, double *gd, 
double *sigma)
/*
Compute stress centrality for the graph in g.  Geodesic distances are assumed to
have been stored in gd, and the path counts in sigma (both being nxn matrices).
It is also assumed that stress has been initialized to 0.  (Note that this is
now a legacy routine, and is no longer used -- see betweenness_R above.)
*/
{
  long int n,i,j,k;

  /*Set up stuff*/
  n=*pn;
  /*Cycle through each triad, accumulating paths*/
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      for(k=0;k<n;k++){
        if((j!=i)&&(k!=i)&&(j!=k)&&(gd[j+k*n]>=gd[j+i*n]+gd[i+k*n]))
          stress[i]+=sigma[j+i*n]*sigma[i+k*n];
      }
    }
  }
}

