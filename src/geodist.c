/*
######################################################################
#
# geodist.c
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 7/15/10
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
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


void spsp(int ego, snaNet *g, double *gd, double *sigma, element **pred, int *npred, int checkna)
/*
Solve the single-pair shortest-path problem for ego the graph in g.  The
geodesic distances are stored in gd, and the path counts in sigma (both being
n-vectors).  Predecessors are contained as a vector of lists in pred, with 
ego's entry instead containing the vertices encountered (as a stack).  Treatment
of missing data is determined by checkna; 0 implies no NA checking (missing 
edges treated as present), and 1 or 2 results in omission of missing edges.
*/
{
  element *tovisit,*v,*last,*ep;
  slelement *w;
  int n,i,vv,wv;

  /*Set up stuff*/
  n=g->n;
  for(i=0;i<n;i++){
    gd[i]=R_PosInf;
    sigma[i]=0.0;
    pred[i]=NULL;
    npred[i]=0;
  }
  
  /*Solve the shortest path problem*/
  tovisit=enqueue(NULL,(double)ego,NULL);
  last=tovisit;
  sigma[ego]=1.0;
  gd[ego]=0.0;
  while(tovisit!=NULL){
    /*Pull the first element from the queue*/
    v=tovisit;
    vv=(int)(v->val);
    if(last==tovisit)
      last=NULL;
    tovisit=v->next;
    npred[ego]++;
    pred[ego]=push(pred[ego],(double)vv,NULL);
    /*Walk the out-neighborhood of v*/
    for(w=snaFirstEdge(g,vv,1);w!=NULL;w=w->next[0])
      if((!checkna)||((w->dp!=NULL)&&(!ISNAN(*(double *)(w->dp))))){
        wv=(int)(w->val);
        if(gd[wv]==R_PosInf){
          gd[wv]=gd[vv]+1.0;
          /*Insert at the end using a custom adjustment*/
          ep=(element *)R_alloc(1,sizeof(element));
          ep->val=w->val;
          ep->dp=NULL;
          ep->next=NULL;
          if(last!=NULL)
            last->next=ep;
          else
            tovisit=ep;
          last=ep;
        }
        if(gd[wv]==gd[vv]+1.0){
          sigma[wv]+=sigma[vv];
          pred[wv]=push(pred[wv],(double)vv, NULL);
          npred[wv]++;
        }
      }
  }
}


void spsp_val(int ego, snaNet *g, double *gd, double *sigma, element **pred, int *npred, int checkna)
/*
Compute a single-pair shortest-path solution for the valued graph in g.  The
geodesic distances are stored in gd, and the path counts in sigma (both being
n-vectors).  Predecessors are stored as an n-vector of lists (pred), with 
ego's entry instead containing the vertices encountered (as a stack).  
Treatment of missing data is determined by checkna; 0 implies no NA checking 
(missing edges treated as present), and 1 or 2 results in omission of missing 
edges.
*/
{
  element *tovisit,*v,*ep,*ep2;
  slelement *w;
  int n,i,*x,vv,wv;
  double ev;
  
  n=g->n; 
  for(i=0;i<n;i++){
    gd[i]=R_PosInf;
    sigma[i]=0.0;
    npred[i]=0;
    pred[i]=NULL;
  }
  
  /*Solve the shortest path problem for ego*/
  x=(int *)R_alloc(1,sizeof(int));
  x[0]=ego;
  tovisit=listInsert(NULL,0.0,(void *)x);
  gd[ego]=0.0;
  sigma[ego]=1.0;
  while(tovisit!=NULL){
    /*Pull the first element from the queue*/
    v=tovisit;
    tovisit=v->next;
    vv=*((int *)(v->dp));
    pred[ego]=push(pred[ego],(double)vv,NULL);
    npred[ego]++;
    /*Walk the out-neighborhood of v*/
    for(w=snaFirstEdge(g,vv,1);w!=NULL;w=w->next[0])
      if((!checkna)||((w->dp!=NULL)&&(!ISNAN(*(double *)(w->dp))))){
        ev=*((double *)(w->dp));
        wv=(int)(w->val);
        if(gd[wv]>gd[vv]+ev){
          /*Set new shortest distance*/
          gd[wv]=gd[vv]+ev;
          /*Reset sigma and pred*/
          sigma[wv]=0.0;
          npred[wv]=0;
          pred[wv]=NULL;
          /*If w not in queue, add it; else, update its position*/
          if(tovisit==NULL){
            x=(int *)R_alloc(1,sizeof(int));
            x[0]=wv;
            tovisit=listInsert(tovisit,gd[wv],(void *)x);
          }else if(tovisit->next==NULL){
            if(*(int *)(tovisit->dp)==wv){
              tovisit->val=gd[wv];
            }else{
              x=(int *)R_alloc(1,sizeof(int));
              x[0]=wv;
              tovisit=listInsert(tovisit,gd[wv],(void *)x);
            }
          }else{
            /*Look for w*/
            for(ep=tovisit; (ep->next!=NULL) && ((*(int *)(ep->next->dp))!=wv); ep=ep->next);
            /*If w not in queue, add it in order*/
            if(ep->next==NULL){
              x=(int *)R_alloc(1,sizeof(int));
              x[0]=wv;
              tovisit=listInsert(tovisit,gd[wv],(void *)x);
            }else{            /*Else, update and re-insert*/
              ep2=ep->next;
              ep->next=ep2->next;
              ep2->val=gd[wv];
              if(ep2->val<=tovisit->val){
                ep2->next=tovisit;
                tovisit=ep2;
              }else{
                for(ep=tovisit;(ep->next!=NULL)&&(ep->next->val<ep2->val); ep=ep->next);
                if(ep->next==NULL){
                  ep2->next=NULL;
                  ep->next=ep2;
                }else{
                  ep2->next=ep->next;
                  ep->next=ep2;
                }
              }
            }
          }
        }
        /*Increment the number of shortest paths, if needed*/
        if(gd[wv]==gd[vv]+ev){
          sigma[wv]+=sigma[vv];
          npred[wv]++;
          pred[wv]=push(pred[wv],(double)vv,NULL);
        }
    }
  }
}


void geodist_adj_R(double *g, double *pn, double *gd, double *sigma)
/*
Compute geodesics for the adjacency matrix in g.  The geodesic distances are
stored in gd, and the path counts in sigma (both being nxn matrices).  Note that
these should be initialized to all infs and all 0s, respectively.  (This 
function is a holdover, as it processes the network in adjacency form.  I'm
leaving it for now, however, in case it is needed for something.)
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


SEXP geodist_R(SEXP mat, SEXP sn, SEXP sm, SEXP scheckna, SEXP scalcsig, SEXP scalcpred)
/*
Compute geodesics for the graph in mat.  The results are returned as a list of objects: an nxn distance matrix (gd), an nxn matrix of path counts (sigma), and a list of lists of predecessor vectors.  Calculation of the latter two can be suppressed via use of the scalcsig and scalcpred arguments (respectively).  Treatment of missing data is determined by checkna; 0 implies no NA checking (missing edges treated as present), and 1 or 2 results in omission of missing edges.

By the way, the various functions in this file should be rationalized -- geodist_R and geodist_val_R could be consolidated, and could serve as wrappers for spsp and spsp_val.  Expect this to happen in a future release.
*/
{
  snaNet *g;
  element *tovisit,*v,*last,*ep,**pred=NULL;
  slelement *w;
  int n,i,j,k,checkna,calcpred,calcsig,*npred=NULL,vv,wv,pc=0;
  double *gd,*sigma=NULL;
  SEXP sgd,ssigma=R_NilValue,allpredl=R_NilValue,predl=R_NilValue,outlist,pl;
  const void *vmax;

  /*Coerce inputs*/
  PROTECT(mat=coerceVector(mat,REALSXP)); pc++;
  PROTECT(sn=coerceVector(sn,INTSXP)); pc++;
  PROTECT(sm=coerceVector(sm,INTSXP)); pc++;
  PROTECT(scheckna=coerceVector(scheckna,INTSXP)); pc++;
  PROTECT(scalcpred=coerceVector(scalcpred,INTSXP)); pc++;
  PROTECT(scalcsig=coerceVector(scalcsig,INTSXP)); pc++;
  checkna=INTEGER(scheckna)[0];
  calcpred=INTEGER(scalcpred)[0];
  calcsig=INTEGER(scalcsig)[0];
  n=INTEGER(sn)[0];

  /*Allocate memory for outputs*/
  PROTECT(sgd=allocVector(REALSXP,n*n)); pc++;
  gd=REAL(sgd);
  if(calcsig){
    PROTECT(ssigma=allocVector(REALSXP,n*n)); pc++;
    sigma=REAL(ssigma);
  }
  if(calcpred){
    PROTECT(allpredl=allocVector(VECSXP,n)); pc++;
    pred=(element **)R_alloc(n,sizeof(element *));
    npred=(int *)R_alloc(n,sizeof(int));
  }
  
  /*Set up stuff*/
  GetRNGstate();
  g=elMatTosnaNet(REAL(mat),INTEGER(sn),INTEGER(sm));
  PutRNGstate();
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      gd[i+n*j]=R_PosInf;
      if(calcsig)
        sigma[i+n*j]=0.0;
    }
  }
  /*Solve the shortest path problem for each vertex*/
  for(i=0;i<n;i++){
    R_CheckUserInterrupt();
    if(calcpred){
      for(j=0;j<n;j++){
        pred[j]=NULL;
        npred[j]=0;
      }
    }
    vmax=vmaxget();   /*Set memory threshold, so that we can recycle later*/
    tovisit=enqueue(NULL,(double)i,NULL);
    last=tovisit;
    gd[i+i*n]=0.0;
    if(calcsig)
      sigma[i+i*n]=1.0;
    while(tovisit!=NULL){
      /*Pull the first element from the queue*/
      v=tovisit;
      vv=(int)(v->val);
      if(last==tovisit)
        last=NULL;
      tovisit=v->next;
      if(calcpred){
        npred[i]++;
        pred[i]=push(pred[i],(double)vv,NULL);
      }
      /*Walk the out-neighborhood of v*/
      for(w=snaFirstEdge(g,vv,1);w!=NULL;w=w->next[0])
        if((!checkna)||((w->dp!=NULL)&&(!ISNAN(*(double *)(w->dp))))){
          wv=(int)(w->val);
          if(gd[i+n*wv]==R_PosInf){
            gd[i+n*wv]=gd[i+n*vv]+1.0;
            /*Insert at the end using a custom adjustment*/
            ep=(element *)R_alloc(1,sizeof(element));
            ep->val=w->val;
            ep->dp=NULL;
            ep->next=NULL;
            if(last!=NULL)
              last->next=ep;
            else
              tovisit=ep;
            last=ep;
          }
          if(gd[i+n*wv]==gd[i+n*vv]+1.0){
            if(calcsig)
              sigma[i+n*wv]+=sigma[i+n*vv];
            if(calcpred){
              pred[wv]=enqueue(pred[wv],(double)vv, NULL);
              npred[wv]++;
            }
          }
        }
    }
    /*Store predecessory lists if collecting*/
    if(calcpred){
      PROTECT(predl=allocVector(VECSXP,n));
      for(j=0;j<n;j++){
        if(npred[j]>0){
          PROTECT(pl=allocVector(INTSXP,npred[j]));
          for(k=0,ep=pred[j];ep!=NULL;ep=ep->next)
            INTEGER(pl)[k++]=(int)(ep->val)+1;
          SET_VECTOR_ELT(predl,j,pl);
          UNPROTECT(1);
        }else
          SET_VECTOR_ELT(predl,j,R_NilValue);
      }
      SET_VECTOR_ELT(allpredl,i,predl);
      UNPROTECT(1);
    }
    /*Unprotect locally allocated memory*/
    vmaxset(vmax);
  }
  
  /*Prepare and return the results*/
  PROTECT(outlist=allocVector(VECSXP,3)); pc++;
  SET_VECTOR_ELT(outlist,0,sgd);
  if(calcsig)
    SET_VECTOR_ELT(outlist,1,ssigma);
  else
    SET_VECTOR_ELT(outlist,1,R_NilValue);
  if(calcsig)
    SET_VECTOR_ELT(outlist,2,allpredl);
  else
    SET_VECTOR_ELT(outlist,2,R_NilValue);
  UNPROTECT(pc);
  return outlist;
}


SEXP geodist_val_R(SEXP mat, SEXP sn, SEXP sm, SEXP scheckna, SEXP scalcsig, SEXP scalcpred)
/*
Compute geodesics for the valued graph in mat.  The results are returned as a list of objects: an nxn distance matrix (gd), an nxn matrix of path counts (sigma), and a list of lists of predecessor vectors.  Calculation of the latter two can be suppressed via use of the scalcsig and scalcpred arguments (respectively).  Treatment of missing data is determined by checkna; 0 implies no NA checking (missing edges treated as present), and 1 or 2 results in omission of missing edges.
*/
{
  snaNet *g;
  element *tovisit,*v,*ep,*ep2,**pred=NULL;
  slelement *w;
  int n,i,j,k,*x,vv,wv,*npred=NULL,pc=0,checkna,calcpred,calcsig;
  double ev,*gd,*sigma=NULL;
  SEXP sgd,ssigma=R_NilValue,allpredl=R_NilValue,predl=R_NilValue,outlist,pl;
  const void *vmax;
  
  /*Coerce inputs*/
  PROTECT(mat=coerceVector(mat,REALSXP)); pc++;
  PROTECT(sn=coerceVector(sn,INTSXP)); pc++;
  PROTECT(sm=coerceVector(sm,INTSXP)); pc++;
  PROTECT(scheckna=coerceVector(scheckna,INTSXP)); pc++;
  PROTECT(scalcpred=coerceVector(scalcpred,INTSXP)); pc++;
  PROTECT(scalcsig=coerceVector(scalcsig,INTSXP)); pc++;
  checkna=INTEGER(scheckna)[0];
  calcpred=INTEGER(scalcpred)[0];
  calcsig=INTEGER(scalcsig)[0];
  n=INTEGER(sn)[0];

  /*Allocate memory for outputs*/
  PROTECT(sgd=allocVector(REALSXP,n*n)); pc++;
  gd=REAL(sgd);
  if(calcsig){
    PROTECT(ssigma=allocVector(REALSXP,n*n)); pc++;
    sigma=REAL(ssigma);
  }
  if(calcpred){
    PROTECT(allpredl=allocVector(VECSXP,n)); pc++;
    pred=(element **)R_alloc(n,sizeof(element *));
    npred=(int *)R_alloc(n,sizeof(int));
  }
  
  /*Set up stuff*/
  GetRNGstate();
  g=elMatTosnaNet(REAL(mat),INTEGER(sn),INTEGER(sm));
  PutRNGstate();
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      gd[i+n*j]=R_PosInf;
      if(calcsig)
        sigma[i+n*j]=0.0;
    }
  }
  
  /*Solve the shortest path problem for each vertex*/
  for(i=0;i<n;i++){
    R_CheckUserInterrupt();
    vmax=vmaxget();   /*Set memory threshold, so that we can recycle later*/
    x=(int *)R_alloc(1,sizeof(int));
    x[0]=i;
    tovisit=listInsert(NULL,0.0,(void *)x);
    gd[i+i*n]=0.0;
    if(calcsig)
      sigma[i+i*n]=1.0;
    while(tovisit!=NULL){
      /*Pull the first element from the queue*/
      v=tovisit;
      tovisit=v->next;
      vv=*((int *)(v->dp));
      if(calcpred){
        npred[i]++;
        pred[i]=push(pred[i],(double)vv,NULL);
      }
      /*Walk the out-neighborhood of v*/
      for(w=snaFirstEdge(g,vv,1);w!=NULL;w=w->next[0])
        if((!checkna)||((w->dp!=NULL)&&(!ISNAN(*(double *)(w->dp))))){
          ev=*((double *)(w->dp));
          wv=(int)(w->val);
          if(gd[i+n*wv]>gd[i+n*vv]+ev){
            /*Set new shortest distance*/
            gd[i+n*wv]=gd[i+n*vv]+ev;
            /*Reset sigma and predecessor lists, if needed*/
            if(calcsig)
              sigma[i+n*wv]=0.0;
            if(calcpred){
              npred[wv]=0;
              pred[wv]=NULL;
            }
            /*If w not in queue, add it; else, update its position*/
            if(tovisit==NULL){
              x=(int *)R_alloc(1,sizeof(int));
              x[0]=wv;
              tovisit=listInsert(tovisit,gd[i+n*wv],(void *)x);
            }else if(tovisit->next==NULL){
              if(*(int *)(tovisit->dp)==wv){
                tovisit->val=gd[i+n*wv];
              }else{
                x=(int *)R_alloc(1,sizeof(int));
                x[0]=wv;
                tovisit=listInsert(tovisit,gd[i+n*wv],(void *)x);
              }
            }else{
              /*Look for w*/
              for(ep=tovisit; (ep->next!=NULL) && ((*(int *)(ep->next->dp))!=wv); ep=ep->next);
              /*If w not in queue, add it in order*/
              if(ep->next==NULL){
                x=(int *)R_alloc(1,sizeof(int));
                x[0]=wv;
                tovisit=listInsert(tovisit,gd[i+n*wv],(void *)x);
              }else{            /*Else, update and re-insert*/
                ep2=ep->next;
                ep->next=ep2->next;
                ep2->val=gd[i+n*wv];
                if(ep2->val<=tovisit->val){
                  ep2->next=tovisit;
                  tovisit=ep2;
                }else{
                  for(ep=tovisit;(ep->next!=NULL)&&(ep->next->val<ep2->val); ep=ep->next);
                  if(ep->next==NULL){
                    ep2->next=NULL;
                    ep->next=ep2;
                  }else{
                    ep2->next=ep->next;
                    ep->next=ep2;
                  }
                }
              }
            }
          }
          /*Increment the number of shortest paths, if needed*/
          if(gd[i+n*wv]==gd[i+n*vv]+ev){
            if(calcsig)
              sigma[i+n*wv]+=sigma[i+n*vv];
            if(calcpred){
              npred[wv]++;
              pred[wv]=push(pred[wv],(double)vv,NULL);
            }
          }
        }
    }
    /*Store predecessory lists if collecting*/
    if(calcpred){
      PROTECT(predl=allocVector(VECSXP,n));
      for(j=0;j<n;j++){
        if(npred[j]>0){
          PROTECT(pl=allocVector(INTSXP,npred[j]));
          for(k=0,ep=pred[j];ep!=NULL;ep=ep->next)
            INTEGER(pl)[k++]=(int)(ep->val)+1;
          SET_VECTOR_ELT(predl,j,pl);
          UNPROTECT(1);
        }else
          SET_VECTOR_ELT(predl,j,R_NilValue);
      }
      SET_VECTOR_ELT(allpredl,i,predl);
      UNPROTECT(1);
    }
    /*Unprotect locally allocated memory*/
    vmaxset(vmax);
  }
  
  /*Prepare and return the results*/
  PROTECT(outlist=allocVector(VECSXP,3)); pc++;
  SET_VECTOR_ELT(outlist,0,sgd);
  if(calcsig)
    SET_VECTOR_ELT(outlist,1,ssigma);
  else
    SET_VECTOR_ELT(outlist,1,R_NilValue);
  if(calcsig)
    SET_VECTOR_ELT(outlist,2,allpredl);
  else
    SET_VECTOR_ELT(outlist,2,R_NilValue);
  UNPROTECT(pc);
  return outlist;
}


void maxflow_EK_R(double *g,int *pn,int *psource,int *psink,double *flow)
/*
Determine the maxmimum flow from source to sink using the Edmonds-Karp 
algorithm.  (This implementation is an adaptation of one provided on Wikipedia
(entry: "Edmonds-Karp algorithm," 4/18/09), for what it's worth.)
*/
{
  int i,j,t,n,p,q,source,sink,flag,*pre,*que;
  double *fmat,*d,f;
  
  n=*pn;
  source=*psource;
  sink=*psink;
  
  /*Rprintf("Entered with source %d, sink %d\n",source,sink);*/
  if(source==sink)   /*If source==sink, just set to infinity and exit*/
    *flow=R_PosInf;
  else{              /*Alas, looks like we'll have to do some work!*/
    /*Initialize everything*/
    fmat=(double *)R_alloc(n*n,sizeof(double));
    pre=(int *)R_alloc(n,sizeof(int));
    que=(int *)R_alloc(n,sizeof(int));
    d=(double *)R_alloc(n,sizeof(double));
    for(i=0;i<n;i++){
      que[i]=source;
      for(j=0;j<n;j++)
        fmat[i+j*n]=0.0;
    }
    flag=0;
    /*Calculate the flow*/
    while(!flag){
      R_CheckUserInterrupt();
      /*Rprintf("Iterating...\n");
      for(i=0;i<n;i++){
        for(j=0;j<n;j++)
          Rprintf("%.1f ",fmat[i+j*n]);
        Rprintf("\n");
      }*/
      for(i=0;i<n;i++){
        pre[i]=0;
        que[i]=source;
      }
      t=source;
      pre[t]=source+1;
      d[t]=R_PosInf;
      for(p=q=0;(p<=q)&&(!pre[sink]);t=que[p++]){
        /*Rprintf("\tp=%d, q=%d\n\t",p,q);
        for(i=0;i<n;i++)
          Rprintf("%d ",pre[i]);
        Rprintf("\n\t");
        for(i=0;i<n;i++)
          Rprintf("%d ",que[i]);
        Rprintf("\n");*/
        for(i=0;i<n;i++){
          if((!pre[i])&&(f=g[t+i*n]-fmat[t+i*n])){
            que[q++]=i;
            pre[i]=t+1;
            d[i] = d[t]<f ? d[t] : f;
          }else{
            if((!pre[i])&&(f=fmat[i+t*n])){
              que[q++]=i;
              pre[i]=-t-1;
              d[i] = d[t]<f ? d[t] : f;
            }
          }
        }
      }
      if(!pre[sink])
        flag++;
      if(!flag){
        for(i=sink;i!=source;){
          /*Rprintf("i=%d, pre[%d]=%d\n",i,i,pre[i]);*/
          if(pre[i]>0){
            fmat[pre[i]-1+i*n]+=d[sink];
            i=pre[i]-1; 
          }else{
            fmat[i+(-pre[i]-1)*n]-=d[sink];
            i=-pre[i]-1;
          }
        }
      }
    }
    for(f=0.0,i=0;i<n;i++)
      f+=fmat[source+i*n];
    *flow=f;
  }
}
