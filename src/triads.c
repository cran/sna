/*
######################################################################
#
# triads.c
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 2/27/13
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
#
# Part of the R/sna package
#
# This file contains routines related to the classification and
# counting of triads.
#
######################################################################
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "triads.h"

int triad_classify(int *g, int gn, int i, int j, int k, int gm)
/*
If gm=1, compute a Holland and Leinhardt classification for the {i,j,k} triad
in graph g; otherwise, compute the four-state undirected classification.  Note
that this routine assumes dichotomized data.  This function assumes that g is encoded in adjacency matrix form.

This routine is not intended to be called from R.
*/
{
  int m,a,n,di,dj,dk;
  
  /*If we are in the undirected case, take the easy way out*/
  if(!gm)
    return g[i+j*gn]+g[j+k*gn]+g[i+k*gn];
      
  /*Get MAN information*/
  m=g[i+j*gn]*g[j+i*gn]+g[i+k*gn]*g[k+i*gn]+g[j+k*gn]*g[k+j*gn];
  n=(1-g[i+j*gn])*(1-g[j+i*gn])+(1-g[i+k*gn])*(1-g[k+i*gn])+(1-g[j+k*gn])*(1-g[k+j*gn]);
  a=3-m-n;

  /*Now classify, using dyad census as a first cut*/
  if(n==3)                    /*003*/
    return 0;
  else if((a==1)&&(n==2))     /*012*/
    return 1; 
  else if((m==1)&&(n==2))     /*102*/
    return 2; 
  else if((a==2)&&(n==1)){    /*021*/
    di=g[i+j*gn]+g[i+k*gn];
    if(di==2)
      return 3;                 /*021D*/
    dj=g[j+i*gn]+g[j+k*gn];
    if(dj==2)
      return 3;                 /*021D*/
    dk=g[k+i*gn]+g[k+j*gn];
    if(dk==2)
      return 3;                 /*021D*/
    di=g[j+i*gn]+g[k+i*gn];
    if(di==2)
      return 4;                 /*021U*/
    dj=g[i+j*gn]+g[k+j*gn];
    if(dj==2)
      return 4;                 /*021U*/
    dk=g[i+k*gn]+g[j+k*gn];
    if(dk==2)
      return 4;                 /*021U*/
    return 5;                   /*021C*/
  }else if((m==1)&&(n==1)){   /*111*/
    di=g[j+i*gn]+g[k+i*gn];
    if((di==0)||(di==2))
      return 6;                 /*111D*/
    dj=g[i+j*gn]+g[k+j*gn];
    if((dj==0)||(dj==2))
      return 6;                 /*111D*/
    return 7;                   /*111U*/
  }else if(a==3){            /*030*/
    di=g[j+i*gn]+g[k+i*gn];
    if((di==2)||(di==0))
      return 8;                 /*030T*/ 
    dj=g[i+j*gn]+g[k+j*gn];
    if((dj==2)||(dj==0))
      return 8;                 /*030T*/ 
    return 9;                   /*030C*/ 
  }else if((m==2)&&(n==1))    /*201*/
    return 10; 
  else if((m==1)&&(a==2)){    /*120*/
    di=g[j+i*gn]+g[k+i*gn];
    if(di==0)
      return 11;                /*120D*/
    dj=g[i+j*gn]+g[k+j*gn];
    if(dj==0)
      return 11;                /*120D*/
    dk=g[i+k*gn]+g[j+k*gn];
    if(dk==0)
      return 11;                /*120D*/
    di=g[i+j*gn]+g[i+k*gn];
    if(di==0)
      return 12;                /*120U*/
    dj=g[j+i*gn]+g[j+k*gn];
    if(dj==0)
      return 12;                /*120U*/
    dk=g[k+i*gn]+g[k+j*gn];
    if(dk==0)
      return 12;                /*120U*/
    return 13;                  /*120C*/
  }else if((m==2)&&(a==1))    /*210*/
    return 14; 
  else                        /*300*/
    return 15; 
}


int triad_classify_el(snaNet *g, int i, int j, int k, int gm, int checkna)
/*
If gm=1, compute a Holland and Leinhardt classification for the {i,j,k} triad
in graph g; otherwise, compute the four-state undirected classification.  Note
that this routine assumes dichotomized data; if checkna is true, edge
states are checked for missingness.  Specifically, checkna=1 results in a return
value of NA_INTEGER for triads with missing edges, and checkna=2 results in 
such edges being treated as absent (with checkna=0 implying no checks).

This routine is not intended to be called from R.
*/
{
  int sij=0,sji=0,sjk=0,skj=0,sik=0,ski=0,m,a,n,di,dj,dk;

  /*Get the raw edge states*/
  sij=snaIsAdjacent(i,j,g,checkna);
  sjk=snaIsAdjacent(j,k,g,checkna);
  sik=snaIsAdjacent(i,k,g,checkna);
  if(gm){
    sji=snaIsAdjacent(j,i,g,checkna);
    skj=snaIsAdjacent(k,j,g,checkna);
    ski=snaIsAdjacent(k,i,g,checkna);
  }

  /*If necessary, check for missingness*/
  if(checkna==1){
    if(IISNA(sij)||IISNA(sjk)||IISNA(sik))
      return NA_INTEGER;
    if(gm)
      if(IISNA(sji)||IISNA(skj)||IISNA(ski))
        return NA_INTEGER;
  }

  /*If we are in the undirected case, take the easy way out*/
  if(!gm){
    return sij+sjk+sik;
  }
      
  /*Get MAN information*/
  m=(sij*sji)+(sjk*skj)+(sik*ski);
  n=((sij+sji)==0)+((sjk+skj)==0)+((sik+ski)==0);
  a=3-m-n;

  /*Now classify, using dyad census as a first cut*/
  if(n==3)                    /*003*/
    return 0;
  else if((a==1)&&(n==2))     /*012*/
    return 1; 
  else if((m==1)&&(n==2))     /*102*/
    return 2; 
  else if((a==2)&&(n==1)){    /*021*/
    di=sij+sik;
    if(di==2)
      return 3;                 /*021D*/
    dj=sji+sjk;
    if(dj==2)
      return 3;                 /*021D*/
    dk=ski+skj;
    if(dk==2)
      return 3;                 /*021D*/
    di=sji+ski;
    if(di==2)
      return 4;                 /*021U*/
    dj=sij+skj;
    if(dj==2)
      return 4;                 /*021U*/
    dk=sik+sjk;
    if(dk==2)
      return 4;                 /*021U*/
    return 5;                   /*021C*/
  }else if((m==1)&&(n==1)){   /*111*/
    di=sji+ski;
    if((di==0)||(di==2))
      return 6;                 /*111D*/
    dj=sij+skj;
    if((dj==0)||(dj==2))
      return 6;                 /*111D*/
    return 7;                   /*111U*/
  }else if(a==3){            /*030*/
    di=sji+ski;
    if((di==2)||(di==0))
      return 8;                 /*030T*/ 
    dj=sij+skj;
    if((dj==2)||(dj==0))
      return 8;                 /*030T*/ 
    return 9;                   /*030C*/ 
  }else if((m==2)&&(n==1))    /*201*/
    return 10; 
  else if((m==1)&&(a==2)){    /*120*/
    di=sji+ski;
    if(di==0)
      return 11;                /*120D*/
    dj=sij+skj;
    if(dj==0)
      return 11;                /*120D*/
    dk=sik+sjk;
    if(dk==0)
      return 11;                /*120D*/
    di=sij+sik;
    if(di==0)
      return 12;                /*120U*/
    dj=sji+sjk;
    if(dj==0)
      return 12;                /*120U*/
    dk=ski+skj;
    if(dk==0)
      return 12;                /*120U*/
    return 13;                  /*120C*/
  }else if((m==2)&&(a==1))    /*210*/
    return 14; 
  else                        /*300*/
    return 15; 
}


void triad_classify_R(int *g, int *tt, int *gm)
/*
Given a triadic adjacency matrix, classify the triad in question.  Note that this routine assumes dichotomized data.  (This is a wrapper for triad_classify.)

This routine may be called from R using .C.
*/
{

/*Perform the classification*/
*tt=triad_classify(g,3,0,1,2,*gm);
}


void transitivity_R(double *mat, int *n, int *m, double *t, int *meas, int *checkna)
/*
Compute transitivity information for the (edgelist) network in mat.  This is
stored in t, with t[0] being the number of ordered triads at risk for
transitivity, and t[1] being the number satisfying the condition.  The 
definition used is controlled by meas, with meas==1 implying the weak condition
(a->b->c => a->c), meas==0 implying the strong condition (a->b->c <=>a->c), 
meas==2 implying the rank condition (a->c >= min(a->b,b->c)), and meas==3
implying Dekker's correlation measure (cor(a->c,a->b*b->c)).
If checkna==0, the measures are computed without missingness checks (i.e.,
treating NA edges as present).  If checkna==1, any triad containing missing
edges is omitted from the total count.  Finally, if checkna==2, missing edges
are treated as absent by the routine.

This routine may be called from R using .C.
*/
{
  int i,j,k,sij,sjk,sik;
  double ev;
  snaNet *g;
  slelement *jp,*kp,*ikp;

  /*Form the snaNet and initialize t*/
  GetRNGstate();
  //Rprintf("Building network, %d vertices and %d edges\n",*n,*m);
  g=elMatTosnaNet(mat,n,m);
  //Rprintf("Build complete.  Proceeding.\n");
  PutRNGstate();
  t[0]=t[1]=0.0;

  /*Get the transitivity information*/
  switch(*meas){
    case 0:                    /*"Strong" form: i->j->k <=> i->k*/
      for(i=0;i<g->n;i++)
        for(j=0;j<g->n;j++)
          if(i!=j){
            for(k=0;k<g->n;k++)
              if((j!=k)&&(i!=k)){
                sij=snaIsAdjacent(i,j,g,*checkna);
                sjk=snaIsAdjacent(j,k,g,*checkna);
                sik=snaIsAdjacent(i,k,g,*checkna);
                if(!(IISNA(sij)||IISNA(sjk)||IISNA(sik))){
                  t[0]+=sij*sjk*sik+(1-sij*sjk)*(1-sik);
                  t[1]++;
                }
              }
        }
      break;
    case 1:                      /*"Weak" form: i->j->k => i->k*/
      for(i=0;i<g->n;i++){
        for(jp=snaFirstEdge(g,i,1);jp!=NULL;jp=jp->next[0]){
          if((i!=(int)(jp->val))&&((*checkna==0)||(!ISNAN(*((double *)(jp->dp)))))){ /*Case 1 acts like case 2 here*/
            for(kp=snaFirstEdge(g,(int)(jp->val),1);kp!=NULL;kp=kp->next[0]){
              if(((int)(jp->val)!=(int)(kp->val))&&(i!=(int)(kp->val))&& ((*checkna==0)||(!ISNAN(*((double *)(kp->dp)))))){
                sik=snaIsAdjacent(i,(int)(kp->val),g,*checkna);
                if(!IISNA(sik)){  /*Not counting in case 1 (but am in case 2)*/
                  t[0]+=sik;
                  t[1]++;
                }
              }
            }
          }
        }
      }
      break;
    case 2:                    /*"Rank" form: i->k >= min(i->j,j->k)*/
      for(i=0;i<g->n;i++){
        for(jp=snaFirstEdge(g,i,1);jp!=NULL;jp=jp->next[0]){
          if((i!=(int)(jp->val))&&((*checkna==0)||(!ISNAN(*((double *)(jp->dp)))))){ /*Case 1 acts like case 2 here*/
            for(kp=snaFirstEdge(g,(int)(jp->val),1);kp!=NULL;kp=kp->next[0]){
              if(((int)(jp->val)!=(int)(kp->val))&&(i!=(int)(kp->val))&& ((*checkna==0)||(!ISNAN(*((double *)(kp->dp)))))){
                sik=snaIsAdjacent(i,(int)(kp->val),g,*checkna);
                if(!IISNA(sik)){  /*Not counting in case 1 (but am in case 2)*/
                  if(sik){
                    ikp=slistSearch(g->oel[i],kp->val); /*We already verified that it is here*/
                    ev=*((double *)(ikp->dp));
                  }else{
                    ev=0.0;
                  }
                  if((*checkna==0)||(!ISNAN(ev))){
                    t[0]+=(ev>=MIN(*((double *)(kp->dp)),*((double *)(jp->dp))));
                    t[1]++;
                  }
                }
              }
            }
          }
        }
      }
      break;
    case 3:                    /*"Corr" form: corr(i->k, i->j * j->k)*/
      error("Edgelist computation not currently supported for correlation measure in gtrans.\n");
      break;
  }
}


void triad_census_R(double *mat, int *n, int *m, double *t, int *gm, int *checkna)
/*
Compute a Holland and Leinhardt triad census for the graph with edgelist
matrix mat.  It is regrettably non-optimized, although it does at least avoid
having to store the entire adjacency matrix.  If checkna==0, the census is
computed without missingness checks (i.e., treating NA edges as present).  If
checkna==1, any triad containing missing edges is omitted from the total count.
Finally, if checkna==2, missing edges are treated as absent by the routine.

This routine may be called from R using .C.
*/
{
  int i,j,k,tc;
  snaNet *g;

  /*Form the snaNet*/
  GetRNGstate();
  g=elMatTosnaNet(mat,n,m);
  PutRNGstate();

  /*Clear out triad structure*/
  for(i=0;i<4+(*gm)*12;i++)     /*Vector length depends on gm*/
    t[i]=0.0;
  /*Get the triad counts*/
  for(i=0;i<*n;i++)
    for(j=i+1;j<*n;j++)
      for(k=j+1;k<*n;k++){
        tc=triad_classify_el(g,i,j,k,*gm,*checkna);
        if(!IISNA(tc))
          t[tc]++;
      }
}
