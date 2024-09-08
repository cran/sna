/*
######################################################################
#
# randomgraph.c
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 9/07/24
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
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
#include <R_ext/Utils.h>
#include "randomgraph.h"

void bn_cftp_R(int *g, int *pn, double *pi, double *sigma, double *rho, double *d, int *pmaxiter, int *sibdichot)
{
  int *lparents,*uparents,lec,uec;
  double *lne,lnpar,lnsib,lndblr,ep,*coins,*ctemp;
  int ostate,*lb,*ub,converged,mismatch,t,maxiter,n,i,j,k,x,*r,*c,*temp;
  
  /*Initialize various things*/
  n=(int)*pn;
  GetRNGstate();
  lparents=(int *)R_alloc(n*n,sizeof(int));
  uparents=(int *)R_alloc(n*n,sizeof(int));
  lb=(int *)R_alloc(n*n,sizeof(int));
  ub=(int *)R_alloc(n*n,sizeof(int));
  lne=(double *)R_alloc(n*n,sizeof(double));
  for(i=0;i<n;i++)                                    /*Consider refining...*/
    for(j=0;j<n;j++)  
      lne[i+j*n] = ((d[i+j*n]<1.0) ? log(1.0-d[i+j*n]) : -DBL_MAX); 
  lnpar = ((*pi<1.0) ? log(1.0-*pi) : -DBL_MAX);
  lnsib = ((*sigma<1.0) ? log(1.0-*sigma) : -DBL_MAX);
  lndblr = ((*rho<1.0) ? log(1.0-*rho) : -DBL_MAX);
  t=n*(n-1);
  maxiter=(int)(*pmaxiter);
  if(maxiter<=0)
    maxiter=INT_MAX;
  converged=0;
  r=(int *)R_alloc(t,sizeof(int));
  c=(int *)R_alloc(t,sizeof(int));
  coins=(double *)R_alloc(t,sizeof(double));
  for(i=0;i<t;i++){
    r[i]=(int)floor(runif(0.0,1.0)*n);
    do{
      c[i]=(int)floor(runif(0.0,1.0)*n);
    }while(r[i]==c[i]);
    coins[i]=runif(0.0,1.0);
  }
 
  /*Run the CFTP loop*/
  Rprintf("t=%d, maxiter=%d, converged=%d\n",t,maxiter,converged);
  while((!converged)&&(t<maxiter)){
    Rprintf("\tt=%d\n",t);
    /*Initialize bounding processes*/
    for(i=0;i<n;i++){
      for(j=0;j<n;j++){
        lb[i+j*n]=0;
        lparents[i+j*n]=0;
        if(i==j){
          ub[i+j*n]=0;
          uparents[i+j*n]=0;
        }else{
          ub[i+j*n]=1;
          uparents[i+j*n]=n-2;
        }
      } 
    }
    lec=0;
    uec=n*(n-1);
    /*Evolve forward in time*/  
    for(i=-t;i<0;i++){
      R_CheckUserInterrupt();           /*Check for interrupts*/
      /*Update chains*/
      j=r[t+i];
      k=c[t+i];
//      Rprintf("\tj=%d, k=%d\n",j,k);
      if(!converged){
        /*Update lower bound*/
        ostate=lb[j+k*n];
        if(*sibdichot)
          ep=1.0-exp(lne[j+k*n]+lb[k+j*n]*lnpar+(lparents[j+k*n]>0)*lnsib+ lb[k+j*n]*(lparents[j+k*n]>0)*lndblr);
        else
          ep=1.0-exp(lne[j+k*n]+lb[k+j*n]*lnpar+lparents[j+k*n]*lnsib+ lb[k+j*n]*lparents[j+k*n]*lndblr);
        if(coins[t+i]<=ep){
          lb[j+k*n]=1;    /*Set the edge*/
          /*If something has changed update the parent count*/
          if(ostate==0){
            /*Rprintf("\tAdded edge update\n");*/
            for(x=0;x<n;x++)
              if((lb[j+x*n])&&(j!=x)&&(x!=k)){
                lparents[k+x*n]++;
                lparents[x+k*n]++;
              }
            lec++;
          }
        }else{
          lb[j+k*n]=0;   /*Unset the edge*/
          /*If something has changed update the parent count*/
          if(ostate==1){
            /*Rprintf("\tDeleted edge update\n");*/
            for(x=0;x<n;x++)
              if((lb[j+x*n])&&(j!=x)&&(x!=k)){
                lparents[k+x*n]--;
                lparents[x+k*n]--;
              }
            lec--;
          }
        }
        /*Update upper bound*/
        ostate=ub[j+k*n];
        if(*sibdichot)
          ep=1.0-exp(lne[j+k*n]+ub[k+j*n]*lnpar+(uparents[j+k*n]>0)*lnsib+ ub[k+j*n]*(uparents[j+k*n]>0)*lndblr);
        else
          ep=1.0-exp(lne[j+k*n]+ub[k+j*n]*lnpar+uparents[j+k*n]*lnsib+ ub[k+j*n]*uparents[j+k*n]*lndblr);
        if(coins[t+i]<=ep){
          ub[j+k*n]=1;    /*Set the edge*/
          /*If something has changed update the parent count*/
          if(ostate==0){
            /*Rprintf("\tAdded edge update\n");*/
            for(x=0;x<n;x++)
              if((ub[j+x*n])&&(j!=x)&&(x!=k)){
                uparents[k+x*n]++;
                uparents[x+k*n]++;
              }
            uec++;
          }
        }else{
          ub[j+k*n]=0;   /*Unset the edge*/
          /*If something has changed update the parent count*/
          if(ostate==1){
            /*Rprintf("\tDeleted edge update\n");*/
            for(x=0;x<n;x++)
              if((ub[j+x*n])&&(j!=x)&&(x!=k)){
                uparents[k+x*n]--;
                uparents[x+k*n]--;
              }
            uec--;
          }
        }
        /*Test for convergence*/
//        Rprintf("\tlec=%d, uec=%d, diff=%d\n",lec,uec,uec-lec);
        if(lec==uec){  /*Try filtering by edge count first*/
          mismatch=0;
          for(j=0;(j<n)&&(!mismatch);j++)
            for(k=0;(k<n)&&(!mismatch);k++)
              if(lb[j+k*n]!=ub[j+k*n])
                mismatch++;
          if(!mismatch){              /*We have a winner!*/
            Rprintf("Converged!\n");
            converged++;
            for(j=0;j<n;j++)          /*Copy lower bound into g*/
              for(k=0;k<n;k++)
                g[j+k*n]=lb[j+k*n];
          }
        }
      }else{
        /*Update g*/
        ostate=g[j+k*n];
        if(*sibdichot)
          ep=1.0-exp(lne[j+k*n]+g[k+j*n]*lnpar+(lparents[j+k*n]>0)*lnsib+ g[k+j*n]*(lparents[j+k*n]>0)*lndblr);
        else
          ep=1.0-exp(lne[j+k*n]+g[k+j*n]*lnpar+lparents[j+k*n]*lnsib+ g[k+j*n]*lparents[j+k*n]*lndblr);
        if(coins[t+i]<=ep){
          g[j+k*n]=1;    /*Set the edge*/
          /*If something has changed update the parent count*/
          if(ostate==0){
            /*Rprintf("\tAdded edge update\n");*/
            for(x=0;x<n;x++)
              if((g[j+x*n])&&(j!=x)&&(x!=k)){
                lparents[k+x*n]++;
                lparents[x+k*n]++;
              }
          }
        }else{
          g[j+k*n]=0;   /*Unset the edge*/
          /*If something has changed update the parent count*/
          if(ostate==1){
            /*Rprintf("\tDeleted edge update\n");*/
            for(x=0;x<n;x++)
              if((g[j+x*n])&&(j!=x)&&(x!=k)){
                lparents[k+x*n]--;
                lparents[x+k*n]--;
              }
          }
        }
      }
    }
    /*If needed, recede further into the past*/
    if(!converged){
      /*First, save old random inputs*/
      temp=(int *)R_alloc(2*t,sizeof(int));
      for(i=0;i<t;i++)
        temp[i]=r[i];
      r=temp;
      temp=(int *)R_alloc(2*t,sizeof(int));
      for(i=0;i<t;i++)
        temp[i]=c[i];
      c=temp;
      ctemp=(double *)R_alloc(2*t,sizeof(double));
      for(i=0;i<t;i++)
        ctemp[i]=coins[i];
      coins=ctemp;
      /*Add new random inputs*/
      for(i=t;i<2*t;i++){
        r[i]=(int)floor(runif(0.0,1.0)*n);
        do{
          c[i]=(int)floor(runif(0.0,1.0)*n);
        }while(r[i]==c[i]);
        coins[i]=runif(0.0,1.0);
      }
      /*Update t*/
      t*=2;
    }
  }
  if(t>=maxiter){
    warning("Maximum CFTP iterations exceeded; returning missing values.  (Your sample may also be biased.)\n");
    for(i=0;i<n;i++){
      for(j=0;j<n;j++){
        g[i+j*n]=NA_INTEGER;
      }
    }
  }
  /*Save the RNG state*/
  PutRNGstate();
}


void bn_mcmc_R(int *g, double *pn, double *pdraws, double *pburn, int *pthin, double *pi, double *sigma, double *rho, double *d, double *delta, double *epsilon, int *sibdichot, double *maxedge)
{
  long int n,i,j,k,x,draws,burn,bc,*parents,*odeg;
  double *lne,*lni,lnpar,lnsib,lndblr,ep,lnsat,ec;
  int thin,tc,ostate,stopflag;
  
  /*Initialize various things*/
  n=(long int)*pn;
  draws=(long int)*pdraws;
  burn=(long int)*pburn;
  thin=(int)*pthin;
  GetRNGstate();
  parents=(long int *)R_alloc(n*n,sizeof(long int));
  odeg=(long int *)R_alloc(n,sizeof(long int));
  lne=(double *)R_alloc(n*n,sizeof(double));        /*Non-excitation event prob*/
  lni=(double *)R_alloc(n*n,sizeof(double));        /*Non-inhibitory event prob*/
  ec=0.0;                                           /*Edge count*/
  for(i=0;i<n;i++){                                 /*Start w/empty graph*/
    odeg[i]=0;
    for(j=0;j<n;j++){
      parents[i+j*n]=0;
    }    
  }
  for(i=0;i<n;i++)                                    /*Consider refining...*/
    for(j=0;j<n;j++){  
      lne[i+j*n] = ((d[i+j*n]<1.0) ? log(1.0-d[i+j*n]) : -DBL_MAX); 
      lni[i+j*n] = ((epsilon[i+j*n]<1.0) ? log(1.0-epsilon[i+j*n]) : -DBL_MAX); 
    }
  lnpar = ((*pi<1.0) ? log(1.0-*pi) : -DBL_MAX);
  lnsib = ((*sigma<1.0) ? log(1.0-*sigma) : -DBL_MAX);
  lndblr = ((*rho<1.0) ? log(1.0-*rho) : -DBL_MAX);
  lnsat = ((*delta<1.0) ? log(1.0-*delta) : -DBL_MAX);

  /*Add edges from the seed graph*/
  for(j=0;j<n;j++){
    for(i=0;i<n;i++){
      if(g[i*draws+j*n*draws]){
        ec++;                           /*Increment the edge count*/
	odeg[i]++;
        for(x=0;x<n;x++)
          if((g[i*draws+x*n*draws])&&(i!=x)&&(x!=j)){
            parents[j+x*n]++;
            parents[x+j*n]++;
          }
      }
    }
  }

  /*Run the MCMC loop*/
  bc=0;
  tc=0;
  stopflag=0;              /*Early stopping flag*/
  for(i=0;(i<draws)&&(!stopflag);i++){
    /*
    if(bc<burn)
      Rprintf("Burn-in Iteration %ld\n",bc);
    else
      Rprintf("Iteration %ld\n",i);
    */
    /*Rprintf("\tEdge selection\n");*/
    /*Select an edge to update*/
    j=(long int)floor(runif(0.0,1.0)*n);
    do{
      k=(long int)floor(runif(0.0,1.0)*n);
    }while(j==k);
    /*Rprintf("\tDrawing and updating\n");*/
    /*Redraw the edge from the full conditional*/
    ostate=g[i+j*draws+k*n*draws];
    if(*sibdichot)
      ep=1.0-exp(lne[j+k*n]+g[i+k*draws+j*n*draws]*lnpar+(parents[j+k*n]>0)*lnsib+ g[i+k*draws+j*n*draws]*(parents[j+k*n]>0)*lndblr);
    else
      ep=1.0-exp(lne[j+k*n]+g[i+k*draws+j*n*draws]*lnpar+parents[j+k*n]*lnsib+ g[i+k*draws+j*n*draws]*parents[j+k*n]*lndblr);
    ep*=exp(odeg[j]*lnsat+lni[j+k*n]);
    if(runif(0.0,1.0)<=ep){
      ec+=1.0-g[i+j*draws+k*n*draws];  /*Increment the edge count*/
      g[i+j*draws+k*n*draws]=1;    /*Set the edge*/
      /*If something has changed update the parent count and outdegree count*/
      if(ostate==0){
        odeg[j]++;
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
      ec-=g[i+j*draws+k*n*draws];  /*Decrement the edge count*/
      g[i+j*draws+k*n*draws]=0;   /*Unset the edge*/
      /*If something has changed update the parent and outdegree count*/
      if(ostate==1){
        odeg[j]--;
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
    /*Check for early stopping condition*/
    if(ec > *maxedge){     /*Stop if density guard is tripped*/
      stopflag=1;       
      *maxedge = -1.0;     /*Set to negative to show the trip state*/
    }
    /*Burn-in check*/
    if(bc<burn){
      /*If burning, increment burn count and don't move*/
      i--;   
      bc++;
    }else{ 
      /*Rprintf("\tThinning?\n");*/
      if(tc%thin==thin-1){
        if(i<draws-1){
          /*Increment the state*/
          for(j=0;j<n;j++)
            for(k=0;k<n;k++)
              g[i+1+j*draws+k*n*draws]=g[i+j*draws+k*n*draws];
        }
      }else
        i--;   /*Thinning*/
      tc++;
    }
  }
  /*Save the RNG state*/
  PutRNGstate();
}
void bn_mcmc_old_R(int *g, double *pn, double *pdraws, double *pburn, int *pthin, double *pi, double *sigma, double *rho, double *d)
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


SEXP rgbern_R(SEXP sn, SEXP stp, SEXP sdirected, SEXP sloops, SEXP spmode)
/*
Draw a Bernoulli graph using the algorithm of Batagelj and Brandes (2005), with minor modifications to allow for directed/undirected graphs, loops, and row/column heterogeneity.  The Bernoulli parameters are contained in stp, as either a vector or single value (as appropriate, given the mode).  spmode indicates the type of heterogeneity to employ:
  0 (BERNHOM) - Homogeneous
  1 (BERNROW) - Row-heterogeneous
  2 (BERNCOL) - Column heterogeneous
  3 (BERNHET) - Heterogeneous  (the brute force method is used here)
The resulting graph is returned as an sna edgelist, complete with "n" attribute.  Multiple calls to this function may be used to obtain multiple graphs.  Note
*/
{
  int i,j,n,m,directed,loops,pmode,pc=0;
  double *tp,*g,r,c,en,w;
  element *el,*ep;
  SEXP sg,sn2,dim;
  
  /*Initalize stuff*/
  PROTECT(sn=coerceVector(sn,INTSXP)); pc++;
  PROTECT(stp=coerceVector(stp,REALSXP)); pc++;
  PROTECT(sdirected=coerceVector(sdirected,INTSXP)); pc++;
  PROTECT(sloops=coerceVector(sloops,INTSXP)); pc++;
  PROTECT(spmode=coerceVector(spmode,INTSXP)); pc++;
  n=INTEGER(sn)[0];
  tp=REAL(stp);
  directed=INTEGER(sdirected)[0];
  loops=INTEGER(sloops)[0];
  pmode=INTEGER(spmode)[0];
  GetRNGstate();
  m=0;
  el=NULL;
    
  /*Generate the graph, using the appropriate mode*/
  switch(pmode){
    case BERNHOM:  /*For each row, use a waiting time scheme*/
      if(directed)            /*Get maximum number of edges*/
        en=(double)(n*(n-1.0));
      else
        en=(double)(n*(n-1.0)/2.0);
      en+=(double)(loops*n);
      w=-1.0;                 /*Draw edges*/
      while(w<en){
        w+=1.0+rgeom(tp[0]);
        if(w<en){
          if(directed){
            if(loops){
              r=fmod(w,(double)n);
              c=floor(w/(double)n);
            }else{
              c=floor(w/(n-1.0));
              r=fmod(w,n-1.0)+(fmod(w,n-1.0)>=c);
            }
          }else{
            if(loops){
              c=n-floor(sqrt(n*(n+1.0)-2.0*w-1.75)-0.5)-1.0;
              r=w-c*(n-1.0)+c*(c-1.0)/2.0;
            }else{
              c=n-2.0-floor(sqrt(n*(n-1.0)-2.0*w-1.75)-0.5);
              r=w+c*((c+1.0)/2.0-n+1.0)+1.0;
            }
          }
          el=enqueue(el,r+c*n,NULL);
          m++;
          if((!directed)&&(r!=c)){
            el=enqueue(el,c+r*n,NULL);
            m++;
          }
        }
      }
      break;
    case BERNROW:  /*For each row, use a waiting time scheme*/
      for(i=0;i<n;i++){    /*Walk through the rows*/
        if(directed)            /*Get maximum number of edges*/
          en=(double)(n-1.0);
        else
          en=(double)i;
        en+=loops;
        w=-1.0;                 /*Draw edges*/
        while(w<en){
          w+=1.0+rgeom(tp[i]);
          if(w<en){
            el=enqueue(el,i+(w+(!loops)*(w>=(double)i))*n,NULL);
            m++;
            if((!directed)&&((!loops)||(w!=(double)i))){
              el=enqueue(el,(w+(!loops)*(w>=(double)i))+i*n,NULL);
              m++;
            }
          }
        }
      }
      break;
    case BERNCOL:
      for(i=0;i<n;i++){    /*Walk through the cols*/
        if(directed)            /*Get maximum number of edges*/
          en=(double)(n-1.0);
        else
          en=(double)i;
        en+=loops;
        w=-1.0;                 /*Draw edges*/
        while(w<en){
          w+=1.0+rgeom(tp[i]);
          if(w<en){
            el=enqueue(el,(w+(!loops)*(w>=(double)i))+i*n,NULL);
            m++;
            if((!directed)&&((!loops)||(w!=(double)i))){
              el=enqueue(el,i+(w+(!loops)*(w>=(double)i))*n,NULL);
              m++;
            }
          }
        }
      }
      break;
    case BERNHET:  /*No shortcuts, just draw directly*/
      for(i=0;i<n;i++)
        for(j=i*(!directed);j<n;j++)
          if(loops||(i!=j))
            if(runif(0.0,1.0)<tp[i+j*n]){
              el=enqueue(el,(double)(i+j*n),NULL);
              m++;
              if((!directed)&&(i!=j)){
                el=enqueue(el,(double)(j+i*n),NULL);
                m++;
              }
            }
      break;
  }
  
  /*Store the result*/
  PROTECT(sg=allocVector(REALSXP,3*m)); pc++;
  g=REAL(sg);
  for(i=0,ep=el;ep!=NULL;ep=ep->next){
    c=floor((ep->val)/(double)n);
    r=fmod(ep->val,(double)n);
    g[i]=r+1;
    g[i+m]=c+1;
    g[i+2*m]=1.0;
    i++;
  }
  PROTECT(sn2=allocVector(INTSXP,1)); pc++;  /*Set graph size attribute*/
  INTEGER(sn2)[0]=n;
  setAttrib(sg,install("n"), sn2);
  PROTECT(dim=allocVector(INTSXP, 2)); pc++; /*Set dimension attribute*/
  INTEGER(dim)[0] = m; 
  INTEGER(dim)[1] = 3;
  setAttrib(sg,R_DimSymbol,dim);

  /*Unprotect and return*/
  PutRNGstate();
  UNPROTECT(pc);
  return sg;
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
