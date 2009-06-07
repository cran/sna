/*
######################################################################
#
# layout.c
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 12/19/04
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines related to computation of vertex layouts
# for gplot and gplot3d (i.e., the gplot.layout.* and gplot3d.layout.*
# functions).
#
######################################################################
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include "layout.h"


/*LAYOUT UTILITY ROUTINES----------------------------------------------*/

double angdist(double a, double b, double ilen)
/*This little routine calculates the angular distance between a and b,
assuming a periodicity of ilen.*/
{
  double minang,maxang,dis;
  
  /*Order the angles for convenience*/
  minang= ((a<=b) ? a : b);  
  maxang= ((a>b) ? a : b);  

  /*Find the distance, shifting if need be*/
  dis=maxang-minang;
  if(dis>ilen)
    dis-=ilen;

  /*Return the result*/
  return dis;
}

double poldist(double ra,double ta,double rb,double tb)
/*Return the Euclidean distance between a and b, where both points are given
in polar coordinates.*/
{
  return sqrt(ra*ra+rb*rb-2.0*ra*rb*cos(ta-tb));
}

double pollinedist(double ra,double ta,double rb,double tb, double rc, double tc)
/*Return the shortest Euclidean distance between a and the line spanning
b and c, where all points are given in polar coordinates.*/
{
  double A,B,dpol;
  
  A=ra*(rb*sin(ta-tb)-rc*sin(ta-tc))+rb*rc*sin(tb-tc);
  B=(rb*cos(tb)-rc*cos(tc)) * sqrt(1.0+pow(rb*sin(tb)-rc*sin(tc),2.0)/pow(rb*cos(tb)-rc*cos(tc),2.0));
  dpol=fabs(A/B);

  return dpol; 
}

int poledgecross(double ra, double ta, double rb, double tb, double rc, double tc, double rd, double td)
/*Checks for an edge cross between {a,b} and {c,d}, where all points are
specified in polar coordinates.  poledgecross() returns 1 if a cross occurs,
or 0 otherwise.*/
{
  double ax,ay,bx,by,cx,cy,dx,dy,denom,sint,tint,scx,scy,sdx;

  /*Convert to Cartesian coordinates*/
  ax=ra*cos(ta);
  ay=ra*sin(ta);
  bx=rb*cos(tb);
  by=rb*sin(tb);
  cx=rc*cos(tc);
  cy=rc*sin(tc);
  dx=rd*cos(td);
  dy=rd*sin(td);

  /*Compute the denomenators for the intersection*/
  denom=(ax-bx)*(cy-dy)-(ay-by)*(cx-dx);
  if(denom==0.0){  /*Interrupt if parallel lines*/
    /*Check for horizontility/verticality*/
    if(ax==bx){
      if((ax==cx)&&(ININT(cx,ax,bx)||ININT(dx,ax,bx)))
        return 1;
      else
        return 0;        
    }
    if(ay==by){
      if((ay==cy)&&(ININT(cy,ay,by)||ININT(dy,ay,by)))
        return 1;
      else
        return 0;        
    }
    /*Check for collinearity*/
    scx=(cx-ax)/(bx-ax);
    scy=(cy-ay)/(by-ay);
    if(scx!=scy)    /*Parallel&&!collin => not intersect*/
      return 0;
    /*If collinear, try for an intersection*/
    sdx=(dx-ax)/(bx-ax);
    if(((scx>=0.0)&&(scx<=1.0))||((sdx>=0.0)&&(sdx<=1.0)))
      return 1;
    else
      return 0;
  }
  sint=(ax*(cy-dy)+ay*(dx-cx)+cx*dy-cy*dx)/denom;
  tint=(ax*(by-cy)+ay*(cx-bx)+bx*cy-by*cx)/(-denom);
  
  /*Return the result*/
  if((sint>=0.0)&&(sint<=1.0)&&(tint>=0.0)&&(tint<=1.0))
    return 1;
  else
    return 0;
}

/*TWO-DIMENSIONAL LAYOUT ROUTINES--------------------------------------*/

void gplot_layout_fruchtermanreingold_R(double *d, int *pn, int *pm, int *pniter, double *pmaxdelta, double *pvolume, double *pcoolexp, double *prepulserad, double *x, double *y)
/*
Calculate a two-dimensional Fruchterman-Reingold layout for (symmetrized) 
edgelist matrix d.  Positions (stored in (x,y)) should be initialized
prior to calling this routine.
*/
{
  double frk,maxdelta,volume,coolexp,repulserad,t,ded,xd,yd,*dx,*dy;
  double rf,af;
  int n,j,k,niter,i,m,l;
  
  /*Define various things*/
  n=(int)*pn;
  m=(int)*pm;
  niter=*pniter;
  maxdelta=*pmaxdelta;
  volume=*pvolume;
  coolexp=*pcoolexp;
  repulserad=*prepulserad;
  frk=sqrt(volume/(double)n); /*Define the F-R constant*/

  /*Allocate memory for transient structures*/
  dx=(double *)R_alloc(n,sizeof(double));
  dy=(double *)R_alloc(n,sizeof(double));
  /*Run the annealing loop*/
  for(i=niter;i>=0;i--){
    /*Set the temperature (maximum move/iteration)*/
    t=maxdelta*pow(i/(double)niter,coolexp);
    /*Clear the deltas*/
    for(j=0;j<n;j++){
      dx[j]=0.0;
      dy[j]=0.0;
    }
    /*Increment deltas for each undirected pair*/
    for(j=0;j<n;j++)
      for(k=j+1;k<n;k++){
        /*Obtain difference vector*/
        xd=x[j]-x[k];
        yd=y[j]-y[k];
        ded=sqrt(xd*xd+yd*yd);  /*Get dyadic euclidean distance*/
        xd/=ded;                /*Rescale differences to length 1*/
        yd/=ded;
        /*Calculate repulsive "force"*/
        rf=frk*frk*(1.0/ded-ded*ded/repulserad);
        dx[j]+=xd*rf;        /*Add to the position change vector*/
        dx[k]-=xd*rf;
        dy[j]+=yd*rf;
        dy[k]-=yd*rf;
      }
    /*Calculate the attractive "force"*/
    for(j=0;j<m;j++){
      k=(int)d[j]-1;
      l=(int)d[j+m]-1;
      if(k<l){
        xd=x[k]-x[l];
        yd=y[k]-y[l];
        ded=sqrt(xd*xd+yd*yd);  /*Get dyadic euclidean distance*/
        xd/=ded;                /*Rescale differences to length 1*/
        yd/=ded;
        af=d[j+2*m]*ded*ded/frk;
        dx[k]-=xd*af;        /*Add to the position change vector*/
        dx[l]+=xd*af;
        dy[k]-=yd*af;
        dy[l]+=yd*af;
      }
    }
    /*Dampen motion, if needed, and move the points*/
    for(j=0;j<n;j++){
      ded=sqrt(dx[j]*dx[j]+dy[j]*dy[j]);
      if(ded>t){                 /*Dampen to t*/
        ded=t/ded;
        dx[j]*=ded;
        dy[j]*=ded;
      }
      x[j]+=dx[j];               /*Update positions*/
      y[j]+=dy[j];
    }
  }
}


void gplot_layout_kamadakawai_R(int *pn, int *pniter, double *elen, double *pinitemp, double *pcoolexp, double *pkkconst, double *psigma, double *x, double *y)
{
  double initemp,coolexp,sigma,temp,candx,candy;
  double dpot,odis,ndis,osqd,nsqd,kkconst;
  int niter,n,i,j,k;
  
  /*Define various things*/
  n=(int)*pn;
  niter=*pniter;
  initemp=*pinitemp;
  coolexp=*pcoolexp;
  kkconst=*pkkconst;
  sigma=*psigma;
  GetRNGstate();   /*Get the RNG state*/
  
  /*Perform the annealing loop*/
  temp=initemp;
  for(i=0;i<niter;i++){
    /*Update each vertex*/
    for(j=0;j<n;j++){
      /*Draw the candidate via a gaussian perturbation*/
      candx=rnorm(x[j],sigma*temp/initemp);
      candy=rnorm(y[j],sigma*temp/initemp);
      /*Calculate the potential difference for the new position*/
      dpot=0.0;
      for(k=0;k<n;k++)  /*Potential differences for pairwise effects*/
        if(j!=k){
          odis=sqrt((x[j]-x[k])*(x[j]-x[k])+(y[j]-y[k])*(y[j]-y[k]));
          ndis=sqrt((candx-x[k])*(candx-x[k])+(candy-y[k])*(candy-y[k]));
          osqd=(odis-elen[j+k*n])*(odis-elen[j+k*n]);
          nsqd=(ndis-elen[j+k*n])*(ndis-elen[j+k*n]);
          dpot+=kkconst*(osqd-nsqd)/(elen[j+k*n]*elen[j+k*n]);
        }
      /*Make a keep/reject decision*/
      if(log(runif(0.0,1.0))<dpot/temp){
        x[j]=candx;
        y[j]=candy;
      }
    }
    /*Cool the system*/
    temp*=coolexp;
  }
  PutRNGstate();   /*Update the RNG*/
}

void gplot_layout_target_R(int *d, double *pn, int *pniter, double *elen, double *radii, int *core, double *pdisconst, double *pcrossconst, double *prepconst, double *pminpdis, double *pinitemp, double *pcoolexp, double *pmaxdelta, double *theta)
{
  double initemp,coolexp,maxdelta,temp,c,dpot,odis,ndis,osqd,nsqd,repconst;
  double odjk,odjl,odjekl,ndjk,ndjl,ndjekl,opot,npot,disconst,crossconst;
  double minpdis4;
  int niter;
  long int n,i,j,k,l,m;
  
  /*Define various things*/
  n=(long int)*pn;
  niter=*pniter;
  initemp=*pinitemp;
  coolexp=*pcoolexp;
  maxdelta=*pmaxdelta;
  repconst=*prepconst;
  disconst=*pdisconst;
  crossconst=*pcrossconst;
  minpdis4=pow(*pminpdis,4.0);
  GetRNGstate();   /*Get the RNG state*/
  
  /*Perform the first annealing loop: core layout w/mutuals*/
  temp=initemp;
  for(i=0;i<niter;i++){
    /*Update each core vertex*/
    for(j=0;j<n;j++)
      if(core[j]){
        /*Draw the candidate via an angular perturbation*/
        c=theta[j]+runif(-1.0,1.0)*temp/initemp*maxdelta;
        while(c>=2.0*PI)   /*Map to [0,2pi) interval*/
          c-=2.0*PI;
        while(c<0.0)
          c+=2.0*PI;
        /*Calculate the potential difference for the new position*/
        dpot=0.0;
        for(k=0;k<n;k++)  
          if((j!=k)&&(core[k])){
            /*Angular distance potential per core alter*/
            odis=poldist(radii[j],theta[j],radii[k],theta[k]);
            ndis=poldist(radii[j],c,radii[k],theta[k]);
            osqd=(odis-elen[j+k*n])*(odis-elen[j+k*n]);
            nsqd=(ndis-elen[j+k*n])*(ndis-elen[j+k*n]);
            dpot+=disconst*(osqd-nsqd)/(elen[j+k*n]*elen[j+k*n]);
            /*Edge repulsion/line crossing potential*/
            for(l=0;l<n;l++)
              if(temp>1.0){
                /*Penalize edge crossings*/
                for(m=0;m<n;m++)
                  if((j!=l)&&(j!=m)&&(k!=l)&&(k!=m)&&(l!=m) && core[l]&&core[m] && ((d[j+k*n]>0)&&(d[k+j*n]>0)) && ((d[l+m*n]>0)&&(d[m+l*n]>0))){
                    opot=(double)poledgecross(radii[j],theta[j], radii[k],theta[k], radii[l],theta[l], radii[m],theta[m]);
                    npot=(double)poledgecross(radii[j],c, radii[k],theta[k], radii[l],theta[l], radii[m],theta[m]);
                    dpot+=crossconst*(opot-npot);  /*Smaller is better*/
                  }
              }else{
                /*Repel j from reciprocated, core edges*/
                if((j!=l)&&(k!=l)&&(k!=l) && core[l] && ((d[k+l*n]>0)&&(d[l+k*n]>0))){
                  /*Calculate old potential*/
                  odjk=poldist(radii[j],theta[j],radii[k],theta[k]);
                  odjl=poldist(radii[j],theta[j],radii[l],theta[l]);
                  odjekl=pollinedist(radii[j],theta[j], radii[k],theta[k],radii[l],theta[l]);
                  if((odjekl<=odjk)&&(odjekl<=odjl))
                    opot=repconst/(odjekl*odjekl);
                  else
                    opot=0.0;
                  /*Calculate new potential*/
                  ndjk=poldist(radii[j],c,radii[k],theta[k]);
                  ndjl=poldist(radii[j],c,radii[l],theta[l]);
                  ndjekl=pollinedist(radii[j],c, radii[k],theta[k],radii[l],theta[l]);
                  if((ndjekl<=ndjk)&&(ndjekl<=ndjl))
                    npot=repconst/(ndjekl*ndjekl);
                  else
                    npot=0.0;
                  /*Add difference*/
                  dpot+=(opot-npot)/temp;  /*Smaller is better*/
                }
              }
          }
        /*Make a keep/reject decision*/
        if(log(runif(0.0,1.0))<dpot/temp){
          theta[j]=c;
        }
      }
    /*Cool the system*/
    temp*=coolexp;
  }

  /*Perform the second annealing loop: add asymmetric edges*/
  temp=1.0;
  for(i=0;i<niter;i++){
    /*Update each core vertex*/
    for(j=0;j<n;j++)
      if(core[j]){
        /*Draw the candidate via an angular perturbation*/
        c=theta[j]+runif(-1.0,1.0)*temp/initemp*maxdelta;
        while(c>=2.0*PI)   /*Map to [0,2pi) interval*/
          c-=2.0*PI;
        while(c<0.0)
          c+=2.0*PI;
        /*Calculate the potential difference for the new position*/
        dpot=0.0;
        for(k=0;k<n;k++)  
          if((j!=k)&&(core[k])){
            /*Edge repulsion potential*/
            for(l=0;l<n;l++)
                /*Repel j from all core edges*/
                if((j!=l)&&(k!=l)&&(k!=l) && core[l] && ((d[k+l*n]>0)||(d[l+k*n]>0))){
                  /*Calculate old potential*/
                  odjk=poldist(radii[j],theta[j],radii[k],theta[k]);
                  odjl=poldist(radii[j],theta[j],radii[l],theta[l]);
                  odjekl=pollinedist(radii[j],theta[j], radii[k],theta[k],radii[l],theta[l]);
                  if((odjekl<odjk)&&(odjekl<odjl))
                    opot=repconst/(odjekl*odjekl);
                  else
                    opot=0.0;
                  /*Calculate new potential*/
                  ndjk=poldist(radii[j],c,radii[k],theta[k]);
                  ndjl=poldist(radii[j],c,radii[l],theta[l]);
                  ndjekl=pollinedist(radii[j],c, radii[k],theta[k],radii[l],theta[l]);
                  if((ndjekl<=ndjk)&&(ndjekl<=ndjl))
                    npot=repconst/(ndjekl*ndjekl);
                  else
                    npot=0.0;
                  /*Add difference*/
                  dpot+=(opot-npot)/temp;  /*Smaller is better*/
                }
          }
        /*Make a keep/reject decision*/
        if(log(runif(0.0,1.0))<dpot/temp){
          theta[j]=c;
        }
      }
    /*Cool the system*/
    temp*=coolexp;
  }

  /*Perform the third annealing loop: place peripheral vertices*/
  temp=initemp;
  for(i=0;i<niter;i++){
    /*Update each peripheral vertex*/
    for(j=0;j<n;j++)
      if(!core[j]){
        /*Draw the candidate via an angular perturbation*/
        c=theta[j]+runif(-1.0,1.0)*temp/initemp*maxdelta;
        while(c>=2.0*PI)   /*Map to [0,2pi) interval*/
          c-=2.0*PI;
        while(c<0.0)
          c+=2.0*PI;
        /*Calculate the potential difference for the new position*/
        dpot=0.0;
        for(k=0;k<n;k++)  
          if(j!=k){
            /*Repulse j from all other vertices*/
            odis=poldist(radii[j],theta[j],radii[k],theta[k]);
            ndis=poldist(radii[j],c,radii[k],theta[k]);
            osqd=odis*odis;
            nsqd=ndis*ndis;
            dpot+=minpdis4/osqd-minpdis4/ndis;
            /*Attract j towards all adjacent vertices*/
            if((d[j+k*n]>0)||(d[k+j*n]>0)){
              odis=poldist(radii[j],theta[j],radii[k],theta[k]);
              ndis=poldist(radii[j],c,radii[k],theta[k]);
              osqd=odis*odis;
              nsqd=ndis*ndis;
              dpot+=nsqd-odis;
            }
            /*Edge repulsion potential*/
            for(l=0;l<n;l++)
              if(temp<=1.0){
                /*Repel j from all edges*/
                if((j!=l)&&(k!=l)&&(k!=l) && ((d[k+l*n]>0)||(d[l+k*n]>0))){
                  /*Calculate old potential*/
                  odjk=poldist(radii[j],theta[j],radii[k],theta[k]);
                  odjl=poldist(radii[j],theta[j],radii[l],theta[l]);
                  odjekl=pollinedist(radii[j],theta[j], radii[k],theta[k],radii[l],theta[l]);
                  if((odjekl<odjk)&&(odjekl<odjl))
                    opot=repconst/(odjekl*odjekl);
                  else
                    opot=0.0;
                  /*Calculate new potential*/
                  ndjk=poldist(radii[j],c,radii[k],theta[k]);
                  ndjl=poldist(radii[j],c,radii[l],theta[l]);
                  ndjekl=pollinedist(radii[j],theta[j], radii[k],theta[k],radii[l],theta[l]);
                  if((ndjekl<ndjk)&&(ndjekl<ndjl))
                    npot=repconst/(ndjekl*ndjekl);
                  else
                    npot=0.0;
                  /*Add difference*/
                  dpot+=(opot-npot)/temp;  /*Smaller is better*/
                }
              }
          }
        /*Make a keep/reject decision*/
        if(log(runif(0.0,1.0))<dpot/temp){
          theta[j]=c;
        }
      }
    /*Cool the system*/
    temp*=coolexp;
  }

  PutRNGstate();   /*Update the RNG*/
}


/*THREE-DIMENSIONAL LAYOUT ROUTINES------------------------------------*/

void gplot3d_layout_fruchtermanreingold_R(double *d, int *pn, int *pm, int *pniter, double *pmaxdelta, double *pvolume, double *pcoolexp, double *prepulserad, double *x, double *y, double *z)
/*
Calculate a three-dimensional Fruchterman-Reingold layout for (symmetrized) 
edgelist matrix d.  Positions (stored in (x,y,z)) should be initialized
prior to calling this routine.
*/
{
  double frk,maxdelta,volume,coolexp,repulserad,t,ded,xd,yd,zd,*dx,*dy,*dz;
  double rf,af;
  int n,j,k,niter,i,m,l;
  
  /*Define various things*/
  n=(int)*pn;
  m=(int)*pm;
  niter=*pniter;
  maxdelta=*pmaxdelta;
  volume=*pvolume;
  coolexp=*pcoolexp;
  repulserad=*prepulserad;
  frk=pow(volume/(double)n,1.0/3.0); /*Define the F-R constant*/

  /*Allocate memory for transient structures*/
  dx=(double *)R_alloc(n,sizeof(double));
  dy=(double *)R_alloc(n,sizeof(double));
  dz=(double *)R_alloc(n,sizeof(double));
  /*Run the annealing loop*/
  for(i=niter;i>=0;i--){
    /*Set the temperature (maximum move/iteration)*/
    t=maxdelta*pow(i/(double)niter,coolexp);
    /*Clear the deltas*/
    for(j=0;j<n;j++){
      dx[j]=0.0;
      dy[j]=0.0;
      dz[j]=0.0;
    }
    /*Increment deltas for each undirected pair*/
    for(j=0;j<n;j++)
      for(k=j+1;k<n;k++){
        /*Obtain difference vector*/
        xd=x[j]-x[k];
        yd=y[j]-y[k];
        zd=z[j]-z[k];
        ded=sqrt(xd*xd+yd*yd+zd*zd);  /*Get dyadic euclidean distance*/
        xd/=ded;                      /*Rescale differences to length 1*/
        yd/=ded;
        zd/=ded;
        /*Calculate repulsive "force"*/
        rf=frk*frk*(1.0/ded-ded*ded/repulserad);
        dx[j]+=xd*rf;        /*Add to the position change vector*/
        dx[k]-=xd*rf;
        dy[j]+=yd*rf;
        dy[k]-=yd*rf;
        dz[j]+=zd*rf;
        dz[k]-=zd*rf;
      }
    /*Calculate the attractive "force"*/
    for(j=0;j<m;j++){
      k=(int)d[j]-1;
      l=(int)d[j+m]-1;
      if(k<l){
        xd=x[k]-x[l];
        yd=y[k]-y[l];
        zd=z[k]-z[l];
        ded=sqrt(xd*xd+yd*yd+zd*zd);  /*Get dyadic euclidean distance*/
        xd/=ded;                    /*Rescale differences to length 1*/
        yd/=ded;
        zd/=ded;
        af=d[j+2*m]*ded*ded/frk;
        dx[k]-=xd*af;        /*Add to the position change vector*/
        dx[l]+=xd*af;
        dy[k]-=yd*af;
        dy[l]+=yd*af;
        dz[k]-=zd*af;
        dz[l]+=zd*af;
      }
    }
    /*Dampen motion, if needed, and move the points*/
    for(j=0;j<n;j++){
      ded=sqrt(dx[j]*dx[j]+dy[j]*dy[j]+dz[j]*dz[j]);
      if(ded>t){                 /*Dampen to t*/
        ded=t/ded;
        dx[j]*=ded;
        dy[j]*=ded;
        dz[j]*=ded;
      }
      x[j]+=dx[j];               /*Update positions*/
      y[j]+=dy[j];
      z[j]+=dz[j];
    }
  }
}


void gplot3d_layout_kamadakawai_R(double *pn, int *pniter, double *elen, double *pinitemp, double *pcoolexp, double *pkkconst, double *psigma, double *x, double *y, double *z)
{
  double initemp,coolexp,sigma,temp,cx,cy,cz;
  double dpot,odis,ndis,osqd,nsqd,kkconst;
  int niter;
  long int n,i,j,k;
  
  /*Define various things*/
  n=(long int)*pn;
  niter=*pniter;
  initemp=*pinitemp;
  coolexp=*pcoolexp;
  kkconst=*pkkconst;
  sigma=*psigma;
  GetRNGstate();   /*Get the RNG state*/
  
  /*Perform the annealing loop*/
  temp=initemp;
  for(i=0;i<niter;i++){
    /*Update each vertex*/
    for(j=0;j<n;j++){
      /*Draw the candidate via a gaussian perturbation*/
      cx=rnorm(x[j],sigma*temp/initemp);
      cy=rnorm(y[j],sigma*temp/initemp);
      cz=rnorm(z[j],sigma*temp/initemp);
      /*Calculate the potential difference for the new position*/
      dpot=0.0;
      for(k=0;k<n;k++)  /*Potential differences for pairwise effects*/
        if(j!=k){
          odis=sqrt((x[j]-x[k])*(x[j]-x[k])+(y[j]-y[k])*(y[j]-y[k]) +(z[j]-z[k])*(z[j]-z[k]));
          ndis=sqrt((cx-x[k])*(cx-x[k])+(cy-y[k])*(cy-y[k]) +(cz-z[k])*(cz-z[k]));
          osqd=(odis-elen[j+k*n])*(odis-elen[j+k*n]);
          nsqd=(ndis-elen[j+k*n])*(ndis-elen[j+k*n]);
          dpot+=kkconst*(osqd-nsqd)/(elen[j+k*n]*elen[j+k*n]);
        }
      /*Make a keep/reject decision*/
      if(log(runif(0.0,1.0))<dpot/temp){
        x[j]=cx;
        y[j]=cy;
        z[j]=cz;
      }
    }
    /*Cool the system*/
    temp*=coolexp;
  }
  PutRNGstate();   /*Update the RNG*/
}
