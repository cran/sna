/*
######################################################################
#
# likelihood.c
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 3/10/05
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
#
# Part of the R/sna package
#
# This file contains routines related likelihood calculation for
# stochastic network models.
#
######################################################################
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "likelihood.h"

double bn_lpka(long int k,double pi, double sigma, double rho, double d)
/*Return the log conditional probability of an asymetric dyad with k parents, 
under a biased net model.  Parameters are as follows:
  sigma = Pr( x->y | xSy ) (Transitivity or "sibling" bias)
  pi = Pr( x->y | y->x ) (Reciprocity or "parent" bias)
  rho = Pr( x->y | y->x & xSy ) ("Double role" bias, or interaction effect)
  d = Pr( x->y ) (Baseline density effect)
*/
{
  double dk,lnbpar,lnbsib,lnbdblr,lne,lpe1;

  /*Create building blocks*/
  dk=(double)k;
  lnbpar=log(1.0-pi);
  lnbsib=log(1.0-sigma)*dk;
  lnbdblr=log(1.0-rho)*dk;
  lne=log(1-d);

  /*Precompute stuff*/
  lpe1=log(1.0-exp(lnbsib+lne));
  
  /*Return the log-probability of an asymetric dyad*/
  return lpe1+lnbpar+lnbsib+lnbdblr+lne;
}

double bn_lpkm(long int k,double pi, double sigma, double rho, double d)
/*Return the log conditional probability of a mutual dyad with k parents, 
under a biased net model.  Parameters are as follows:
  sigma = Pr( x->y | xSy ) (Transitivity or "sibling" bias)
  pi = Pr( x->y | y->x ) (Reciprocity or "parent" bias)
  rho = Pr( x->y | y->x & xSy ) ("Double role" bias, or interaction effect)
  d = Pr( x->y ) (Baseline density effect)
*/
{
  double dk,lne,lnbsib,lnbpar,lnbdblr,lpe1,lpe2;

  /*Create building blocks*/
  dk=(double)k;
  lnbpar=log(1.0-pi);
  lnbsib=log(1.0-sigma)*dk;
  lnbdblr=log(1.0-rho)*dk;
  lne=log(1-d);

  /*Precompute stuff*/
  lpe1=log(1.0-exp(lnbsib+lne));
  lpe2=log(1.0-exp(lnbpar+lnbsib+lnbdblr+lne));
  
  /*Return the log-probability of a mutual dyad*/
  return lpe1+lpe2;
}

double bn_lpkn(long int k,double pi, double sigma, double rho, double d)
/*Return the log conditional probability of a null dyad with k parents, 
under a biased net model.  Parameters are as follows:
  sigma = Pr( x->y | xSy ) (Transitivity or "sibling" bias)
  pi = Pr( x->y | y->x ) (Reciprocity or "parent" bias)
  rho = Pr( x->y | y->x & xSy ) ("Double role" bias, or interaction effect)
  d = Pr( x->y ) (Baseline density effect)
*/
{
  double dk,lnbsib,lnbdblr,lnbpar,lne,p1,p2;

  /*Create building blocks*/
  dk=(double)k;
  lnbpar=log(1.0-pi);
  lnbsib=dk*log(1.0-sigma);
  lnbdblr=dk*log(1.0-rho);
  lne=log(1.0-d);

  /*Precompute stuff*/
  p1=1.0-exp(lnbsib+lne);
  p2=1.0+exp(lnbpar+lnbsib+lnbdblr+lne);

  /*Return the log-probability of a null dyad*/
  return log(1.0-p1*p2);
}

double bn_lpt_M(long int m,double pi,double sigma,double rho,double d)
{
  return log(1.0 - (1.0-pi)*pow(1.0-rho,(double)m)*pow(1.0-sigma,(double)m)*(1.0-d)) + log(1.0-pow(1.0-sigma,(double)m)*(1.0-d));
}

double bn_lpt_a(long int m,double pi,double sigma,double rho,double d)
{
  return log(1.0-pow(1.0-sigma,(double)m)*(1.0-d)) + log((1.0-pi)*pow(1.0-rho,(double)m)*pow(1.0-sigma,(double)m)*(1.0-d));
}

double bn_lpt_N(long int m,double pi,double sigma,double rho,double d)
{
  double calc;

  calc = 1.0 - exp(bn_lpt_M(m,pi,sigma,rho,d)) - 2.0*exp(bn_lpt_a(m,pi,sigma,rho,d));
  /*Check for numerical pathologies*/
  if(calc<0.0)
    calc=0.0;
    
  return log(calc);
}

double bn_lpt_Mp1(long int m,double pi,double sigma,double rho,double d)
{
  return bn_lpt_M(m+1,pi,sigma,rho,d);
}

double bn_lpt_ap1(long int m,double pi,double sigma,double rho,double d)
{
  return bn_lpt_a(m+1,pi,sigma,rho,d);
}

double bn_lpt_Np1(long int m,double pi,double sigma,double rho,double d)
{
  return bn_lpt_N(m+1,pi,sigma,rho,d);
}

double bn_lpt_M1(double pi,double sigma,double rho,double d)
{
  return log(sigma*(1.0-(1.0-sigma)*(1.0-rho)));
}

double bn_lpt_a1(double pi,double sigma,double rho,double d)
{
  return log(sigma*(1.0-sigma)*(1.0-rho));
}

double bn_lpt_N1(double pi,double sigma,double rho,double d)
{
  return log(1.0-sigma*(1.0+(1.0-sigma)*(1.0-rho)));
}

double bn_lpt_Sr(double pi,double sigma,double rho,double d)
{
  return log(1.0-(1.0-sigma)*(1.0-rho));
}

double bn_lpt_1mSr(double pi,double sigma,double rho,double d)
{
  return log((1.0-sigma)*(1.0-rho));
}

double bn_lpt(int xy, int yx, int yz, int zy, int xz, int zx, long int kxy, long int kyz, long int kxz, double pi, double sigma, double rho, double d)
/*Return the log conditional probability of an x,y,z triad with dyad parent
counts specified by the respective k parameters, under a biased net model.  
Bias parameters are as follows:
  sigma = Pr( x->y | xSy ) (Transitivity or "sibling" bias)
  pi = Pr( x->y | y->x ) (Reciprocity or "parent" bias)
  rho = Pr( x->y | y->x & xSy ) ("Double role" bias, or interaction effect)
  d = Pr( x->y ) (Baseline density effect)
Note that this uses the table in Skvoretz (2003).
*/
{

  if(xy>0){
    if(yx>0){
      if(yz>0){
        if(zy>0){
          if(xz>0){
            if(zx>0){  /*1 1 1 1 1 1*/
              return log((exp(bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_M(kyz,pi,sigma,rho,d) + bn_lpt_Mp1(kxz,pi,sigma,rho,d)) + exp(bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_Mp1(kyz,pi,sigma,rho,d) + bn_lpt_M(kxz,pi,sigma,rho,d))+exp(bn_lpt_Mp1(kxy,pi,sigma,rho,d) + bn_lpt_M(kyz,pi,sigma,rho,d) + bn_lpt_M(kxz,pi,sigma,rho,d)))/3.0
+ (exp(bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d))*(exp(bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_Mp1(kyz,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_M(kyz,pi,sigma,rho,d)) + exp(bn_lpt_Mp1(kxy,pi,sigma,rho,d) + bn_lpt_M(kyz,pi,sigma,rho,d)))
+ exp(bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d))*(exp(bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_Mp1(kxz,pi,sigma,rho,d))+2.0*exp(bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_M(kxz,pi,sigma,rho,d))+exp(bn_lpt_Mp1(kxy,pi,sigma,rho,d) + bn_lpt_M(kxz,pi,sigma,rho,d)))
+ exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d))*(exp(bn_lpt_Mp1(kyz,pi,sigma,rho,d) + bn_lpt_M(kxz,pi,sigma,rho,d))+2.0*exp(bn_lpt_M(kyz,pi,sigma,rho,d) + bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_M(kxz,pi,sigma,rho,d)) + exp(bn_lpt_M(kyz,pi,sigma,rho,d) + bn_lpt_Mp1(kxz,pi,sigma,rho,d))))/3.0
+4.0*(exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_M(kyz,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d))+exp(bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d)) + exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_M(kxz,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d)))/3.0
+ (exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_M(kxz,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d)) + exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_M(kyz,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d)) + exp(bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d)))/3.0);
            }else{     /*1 1 1 1 1 0*/
              return log(exp(bn_lpt_M(kxy,pi,sigma,rho,d)) * (exp(bn_lpt_ap1(kxz,pi,sigma,rho,d) + bn_lpt_M(kyz,pi,sigma,rho,d)) + exp(bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_Mp1(kyz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d)) + exp(bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_M(kyz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d)))/3.0 + 2.0*exp(bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d))/3.0 + exp(bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d))*(exp(bn_lpt_ap1(kxz,pi,sigma,rho,d)) + exp(bn_lpt_a(kxz,pi,sigma,rho,d)+bn_lpt_1mSr(pi,sigma,rho,d)))/3.0);
            }
          }else{
            if(zx>0){  /*1 1 1 1 0 1*/
              return log(exp(bn_lpt_M(kyz,pi,sigma,rho,d)) * (exp(bn_lpt_ap1(kxz,pi,sigma,rho,d) + bn_lpt_M(kxy,pi,sigma,rho,d)) + exp(bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_Mp1(kxy,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d)) + exp(bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d)))/3.0 + 2.0*exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_M(kyz,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d))/3.0 + exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_M(kyz,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d))*(exp(bn_lpt_ap1(kxz,pi,sigma,rho,d)) + exp(bn_lpt_a(kxz,pi,sigma,rho,d)+bn_lpt_1mSr(pi,sigma,rho,d)))/3.0);
            }else{     /*1 1 1 1 0 0*/
              return bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_M(kyz,pi,sigma,rho,d) + log(exp(bn_lpt_Np1(kxz,pi,sigma,rho,d))+2.0*exp(bn_lpt_N(kxz,pi,sigma,rho,d))) - log(3.0);
            }
          }
        }else{
          if(xz>0){ 
            if(zx>0){  /*1 1 1 0 1 1*/
              return log(exp(bn_lpt_M(kxy,pi,sigma,rho,d)) * (exp(bn_lpt_ap1(kyz,pi,sigma,rho,d) + bn_lpt_M(kxz,pi,sigma,rho,d)) + exp(bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_Mp1(kxz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d)) + exp(bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_M(kxz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d)))/3.0 + 2.0*exp(bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d))/3.0 + exp(bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d))*(exp(bn_lpt_ap1(kyz,pi,sigma,rho,d)) + exp(bn_lpt_a(kyz,pi,sigma,rho,d)+bn_lpt_1mSr(pi,sigma,rho,d)))/3.0);
            }else{     /*1 1 1 0 1 0*/
              return log(exp(bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_ap1(kxz,pi,sigma,rho,d)) + exp(bn_lpt_ap1(kxz,pi,sigma,rho,d) + bn_lpt_a(kyz,pi,sigma,rho,d))+exp(bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d))) + bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d) - log(3.0);
            }
          }else{
            if(zx>0){  /*1 1 1 0 0 1*/
              return bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_a(kyz,pi,sigma,rho,d) + log(exp(bn_lpt_ap1(kxz,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d))) - log(3.0);
            }else{     /*1 1 1 0 0 0*/
              return bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_a(kyz,pi,sigma,rho,d) + log(exp(bn_lpt_Np1(kxz,pi,sigma,rho,d))+2.0*exp(bn_lpt_N(kxz,pi,sigma,rho,d))) - log(3.0);
            }
          }
        }
      }else{
        if(zy>0){
          if(xz>0){
            if(zx>0){  /*1 1 0 1 1 1*/
              return log(exp(bn_lpt_M(kxz,pi,sigma,rho,d)) * (exp(bn_lpt_ap1(kyz,pi,sigma,rho,d) + bn_lpt_M(kxy,pi,sigma,rho,d)) + exp(bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_Mp1(kxy,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d)) + exp(bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_M(kxy,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d)))/3.0 + 2.0*exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_M(kxz,pi,sigma,rho,d) + bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d))/3.0 + exp(bn_lpt_M(kxz,pi,sigma,rho,d) + bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d))*(exp(bn_lpt_ap1(kyz,pi,sigma,rho,d)) + exp(bn_lpt_a(kyz,pi,sigma,rho,d)+bn_lpt_1mSr(pi,sigma,rho,d)))/3.0);
            }else{     /*1 1 0 1 1 0*/
              return bn_lpt_M(kxy,pi,sigma,rho,d) + log(exp(bn_lpt_ap1(kyz,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d))) + bn_lpt_a(kxz,pi,sigma,rho,d) - log(3.0);
            }
          }else{
            if(zx>0){  /*1 1 0 1 0 1*/
              return log(exp(bn_lpt_Mp1(kxy,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_M(kxy,pi,sigma,rho,d)) + 4.0*exp(bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_Sr(pi,sigma,rho,d))) + bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d) - log(3.0);
            }else{     /*1 1 0 1 0 0*/
              return bn_lpt_M(kxy,pi,sigma,rho,d)+bn_lpt_a(kyz,pi,sigma,rho,d)+ bn_lpt_N(kxz,pi,sigma,rho,d);
            }
          }
        }else{
          if(xz>0){
            if(zx>0){  /*1 1 0 0 1 1*/
              return bn_lpt_M(kxy,pi,sigma,rho,d) + log(exp(bn_lpt_Np1(kyz,pi,sigma,rho,d))+2.0*exp(bn_lpt_N(kyz,pi,sigma,rho,d))) + bn_lpt_a(kxz,pi,sigma,rho,d) - log(3.0);
            }else{     /*1 1 0 0 1 0*/
              return bn_lpt_M(kxy,pi,sigma,rho,d) + log(exp(bn_lpt_Np1(kyz,pi,sigma,rho,d))+2.0*exp(bn_lpt_N(kyz,pi,sigma,rho,d))) + bn_lpt_a(kxz,pi,sigma,rho,d) - log(3.0);
            }
          }else{
            if(zx>0){  /*1 1 0 0 0 1*/
              return bn_lpt_M(kxy,pi,sigma,rho,d)+bn_lpt_N(kyz,pi,sigma,rho,d)+ bn_lpt_a(kxz,pi,sigma,rho,d);
            }else{     /*1 1 0 0 0 0*/
              return bn_lpt_M(kxy,pi,sigma,rho,d)+bn_lpt_N(kyz,pi,sigma,rho,d)+ bn_lpt_N(kxz,pi,sigma,rho,d);
            }
          }
        }
      }
    }else{
      if(yz>0){
        if(zy>0){
          if(xz>0){
            if(zx>0){  /*1 0 1 1 1 1*/
              return log(exp(bn_lpt_M(kyz,pi,sigma,rho,d)) * (exp(bn_lpt_ap1(kxy,pi,sigma,rho,d) + bn_lpt_M(kxz,pi,sigma,rho,d)) + exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_Mp1(kxz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d)) + exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_M(kxz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d)))/3.0 + 2.0*exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_M(kyz,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d))/3.0 + exp(bn_lpt_M(kyz,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_Sr(pi,sigma,rho,d))*(exp(bn_lpt_ap1(kxy,pi,sigma,rho,d)) + exp(bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_1mSr(pi,sigma,rho,d)))/3.0);
            }else{     /*1 0 1 1 1 0*/
              return bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_M(kyz,pi,sigma,rho,d) + log(exp(bn_lpt_ap1(kxz,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d))) - log(3.0);
            }
          }else{
            if(zx>0){  /*1 0 1 1 0 1*/
              return log(exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_ap1(kxz,pi,sigma,rho,d)) + exp(bn_lpt_ap1(kxy,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d))+exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d))) + bn_lpt_M(kyz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d) - log(3.0);
            }else{     /*1 0 1 1 0 0*/
              return bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_M(kxy,pi,sigma,rho,d) + log(exp(bn_lpt_Np1(kxz,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_N(kxz,pi,sigma,rho,d))) - log(3.0);
            }
          }
        }else{
          if(xz>0){
            if(zx>0){  /*1 0 1 0 1 1*/
              return bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_a(kxy,pi,sigma,rho,d) + log(exp(bn_lpt_Mp1(kxz,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_M(kxz,pi,sigma,rho,d)) + 4.0*exp(bn_lpt_a(kxz,pi,sigma,rho,d)+bn_lpt_Sr(pi,sigma,rho,d))) - log(3.0);
            }else{     /*1 0 1 0 1 0*/
              return bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_a(kyz,pi,sigma,rho,d) + log(exp(bn_lpt_ap1(kxz,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d))) - log(3.0);
            }
          }else{
            if(zx>0){  /*1 0 1 0 0 1*/
              return bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_a(kyz,pi,sigma,rho,d) + log(exp(bn_lpt_ap1(kxz,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d))) - log(3.0);
            }else{     /*1 0 1 0 0 0*/
              return bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_a(kyz,pi,sigma,rho,d) + log(exp(bn_lpt_Np1(kxz,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_N(kxz,pi,sigma,rho,d) + bn_lpt_N1(pi,sigma,rho,d))) - log(3.0);
            }
          }
        }
      }else{
        if(zy>0){
          if(xz>0){
            if(zx>0){  /*1 0 0 1 1 1*/
              return log(exp(bn_lpt_ap1(kxy,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d))) + bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_M(kxz,pi,sigma,rho,d) - log(3.0);
            }else{     /*1 0 0 1 1 0*/
              return bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_a(kyz,pi,sigma,rho,d)+ bn_lpt_a(kxz,pi,sigma,rho,d);
            }
          }else{
            if(zx>0){  /*1 0 0 1 0 1*/
              return log(exp(bn_lpt_ap1(kxy,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d))) + bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d) - log(3.0);
            }else{     /*1 0 0 1 0 0*/
              return bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_a(kyz,pi,sigma,rho,d)+ bn_lpt_N(kxz,pi,sigma,rho,d);
            }
          }
        }else{
          if(xz>0){
            if(zx>0){  /*1 0 0 0 1 1*/
              return bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_N(kyz,pi,sigma,rho,d)+ bn_lpt_M(kxz,pi,sigma,rho,d);
            }else{     /*1 0 0 0 1 0*/
               return bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_N(kyz,pi,sigma,rho,d)+ bn_lpt_a(kxz,pi,sigma,rho,d);
           }
          }else{
            if(zx>0){  /*1 0 0 0 0 1*/
              return bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_N(kyz,pi,sigma,rho,d)+ bn_lpt_a(kxz,pi,sigma,rho,d);
            }else{     /*1 0 0 0 0 0*/
              return bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_N(kyz,pi,sigma,rho,d)+ bn_lpt_N(kxz,pi,sigma,rho,d);
            }
          }
        }
      }
    }
  }else{
    if(yx>0){
      if(yz>0){
        if(zy>0){
          if(xz>0){
            if(zx>0){  /*0 1 1 1 1 1*/
              return log(exp(bn_lpt_M(kyz,pi,sigma,rho,d))*(exp(bn_lpt_ap1(kxy,pi,sigma,rho,d)+bn_lpt_M(kxz,pi,sigma,rho,d))+exp(bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_Mp1(kxz,pi,sigma,rho,d)+bn_lpt_1mSr(pi,sigma,rho,d))+exp(bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_M(kxz,pi,sigma,rho,d)+bn_lpt_1mSr(pi,sigma,rho,d)))/3.0+2.0*exp(bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_M(kyz,pi,sigma,rho,d)+bn_lpt_a(kxz,pi,sigma,rho,d)+bn_lpt_Sr(pi,sigma,rho,d)+bn_lpt_1mSr(pi,sigma,rho,d))/3.0+exp(bn_lpt_M(kyz,pi,sigma,rho,d)+bn_lpt_a(kxz,pi,sigma,rho,d)+bn_lpt_Sr(pi,sigma,rho,d))*(exp(bn_lpt_a(kxy,pi,sigma,rho,d))+1.0+exp(bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_1mSr(pi,sigma,rho,d)))/3.0+2.0*exp(bn_lpt_N(kxy,pi,sigma,rho,d)+bn_lpt_M(kyz,pi,sigma,rho,d)+bn_lpt_M(kxz,pi,sigma,rho,d)+bn_lpt_a1(pi,sigma,rho,d))/3.0+2.0*exp(bn_lpt_M(kyz,pi,sigma,rho,d))*(exp(bn_lpt_N(kxy,pi,sigma,rho,d)+bn_lpt_a(kxz,pi,sigma,rho,d)+bn_lpt_a1(pi,sigma,rho,d)+bn_lpt_Sr(pi,sigma,rho,d))+exp(bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_N(kxz,pi,sigma,rho,d)+bn_lpt_M1(pi,sigma,rho,d)+bn_lpt_1mSr(pi,sigma,rho,d)))/3.0);
            }else{     /*0 1 1 1 1 0*/
              return bn_lpt_a(kxy,pi,sigma,rho,d) + log(exp(bn_lpt_Mp1(kyz,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_M(kyz,pi,sigma,rho,d)) + 4.0*exp(bn_lpt_a(kyz,pi,sigma,rho,d)+bn_lpt_Sr(pi,sigma,rho,d))) + bn_lpt_a(kxz,pi,sigma,rho,d) - log(3.0);
            }
          }else{
            if(zx>0){  /*0 1 1 1 0 1*/
              return bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_M(kyz,pi,sigma,rho,d) + log(exp(bn_lpt_ap1(kxz,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d))) - log(3.0);
            }else{     /*0 1 1 1 0 0*/
              return bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_M(kyz,pi,sigma,rho,d)+ bn_lpt_N(kxz,pi,sigma,rho,d);
            }
          }
        }else{
          if(xz>0){  
            if(zx>0){  /*0 1 1 0 1 1*/
              return bn_lpt_a(kxy,pi,sigma,rho,d) + log(exp(bn_lpt_ap1(kyz,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d))) + bn_lpt_M(kxz,pi,sigma,rho,d) - log(3.0);
            }else{     /*0 1 1 0 1 0*/
              return bn_lpt_a(kxy,pi,sigma,rho,d) + log(exp(bn_lpt_ap1(kyz,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d))) + bn_lpt_a(kxz,pi,sigma,rho,d) - log(3.0);
            }
          }else{
            if(zx>0){  /*0 1 1 0 0 1*/
              return bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_a(kyz,pi,sigma,rho,d)+ bn_lpt_a(kxz,pi,sigma,rho,d);
            }else{     /*0 1 1 0 0 0*/
              return bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_a(kyz,pi,sigma,rho,d)+ bn_lpt_N(kxz,pi,sigma,rho,d);
            }
          }
        }
      }else{
        if(zy>0){
          if(xz>0){
            if(zx>0){  /*0 1 0 1 1 1*/
              return log(exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_ap1(kyz,pi,sigma,rho,d)) + exp(bn_lpt_ap1(kxy,pi,sigma,rho,d) + bn_lpt_a(kyz,pi,sigma,rho,d)) + exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_a(kyz,pi,sigma,rho,d)+bn_lpt_1mSr(pi,sigma,rho,d)) + 2.0*exp(bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_N(kyz,pi,sigma,rho,d)+bn_lpt_a1(pi,sigma,rho,d)) + 2.0*exp(bn_lpt_N(kxy,pi,sigma,rho,d)+bn_lpt_a1(pi,sigma,rho,d)+bn_lpt_a(kyz,pi,sigma,rho,d))) + bn_lpt_M(kxz,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d) - log(3.0);
            }else{     /*0 1 0 1 1 0*/
              return log(exp(bn_lpt_ap1(kxy,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_a(kxy,pi,sigma,rho,d) + bn_lpt_1mSr(pi,sigma,rho,d)) + 2.0*exp(bn_lpt_N(kxy,pi,sigma,rho,d) + bn_lpt_a1(pi,sigma,rho,d))) + bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d) - log(3.0);
            }
          }else{
            if(zx>0){  /*0 1 0 1 0 1*/
              return log(exp(bn_lpt_ap1(kxy,pi,sigma,rho,d))+2.0*exp(bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_1mSr(pi,sigma,rho,d))+2.0*exp(bn_lpt_N(kxy,pi,sigma,rho,d)+bn_lpt_a1(pi,sigma,rho,d))) + bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d) - log(3.0);
            }else{     /*0 1 0 1 0 0*/
              return bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_a(kyz,pi,sigma,rho,d)+ bn_lpt_N(kxz,pi,sigma,rho,d);
            }
          }
        }else{
          if(xz>0){
            if(zx>0){  /*0 1 0 0 1 1*/
              return bn_lpt_a(kxy,pi,sigma,rho,d) + log(exp(bn_lpt_Np1(kyz,pi,sigma,rho,d))+2.0*exp(bn_lpt_N(kyz,pi,sigma,rho,d)+bn_lpt_N1(pi,sigma,rho,d))) + bn_lpt_M(kxz,pi,sigma,rho,d) - log(3.0);
            }else{     /*0 1 0 0 1 0*/
              return bn_lpt_a(kxy,pi,sigma,rho,d) + log(exp(bn_lpt_Np1(kyz,pi,sigma,rho,d))+2.0*exp(bn_lpt_N(kyz,pi,sigma,rho,d)+bn_lpt_N1(pi,sigma,rho,d))) + bn_lpt_a(kxz,pi,sigma,rho,d) - log(3.0);
            }
          }else{
            if(zx>0){  /*0 1 0 0 0 1*/
              return bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_N(kyz,pi,sigma,rho,d)+ bn_lpt_a(kxz,pi,sigma,rho,d);
            }else{     /*0 1 0 0 0 0*/
              return bn_lpt_a(kxy,pi,sigma,rho,d)+bn_lpt_N(kyz,pi,sigma,rho,d)+ bn_lpt_N(kxz,pi,sigma,rho,d);
            }
          }
        }
      }
    }else{
      if(yz>0){
        if(zy>0){
          if(xz>0){
            if(zx>0){  /*0 0 1 1 1 1*/
              return log(exp(bn_lpt_Np1(kxy,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_N(kxy,pi,sigma,rho,d) + bn_lpt_N1(pi,sigma,rho,d))) + bn_lpt_M(kyz,pi,sigma,rho,d) + bn_lpt_M(kxz,pi,sigma,rho,d) - log(3.0);
            }else{     /*0 0 1 1 1 0*/
              return bn_lpt_N(kxy,pi,sigma,rho,d)+bn_lpt_M(kyz,pi,sigma,rho,d)+ bn_lpt_a(kxz,pi,sigma,rho,d);
            }
          }else{
            if(zx>0){  /*0 0 1 1 0 1*/
              return log(exp(bn_lpt_Np1(kxy,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_N(kxy,pi,sigma,rho,d) + bn_lpt_N1(pi,sigma,rho,d))) + bn_lpt_M(kyz,pi,sigma,rho,d) + bn_lpt_a(kxz,pi,sigma,rho,d) - log(3.0);
            }else{     /*0 0 1 1 0 0*/
              return bn_lpt_N(kxy,pi,sigma,rho,d)+bn_lpt_M(kyz,pi,sigma,rho,d)+ bn_lpt_N(kxz,pi,sigma,rho,d);
            }
          }
        }else{
          if(xz>0){
            if(zx>0){  /*0 0 1 0 1 1*/
              return bn_lpt_N(kxy,pi,sigma,rho,d)+bn_lpt_a(kyz,pi,sigma,rho,d)+ bn_lpt_M(kxz,pi,sigma,rho,d);
            }else{     /*0 0 1 0 1 0*/
              return bn_lpt_N(kxy,pi,sigma,rho,d)+bn_lpt_a(kyz,pi,sigma,rho,d)+ bn_lpt_a(kxz,pi,sigma,rho,d);
            }
          }else{
            if(zx>0){  /*0 0 1 0 0 1*/
              return bn_lpt_N(kxy,pi,sigma,rho,d)+bn_lpt_a(kyz,pi,sigma,rho,d)+ bn_lpt_a(kxz,pi,sigma,rho,d);
            }else{     /*0 0 1 0 0 0*/
              return bn_lpt_N(kxy,pi,sigma,rho,d)+bn_lpt_a(kyz,pi,sigma,rho,d)+ bn_lpt_N(kxz,pi,sigma,rho,d);
            }
          }
        }
      }else{
        if(zy>0){
          if(xz>0){
            if(zx>0){  /*0 0 0 1 1 1*/
              return log(exp(bn_lpt_Np1(kxy,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_N(kxy,pi,sigma,rho,d) + bn_lpt_N1(pi,sigma,rho,d))) + bn_lpt_a(kyz,pi,sigma,rho,d) + bn_lpt_M(kxz,pi,sigma,rho,d) - log(3.0);
            }else{     /*0 0 0 1 1 0*/
              return bn_lpt_N(kxy,pi,sigma,rho,d)+bn_lpt_a(kyz,pi,sigma,rho,d)+ bn_lpt_a(kxz,pi,sigma,rho,d);
            }
          }else{
            if(zx>0){  /*0 0 0 1 0 1*/
              return log(exp(bn_lpt_Np1(kxy,pi,sigma,rho,d)) + 2.0*exp(bn_lpt_N(kxy,pi,sigma,rho,d) + bn_lpt_N1(pi,sigma,rho,d))) + bn_lpt_a(kxz,pi,sigma,rho,d) + bn_lpt_a(kyz,pi,sigma,rho,d) - log(3.0);
            }else{     /*0 0 0 1 0 0*/
              return bn_lpt_N(kxy,pi,sigma,rho,d)+bn_lpt_a(kyz,pi,sigma,rho,d)+ bn_lpt_N(kxz,pi,sigma,rho,d);
            }
          }
        }else{
          if(xz>0){
            if(zx>0){  /*0 0 0 0 1 1*/
              return bn_lpt_N(kxy,pi,sigma,rho,d)+bn_lpt_N(kyz,pi,sigma,rho,d)+ bn_lpt_M(kxz,pi,sigma,rho,d);
            }else{     /*0 0 0 0 1 0*/
              return bn_lpt_N(kxy,pi,sigma,rho,d)+bn_lpt_N(kyz,pi,sigma,rho,d)+ bn_lpt_a(kxz,pi,sigma,rho,d);
            }
          }else{
            if(zx>0){  /*0 0 0 0 0 1*/
              return bn_lpt_N(kxy,pi,sigma,rho,d)+bn_lpt_N(kyz,pi,sigma,rho,d)+ bn_lpt_a(kxz,pi,sigma,rho,d);
            }else{     /*0 0 0 0 0 0*/
              return bn_lpt_N(kxy,pi,sigma,rho,d)+bn_lpt_N(kyz,pi,sigma,rho,d)+ bn_lpt_N(kxz,pi,sigma,rho,d);
            }
          }
        }
      }
    }
  }
}

void bn_dyadstats_R(int *g, double *pn, double *stats)
/*Return the matrix of sufficient dyad statistics for a biased net model,
using the method of Skvoretz et al.*/
{
  long int n,i,j,k,parents;

  /*Initialize the stats structure*/
  n=(long int)*pn;  
  for(i=0;i<n-1;i++){
    stats[i]=(double)i;
    for(j=1;j<4;j++)
      stats[i+(n-1)*j]=0.0;
  }
  /*Accumulate the stats (note -- this could be made much faster)*/
  for(i=0;i<n;i++)
    for(j=i+1;j<n;j++){
      /*Count the parents of i and j*/
      parents=0;
      for(k=0;k<n;k++)
        if((g[k+i*n]>0)&&(g[k+j*n]>0))
          parents++;
      /*Increment the appropriate statistic*/
      if(g[i+j*n]>0){
        if(g[j+i*n]>0){              /*Mutual*/
          stats[parents+(n-1)]++;
        }else{                        /*Asym*/
          stats[parents+(n-1)*2]++;
        }
      }else{
        if(g[j+i*n]>0){              /*Asym*/
          stats[parents+(n-1)*2]++;
        }else{                        /*Null*/
          stats[parents+(n-1)*3]++;
        }
      }
    }
}

void bn_triadstats_R(int *g, double *pn, double *stats)
/*Compute the matrix of parent statistics for a biased net model,
using the method of Skvoretz (2003). In particular, stats will be filled
with an adjacency matrix containing the dyadic parent counts.*/
{
  long int n,i,j,k;

  n=(long int)*pn;
  /*Accumulate parent counts*/
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      if(i<j){
        /*Accumulate parent count*/
        for(k=0;k<n;k++)
          if((g[k+i*n]>0)&&(g[k+j*n]>0))
            stats[i+j*n]++;
      }else if(i==j){
        stats[i+j*n]=0.0;       /*Treat diagonal as zero*/
      }else{
        stats[i+j*n]=stats[j+i*n];  /*Use what we already have*/
      }
}

void bn_lpl_dyad_R(double *stats, double *psr, double *pi, double *sigma, double *rho, double *d, double *lpl)
/*Return the dyadic log pseudolikelihood for a given graph under a biased net 
model.  Parameters are as follows:
  stats = a *psr x 4 matrix, whose rows contain # parents, # mutuals, 
    # asymmetrics, and # nulls (in order)
  sigma = Pr( x->y | xSy ) (Transitivity or "sibling" bias)
  pi = Pr( x->y | y->x ) (Reciprocity or "parent" bias)
  rho = Pr( x->y | y->x & xSy ) ("Double role" bias, or interaction effect)
  d = Pr( x->y ) (Baseline density effect)
  lpl = a pointer to the likelihood
*/
{
  long int sr,i;

  /*Calculate the log pseudolikelihood*/
  *lpl=0.0;
  sr=(long int)(*psr);
  for(i=0;i<sr;i++){
    *lpl+=bn_lpkm((long int)stats[i],*pi,*sigma,*rho,*d)*stats[i+sr];
    *lpl+=bn_lpka((long int)stats[i],*pi,*sigma,*rho,*d)*stats[i+2*sr];
    *lpl+=bn_lpkn((long int)stats[i],*pi,*sigma,*rho,*d)*stats[i+3*sr];
  }
}

void bn_lpl_triad_R(int *g, double *stats, double *pn, double *pi, double *sigma, double *rho, double *d, double *lpl)
/*Return the triadic log pseudolikelihood for a given graph under a biased net 
model, using the method of Skvoretz (2003).  Parameters are as follows:
  g = a *pn x *pn adjacency matrix
  stats = a *pn x *pn matrix, whose i,jth cell contains the number of parents
    of the {i,j} dyad
  sigma = Pr( x->y | xSy ) (Transitivity or "sibling" bias)
  pi = Pr( x->y | y->x ) (Reciprocity or "parent" bias)
  rho = Pr( x->y | y->x & xSy ) ("Double role" bias, or interaction effect)
  d = Pr( x->y ) (Baseline density effect)
  lpl = a pointer to the likelihood
*/
{
  long int i,j,k,n;

  n=(long int)*pn;
  /*Calculate the log pseudolikelihood*/
  *lpl=0.0;
  for(i=0;i<n;i++)
    for(j=i+1;j<n;j++)
      for(k=j+1;k<n;k++){
        *lpl+=bn_lpt(g[i+j*n],g[j+i*n],g[j+k*n],g[k+j*n],g[i+k*n],g[k+i*n],(long int)stats[i+j*n],(long int)stats[j+k*n],(long int)stats[i+k*n],*pi,*sigma,*rho, *d);
      }
}

void bn_ptriad_R(double *ppi, double *psigma, double *prho, double *pd, double *pt)
/*Return the probability distribution for the triad census under a biased net 
model, using the "lone triad" approximation of Skvoretz et al., 2004.
Parameters are as follows:
  stats = a vector containing the triad census
  sigma = Pr( x->y | xSy ) (Transitivity or "sibling" bias)
  pi = Pr( x->y | y->x ) (Reciprocity or "parent" bias)
  rho = Pr( x->y | y->x & xSy ) ("Double role" bias, or interaction effect)
  d = Pr( x->y ) (Baseline density effect)
  pt = a pointer to the triad probabilities
*/
{
  double pi,sigma,rho,d,M0,a0,N0,M1,a1,N1,Mp1,ap1,Np1,Sr;

  /*Initialize for convenience*/
  pi=*ppi;
  sigma=*psigma;
  rho=*prho;
  d=*pd;

  /*Compute the simplifying factors*/
  M0=d*(pi+(1.0-pi)*d);
  a0=d*(1.0-d)*(1.0-pi);
  N0=(1.0-d)*(1.0-d*(1.0-pi));
  M1=(sigma+(1.0-sigma)*d)*(1.0-(1.0-pi)*(1.0-sigma)*(1.0-rho)*(1.0-d));
  a1=(sigma+(1.0-sigma)*d)*(1.0-pi)*(1.0-sigma)*(1.0-rho)*(1.0-d);
  N1=1.0-(sigma+(1.0-sigma)*d)*(1.0+(1.0-pi)*(1.0-sigma)*(1.0-rho)*(1.0-d));
  Mp1=sigma*(1.0-(1.0-sigma)*(1.0-rho));
  ap1=sigma*(1.0-sigma)*(1.0-rho);
  Np1=1.0-sigma*(1.0-(1.0-sigma)*(1.0-rho)+2.0*(1.0-sigma)*(1.0-rho));
  Sr=1.0-(1.0-sigma)*(1.0-rho);
  
  /*Compute the actual triad probabilities*/
  pt[0]=N0*N0*N0;   /*003*/
  pt[1]=6.0*a0*N0*N0;   /*012*/
  pt[2]=3.0*M0*N0*N0;   /*102*/
  pt[3]=a0*a0*(N1+2.0*N0*Np1);   /*021D*/  /*This fixes a typo in SN article*/
  pt[4]=3.0*a0*a0*N0;   /*021U*/
  pt[5]=6.0*a0*a0*N0;   /*021C*/
  pt[6]=6.0*M0*a0*N0;   /*111D*/
  pt[7]=2.0*M0*a0*(N1+2.0*N0*Np1);   /*111U*/
  pt[8]=2.0*a0*a0*(a1+2.0*a0*(1-Sr)+2.0*N0*ap1);   /*030T*/
  pt[9]=2.0*a0*a0*a0;   /*030C*/
  pt[10]=M0*M0*(N1+2*N0*Np1);   /*201*/
  pt[11]=a0*a0*(M1+2.0*M0+2.0*N0*Mp1+4.0*a0*Sr);   /*120D*/
  pt[12]=M0*a0*(1.0-Sr)*(2.0*a1+a0*(1.0-Sr)+4.0*N0*ap1);   /*120U*/
  pt[13]=2.0*M0*a0*(a1+2.0*a0*(1-Sr)+2*N0*ap1);   /*120C*/
  pt[14]=M0*(2.0*M0*a1+2.0*a0*M1*(1.0-Sr)+2.0*M0*a0*(1.0-Sr)+ 4.0*a0*N0*ap1*Sr+4.0*a0*N0*Mp1*(1.0-Sr)+4.0*M0*N0*ap1+2.0*a0*a1*Sr+6.0*a0*a0*Sr*(1.0-Sr));   /*210*/
  pt[15]=M0*(M0*M1+4.0*a0*N0*Mp1*Sr+2.0*M0*N0*Mp1+5*a0*a0*Sr*Sr+ 2.0*a0*M1*Sr+2.0*M0*a0*Sr);   /*300*/
}

