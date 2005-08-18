######################################################################
#
# sna-operators.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 4/23/05
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains operators associated with the sna package.
#
# Contents:
#   %c%
#   gapply
#   gliop
#
######################################################################


#%c% - Composition of two adjacancy matrices
"%c%"<-function(x,y){
  #Pre-process the raw input
  x<-as.sociomatrix.sna(x)
  y<-as.sociomatrix.sna(y)
  if(!(is.matrix(x)&&is.matrix(y)))
    stop("Single graphs required for composition.")
  #End pre-processing
  round((x%*%y)>0)
}

#gapply - Apply a function to vertex neighborhoods within a graph
gapply<-function(X,MARGIN,STATS,FUN,...,mode="digraph",diag=FALSE,distance=1,thresh=0,simplify=TRUE){
  #Pre-process the raw input
  X<-as.sociomatrix.sna(X)
  if(is.list(X))
    return(lapply(X,gapply,MARGIN,STATS,FUN,...,mode=mode, diag=diag,distance=distance,thresh=thresh,simplify=simplify))
  else if(length(dim(X))>2){
    return(apply(X,1,gapply,MARGIN,STATS,FUN,...,mode=mode, diag=diag,distance=distance,thresh=thresh,simplify=simplify))
  }
  #End pre-processing
  #Match the input function
  fun<-match.fun(FUN)
  #Dichotomize, if needed
  X<-X>thresh
  #If needed, calculate the reachability graph
  if(distance>1)
    X<-geodist(X,inf.replace=Inf)$gdist<=distance
  #Remove unwanted elements
  if(!diag)
    diag(X)<-FALSE
  if(mode=="graph")
    X[lower.tri(X)]<-FALSE
  #Extract the relevant stats
  if(!is.matrix(STATS))
    STATS<-matrix(STATS,nc=1)
  if(length(MARGIN)==1){
    if(MARGIN==1)
      stats<-apply(X,1,function(x){STATS[x,]})
    else if(MARGIN==2)
      stats<-apply(X,2,function(x){STATS[x,]})
  }else if(all(c(1,2)%in%MARGIN))
    stats<-apply(symmetrize(X,rule="weak")>0,1,function(x){STATS[x,]})
  else
    stop("MARGIN must be one of 1, 2, or c(1,2) in gapply.  Exiting.\n")
  #Apply the function and return the result
  if(is.matrix(stats))
    apply(stats,2,fun,...)
  else
    sapply(stats,fun,...,simplify=simplify)
}


#gliop - Return a binary operation on GLI values computed on two graphs (for 
#test routines).
gliop<-function(dat,GFUN,OP="-",g1=1,g2=2,...){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if((!is.list(dat))&&(length(dim(dat))==2))
     dat<-array(dat,dim=c(1,dim(dat)))
   #End pre-processing
   fun<-match.fun(GFUN)
   op<-match.fun(OP)
   if(is.list(dat))
     op(fun(dat[[g1]],...),fun(dat[[g2]],...))
   else
     op(fun(dat[g1,,],...),fun(dat[g2,,],...))
}
