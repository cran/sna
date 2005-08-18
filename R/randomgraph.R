######################################################################
#
# randomgraph.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 8/8/05
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains various routines for random graph generation in
# R.
#
# Contents:
#   rewire.ud
#   rewire.ws
#   rgbn
#   rgmn
#   rgraph
#   rguman
#   rgws
#
######################################################################


#rewire.ud - Perform a uniform dyadic rewiring of a graph or graph stack
rewire.ud<-function(g,p){
  #Pre-process the raw input
  g<-as.sociomatrix.sna(g)
  if(is.list(g))
    return(lapply(g,rewire.ud,p=p))
  #End pre-processing
  #Coerce g to an array
  if(length(dim(g))==2)
    g<-array(g,dim=c(1,NROW(g),NCOL(g)))
  n<-dim(g)[1]
  nv<-dim(g)[2]
  #Perform the rewiring, and return the result
  rewired<-.C("udrewire_R",g=as.double(g),as.double(n),as.double(nv), as.double(p),PACKAGE="sna")
  array(rewired$g,dim=c(n,nv,nv))
}


#rewire.ws - Perform a Watts-Strogatz rewiring of a graph or graph stack
rewire.ws<-function(g,p){
  #Pre-process the raw input
  g<-as.sociomatrix.sna(g)
  if(is.list(g))
    return(lapply(g,rewire.ud,p=p))
  #End pre-processing
  #Coerce g to an array
  if(length(dim(g))==2)
    gi<-array(g,dim=c(1,NROW(g),NCOL(g)))
  go<-gi
  n<-dim(gi)[1]
  nv<-dim(gi)[2]
  #Perform the rewiring, and return the result
  rewired<-.C("wsrewire_R",as.double(gi),go=as.double(go),as.double(n), as.double(nv),as.double(p),PACKAGE="sna")
  array(rewired$go,dim=c(n,nv,nv))
}


#rgbn - Draw from a biased net model
rgbn<-function(n,nv,param=list(pi=0,sigma=0,rho=0,d=0.5),burn=nv*nv*1e3,thin=nv){
  #Allocate memory for the graphs
  g<-array(0,dim=c(n,nv,nv))
  #Get the parameter vector
  p<-rep(0,4)
  if(!is.null(param$pi))
    p[1]<-param$pi
  if(!is.null(param$sigma))
    p[2]<-param$sigma
  if(!is.null(param$rho))
    p[3]<-param$rho
  if(!is.null(param$d))
    p[4]<-param$d
  #Take the draws
  g<-array(.C("bn_mcmc_R",g=as.integer(g),as.double(nv),as.double(n), as.double(burn),as.integer(thin),as.double(p[1]),as.double(p[2]),as.double(p[3]),as.double(p[4]), PACKAGE="sna")$g,dim=c(n,nv,nv))
  #Return the result
  if(dim(g)[1]==1)
    g[1,,]
  else
    g
}


#rgmn - Draw a density-conditioned graph
rgnm<-function(n,nv,m,mode="digraph",diag=FALSE){
  #Allocate the graph stack
  g<-array(0,dim=c(n,nv,nv))
  #Draw numbers for edge placement
  if(mode=="graph"){
    enum<-matrix(1:nv^2,nv,nv)
    enum<-enum[lower.tri(enum,diag=diag)]
  }else if(mode=="digraph"){
    enum<-matrix(1:nv^2,nv,nv)
    diag(enum)<-NA
    enum<-enum[!is.na(enum)]
  }else
    stop("Unsupported mode in rgnm.")
  #Place the edges
  for(i in 1:n){
    g[i,,][sample(enum,m)]<-1
    if(mode=="graph")
      g[i,,][upper.tri(g[i,,])]<-t(g[i,,])[upper.tri(g[i,,])]
  }
  #Return the results
  if(n>1)
    g
  else
    g[1,,]
}


#rgraph - Draw a Bernoulli graph.
rgraph<-function(n,m=1,tprob=0.5,mode="digraph",diag=FALSE,replace=FALSE,tielist=NULL){
   if(m==1){
      if(is.null(tielist)){
         if(length(dim(tprob))==2){
            g<-matrix(ncol=n,nrow=n)
            for(i in 1:n)
               for(j in 1:n)
                  g[i,j]<-sample(0:1,1,replace=TRUE,prob=c(1-tprob[i,j],tprob[i,j]))
         }else{
            g<-matrix(data=sample(0:1,n*n,replace=TRUE,prob=c(1-tprob,tprob)),ncol=n,nrow=n)
         }
      }else{
         g<-matrix(data=sample(as.vector(tielist),n*n,replace=replace),ncol=n,nrow=n)
      }
      if(!diag)
         diag(g)<-0
      if(mode!="digraph")
#         for(i in 1:n)
#            g[i,c(1:i)[-i]]<-g[c(1:i)[-i],i]
         g[upper.tri(g)]<-t(g)[upper.tri(g)]
   }else{
      g<-array(dim=c(m,n,n))
      if(is.null(tielist)){
         if(length(dim(tprob))==2){
            for(i in 1:m)
               for(j in 1:n)
                  for(k in 1:n)
                     g[i,j,k]<-sample(0:1,1,replace=TRUE,prob=c(1-tprob[j,k],tprob[j,k]))
         }else if(length(tprob)==m){
            for(i in 1:m){
               g[i,,]<-array(sample(0:1,n*n,replace=TRUE,prob=c(1-tprob[i],tprob[i])),dim=c(n,n))
            }
         }else{
            for(i in 1:m){
               g[i,,]<-array(sample(0:1,n*n,replace=TRUE,prob=c(1-tprob,tprob)),dim=c(n,n))
            }
         }
      }else{
         if(length(dim(tielist))==3){
            for(i in 1:m)
               g[i,,]<-array(sample(as.vector(tielist[i,,]),n*n,replace=replace),dim=c(n,n))
         }else{
            for(i in 1:m)
               g[i,,]<-array(sample(as.vector(tielist),n*n,replace=replace),dim=c(n,n))
         }
      }
      if(!diag)
         for(i in 1:m)
            diag(g[i,,])<-0
      if(mode!="digraph")
         for(i in 1:m){
#            for(j in 1:n)
#               g[i,j,c(1:j)[-j]]<-g[i,c(1:j)[-j],j]
            temp<-g[i,,]
            temp[upper.tri(temp)]<-t(temp)[upper.tri(temp)]
            g[i,,]<-temp
         }
   }
   g
}


#rguman - Draw from the U|MAN graph distribution
rguman<-function(n,nv,mut=0.25,asym=0.5,null=0.25,method=c("probability","exact")){
  #Create the output structure
  g<-array(0,dim=c(n,nv,nv))
  #Create the dyad list
  dl<-matrix(1:(nv^2),nv,nv)
  dlu<-dl[upper.tri(dl)]
  dll<-t(dl)[upper.tri(dl)]
  ndl<-length(dlu)      #Number of dyads
  #Perform a reality check
  if((match.arg(method)=="exact")&&(mut+asym+null!=ndl))
    stop("Sum of dyad counts must equal number of dyads for method==exact.\n")
  else if((match.arg(method)=="probability")&&(mut+asym+null!=1)){
      s<-mut+asym+null
      mut<-mut/s; asym<-asym/s; null<-null/s
    }    
  #Draw the graphs
  for(i in 1:n){
    #Determine the number of dyads in each class
    if(match.arg(method)=="probability"){
      mc<-rbinom(1,ndl,mut)
      ac<-rbinom(1,ndl-mc,asym/(asym+null))
      nc<-ndl-mc-ac
    }else{
      mc<-mut
      ac<-asym
      nc<-null
    }
    #Draw the dyad states 
    ds<-sample(rep(1:3,times=c(mc,ac,nc)))
    #Place edges accordingly
    if(mc>0){
      g[i,,][dlu[ds==1]]<-1                      #Mutuals
      g[i,,][dll[ds==1]]<-1
    }
    if(ac>0){
      g[i,,][dlu[ds==2]]<-rbinom(ac,1,0.5)       #Asymetrics
      g[i,,][dll[ds==2]]<-1-g[i,,][dlu[ds==2]]
    }
  }
  #Return the result
  if(n>1)
    g
  else
    g[1,,]
}


#rgws - Draw a graph from the Watts-Strogatz model
rgws<-function(n,nv,d,z,p){
  #Begin by creating the lattice
  tnv<-nv^d
  temp<-vector()
  nums<-1:nv
  count<-tnv/nv
  for(i in 1:d){
    temp<-cbind(temp,rep(nums,count))
    nums<-rep(nums,each=nv)
    count<-count/nv
  }
  lat<-as.matrix(dist(temp,method="manhattan"))<=z  #Identify nearest neighbors
  diag(lat)<-0
  #Create n copies of the lattice
  if(n>1)
    lat<-apply(lat,c(1,2),rep,n)
  else
    lat<-array(lat,dim=c(1,tnv,tnv))
  #Rewire the copies
  g<-lat
  lat<-array(.C("wsrewire_R",as.double(lat),g=as.double(g),as.double(n), as.double(tnv),as.double(p),PACKAGE="sna")$g,dim=c(n,tnv,tnv))
  #Return the result
  if(n>1)
    lat
  else
    lat[1,,]
}

