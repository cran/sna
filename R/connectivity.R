######################################################################
#
# connectivity.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 12/26/04
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains various routines associated with connectivity
# properties (including geodesic distance and friends).
#
# Contents:
#  component
#  component.dist
#  geodist
#  isolates
#  is.connected
#  is.isolate
#  reachability
#
######################################################################


#component.dist - Returns a data frame containing a vector of length n such that
#the ith element contains the number of components of G having size i, and a 
#vector of length n giving component membership.  Component strength is 
#determined by the rule which is used to symmetrize the matrix; this controlled 
#by the eponymous parameter given to the symmetrize command.
component.dist<-function(dat,connected=c("strong","weak","unilateral","recursive")){
   n<-dim(dat)[2]
   #Symmetrize dat based on the connectedness rule
   if(any(dat!=t(dat)))  #Don't bother with this unless we need to do so
      dat<-switch(match.arg(connected),
         "weak"=symmetrize(dat,rule="weak"),
         "unilateral"=reachability(dat),
         "strong"=symmetrize(reachability(dat),rule="strong"),
         "recursive"=symmetrize(dat,rule="strong")
      )
   #Warn of non-uniqueness in the unilateral case, if need be
   if(match.arg(connected)=="unilateral")
      if(any(dat!=t(dat)))
         warning("Nonunique unilateral component partition detected in component.dist.  Problem vertices will be arbitrarily assigned to one of their components.\n")
   #Perform initial setup
   membership<-rep(0,n)
   #Call the C routine, which performs a fast BFS
   membership<-.C("component_dist_R",as.double(dat),as.double(n), membership=as.double(membership),PACKAGE="sna")$membership
   #Return the results
   o<-list()
   o$membership<-membership          #Copy memberships
   o$csize<-vector()
   for(i in 1:max(membership))           #Extract component sizes
      o$csize[i]<-length(membership[membership==i])
   o$cdist<-vector()
   for(i in 1:n)                                     #Find component size distribution
      o$cdist[i]<-length(o$csize[o$csize==i])
   o
}


#components - Find the number of (maximal) components within a given graph
components<-function(dat,connected="strong",comp.dist.precomp=NULL){
   #Use component.dist to get the distribution
   if(!is.null(comp.dist.precomp))
      cd<-comp.dist.precomp
   else
      cd<-component.dist(dat,connected=connected)
   #Return the result
   length(unique(cd$membership))
}


#geodist - Find the numbers and lengths of geodesics among nodes in a graph 
#using a BFS, a la Brandes (2000).  (Thanks, Ulrik!)
geodist<-function(dat,inf.replace=Inf){
   n<-dim(dat)[2]
   #Initialize the matrices
   sigma<-matrix(0,nrow=n,ncol=n)
   gd<-matrix(Inf,nrow=n,ncol=n)
   #Perform the calculation
   geo<-.C("geodist_R",as.double(dat),as.double(n),gd=as.double(gd), sigma=as.double(sigma),NAOK=TRUE,PACKAGE="sna")
   #Return the results
   o<-list()
   o$counts<-matrix(geo$sigma,n,n)
   o$gdist<-matrix(geo$gd,n,n)
   o$gdist[o$gdist==Inf]<-inf.replace  #Patch Infs, if desired
   o
}


#isolates - Returns a list of the isolates in a given graph or stack
isolates<-function(dat,diag=FALSE){
   if(length(dim(dat))>2){
      o<-vector()
      for(g in 1:dim(dat)[1])
         o<-c(o,list(seq(1:dim(dat)[2])[is.isolate(dat,g=g,ego=1:dim(dat)[2],diag=diag)]))
   }else
      o<-seq(1:dim(dat)[2])[is.isolate(dat,ego=1:dim(dat)[2],diag=diag)]
   o
}


#is.connected - Determine whether or not one or more graphs are connected
is.connected<-function(g,connected="strong",comp.dist.precomp=NULL){
  #Calculate numbers of components
  if(is.matrix(g)){
    comp<-components(g,connected=connected,comp.dist.precomp=comp.dist.precomp)
  }else{
    comp<-apply(g,1,components,connected=connected)
  }
  #Return the result
  comp==1
}


#is.isolate - Returns TRUE iff ego is an isolate
is.isolate<-function(dat,ego,g=1,diag=FALSE){
   if(length(dim(dat))>2)
      d<-dat[g,,]
   else
      d<-dat
   if(!diag)
      diag(d)<-NA
   o<-vector()
   for(i in 1:length(ego))
      o<-c(o,all(is.na(d[ego[i],])|(d[ego[i],]==0))&all(is.na(d[,ego[i]])|(d[,ego[i]]==0)))
   o   
}


#reachability - Find the reachability matrix of a graph.
reachability<-function(dat,geodist.precomp=NULL){
   #Get the counts matrix
   if(is.null(geodist.precomp))
      cnt<-geodist(dat)$counts
   else
      cnt<-geodist.precomp$counts
   #Dichotomize and return
   apply(cnt>0,c(1,2),as.numeric)
}
