######################################################################
#
# gli.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 8/8/05
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains various routines for the calculation of 
# graph-level indices.
#
# Contents:
#   centralization
#   connectedness
#   dyad.census
#   efficiency
#   gden
#   grecip
#   gtrans
#   hierarchy
#   lubness
#   mutuality
#   triad.census
#   triad.classify
#
######################################################################


#centralization - Find the centralization of a graph (for some arbitrary 
#centrality measure)
centralization<-function(dat,FUN,g=1,mode="digraph",diag=FALSE,normalize=TRUE,...){
   #Find the centrality function
   fun<-match.fun(FUN)
   #Grab the vector of centralities
   cv<-fun(dat,g=g,gmode=mode,diag=diag,...)
   #Find the empirical maximum
   cmax<-max(cv)
   #Now, for the absolute deviations....
   cent<-sum(cmax-cv)
   #If we're normalizing, we'll need to get the theoretical max from our centrality function
   if(normalize)
      cent<-cent/fun(dat,g=g,gmode=mode,diag=diag,tmaxdev=TRUE,...)
   #Return the centralization
   cent
}


#connectedness - Find the Krackhardt connectedness of a graph or graph stack
connectedness<-function(dat,g=NULL){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat)){
     if(is.null(g))
       g<-1:length(dat)
     return(sapply(dat[g],connectedness))
   }
   #End pre-processing
   if((!is.null(g))&&(length(dim(dat))>2))
     dat<-dat[g,,]
   if(length(dim(dat))>2)
     con<-apply(dat,1,function(x){r<-reachability(symmetrize(x,rule="weak")); gden(r,diag=FALSE)})
   else
     con<-gden(reachability(symmetrize(dat,rule="weak")),diag=FALSE)
   #Return the result
   con
}


#dyad.census - Return the Holland and Leinhardt MAN dyad census for a given 
#graph or graph stack
dyad.census<-function(dat,g=NULL){
   #Define an internal function
   intcalc<-function(m,meas){
      switch(meas,
         mut=sum(m[upper.tri(m)]&t(m)[upper.tri(m)],na.rm=TRUE),
         asym=sum(xor(m[upper.tri(m)],t(m)[upper.tri(m)]),na.rm=TRUE),
         null=sum(!m[upper.tri(m)]&!t(m)[upper.tri(m)],na.rm=TRUE)
      )
   }
   #Organize the data
   dat<-as.sociomatrix.sna(dat)
   #Perform the census
   if(is.list(dat)){
     if(is.null(g))
       g<-1:length(dat)
     man<-cbind(sapply(dat,intcalc,"mut"),sapply(dat,intcalc,"asym"), sapply(dat,intcalc,"null"))
   }else if(length(dim(dat))>2){
     if(is.null(g))
       g<-1:dim(dat)[1]
     dat<-dat[g,,,drop=FALSE]
     man<-cbind(apply(dat,1,intcalc,"mut"),apply(dat,1,intcalc,"asym"), apply(dat,1,intcalc,"null"))
   }else{
     man<-cbind(intcalc(dat,"mut"),intcalc(dat,"asym"),intcalc(dat,"null"))
   }
   colnames(man)<-c("Mut","Asym","Null")
   #Return the result
   man
}


#efficiency - Find the Krackhardt efficiency of a graph or graph stack
efficiency<-function(dat,g=NULL,diag=FALSE){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat)){
     if(is.null(g))
       g<-1:length(dat)
     return(sapply(dat[g],efficiency,diag=diag))
   }
   #End pre-processing
   #Define an internal function, for convenience
   inteff<-function(g,diag){
      comsz<-component.dist(g,connected="weak")$csize
      reqedge<-sum(comsz-1)     #Get required edges
      maxv<-sum(comsz*(comsz-(!diag))-(comsz-1)) #Get maximum violations
      if(!diag)
         g<-diag.remove(g)
      edgec<-sum(g,na.rm=TRUE)  #Get count of actual edges
      1-(edgec-reqedge)/maxv
   }
   #Perform the actual calculation
   if(length(dim(dat))>2){
      if(is.null(g))
        g<-1:dim(dat)[1]
      eff<-apply(dat[g,,,drop=FALSE],1,inteff,diag=diag)
   }else
      eff<-inteff(dat,diag=diag)
   #Return the result
   eff
}


#gden - Compute the density of an input graph or graph stack.
gden<-function(dat,g=NULL,diag=FALSE,mode="digraph"){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat)){
     if(is.null(g))
       g<-1:length(dat)
     return(sapply(dat[g],gden,diag=diag,mode=mode))
   }
   #End pre-processing
   n<-dim(dat)[2]
   if(length(dim(dat))>2){     #Is this a stack?
      if(!is.null(g)){                 #Were individual graphs selected?
         gn<-length(g)
         d<-dat[g,,]
      }else{
         d<-dat
         gn<-dim(dat)[1]
      }
   }else{
      d<-dat
      gn<-1
   }
   if(gn==1){     #Only one graph - convert to stack format
      temp<-array(dim=c(1,n,n))
      temp[1,,]<-d
      d<-temp
   }
   if(!diag)           #If not using the diagonal, remove it
      d<-diag.remove(d)
   if(mode=="graph")    #If this is a simple graph, remove one triangle of each matrix
      d<-upper.tri.remove(d)
   #Find number of ties, and counts
   count<-apply(d,1,sum,na.rm=TRUE)
   count/nties(d[1,,],mode=mode,diag=diag)
}


#grecip - Compute the reciprocity of an input graph or graph stack.
grecip<-function(dat,g=NULL,measure=c("dyadic","dyadic.nonnull","edgewise")){
  #Pre-process the raw input
  dat<-as.sociomatrix.sna(dat)
  if(is.list(dat)){
    if(is.null(g))
      g<-1:length(dat)
    return(sapply(dat[g],grecip,measure=measure))
  }
  #End pre-processing
  m<-stackcount(dat)         #How many graphs are there?
  if(is.null(g))             #Create g, if needed
    g<-1:m
  dc<-dyad.census(dat,g=g)   #Obtain the dyad census for all specified graphs
  #Return the appropriate measure
  switch(match.arg(measure),
    dyadic=(dc[,1]+dc[,3])/(dc[,1]+dc[,2]+dc[,3]),
    dyadic.nonnull=dc[,1]/(dc[,1]+dc[,2]),
    edgewise=2*dc[,1]/(2*dc[,1]+dc[,2])
  )
}


#gtrans - Compute the transitivity of an input graph or graph stack.
gtrans<-function(dat,g=NULL,diag=FALSE,mode="digraph",measure=c("weak","strong","weakcensus","strongcensus")){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat)){
     if(is.null(g))
       g<-1:length(dat)
     return(sapply(dat[g],gtrans,diag=diag,mode=mode,measure=measure))
   }
   #End pre-processing
   n<-dim(dat)[2]
   if(length(dim(dat))>2){     #Is this a stack?
      if(!is.null(g)){                 #Were individual graphs selected?
         gn<-length(g)
         d<-dat[g,,]
      }else{
         d<-dat
         gn<-dim(dat)[1]
      }
   }else{
      d<-dat
      gn<-1
   }
   if(gn==1){     #Only one graph - convert to stack format
      temp<-array(dim=c(1,n,n))
      temp[1,,]<-d
      d<-temp
   }
   if(!diag)           #If not using the diagonal, remove it
      d<-diag.remove(d,remove.val=0)
   #Compute the appropriate transitivity indices
   t<-vector()
   for(i in 1:gn){
      #Prepare the transitivity test matrices
      dsqt<-(d[i,,]%*%d[i,,])
      dt<-d[i,,]>0
      #NA the diagonal, if needed
      if(!diag){
         diag(dt)<-NA
         diag(dsqt)<-NA
       }
      #Compute the transitivity
      t[i]<-switch(match.arg(measure),
         strong=sum(dt*dsqt+(!dt)*(NCOL(d[i,,])-2-dsqt),na.rm=TRUE) / (choose(NCOL(d[i,,]),3)*6),
         strongcensus=sum(dt*dsqt+(!dt)*(NCOL(d[i,,])-2-dsqt),na.rm=TRUE),
         weak=sum(dt*dsqt,na.rm=TRUE)/sum(dsqt,na.rm=TRUE),
         weakcensus=sum(dt*dsqt,na.rm=TRUE)
      )
      if(is.nan(t[i]))  #By convention, map undefined case to 1
        t[i]<-1
   }
   #Return the result
   t
}


#hierarchy - Find the hierarchy score of a graph or graph stack
hierarchy<-function(dat,g=NULL,measure=c("reciprocity","krackhardt")){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat)){
     if(is.null(g))
       g<-1:length(dat)
     return(sapply(dat[g],hierarchy,measure=measure))
   }
   #End pre-processing
   if(is.null(g))
     g<-1:stackcount(dat)
   if(match.arg(measure)=="reciprocity")  #Use reciprocity scores
         h<-1-grecip(dat,g)
   else if(match.arg(measure)=="krackhardt"){ #Calculate the Krackhardt reciprocity
      d<-array(dim=c(length(g),dim(dat)[2],dim(dat)[2]))
      if(length(dim(dat))>2)
         d<-dat[g,,,drop=FALSE]
      else
         d[1,,]<-dat
      h<-1-apply(d,1,function(x){r<-reachability(x); grecip(r,measure="dyadic.nonnull")})
   }
   #Return the result
   h
}


#lubness - Find Krackhardt's Least Upper Boundedness of a graph or graph stack
lubness<-function(dat,g=NULL){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat)){
     if(is.null(g))
       g<-1:length(dat)
     return(sapply(dat[g],lubness))
   }
   #End pre-processing
   #Define an internal function, for convenience
   intlub<-function(g){
      r<-reachability(g)    #Get reachability (in directed paths) of g
      cd<-component.dist(g,connected="weak")  #Get weak components of g
      nolub<-0
      maxnolub<-0
      for(i in 1:max(cd$membership)){   #Walk through the components
         vi<-(1:dim(g)[1])[cd$membership==i]  #Get the vertices of component i
         if(length(vi)>2){  #Components must be of size 3
           #Accumulate violations
           viol<-as.double(0)
           viol<-.C("lubness_con_R",as.double(g[vi,vi]), as.double(length(vi)),as.integer(r[vi,vi]),viol=viol,PACKAGE="sna")$viol
           nolub<-nolub+viol
           #Also accumulate maximum violations
           maxnolub<-maxnolub+(length(vi)-1)*(length(vi)-2)/2 
         }
      }
      #Return 1-violations/max(violations)
      1-nolub/maxnolub
   }
   #Perform the actual calculation
   if(length(dim(dat))>2){
      if(!is.null(g))
        dat<-dat[g,,,drop=FALSE]
      lub<-apply(dat,1,intlub)
   }
   else
      lub<-intlub(dat)
   #Return the result
   lub
}


#mutuality - Find the number of mutual (i.e., reciprocated) edges in a graph
mutuality<-function(dat,g=NULL){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat)){
     if(is.null(g))
       g<-1:length(dat)
     return(sapply(dat[g],mutuality))
   }
   #End pre-processing
   n<-dim(dat)[2]
   if(length(dim(dat))>2){     #Is this a stack?
      if(!is.null(g)){                 #Were individual graphs selected?
         gn<-length(g)
         d<-dat[g,,]
      }else{
         d<-dat
         gn<-dim(dat)[1]
      }
   }else{
      d<-dat
      gn<-1
   }
   if(gn==1){     #Only one graph - convert to stack format
      temp<-array(dim=c(1,n,n))
      temp[1,,]<-d
      d<-temp
   }
   #Find numbers of mutuals
   m<-apply(d,1,function(a){sum(a[upper.tri(a)]*t(a)[upper.tri(a)],na.rm=TRUE)})
   m
}


#triad.census - Conduct a Davis and Leinhardt triad census for a graph or graph stack
triad.census<-function(dat,g=NULL,mode=c("digraph","graph")){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat)){
     if(is.null(g))
       g<-1:length(dat)
     return(sapply(dat[g],triad.census,mode=mode))
   }
   #End pre-processing
   #First, define the triad class vector
   tc<-switch(match.arg(mode),
     graph=0:3,
     digraph=c("003","012","102","021D","021U","021C","111D","111U","030T", "030C","201","120D","120U","120C","210","300")
   )
   #Organize the data
   if(length(dim(dat))>2){
      if(is.null(g))
        d<-dat
      else
        d<-dat[g,,,drop=FALSE]
      rnam<-dimnames(d)[[1]]
   }else{
      d<-array(dim=c(1,dim(dat)[1],dim(dat)[2]))
      d[1,,]<-dat
      rnam<-NULL
   }
   d<-diag.remove(d,remove.val=0)  #Remove any diagonals
   #Obtain triad census scores
   tcm<-vector()
   n<-as.integer(dim(d)[2])
   for(i in 1:dim(d)[1]){
     tcv<-as.double(rep(0,length(tc)))
     gm<-as.integer(switch(match.arg(mode),graph=0,digraph=1))
     tcv<-.C("triad_census_R",as.integer(d[i,,]),n,tcv=tcv,gm,PACKAGE="sna")$tcv
     tcm<-rbind(tcm,tcv)
   }
   colnames(tcm)<-tc
   rownames(tcm)<-rnam
   #Return the result
   tcm
}


#triad.classify - Return the Davis and Leinhardt classification of a given triad
triad.classify<-function(dat,g=1,tri=c(1,2,3),mode=c("digraph","graph")){
   #Zeroth step: extract the triad
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
      d<-dat[[g]][tri,tri]
   else if(length(dim(dat))==2)
      d<-dat[tri,tri]
   else
      d<-dat[g,tri,tri]
   #First, classify as NA if any entries are missing
   if(any(is.na(d[upper.tri(d)|lower.tri(d)])))
      return(NA)
   #Next, define the triad class vector
   tc<-switch(match.arg(mode),
     graph=0:3,
     digraph=c("003","012","102","021D","021U","021C","111D","111U","030T", "030C","201","120D","120U","120C","210","300")
   )
   #Classify the triad
   tt<-as.integer(0)
   gm<-as.integer(switch(match.arg(mode),graph=0,digraph=1))
   tt<-.C("triad_classify_R",as.integer(d),tt=tt,gm,PACKAGE="sna")$tt
   tc[tt]
}
