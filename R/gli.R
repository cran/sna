######################################################################
#
# gli.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 1/09/04
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
connectedness<-function(dat,g=1:stackcount(dat)){
   d<-array(dim=c(length(g),dim(dat)[2],dim(dat)[2]))
   if(length(dim(dat))>2)
      d<-dat[g,,,drop=FALSE]
   else
      d[1,,]<-dat
   con<-apply(d,1,function(x){r<-reachability(symmetrize(x,rule="weak")); gden(r,diag=FALSE)})
   #Return the result
   con
}


#dyad.census - Return the Holland and Leinhardt MAN dyad census for a given 
#graph or graph stack
dyad.census<-function(dat,g=1:stackcount(dat)){
   #Define an internal function
   intcalc<-function(m,meas){
      switch(meas,
         mut=sum(m[upper.tri(m)]&t(m)[upper.tri(m)],na.rm=TRUE),
         asym=sum(xor(m[upper.tri(m)],t(m)[upper.tri(m)]),na.rm=TRUE),
         null=sum(!m[upper.tri(m)]&!t(m)[upper.tri(m)],na.rm=TRUE)
      )
   }
   #Organize the data
   if(length(dim(dat))>2)
      d<-dat[g,,,drop=FALSE]
   else{
      d<-array(dim=c(1,dim(dat)[1],dim(dat)[2]))
      d[1,,]<-dat
   }
   #Perform the census
   man<-cbind(apply(d,1,intcalc,"mut"),apply(d,1,intcalc,"asym"),apply(d,1,intcalc,"null"))
   colnames(man)<-c("Mut","Asym","Null")
   #Return the result
   man
}


#efficiency - Find the Krackhardt efficiency of a graph or graph stack
efficiency<-function(dat,g=1:stackcount(dat),diag=FALSE){
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
   d<-array(dim=c(length(g),dim(dat)[2],dim(dat)[2]))
   if(length(dim(dat))>2)
      d<-dat[g,,,drop=FALSE]
   else
      d[1,,]<-dat
   eff<-apply(d,1,inteff,diag=diag)
   #Return the result
   eff
}


#gden - Compute the density of an input graph or graph stack.
gden<-function(dat,g=NULL,diag=FALSE,mode="digraph"){
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
grecip<-function(dat,g=NULL,measure=c("dyadic","edgewise")){
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
   #Find mean lack of reciprocation
   gr<-vector()
   for(i in 1:gn){
      temp<-d[i,,]
      gr[i]<-switch(match.arg(measure),
         dyadic=1-mean(abs(temp-t(temp))[upper.tri(temp)],na.rm=TRUE),
         edgewise=1-mean(abs(temp-t(temp))[upper.tri(temp)&(temp|t(temp))],na.rm=TRUE)
      )
   }
   gr
}


#gtrans - Compute the transitivity of an input graph or graph stack.
gtrans<-function(dat,g=NULL,diag=FALSE,mode="digraph",measure=c("weak","strong","weakcensus","strongcensus")){
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
   if(mode=="graph")    #If this is a simple graph, remove one triangle of each matrix
      d<-upper.tri.remove(d,remove.val=0)
   #Find the fraction of transitive triads, removing any NAs 
   t<-vector()
   for(i in 1:gn){
      #Prepare the transitivity test matrices
      dsqt<-(d[i,,]%*%d[i,,])>0
      dt<-d[i,,]>0
      #NA the appropriate portions of the graph, if needed (this to prevent their being counted)
      if(!diag)
         diag(dt)<-NA
      if(mode=="graph")
         dt[upper.tri(dt)]<-NA              
      #Compute the transitivity
      t[i]<-switch(match.arg(measure),
         strong=mean(dsqt==dt,na.rm=TRUE),
         strongcensus=sum(dsqt==dt,na.rm=TRUE),
         weak=mean(dsqt&dt,na.rm=TRUE),
         weakcensus=sum(dsqt&dt,na.rm=TRUE)
      )
   }
   #Return the result
   t
}


#hierarchy - Find the hierarchy score of a graph or graph stack
hierarchy<-function(dat,g=1:stackcount(dat),measure=c("reciprocity","krackhardt")){
   if(match.arg(measure)=="reciprocity")  #Use reciprocity scores
         h<-1-grecip(dat,g)
   else if(match.arg(measure)=="krackhardt"){ #Calculate the Krackhardt reciprocity
      d<-array(dim=c(length(g),dim(dat)[2],dim(dat)[2]))
      if(length(dim(dat))>2)
         d<-dat[g,,,drop=FALSE]
      else
         d[1,,]<-dat
      h<-1-apply(d,1,function(x){r<-reachability(x); grecip(r,measure="edgewise")})
   }
   #Return the result
   h
}


#lubness - Find Krackhardt's Least Upper Boundedness of a graph or graph stack
lubness<-function(dat,g=1:stackcount(dat)){
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
           maxnolub<-maxnolub+(length(vi))*(length(vi)-1)/2 
         }
      }
      #Return 1-violations/max(violations)
      1-nolub/maxnolub
   }
   #Perform the actual calculation
   d<-array(dim=c(length(g),dim(dat)[2],dim(dat)[2]))
   if(length(dim(dat))>2)
      d<-dat[g,,,drop=FALSE]
   else
      d[1,,]<-dat
   lub<-apply(d,1,intlub)
   #Return the result
   lub
}


#mutuality - Find the number of mutual (i.e., reciprocated) edges in a graph
mutuality<-function(dat,g=NULL){
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
triad.census<-function(dat,g=1:stackcount(dat)){
   #First, define the triad class vector
   tc<-c("003","012","102","021D","021U","021C","111D","111U","030T","030C","201","120D","120U","120C","210","300")
   #Organize the data
   if(length(dim(dat))>2){
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
     tcv<-as.double(rep(0,16))
     tcv<-.C("triad_census_R",as.integer(d[i,,]),n,tcv=tcv,PACKAGE="sna")$tcv
     tcm<-rbind(tcm,tcv)
   }
   colnames(tcm)<-tc
   rownames(tcm)<-rnam
   #Return the result
   tcm
}
#triad.census<-function(dat,g=1:stackcount(dat)){
#   #First, define the triad class vector
#   tc<-c("003","012","102","021D","021U","021C","111D","111U","030T","030C","201","120D","120U","120C","210","300")
#   #Define an internal census function
#   intcalc<-function(m,tcv){
#      n<-dim(m)[1]
#      tricent<-rep(0,length(tcv))
#      for(i in 1:n)
#         for(j in i:n)
#            for(k in j:n)
#               if((i!=j)&&(j!=k)&&(i!=k)){
#                  tric<-triad.classify(m,tri=c(i,j,k))
#                  tricent[tcv==tric]<-tricent[tcv==tric]+1
#               }
#      tricent
#   }
#   #Organize the data
#   if(length(dim(dat))>2)
#      d<-dat[g,,,drop=FALSE]
#   else{
#      d<-array(dim=c(1,dim(dat)[1],dim(dat)[2]))
#      d[1,,]<-dat
#   }
#   d<-diag.remove(d,remove.val=0)  #Remove any diagonals
#   #Perform the census
#   census<-t(apply(d,1,intcalc,tc))
#   colnames(census)<-tc
#   #Return the result
#   census
#}


#triad.classify - Return the Davis and Leinhardt classification of a given triad
triad.classify<-function(dat,g=1,tri=c(1,2,3)){
   #Zeroth step: extract the triad
   if(length(dim(dat))==2)
      d<-dat[tri,tri]
   else
      d<-dat[g,tri,tri]
   #First, classify as NA if any entries are missing
   if(any(is.na(d[upper.tri(d)|lower.tri(d)])))
      return(NA)
   #Next, define the triad class vector
   tc<-c("003","012","102","021D","021U","021C","111D","111U","030T","030C","201","120D","120U","120C","210","300")
   #Classify the triad
   tt<-as.integer(0)
   tt<-.C("triad_classify_R",as.integer(d),tt=tt,PACKAGE="sna")$tt
   tc[tt]
}
#triad.classify<-function(dat,g=1,tri=c(1,2,3)){
#   #Zeroth step: extract the triad
#   if(length(dim(dat))==2)
#      d<-dat[tri,tri]
#   else
#      d<-dat[g,tri,tri]
#   #First, classify as NA if any entries are missing
#   if(any(is.na(d[upper.tri(d)|lower.tri(d)])))
#      man<-NA
#   else{
#      man<-dyad.census(d)  #Start with the dyad census
#      #Refine the classification using configural properties
#      if(all(man==c(0,2,1))){   #The two asym/one null triad
#         ind<-apply(d,2,sum)  
#         outd<-apply(d,1,sum)
#         if(any(ind==2))
#            man<-c(man,"U")   #"Up" variant
#         else if(any(outd==2))
#            man<-c(man,"D")   #"Down" variant
#         else
#            man<-c(man,"C")   #"Cyclic" variant
#      }else if(all(man==c(1,1,1))){   #The one mut/one asym/one null triad
#         ind<-apply(d,2,sum)  
#         if(any(ind==2))
#            man<-c(man,"D")   #"Down" variant
#         else
#            man<-c(man,"U")   #"Up" variant
#      }else if(all(man==c(0,3,0))){   #The three asym triad
#         ind<-apply(d,2,sum)  
#         if(any(ind==2))
#            man<-c(man,"T")   #"Transitive" variant
#         else
#            man<-c(man,"C")   #"Cyclic" variant
#      }else if(all(man==c(1,2,0))){   #The one mut/two asym triad
#         ind<-apply(d,2,sum)  
#         outd<-apply(d,1,sum)
#         if(any(ind==0))
#            man<-c(man,"D")   #"Down" variant
#         else if(any(outd==0))
#            man<-c(man,"U")   #"Up" variant
#         else
#            man<-c(man,"C")   #"Cyclic" variant
#      }         
#   }
#   #Return the classification
#   paste(man,collapse="")
#}
