######################################################################
#
# connectivity.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 8/14/20
# Licensed under the GNU General Public License version 2 (June, 1991)
# or later.
#
# Part of the R/sna package
#
# This file contains various routines associated with connectivity
# properties (including geodesic distance and friends).
#
# Contents:
#  bicomponent.dist
#  clique.census
#  component.dist
#  component.largest
#  component.size.byvertex
#  components
#  cutpoints
#  geodist
#  isolates
#  is.connected
#  is.isolate
#  kcores
#  kcycle.census
#  kpath.census
#  maxflow
#  neighborhood
#  reachability
#  simmelian
#  structure.statistics
#
######################################################################


#bicomponent.dist - Returns a list containing a vector of length n such that
#the ith element contains the number of components of G having size i, and a 
#vector of length n giving component membership.  Component strength is 
#determined by the rule which is used to symmetrize the matrix; this controlled 
#by the eponymous parameter given to the symmetrize command.
bicomponent.dist<-function(dat,symmetrize=c("strong","weak")){
   #Pre-process the raw input
   dat<-as.edgelist.sna(dat,suppress.diag=TRUE)
   if(is.list(dat))
     return(lapply(dat,bicomponent.dist,symmetrize=symmetrize))
   #End pre-processing
   #Begin routine
   n<-attr(dat,"n")
   #Symmetrize dat based on the connectedness rule
   dat<-symmetrize(dat,rule=match.arg(symmetrize),return.as.edgelist=TRUE)
   #Compute the bicomponents
   bc<-.Call("bicomponents_R",dat,n,NROW(dat),PACKAGE="sna")
   if(length(bc[[1]])>1){                             #Sort by size
     ord<-order(sapply(bc[[1]],length),decreasing=TRUE)
     bc[[1]]<-bc[[1]][ord]
     bc[[2]][bc[[2]]>0]<-match(bc[[2]][bc[[2]]>0],ord)
   }
   bc[[2]][bc[[2]]<0]<-NA
   bc[[1]]<-bc[[1]][sapply(bc[[1]],length)>0]
   #Return the results
   o<-list()
   if(length(bc[[1]])>0){
     o$members<-bc[[1]]                    #Copy membership lists
     names(o$members)<-1:length(o$members)
     o$membership<-bc[[2]]                 #Copy memberships
     o$csize<-sapply(o$members,length)     #Extract component sizes
     names(o$csize)<-1:length(o$csize)
     o$cdist<-tabulate(o$csize,nbins=n)    #Find component size distribution
     names(o$cdist)<-1:n
   }else{
     o$members<-list()
     o$membership<-bc[[2]]
     o$csize<-vector(mode="numeric")
     o$cdist<-rep(0,n)
     names(o$cdist)<-1:n
   }
   o
}


#clique.census - Enumerate all maximal cliques
clique.census<-function(dat,mode="digraph",tabulate.by.vertex=TRUE,clique.comembership=c("none","sum","bysize"),enumerate=TRUE, na.omit=TRUE){
  #Pre-process the raw input
  dat<-as.edgelist.sna(dat)
  if(is.list(dat))
    return(lapply(dat,clique.census,mode=mode, tabulate.by.vertex=tabulate.by.vertex,clique.comembership=clique.comembership,enumerate=enumerate, na.omit=na.omit))
  #End pre-processing
  n<-attr(dat,"n")
  if(is.null(attr(dat,"vnames")))
    vnam<-paste("v",1:n,sep="")
  else
    vnam<-attr(dat,"vnames")
  if(na.omit)
    dat<-dat[!is.na(dat[,3]),,drop=FALSE]  #Drop any edges with NAs
  else
    dat[is.na(dat[,3]),3]<-1               #Else, recode to safe values
  dat<-dat[dat[,1]!=dat[,2],]              #Remove loops
  attr(dat,"n")<-n
  #If called with a digraph, symmetrize
  if(mode=="digraph")
    dat<-symmetrize(dat,rule="strong",return.as.edgelist=TRUE)
  #Compute the census
  clique.comembership<-switch(match.arg(clique.comembership),
    none=0,
    sum=1,
    bysize=2
  )
  census<-.Call("cliques_R",dat,n,NROW(dat),tabulate.by.vertex, clique.comembership,enumerate,PACKAGE="sna")
  #Assemble the results
  maxsize<-census[[1]]
  census<-census[-1]
  names(census)<-c("clique.count","clique.comemb","cliques")
  if(tabulate.by.vertex){
    census[[1]]<-matrix(census[[1]],maxsize,n+1)
    census[[1]]<-census[[1]][,c(n+1,1:n),drop=FALSE]
    rownames(census[[1]])<-1:maxsize
    colnames(census[[1]])<-c("Agg",vnam)
  }else{
    names(census[[1]])<-1:length(census[[1]])
  }
  if(clique.comembership==1){
    census[[2]]<-matrix(census[[2]],n,n)
    rownames(census[[2]])<-vnam
    colnames(census[[2]])<-vnam
  }else if(clique.comembership==2){
    census[[2]]<-array(census[[2]],dim=c(maxsize,n,n))
    dimnames(census[[2]])<-list(1:maxsize,vnam,vnam)
  }
  #Return the non-null components
  pres<-c(TRUE,clique.comembership>0,enumerate>0)
  census[pres]       
}


#component.dist - Returns a data frame containing a vector of length n such that
#the ith element contains the number of components of G having size i, and a 
#vector of length n giving component membership.  Component strength is 
#determined by the rule which is used to symmetrize the matrix; this controlled 
#by the eponymous parameter given to the symmetrize command.
component.dist<-function(dat,connected=c("strong","weak","unilateral","recursive")){
   #Pre-process the raw input
   if(match.arg(connected)%in%c("strong","weak","recursive"))
     dat<-as.edgelist.sna(dat)
   else
     dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(lapply(dat,component.dist,connected=connected))
   else if(length(dim(dat))>2)
     return(apply(dat,1,component.dist,connected=connected))
   #End pre-processing
   #Begin routine
   #Proceed depending on the rule being used
   if(match.arg(connected)%in%c("strong","weak","recursive")){ #Strong, weak, recursive
     n<-attr(dat,"n")
     #Preprocess as needed
     dat<-switch(match.arg(connected),
       "weak"=symmetrize(dat,rule="weak",return.as.edgelist=TRUE),
       "strong"=symmetrize(reachability(dat,return.as.edgelist=TRUE),rule="strong", return.as.edgelist=TRUE),
       "recursive"=symmetrize(dat,rule="strong",return.as.edgelist=TRUE)
     )
     #Find the component information using the leanest available method
     memb<-.C("undirComponents_R",as.double(dat),as.integer(n),as.integer(NROW(dat)), memb=integer(n+1),PACKAGE="sna",NAOK=TRUE)$memb
     csize<-tabulate(memb[-1],memb[1])
     cdist<-rep(0,n)
     cdist[1:max(csize)]<-tabulate(csize,max(csize))
     memb<-memb[-1]
   }else{                                          #Unilateral
     n<-dim(dat)[2]
     dat<-reachability(dat)
     #Warn of non-uniqueness in the unilateral case, if need be
     if(any(dat!=t(dat)))
       warning("Nonunique unilateral component partition detected in component.dist.  Problem vertices will be arbitrarily assigned to one of their components.\n")
     #Find the membership information using a not-too-shabby method
     memb<-.C("component_dist_R",as.double(dat),as.double(n), memb=as.double(rep(0,n)),PACKAGE="sna",NAOK=TRUE)$memb
     csize<-tabulate(memb,max(memb))
     cdist<-rep(0,n)
     cdist[1:max(csize)]<-tabulate(csize,max(csize))
   }
   #Return the results
   o<-list(membership=memb,csize=csize,cdist=cdist)
   o
}


#component.largest - Extract the largest component from a graph
component.largest<-function(dat,connected=c("strong","weak","unilateral", "recursive"), result=c("membership","graph"),return.as.edgelist=FALSE){
    #Deal with network, array, or list data
    dat <- as.edgelist.sna(dat)
    if (is.list(dat))
        return(lapply(dat, component.largest, connected = connected, result = result))
    #We now have a single graph.  Proceed accordingly.
    if(attr(dat,"n")==1){
      if(match.arg(result)=="membership"){
        return(TRUE)
      }else{
        if(return.as.edgelist)
          return(dat)
        else
          return(as.sociomatrix.sna(dat))
      }
    }
    cd<-component.dist(dat,connected=connected)
    lgcmp<-which(cd$csize==max(cd$csize))  #Get largest component(s)
    #Return the appropriate result
    if(match.arg(result)=="membership"){
      cd$membership%in%lgcmp
    }else{
      tokeep<-which(cd$membership%in%lgcmp)
      ovn<-attr(dat,"vnames")
      if(is.null(ovn))
        ovn<-1:attr(dat,"n")
      if(return.as.edgelist){
        sel<-rowSums(apply(dat,1:2,function(z){z%in%tokeep}))==2
        dat<-dat[sel,,drop=FALSE]
        if(NROW(dat)>0){
          dat[,1:2]<-apply(dat,1:2,function(z){match(z,tokeep)})
        }
        attr(dat,"n")<-length(tokeep)
        attr(dat,"vnames")<-ovn[tokeep]
        dat
      }else{
        as.sociomatrix.sna(dat)[tokeep,tokeep,drop=FALSE]
      }
    }
}


#component.size.byvertex
component.size.byvertex<-function(dat, connected=c("strong","weak","unilateral","recursive")){
  #Pre-process the input
  g<-as.edgelist.sna(dat)
  if(is.list(g)){
    return(lapply(g,component.size.byvertex,connected=connected))
  }
  #End pre-processing
  if(match.arg(connected)%in%c("weak","recursive")){ #We have a shortcut for these cases! 
    if(match.arg(connected)=="weak")
      rule<-"weak"
    else
      rule<-"strong"
    g<-symmetrize(g,rule=rule, return.as.edgelist=TRUE) #Must symmetrize!
    cs<-.C("compsizes_R",as.double(g),as.integer(attr(g,"n")),as.integer(NROW(g)), csizes=integer(attr(g,"n")),PACKAGE="sna",NAOK=TRUE)$csizes
  }else{                                            #No shortcut.  Sad!
    cd<-component.dist(dat,connected=match.arg(connected))
    cs<-cd$csize[cd$membership]
  }
  #Return the results
  cs
}



#components - Find the number of (maximal) components within a given graph
components<-function(dat,connected="strong",comp.dist.precomp=NULL){
   #Pre-process the raw input
   dat<-as.edgelist.sna(dat)
   if(is.list(dat))
     return(lapply(dat,components,connected=connected, comp.dist.precomp=comp.dist.precomp))
   #End pre-processing
   #Use component.dist to get the distribution
   if(!is.null(comp.dist.precomp))
      cd<-comp.dist.precomp
   else
      cd<-component.dist(dat,connected=connected)
   #Return the result
   length(unique(cd$membership))
}


#cutpoints - Find the cutpoints of an input graph
cutpoints<-function(dat,mode="digraph",connected=c("strong","weak","recursive"),return.indicator=FALSE){
   #Pre-process the raw input
   dat<-as.edgelist.sna(dat)
   if(is.list(dat))
     return(lapply(dat,cutpoints,mode=mode,connected=connected, return.indicator=return.indicator))
   #End pre-processing
   n<-attr(dat,"n")
   dat<-dat[dat[,1]!=dat[,2],]   #Remove any loops, lest they break things
   attr(dat,"n")<-n
   cp<-rep(0,n)
   if(mode=="graph")
     cp<-.C("cutpointsUndir_R",as.double(dat),as.integer(n), as.integer(NROW(dat)),cp=as.integer(cp),NAOK=TRUE,PACKAGE="sna")$cp
   else{
     dat<-switch(match.arg(connected),
       strong=dat,
       weak=symmetrize(dat,rule="weak",return.as.edgelist=TRUE),
       recursive=symmetrize(dat,rule="strong",return.as.edgelist=TRUE)
     )
     if(match.arg(connected)=="strong")
       cp<-.C("cutpointsDir_R",as.double(dat),as.integer(n), as.integer(NROW(dat)),cp=as.integer(cp),NAOK=TRUE,PACKAGE="sna")$cp
     else
       cp<-.C("cutpointsUndir_R",as.double(dat),as.integer(n), as.integer(NROW(dat)),cp=as.integer(cp),NAOK=TRUE,PACKAGE="sna")$cp
   }
   if(!return.indicator)
     return(which(cp>0))
   else{
     if(is.null(attr(dat,"vnames")))
       names(cp)<-1:n
     else
       names(cp)<-attr(dat,"vnames")
     return(cp>0)
   }   
}


#geodist - Find the numbers and lengths of geodesics among nodes in a graph 
#using a BFS, a la Brandes (2008).  Note that we still need N^2 storage,
#although calculations are done on the edgelist (which should save some time).
#Both valued and unvalued variants are possible -- don't use the valued 
#version unless you need to, since it can be considerably slower.
geodist<-function(dat,inf.replace=Inf,count.paths=TRUE,predecessors=FALSE,ignore.eval=TRUE, na.omit=TRUE){
   #Pre-process the raw input
   dat<-as.edgelist.sna(dat)
   if(is.list(dat))
     return(lapply(dat,geodist,inf.replace=inf.replace,ignore.eval=ignore.eval))
   #End pre-processing
   n<-attr(dat,"n")
   if(na.omit)
     sel<-!is.na(dat[,3])
   else
     sel<-rep(TRUE,NROW(dat))
   dat<-dat[(dat[,1]!=dat[,2])&sel,,drop=FALSE]
   m<-NROW(dat)
   #Initialize the matrices
   #Perform the calculation
   if(ignore.eval)
     geo<-.Call("geodist_R",dat,n,m,as.integer(1),count.paths,predecessors, PACKAGE="sna")
   else{
     if(any(dat[!is.na(dat[,3]),3]<0))
       stop("Negative edge values not currently supported in geodist; transform or otherwise alter them to ensure that they are nonnegative.")
     geo<-.Call("geodist_val_R",dat,n,m,as.integer(1),count.paths,predecessors, PACKAGE="sna")
   }
   #Return the results
   o<-list()
   if(count.paths)
     o$counts<-matrix(geo[[2]],n,n)
   o$gdist<-matrix(geo[[1]],n,n)
   o$gdist[o$gdist==Inf]<-inf.replace  #Patch Infs, if desired
   if(predecessors)
     o$predecessors<-geo[[2+count.paths]]
   o
}


#isolates - Returns a list of the isolates in a given graph or stack
isolates<-function(dat,diag=FALSE){
   #Pre-process the raw input
   dat<-as.edgelist.sna(dat)
   if(is.list(dat))
     return(lapply(dat,isolates,diag))
   #End pre-processing
   n<-attr(dat,"n")
   if(!diag){
     dat<-dat[dat[,1]!=dat[,2],,drop=FALSE]
   }
   which(tabulate(as.vector(dat[,1:2]),n)==0)
}


#is.connected - Determine whether or not one or more graphs are connected
is.connected<-function(g,connected="strong",comp.dist.precomp=NULL){
  #Pre-process the raw input
  g<-as.edgelist.sna(g)
  if(is.list(g))
    return(lapply(g,is.connected,connected=connected, comp.dist.precomp=comp.dist.precomp))
  #End pre-processing
  #Calculate numbers of components
  components(g,connected=connected,comp.dist.precomp=comp.dist.precomp)==1
}


#is.isolate - Returns TRUE iff ego is an isolate
is.isolate<-function(dat,ego,g=1,diag=FALSE){
   #Pre-process the raw input
   dat<-as.edgelist.sna(dat)
   if(is.list(dat))
     return(lapply(dat[g],is.isolate,ego=ego,g=1,diag=diag))
   #End pre-processing
   if(!diag)
      dat<-dat[dat[,1]!=dat[,2],,drop=FALSE]
   dat<-dat[!is.na(dat[,3]),,drop=FALSE]
   noniso<-unique(c(dat[,1],dat[,2]))
   !(ego%in%noniso)
}


#kcores - Perform k-core decomposition of one or more input graphs
kcores<-function(dat,mode="digraph",diag=FALSE,cmode="freeman",ignore.eval=FALSE){
  #Pre-process the raw input
  dat<-as.edgelist.sna(dat,as.digraph=TRUE,suppress.diag=TRUE)
  if(is.list(dat))
    return(lapply(dat,kcores,dat=dat,mode=mode,diag=diag,cmode=cmode, ignore.eval=ignore.eval))
  #End pre-processing
  if(mode=="graph")             #If undirected, force to "indegree"
    cmode<-"indegree"
  n<-attr(dat,"n")
  m<-NROW(dat)
  corevec<-1:n
  dtype<-switch(cmode,
    indegree=0,
    outdegree=1,
    freeman=2
  )
  if(!(cmode%in%c("indegree","outdegree","freeman")))
    stop("Illegal cmode in kcores.\n")
  solve<-.C("kcores_R",as.double(dat),as.integer(n),as.integer(m), cv=as.double(corevec), as.integer(dtype), as.integer(diag), as.integer(ignore.eval), NAOK=TRUE,PACKAGE="sna")
  if(is.null(attr(dat,"vnames")))
    names(solve$cv)<-1:n
  else
    names(solve$cv)<-attr(dat,"vnames")
  solve$cv
}


#kcycle.census - Compute the cycle census of a graph, possibly along with 
#additional information on the inidence of cycles.
kcycle.census<-function(dat,maxlen=3,mode="digraph",tabulate.by.vertex=TRUE,cycle.comembership=c("none","sum","bylength")){
  #Pre-process the raw input
  dat<-as.edgelist.sna(dat)
  if(is.list(dat))
    return(lapply(dat,kcycle.census,maxlen=maxlen,mode=mode, tabulate.by.vertex=tabulate.by.vertex,cycle.comembership=cycle.comembership))
  #End pre-processing
  n<-attr(dat,"n")
  if(is.null(maxlen))
    maxlen<-n
  if(maxlen<2)
    stop("maxlen must be >=2")
  if(is.null(attr(dat,"vnames")))
    vnam<-paste("v",1:n,sep="")
  else
    vnam<-attr(dat,"vnames")
  if(mode=="digraph")
    directed<-TRUE
  else
    directed<-FALSE
  cocycles<-switch(match.arg(cycle.comembership),
    "none"=0,
    "sum"=1,
    "bylength"=2
  )
  #Generate the data structures for the counts
  if(!tabulate.by.vertex)
    count<-rep(0,maxlen-1)
  else
    count<-matrix(0,maxlen-1,n+1)
  if(!cocycles)
    cccount<-NULL
  else if(cocycles==1)
    cccount<-matrix(0,n,n)
  else
    cccount<-array(0,dim=c(maxlen-1,n,n))
  if(is.null(maxlen))
    maxlen<-n
  #Calculate the cycle information
  ccen<-.C("cycleCensus_R",as.integer(dat), as.integer(n), as.integer(NROW(dat)), count=as.double(count), cccount=as.double(cccount), as.integer(maxlen), as.integer(directed), as.integer(tabulate.by.vertex), as.integer(cocycles),PACKAGE="sna")
  #Coerce the cycle counts into the right form
  if(!tabulate.by.vertex){
    count<-ccen$count
    names(count)<-2:maxlen
  }else{
    count<-matrix(ccen$count,maxlen-1,n+1)
    rownames(count)<-2:maxlen
    colnames(count)<-c("Agg",vnam)
  }  
  if(cocycles==1){
    cccount<-matrix(ccen$cccount,n,n)
    rownames(cccount)<-vnam
    colnames(cccount)<-vnam
  }else if(cocycles==2){
    cccount<-array(ccen$cccount,dim=c(maxlen-1,n,n))
    dimnames(cccount)<-list(2:maxlen,vnam,vnam)
  }
  #Return the result
  out<-list(cycle.count=count)
  if(cocycles>0)
    out$cycle.comemb<-cccount
  out
}


#kpath.census - Compute the path census of a graph, possibly along with 
#additional information on the inidence of paths.
kpath.census<-function(dat,maxlen=3,mode="digraph",tabulate.by.vertex=TRUE,path.comembership=c("none","sum","bylength"),dyadic.tabulation=c("none","sum","bylength")){
  #Pre-process the raw input
  dat<-as.edgelist.sna(dat)
  if(is.list(dat))
    return(lapply(dat,kpath.census,maxlen=maxlen,mode=mode, tabulate.by.vertex=tabulate.by.vertex,path.comembership=path.comembership, dyadic.tabulation=dyadic.tabulation))
  #End pre-processing
  n<-attr(dat,"n")
  if(is.null(maxlen))
    maxlen<-n-1
  if(maxlen<1)
    stop("maxlen must be >=1")
  if(is.null(attr(dat,"vnames")))
    vnam<-paste("v",1:n,sep="")
  else
    vnam<-attr(dat,"vnames")
  if(mode=="digraph")
    directed<-TRUE
  else
    directed<-FALSE
  copaths<-switch(match.arg(path.comembership),
    "none"=0,
    "sum"=1,
    "bylength"=2
  )
  dyadpaths<-switch(match.arg(dyadic.tabulation),
    "none"=0,
    "sum"=1,
    "bylength"=2
  )
  #Generate the data structures for the counts
  if(!tabulate.by.vertex)
    count<-rep(0,maxlen)
  else
    count<-matrix(0,maxlen,n+1)
  if(!copaths)
    cpcount<-NULL
  else if(copaths==1)
    cpcount<-matrix(0,n,n)
  else
    cpcount<-array(0,dim=c(maxlen,n,n))
  if(!dyadpaths)
    dpcount<-NULL
  else if(dyadpaths==1)
    dpcount<-matrix(0,n,n)
  else
    dpcount<-array(0,dim=c(maxlen,n,n))
  #Calculate the path information
  pcen<-.C("pathCensus_R",as.double(dat), as.integer(n), as.integer(NROW(dat)), count=as.double(count), cpcount=as.double(cpcount), dpcount=as.double(dpcount), as.integer(maxlen), as.integer(directed), as.integer(tabulate.by.vertex), as.integer(copaths), as.integer(dyadpaths),PACKAGE="sna")
  #Coerce the path counts into the right form
  if(!tabulate.by.vertex){
    count<-pcen$count
    names(count)<-1:maxlen
  }else{
    count<-matrix(pcen$count,maxlen,n+1)
    rownames(count)<-1:maxlen
    colnames(count)<-c("Agg",vnam)
  }  
  if(copaths==1){
    cpcount<-matrix(pcen$cpcount,n,n)
    rownames(cpcount)<-vnam
    colnames(cpcount)<-vnam
  }else if(copaths==2){
    cpcount<-array(pcen$cpcount,dim=c(maxlen,n,n))
    dimnames(cpcount)<-list(1:maxlen,vnam,vnam)
  }
  if(dyadpaths==1){
    dpcount<-matrix(pcen$dpcount,n,n)
    rownames(dpcount)<-vnam
    colnames(dpcount)<-vnam
  }else if(dyadpaths==2){
    dpcount<-array(pcen$dpcount,dim=c(maxlen,n,n))
    dimnames(dpcount)<-list(1:maxlen,vnam,vnam)
  }
  #Return the result
  out<-list(path.count=count)
  if(copaths>0)
    out$path.comemb<-cpcount
  if(dyadpaths>0)
    out$paths.bydyad<-dpcount
  out
}


#maxflow - Return the matrix of maximum flows between positions
maxflow<-function(dat,src=NULL,sink=NULL,ignore.eval=FALSE){
  #Pre-process the raw input
  dat<-as.sociomatrix.sna(dat)
  if(is.list(dat))
    return(lapply(dat,maxflow,src=src,sink=sink,ignore.eval=ignore.eval))
  else if(length(dim(dat))>2)
    return(apply(dat,1,maxflow,src=src,sink=sink,ignore.eval=ignore.eval))
  #End pre-processing
  n<-NROW(dat)
  dat[is.na(dat)]<-0                       #Deal with values and missingness
  if(ignore.eval)
    dat[dat!=0]<-1
  if(length(src)==0)                       #Define sources and sinks
    src<-1:n
  else
    src<-src[(src>0)&(src<=n)]
  if(length(sink)==0)
    sink<-1:n
  else
    sink<-sink[(sink>0)&(sink<=n)]
  fmat<-matrix(nrow=length(src),ncol=length(sink))
  for(i in 1:length(src))
    for(j in 1:length(sink))
      fmat[i,j]<-.C("maxflow_EK_R",as.double(dat),as.integer(NROW(dat)), as.integer(src[i]-1),as.integer(sink[j]-1),flow=as.double(0),NAOK=TRUE,PACKAGE="sna")$flo
  #Return the result
  if(length(src)*length(sink)>1){
    if(is.null(rownames(dat)))
      rownames(fmat)<-src
    else
      rownames(fmat)<-rownames(dat)[src]
    if(is.null(colnames(dat)))
      colnames(fmat)<-sink
    else
      colnames(fmat)<-colnames(dat)[sink]
  }else
    fmat<-as.numeric(fmat)
  fmat
}


#neighborhood - Return the matrix of n-th order neighbors for an input graph
neighborhood<-function(dat,order,neighborhood.type=c("in","out","total"),mode="digraph",diag=FALSE,thresh=0,return.all=FALSE,partial=TRUE){
  #Pre-process the raw input
  dat<-as.sociomatrix.sna(dat)
  if(is.list(dat))
    return(lapply(dat,neighborhood,order=order, neighborhood.type=neighborhood.type,mode=mode,diag=diag,thresh=thresh,return.all=return.all,partial=partial))
  else if(length(dim(dat))>2)
    return(apply(dat,1,neighborhood,order=order, neighborhood.type=neighborhood.type,mode=mode,diag=diag,thresh=thresh,return.all=return.all,partial=partial))
  #End pre-processing
  dat<-dat>thresh           #Dichotomize at threshold
  #Adjust the graph to take care of symmetry or neighborhood type issues
  if((mode=="graph")||(match.arg(neighborhood.type)=="total"))
    dat<-dat|t(dat)
  if(match.arg(neighborhood.type)=="in")
    dat<-t(dat)
  #Extract the neighborhood graphs
  geo<-geodist(dat)
  if(return.all){                     #Return all orders?
    neigh<-array(dim=c(order,NROW(dat),NROW(dat)))
    for(i in 1:order){
      neigh[i,,]<-switch(partial+1,
        geo$gdist<=i,                       #!partial -> order i or less
        geo$gdist==i                        #partial -> exactly order i
      )
      if(!diag)
        diag(neigh[i,,])<-0
    }
  }else{                              #Don't return all orders
    neigh<-switch(partial+1,
      geo$gdist<=order,
      geo$gdist==order
    )
    if(!diag)
      diag(neigh)<-0
  }
  #Return the result
  neigh
}


#reachability - Find the reachability matrix of a graph.
reachability<-function(dat,geodist.precomp=NULL,return.as.edgelist=FALSE,na.omit=TRUE){
   #Pre-process the raw input
   if(!is.null(geodist.precomp)){ #Might as well use a matrix, and not repeat the BFS!
     dat<-as.sociomatrix.sna(dat)
     if(is.list(dat))
       return(lapply(dat,reachability,geodist.precomp=geodist.precomp, return.as.edgelist=return.as.edgelist,na.omit=na.omit))
     else if(length(dim(dat))>2)
       return(unlist(apply(dat,1,function(x,geodist.precomp,return.as.edgelist,na.omit){list(reachability(x, geodist.precomp=geodist.precomp, return.as.edgelist=return.as.edgelist, na.omit=na.omit))}, geodist.precomp=geodist.precomp, return.as.edgelist=return.as.edgelist, na.omit=na.omit),recursive=FALSE))
   }else{                         #Starting from scratch - use the sparse version
     dat<-as.edgelist.sna(dat)
     if(is.list(dat))
       return(lapply(dat,reachability,geodist.precomp=geodist.precomp, return.as.edgelist=return.as.edgelist,na.omit=na.omit))
   }
   #End pre-processing
   if(!is.null(geodist.precomp)){
     #Get the counts matrix
     cnt<-geodist.precomp$counts
     #Dichotomize and return
     if(!return.as.edgelist)
       apply(cnt>0,c(1,2),as.numeric)
     else
       as.edgelist.sna(apply(cnt>0,c(1,2),as.numeric))
   }else{
     n<-attr(dat,"n")
     if(na.omit)
       sel<-!is.na(dat[,3])
     else
       sel<-rep(TRUE,NROW(dat))
     dat<-dat[(dat[,1]!=dat[,2])&sel,,drop=FALSE]
     m<-NROW(dat)
     rg<-.Call("reachability_R",dat,n,m,PACKAGE="sna")
     if(return.as.edgelist)
       rg
     else
       as.sociomatrix.sna(rg)
   }
}


#Function to compute the Simmelian ties for one or more input networks.
#  Arguments:
#    dat - one or more input networks, in any form recognized by as.edgelist.sna
#    dichotomize - logical; should we report whether each dyad has a Simmelian tie?
#        Otherwise, the count of three-clique co-memberships for each dyad is used as
#        an edge value.
#    return.as.edgelist - logical; return the result as an sna edgelist?
#
#  Return value:
#    Either an adjacency matrix or sna edgelist containing the Simmelian tie structure;
#    if multiple networks are supplied, a list of results is returned.
#
simmelian<-function(dat, dichotomize=TRUE, return.as.edgelist=FALSE){
  #Regularize the inputs
  dat <- as.edgelist.sna(dat)
  if(is.list(dat))
    return(lapply(dat, simmelian, dichotomize = dichotomize))
  #Symmetrize the network
  g <- symmetrize(dat, rule="strong", return.as.edgelist=TRUE)
  #Compute the 3-cycle co-memberships
  ccen <- kcycle.census(g, mode="graph", tabulate.by.vertex=FALSE, cycle.comembership="bylength")
  if(dichotomize)
    comemb <- ccen$cycle.comemb[2,,]>0
  else
    comemb <- ccen$cycle.comemb[2,,]
  diag(comemb) <- 0
  #Return the result
  if(return.as.edgelist)
    as.edgelist.sna(comemb)
  else
    comemb
}


#structure.statistics - Return the structure statistics for a given graph
structure.statistics<-function(dat,geodist.precomp=NULL){
  #Pre-process the raw input
  dat<-as.sociomatrix.sna(dat)
  if(is.list(dat))
    return(lapply(dat,structure.statistics,geodist.precomp=geodist.precomp))
  else if(length(dim(dat))>2)
    return(apply(dat,1,structure.statistics,geodist.precomp=geodist.precomp))
  #End pre-processing
  #Get the geodesic distance matrix
  if(is.null(geodist.precomp))
    gd<-geodist(dat)$gdist
  else
    gd<-geodist.precomp$gdist
  #Compute the reachability proportions for each vertex
  ss<-vector()
  for(i in 1:NROW(dat))
    ss[i]<-mean(apply(gd<=i-1,1,mean))
  names(ss)<-0:(NROW(dat)-1)
  ss
}

