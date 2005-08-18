######################################################################
#
# nli.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 4/23/05
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines for calculating node-level indices.
# 
# Contents:
#   betweenness
#   bonpow
#   closeness
#   degree
#   evcent
#   graphcent
#   infocent
#   stresscent
#
######################################################################


#betweenness - Find the betweenness centralities of network positions
betweenness<-function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,tmaxdev=FALSE,cmode="directed",geodist.precomp=NULL,rescale=FALSE){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],betweenness,g=1,nodes=nodes,gmode=gmode, diag=diag,tmaxdev=tmaxdev,cmode=cmode,geodist.precomp=geodist.precomp, rescale=rescale))
   else if((length(g)>1)&&(length(dim(dat))>2))
     return(apply(dat[g,,],1,betweenness,g=1,nodes=nodes,gmode=gmode, diag=diag,tmaxdev=tmaxdev,cmode=cmode,geodist.precomp=geodist.precomp, rescale=rescale))
   #End pre-processing
   if(gmode=="graph")   #If the data is symmetric, treat it as such
      cmode<-"undirected"
   if(tmaxdev){
      #We got off easy: just return the theoretical maximum deviation for the centralization routine
      bet<-switch(cmode,
         directed = (dim(dat)[2]-1)^2*(dim(dat)[2]-2),
         undirected = (dim(dat)[2]-1)^2*(dim(dat)[2]-2)/2
      )
   }else{
      #First, prepare the data
      if(length(dim(dat))>2)
         d<-dat[g,,]
      else
         d<-dat
      n<-dim(d)[1]
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:n
      if(cmode=="undirected")   #Symmetrize if need be
         for(i in 1:n)
            for(j in 1:n)
               if(i!=j)
                  d[i,j]<-max(d[i,j],d[j,i])
      #Do the computation
      if(is.null(geodist.precomp))
         gd<-geodist(d)
      else
         gd<-geodist.precomp
      bet<-rep(0,n)
      bet<-.C("betweenness_R",as.double(d),as.double(n),bet=as.double(bet), as.double(gd$gdist),as.double(gd$counts),NAOK=TRUE,PACKAGE="sna")$bet
      if(cmode=="undirected")
         bet<-bet/2
      #Return the results
      if(rescale)
         bet<-bet/sum(bet)
      bet<-bet[nodes]
   }
   bet
}


#bonpow - Find the Bonacich power centrality scores of network positions
bonpow<-function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,tmaxdev=FALSE,exponent=1,rescale=FALSE,tol=1e-7){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],bonpow,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,exponent=exponent,rescale=rescale,tol=tol))
   else if((length(g)>1)&&(length(dim(dat))>2))
     return(apply(dat[g,,],1,bonpow,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,exponent=exponent,rescale=rescale,tol=tol))
   #End pre-processing
   if(tmaxdev){
      #We got off easy: just return the theoretical maximum deviation for the centralization routine
      if(gmode=="graph")
         ev<-(dim(dat)[2]-2)*sqrt(dim(dat)[2]/2)
      else
         ev<-sqrt(dim(dat)[2])*(dim(dat)[2]-1)
   }else{
      #First, prepare the data
      if(length(dim(dat))>2)
         d<-dat[g,,]
      else
         d<-dat
      n<-dim(d)[1]
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:n
      if(!diag)
         diag(d)<-0
      #Make an identity matrix
      id<-matrix(rep(0,n*n),nrow=n)
      diag(id)<-1
      #Do the computation
      ev<-apply(solve(id-exponent*d,tol=tol)%*%d,1,sum)  #This works, when it works.
      #Apply the Bonacich scaling, by default (sum of squared ev=n)
      ev<-ev*sqrt(n/sum((ev)^2))
      if(rescale)
         ev<-ev/sum(ev)
      ev[nodes]
   }
   ev
}


#closeness - Find the closeness centralities of network positions
closeness<-function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,tmaxdev=FALSE,cmode="directed",geodist.precomp=NULL,rescale=FALSE){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],closeness,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,cmode=cmode,geodist.precomp=geodist.precomp,rescale=rescale))
   else if((length(g)>1)&&(length(dim(dat))>2))
     return(apply(dat[g,,],1,closeness,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,cmode=cmode,geodist.precomp=geodist.precomp,rescale=rescale))
   #End pre-processing
   if(gmode=="graph")   #If the data is symmetric, treat it as such
      cmode<-"undirected"
   if(tmaxdev){
      #We got off easy: just return the theoretical maximum deviation for the centralization routine
      n<-dim(dat)[2]
      clo<-switch(cmode,
         directed = (n-1)*(1-1/n),    #Depends on n subst for max distance
         undirected = (n-2)*(n-1)/(2*n-3)
      )
   }else{
      #First, prepare the data
      if(length(dim(dat))>2)
         d<-dat[g,,]
      else
         d<-dat
      n<-dim(d)[1]
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:n
      if(cmode=="undirected")   #Symmetrize if need be
         d<-symmetrize(d,"weak")
      #Do the computation
      if(is.null(geodist.precomp))
         gd<-geodist(d)
      else
         gd<-geodist.precomp
      clo<-rep(0,n)
      for(i in 1:n)
         clo[i]<-sum(gd$gdist[i,-i])
      clo<-(n-1)/clo
      if(rescale)
         clo<-clo/sum(clo)
      clo<-clo[nodes]
   }
   #Return the results
   clo
}


#degree - Find the degree centralities of network positions
degree<-function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,tmaxdev=FALSE,cmode="freeman",rescale=FALSE){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],degree,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,cmode=cmode,rescale=rescale))
   else if((length(g)>1)&&(length(dim(dat))>2))
     return(apply(dat[g,,],1,degree,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,cmode=cmode,rescale=rescale))
   #End pre-processing
   if(tmaxdev){
      #We got off easy: just return the theoretical maximum deviation for the centralization routine
      if(gmode=="digraph")
        deg<-switch(cmode,
           indegree = (dim(dat)[2]-1)*(dim(dat)[2]-1+as.numeric(diag)),
           outdegree = (dim(dat)[2]-1)*(dim(dat)[2]-1+as.numeric(diag)),
           freeman = (dim(dat)[2]-1)*(2*(dim(dat)[2]-1)-2+as.numeric(diag))
        )
      else
        deg<-switch(cmode,
           indegree = (dim(dat)[2]-1)*(dim(dat)[2]-2+as.numeric(diag)),
           outdegree = (dim(dat)[2]-1)*(dim(dat)[2]-2+as.numeric(diag)),
           freeman = (dim(dat)[2]-1)*(2*(dim(dat)[2]-1)-2+as.numeric(diag))
        )
   }else{
      #First, prepare the data
      if(length(dim(dat))>2)
         d<-dat[g,,]
      else
         d<-dat
      n<-dim(d)[1]
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:n
      if(!diag)
         diag(d)<-NA
      #Do the computation
      deg<-switch(cmode,
         indegree = apply(d,2,sum,na.rm=TRUE),
         outdegree = apply(d,1,sum,na.rm=TRUE),
         freeman = apply(d,2,sum,na.rm=TRUE) + apply(d,1,sum,na.rm=TRUE)
      )
      if(rescale)
         deg<-deg/sum(deg)
      deg<-deg[nodes]
   }
   deg
}


#evcent - Find the eigenvector centralities of network positions
evcent<-function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,tmaxdev=FALSE,rescale=FALSE){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],evcent,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,rescale=rescale))
   else if((length(g)>1)&&(length(dim(dat))>2))
     return(apply(dat[g,,],1,evcent,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,rescale=rescale))
   #End pre-processing
   if(tmaxdev){
      #We got off easy: just return the theoretical maximum deviation for the centralization routine
      if(gmode=="graph"){
         temp<-matrix(0,dim(dat)[2],dim(dat)[2]) #Construct the max
         temp[1,2]<-1                            #deviation structure
         temp[2,1]<-1
         ev<-eigen(temp)$vectors[,1]
         ev<-sum(max(ev)-ev)
      }else
         ev<-dim(dat)[2]-1
   }else{
      #First, prepare the data
      if(length(dim(dat))>2)
         d<-dat[g,,]
      else
         d<-dat
      n<-dim(d)[1]
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:n
      if(!diag)
         diag(d)<-0
      #Do the computation
      ev<-abs(eigen(d)$vectors[,1])
      if(rescale)
         ev<-ev/sum(ev)
      ev<-ev[nodes]
   }
   ev
}


#graphcent - Find the graph centralities of network positions
graphcent<-function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,tmaxdev=FALSE,cmode="directed",geodist.precomp=NULL,rescale=FALSE){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],graphcent,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,cmode=cmode,geodist.precomp=geodist.precomp,rescale=rescale))
   else if((length(g)>1)&&(length(dim(dat))>2))
     return(apply(dat[g,,],1,graphcent,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,cmode=cmode,geodist.precomp=geodist.precomp,rescale=rescale))
   #End pre-processing
   if(gmode=="graph")   #If the data is symmetric, treat it as such
      cmode<-"undirected"
   if(tmaxdev){
      #We got off easy: just return the theoretical maximum deviation for the centralization routine
      n<-dim(dat)[2]
      gc<-switch(cmode,
         directed = (n-1)*(1-1/n),  #Depends on n subst for infinite distance
         undirected = (n-1)/2
      )
   }else{
      #First, prepare the data
      if(length(dim(dat))>2)
         d<-dat[g,,]
      else
         d<-dat
      n<-dim(d)[1]
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:n
      if(cmode=="undirected")   #Symmetrize if need be
         for(i in 1:n)
            for(j in 1:n)
               if(i!=j)
                  d[i,j]<-max(d[i,j],d[j,i])
      #Do the computation
      if(is.null(geodist.precomp))
         gd<-geodist(d)
      else
         gd<-geodist.precomp
      gc<-rep(0,n)
      for(i in 1:n){
         for(j in 1:n)
            if(j!=i)
               gc[i]<-max(gc[i],gd$gdist[i,j])
      }
      gc<-1/gc
      if(rescale)
         gc<-gc/sum(gc)
      gc<-gc[nodes]
   }
   #Return the results
   gc
}


# infocent - Find actor information centrality scores
# Wasserman & Faust pp. 192-197; based on code generously submitted by David
# Barron (thanks!) and tweaked by myself to enable compatibility with the
# centralization() routine.
infocent <- function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,cmode="weak",tmaxdev=FALSE,rescale=FALSE,tol=1e-20){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],infocent,g=1,nodes=nodes,gmode=gmode,diag=diag, cmode=cmode,tmaxdev=tmaxdev,rescale=rescale,tol=tol))
   else if((length(g)>1)&&(length(dim(dat))>2))
     return(apply(dat[g,,],1,infocent,g=1,nodes=nodes,gmode=gmode,diag=diag, cmode=cmode,tmaxdev=tmaxdev,rescale=rescale,tol=tol))
   #End pre-processing
   if(tmaxdev){  #If necessary, return the theoretical maximum deviation
      #We don't know the real maximum value...return the lone dyad instead
      m<-matrix(0,nr=dim(dat)[2],nc=dim(dat)[2])
      m[1,2]<-1
      m[2,1]<-1
      IC<-infocent(m,1,rescale=rescale)  #Get ICs for dyad
      cent<-sum(max(IC)-IC,na.rm=TRUE)    #Return the theoretical max deviation 
   }else{
      #First, prepare the data
      if(length(dim(dat))>2)
         m<-dat[g,,]
      else
         m<-dat
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:dim(dat)[2]
      if(sum(m != t(m),na.rm=TRUE) > 0)   #test to see if directed
         m <- symmetrize(m,rule=cmode)    #if not, we have to symmetrize...
      n <- dim(m)[1]
      if(!diag) 
         diag(m)<-NA   # if diag=F set diagonal to NA
      iso <- is.isolate(m,1:n,diag=diag) # check for isolates
      ix <- which(!iso)
      m <- m[ix,ix]           # remove any isolates (can't invert A otherwise)
      A<-1-m
      A[m==0] <- 1
      diag(A) <- 1 + apply(m, 1, sum, na.rm=TRUE)
      Cn <- solve(A,tol=tol)
      Tr <- sum(diag(Cn))
      R <- apply(Cn, 1, sum)
      IC <- 1/(diag(Cn) + (Tr - 2*R)/n)   # Actor information centrality
      #Add back the isolates
      cent<-rep(0,n)
      cent[ix]<-IC
      #Rescale if needed
      if(rescale)
         cent<-cent/sum(cent)
      #Subset as requested
      cent<-cent[nodes]
   }
   #Return the result
   cent
}


#prestige - Find actor prestige scores from one of several measures
prestige<-function(dat,g=1,nodes=NULL,gmode="digraph",diag=FALSE,cmode="indegree",tmaxdev=FALSE,rescale=FALSE,tol=1e-7){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],prestige,g=1,nodes=nodes,gmode=gmode,diag=diag, cmode=cmode,tmaxdev=tmaxdev,rescale=rescale,tol=tol))
   else if((length(g)>1)&&(length(dim(dat))>2))
     return(apply(dat[g,,],1,prestige,g=1,nodes=nodes,gmode=gmode,diag=diag, cmode=cmode,tmaxdev=tmaxdev,rescale=rescale,tol=tol))
   #End pre-processing
   if(tmaxdev){
      #We got off easy: just return the theoretical maximum deviation for the centralization routine
      n<-dim(dat)[2]
      if(cmode=="indegree")
         p<-degree(dat=matrix(nrow=n,ncol=n),g=1,tmaxdev=TRUE,gmode=gmode,diag=diag,cmode="indegree",rescale=FALSE)
      else if(cmode=="indegree.rownorm")
         p<-degree(dat=matrix(nrow=n,ncol=n),g=1,tmaxdev=TRUE,gmode=gmode,diag=diag,cmode="indegree",rescale=FALSE)
      else if(cmode=="indegree.rowcolnorm")
         p<-degree(dat=matrix(nrow=n,ncol=n),g=1,tmaxdev=TRUE,gmode=gmode,diag=diag,cmode="indegree",rescale=FALSE)
      else if(cmode=="eigenvector")
         p<-evcent(dat=matrix(nrow=n,ncol=n),g=1,tmaxdev=TRUE,gmode=gmode,diag=diag)
      else if(cmode=="eigenvector.rownorm")
         p<-evcent(dat=matrix(nrow=n,ncol=n),g=1,tmaxdev=TRUE,gmode=gmode,diag=diag)
      else if(cmode=="eigenvector.colnorm")
         p<-evcent(dat=matrix(nrow=n,ncol=n),g=1,tmaxdev=TRUE,gmode=gmode,diag=diag)
      else if(cmode=="eigenvector.rowcolnorm")
         p<-evcent(dat=matrix(nrow=n,ncol=n),g=1,tmaxdev=TRUE,gmode=gmode,diag=diag)
      else if(cmode=="domain"){
         p<-(n-1)^2
      }else if(cmode=="domain.proximity"){
         p<-(n-1)^2
      }else
         stop(paste("Cmode",cmode,"unknown.\n"))      
   }else{
      #First, prepare the data
      if(length(dim(dat))>2)
         d<-dat[g,,]
      else
         d<-dat
      n<-dim(d)[1]
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:n
      if(!diag)
         diag(d)<-0
      #Now, perform the computation
      if(cmode=="indegree")
         p<-degree(dat=dat,g=g,nodes=nodes,gmode=gmode,diag=diag,cmode="indegree",rescale=FALSE)
      else if(cmode=="indegree.rownorm")
         p<-degree(dat=make.stochastic(d,mode="row"),g=1,nodes=nodes,gmode=gmode,diag=diag,cmode="indegree",rescale=FALSE)
      else if(cmode=="indegree.rowcolnorm")
         p<-degree(dat=make.stochastic(d,mode="rowcol"),g=1,nodes=nodes,gmode=gmode,diag=diag,cmode="indegree",rescale=FALSE)
      else if(cmode=="eigenvector")
         p<-eigen(t(d))$vector[,1]      
      else if(cmode=="eigenvector.rownorm")
         p<-eigen(t(make.stochastic(d,mode="row")))$vector[,1]      
      else if(cmode=="eigenvector.colnorm")
         p<-eigen(t(make.stochastic(d,mode="col")))$vector[,1]      
      else if(cmode=="eigenvector.rowcolnorm")
         p<-eigen(t(make.stochastic(d,mode="rowcol")))$vector[,1]
      else if(cmode=="domain"){
         r<-reachability(d)
         p<-apply(r,2,sum)-1
      }else if(cmode=="domain.proximity"){
         g<-geodist(d)
         p<-(apply(g$counts>0,2,sum)-1)^2/(apply((g$counts>0)*(g$gdist),2,sum)*(n-1))
         p[is.nan(p)]<-0
      }else
         stop(paste("Cmode",cmode,"unknown.\n"))      
      if(rescale)
         p<-p/sum(p)
      p<-p[nodes]
   }  
   p
}


#stresscent - Find the stress centralities of network positions
stresscent<-function(dat,g=1,nodes=c(1:dim(dat)[2]),gmode="digraph",diag=FALSE,tmaxdev=FALSE,cmode="directed",geodist.precomp=NULL,rescale=FALSE){
   #Pre-process the raw input
   dat<-as.sociomatrix.sna(dat)
   if(is.list(dat))
     return(sapply(dat[g],stresscent,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,cmode=cmode,geodist.precomp=geodist.precomp,rescale=rescale))
   else if((length(g)>1)&&(length(dim(dat))>2))
     return(apply(dat[g,,],1,stresscent,g=1,nodes=nodes,gmode=gmode,diag=diag, tmaxdev=tmaxdev,cmode=cmode,geodist.precomp=geodist.precomp,rescale=rescale))
   #End pre-processing
   if(gmode=="graph")   #If the data is symmetric, treat it as such
      cmode<-"undirected"
   if(tmaxdev){
      #We got off easy: just return the theoretical maximum deviation for the centralization routine
      str<-switch(cmode,
         directed = (dim(dat)[2]-1)^2*(dim(dat)[2]-2),
         undirected = (dim(dat)[2]-1)^2*(dim(dat)[2]-2)/2
      )
   }else{
      #First, prepare the data
      if(length(dim(dat))>2)
         d<-dat[g,,]
      else
         d<-dat
      n<-dim(d)[1]
      if(is.null(nodes))        #Set up node list, if needed
        nodes<-1:n
      if(cmode=="undirected")   #Symmetrize if need be
         for(i in 1:n)
            for(j in 1:n)
               if(i!=j)
                  d[i,j]<-max(d[i,j],d[j,i])
      #Do the computation
      if(is.null(geodist.precomp))
         gd<-geodist(d)
      else
         gd<-geodist.precomp
      str<-rep(0,n)
      str<-.C("stresscent_R",as.double(d),as.double(n),str=as.double(str), as.double(gd$gdist),as.double(gd$counts),NAOK=TRUE,PACKAGE="sna")$str
      if(cmode=="undirected")
         str<-str/2
      #Return the results
      if(rescale)
         str<-str/sum(str)
      str<-str[nodes]
   }
   str
}

