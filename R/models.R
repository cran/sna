######################################################################
#
# models.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 1/05/05
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines related to stochastic models.
#
# Contents:
#   bbnam
#   bbnam.actor
#   bbnam.bf
#   bbnam.fixed
#   bbnam.jntlik
#   bbnam.jntlik.slice
#   bbnam.pooled
#   coef.lnam
#   consensus
#   eval.edgeperturbation
#   lnam
#   netcancor
#   netlm
#   netlogit
#   npostpred
#   potscalered.mcmc
#   print.bbnam
#   print.bbnam.actor
#   print.bbnam.fixed
#   print.bbnam.pooled
#   print.lnam
#   print.netcancor
#   print.netlm
#   print.netlogit
#   print.summary.bbnam
#   print.summary.bbnam.actor
#   print.summary.bbnam.fixed
#   print.summary.bbnam.pooled
#   print.summary.lnam
#   print.summary.netcancor
#   print.summary.netlm
#   print.summary.netlogit
#   pstar
#   se.lnam
#   summary.bbnam
#   summary.bbnam.actor
#   summary.bbnam.fixed
#   summary.bbnam.pooled
#   summary.lnam
#   summary.netcancor
#   summary.netlm
#   summary.netlogit 
#
######################################################################


#bbnam - Draw from Butts' Bayesian Network Accuracy Model.  This version uses a 
#Gibbs' Sampler, and assumes error rates to be drawn from conditionally 
#independent betas for each actor.  Note that dat MUST be an n x n x n array, 
#and that the data in question MUST be dichotomous.  Priors are also assumed to 
#be in the right form (n x 2 matrices of alpha, beta pairs for em and ep, and
#an n x n probability matrix for the network itself), and are not checked; 
#default behavior if no priors are provided is the uninformative case.
#Wrapper function for the various bbnam models
bbnam<-function(dat,model="actor",...){
   if(model=="actor")
      bbnam.actor(dat,...)
   else if(model=="pooled")
      bbnam.pooled(dat,...)
   else if(model=="fixed")
      bbnam.fixed(dat,...)
}


#bbnam.actor - Draw from the error-prob-by-actor model
bbnam.actor<-function(dat,nprior=matrix(rep(0.5,dim(dat)[2]*dim(dat)[3]),nrow=dim(dat)[2],ncol=dim(dat)[3]),emprior=cbind(rep(1,dim(dat)[1]),rep(1,dim(dat)[1])),epprior=cbind(rep(1,dim(dat)[1]),rep(1,dim(dat)[1])),diag=FALSE, mode="digraph",reps=5,draws=1500,burntime=500,quiet=TRUE,anames=paste("a",1:dim(dat)[2],sep=""),onames=paste("o",1:dim(dat)[1],sep=""),compute.sqrtrhat=TRUE){
   #First, collect some basic model parameters and do other "setup" stuff
   m<-dim(dat)[1]
   n<-dim(dat)[2]
   d<-dat
   slen<-burntime+floor(draws/reps)
   out<-list()
   #Remove any data which doesn't count...
   if(mode=="graph")
      d<-upper.tri.remove(d)
   if(!diag)
      d<-diag.remove(d)
   #OK, let's get started.  First, create temp variables to hold draws, and draw
   #initial conditions for the Markov chain
   if(!quiet)
      cat("Creating temporary variables and drawing initial conditions....\n")
   a<-array(dim=c(reps,slen,n,n))
   em<-array(dim=c(reps,slen,m))
   ep<-array(dim=c(reps,slen,m))
   for(k in 1:reps){
      a[k,1,,]<-rgraph(n,1,diag=diag,mode=mode)
      em[k,1,]<-runif(m,0,0.5)
      ep[k,1,]<-runif(m,0,0.5)
   }
   #Let the games begin: draw from the Gibbs' sampler
   for(i in 1:reps){
      for(j in 2:slen){
         if(!quiet)
            cat("Repetition",i,", draw",j,":\n\tDrawing adjacency matrix\n")
         #Create tie probability matrix
         ep.a<-aperm(array(sapply(ep[i,j-1,],rep,n^2),dim=c(n,n,m)),c(3,2,1))
         em.a<-aperm(array(sapply(em[i,j-1,],rep,n^2),dim=c(n,n,m)),c(3,2,1))
         pygt<-apply(d*(1-em.a)+(1-d)*em.a,c(2,3),prod,na.rm=TRUE)
         pygnt<-apply(d*ep.a+(1-d)*(1-ep.a),c(2,3),prod,na.rm=TRUE)
         tieprob<-(nprior*pygt)/(nprior*pygt+(1-nprior)*pygnt)
         #Draw Bernoulli graph
         a[i,j,,]<-rgraph(n,1,tprob=tieprob,mode=mode,diag=diag)
         if(!quiet)
            cat("\tAggregating binomial counts\n")
         cem<-matrix(nrow=m,ncol=2)
         cep<-matrix(nrow=m,ncol=2)
         for(x in 1:m){
               cem[x,1]<-sum((1-d[x,,])*a[i,j,,],na.rm=TRUE)
               cem[x,2]<-sum(d[x,,]*a[i,j,,],na.rm=TRUE)
               cep[x,1]<-sum(d[x,,]*(1-a[i,j,,]),na.rm=TRUE)
               cep[x,2]<-sum((1-d[x,,])*(1-a[i,j,,]),na.rm=TRUE)
         }
         if(!quiet)
            cat("\tDrawing error parameters\n")
         em[i,j,]<-rbeta(m,emprior[,1]+cem[,1],emprior[,2]+cem[,2])
         ep[i,j,]<-rbeta(m,epprior[,1]+cep[,1],epprior[,2]+cep[,2])
      }
   }
   if(!quiet)
      cat("Finished drawing from Markov chain.  Now computing potential scale reduction statistics.\n")
   if(compute.sqrtrhat){
      out$sqrtrhat<-vector()
      for(i in 1:n)
         for(j in 1:n)
            out$sqrtrhat<-c(out$sqrtrhat,potscalered.mcmc(aperm(a,c(2,1,3,4))[,,i,j]))
      for(i in 1:m)
         out$sqrtrhat<-c(out$sqrtrhat,potscalered.mcmc(aperm(em,c(2,1,3))[,,i]),potscalered.mcmc(aperm(ep,c(2,1,3))[,,i]))
      if(!quiet)
         cat("\tMax potential scale reduction (Gelman et al.'s sqrt(Rhat)) for all scalar estimands:",max(out$sqrtrhat[!is.nan(out$sqrtrhat)],na.rm=TRUE),"\n")
   }
   if(!quiet)
      cat("Preparing output.\n")
   #Whew, we're done with the MCMC.  Now, let's get that data together.
   out$net<-array(dim=c(reps*(slen-burntime),n,n))
   for(i in 1:reps)
      for(j in burntime:slen){
         out$net[(i-1)*(slen-burntime)+(j-burntime),,]<-a[i,j,,]
      }
   if(!quiet)
      cat("\tAggregated network variable draws\n")
   out$em<-em[1,(burntime+1):slen,]
   out$ep<-ep[1,(burntime+1):slen,]
   if(reps>=2)
      for(i in 2:reps){
         out$em<-rbind(out$em,em[i,(burntime+1):slen,])
         out$ep<-rbind(out$ep,ep[i,(burntime+1):slen,])
      }
   if(!quiet)
      cat("\tAggregated error parameters\n")
   #Mix up draws (keeping components together, of course!) to reduce dependence
   o<-sample(1:dim(out$em)[1])
   out$net<-out$net[o,,]
   out$em<-out$em[o,]
   out$ep<-out$ep[o,]
   if(!quiet)
      cat("Remixing draws\n")
   #Finish off the output and return it.
   out$anames<-anames
   out$onames<-onames
   out$nactors<-n
   out$nobservers<-m
   out$reps<-reps
   out$draws<-dim(out$em)[1]
   out$burntime<-burntime
   out$model<-"actor"
   class(out)<-c("bbnam.actor","bbnam")
   out
}


#bbnam.bf - Estimate Bayes Factors for the Butts Bayesian Network Accuracy 
#Model.  This implementation relies on monte carlo integration to estimate the 
#BFs, and tests the fixed probability, pooled, and pooled by actor models.
bbnam.bf<-function(dat,nprior=matrix(rep(0.5,dim(dat)[1]^2),nrow=dim(dat)[1],ncol=dim(dat)[1]),em.fp=0.5,ep.fp=0.5,emprior.pooled=c(1,1),epprior.pooled=c(1,1),emprior.actor=cbind(rep(1,dim(dat)[1]),rep(1,dim(dat)[1])),epprior.actor=cbind(rep(1,dim(dat)[1]),rep(1,dim(dat)[1])),diag=FALSE, mode="digraph",reps=1000){
   n<-dim(dat)[1]
   d<-dat
   if(!diag)
      d<-diag.remove(d)
   if(mode=="graph")
      d<-lower.tri.remove(d)
   pfpv<-vector()
   ppov<-vector()
   pacv<-vector()
   #Draw em, ep, and a values for the various models
   for(i in 1:reps){
      a<-rgraph(n,1,tprob=nprior)
      em.pooled<-eval(call("rbeta",1,emprior.pooled[1],emprior.pooled[2]))
      ep.pooled<-eval(call("rbeta",1,epprior.pooled[1],epprior.pooled[2]))
      em.actor<-eval(call("rbeta",n,emprior.actor[,1],emprior.actor[,2]))
      ep.actor<-eval(call("rbeta",n,epprior.actor[,1],epprior.actor[,2]))
      pfpv[i]<-bbnam.jntlik(d,a=a,em=em.fp,ep=ep.fp)
      ppov[i]<-bbnam.jntlik(d,a=a,em=em.pooled,ep=ep.pooled)
      pacv[i]<-bbnam.jntlik(d,a=a,em=em.actor,ep=ep.actor)
   }
   int.lik<-c(mean(pfpv),mean(ppov),mean(pacv))
   int.lik.std<-sqrt(c(var(pfpv),var(ppov),var(pacv)))
   #Find the Bayes Factors
   o<-list()
   o$int.lik<-matrix(nrow=3,ncol=3)
   for(i in 1:3)
      for(j in 1:3){
         if(i!=j)
            o$int.lik[i,j]<-int.lik[i]/int.lik[j]
         else
            o$int.lik[i,i]<-int.lik[i]
      }
   o$int.lik.std<-int.lik.std
   o$reps<-reps
   o$prior.param<-list(nprior,em.fp,ep.fp,emprior.pooled,epprior.pooled,emprior.actor,epprior.actor)
   o$prior.param.names<-c("nprior","em.fp","ep.fp","emprior.pooled","epprior.pooled","emprior.actor","epprior.actor")
   o$model.names<-c("Fixed Error Prob","Pooled Error Prob","Actor Error Prob")
   class(o)<-c("bbnam.bf","bayes.factor")
   o
}


#bbnam.fixed - Draw from the fixed probability error model
bbnam.fixed<-function(dat,nprior=matrix(rep(0.5,dim(dat)[2]^2),nrow=dim(dat)[2],ncol=dim(dat)[2]),em=0.25,ep=0.25,diag=FALSE,mode="digraph",draws=1500,outmode="draws",anames=paste("a",1:dim(dat)[2],sep=""),onames=paste("o",1:dim(dat)[1],sep="")){
   #How many actors are involved?
   m<-dim(dat)[1]
   n<-dim(dat)[2]
   #Check to see if we've been given full matrices (or vectors) of error probs...
   if(length(em)==m*n^2)
      em.a<-em
   else if(length(em)==n^2)
      em.a<-apply(em,c(1,2),rep,m)
   else if(length(em)==m)
      em.a<-aperm(array(sapply(em,rep,n^2),dim=c(n,n,m)),c(3,2,1))
   else if(length(em)==1)
      em.a<-array(rep(em,m*n^2),dim=c(m,n,n))
   if(length(ep)==m*n^2)
      ep.a<-ep
   else if(length(ep)==n^2)
      ep.a<-apply(ep,c(1,2),rep,m)
   else if(length(ep)==m)
      ep.a<-aperm(array(sapply(ep,rep,n^2),dim=c(n,n,m)),c(3,2,1))
   else if(length(ep)==1)
      ep.a<-array(rep(ep,m*n^2),dim=c(m,n,n))
   #Find the network posterior
   pygt<-apply(dat*(1-em.a)+(1-dat)*em.a,c(2,3),prod,na.rm=TRUE)
   pygnt<-apply(dat*ep.a+(1-dat)*(1-ep.a),c(2,3),prod,na.rm=TRUE)
   npost<-(nprior*pygt)/(nprior*pygt+(1-nprior)*pygnt)
   #Send the needed output
   if(outmode=="posterior")
      npost
   else{
      o<-list()
      o$net<-rgraph(n,draws,tprob=npost,diag=diag,mode=mode)
      o$anames<-anames
      o$onames<-onames
      o$nactors<-n
      o$nobservers<-m
      o$draws<-draws
      o$model<-"fixed"
      class(o)<-c("bbnam.fixed","bbnam")
      o
   }
}


#bbnam.jntlik - An internal function for bbnam
bbnam.jntlik<-function(dat,log=FALSE,...){
   p<-sum(sapply(1:dim(dat)[1],bbnam.jntlik.slice,dat=dat,log=TRUE,...))
   if(!log)
      exp(p)
   else
      p
}


#bbnam.jntlik.slice - An internal function for bbnam
bbnam.jntlik.slice<-function(s,dat,a,em,ep,log=FALSE){
   if(length(em)>1)
      em.l<-em[s]
   else
      em.l<-em
   if(length(ep)>1)
      ep.l<-ep[s]
   else
      ep.l<-ep
   p<-sum(log((1-a)*(dat[s,,]*ep.l+(1-dat[s,,])*(1-ep.l))+a*(dat[s,,]*(1-em.l)+(1-dat[s,,])*em.l)),na.rm=TRUE)
   if(!log)
      exp(p)
   else
      p
}


#bbnam.pooled - Draw from the pooled error model
bbnam.pooled<-function(dat,nprior=matrix(rep(0.5,dim(dat)[2]*dim(dat)[3]),nrow=dim(dat)[2],ncol=dim(dat)[3]),emprior=c(1,1),epprior=c(1,1),diag=FALSE, mode="digraph",reps=5,draws=1500,burntime=500,quiet=TRUE,anames=paste("a",1:dim(dat)[2],sep=""),onames=paste("o",1:dim(dat)[1],sep=""),compute.sqrtrhat=TRUE){
   #First, collect some basic model parameters and do other "setup" stuff
   m<-dim(dat)[1]
   n<-dim(dat)[2]
   d<-dat
   slen<-burntime+floor(draws/reps)
   out<-list()
   #Remove any data which doesn't count...
   if(mode=="graph")
      d<-upper.tri.remove(d)
   if(!diag)
      d<-diag.remove(d)
   #OK, let's get started.  First, create temp variables to hold draws, and draw
   #initial conditions for the Markov chain
   if(!quiet)
      cat("Creating temporary variables and drawing initial conditions....\n")
   a<-array(dim=c(reps,slen,n,n))
   em<-array(dim=c(reps,slen))
   ep<-array(dim=c(reps,slen))
   for(k in 1:reps){
      a[k,1,,]<-rgraph(n,1,diag=diag,mode=mode)
      em[k,1]<-runif(1,0,0.5)
      ep[k,1]<-runif(1,0,0.5)
   }
   #Let the games begin: draw from the Gibbs' sampler
   for(i in 1:reps){
      for(j in 2:slen){
         if(!quiet)
            cat("Repetition",i,", draw",j,":\n\tDrawing adjacency matrix\n")
         #Create tie probability matrix
         ep.a<-array(rep(ep[i,j-1],m*n^2),dim=c(m,n,n))
         em.a<-array(rep(em[i,j-1],m*n^2),dim=c(m,n,n))
         pygt<-apply(d*(1-em.a)+(1-d)*em.a,c(2,3),prod,na.rm=TRUE)
         pygnt<-apply(d*ep.a+(1-d)*(1-ep.a),c(2,3),prod,na.rm=TRUE)
         tieprob<-(nprior*pygt)/(nprior*pygt+(1-nprior)*pygnt)
         #Draw Bernoulli graph
         a[i,j,,]<-rgraph(n,1,tprob=tieprob,mode=mode,diag=diag)
         if(!quiet)
            cat("\tAggregating binomial counts\n")
         cem<-vector(length=2)
         cep<-vector(length=2)
         a.a<-apply(a[i,j,,],c(1,2),rep,m)
         cem[1]<-sum((1-d)*a.a,na.rm=TRUE)
         cem[2]<-sum(d*a.a,na.rm=TRUE)
         cep[1]<-sum(d*(1-a.a),na.rm=TRUE)
         cep[2]<-sum((1-d)*(1-a.a),na.rm=TRUE)
         #cat("em - alpha",cem[1],"beta",cem[2]," ep - alpha",cep[1],"beta",cep[2],"\n")
         if(!quiet)
            cat("\tDrawing error parameters\n")
         em[i,j]<-rbeta(1,emprior[1]+cem[1],emprior[2]+cem[2])
         ep[i,j]<-rbeta(1,epprior[1]+cep[1],epprior[2]+cep[2])
      }
   }
   if(!quiet)
      cat("Finished drawing from Markov chain.  Now computing potential scale reduction statistics.\n")
   if(compute.sqrtrhat){
      out$sqrtrhat<-vector()
      for(i in 1:n)
         for(j in 1:n)
            out$sqrtrhat<-c(out$sqrtrhat,potscalered.mcmc(aperm(a,c(2,1,3,4))[,,i,j]))
      out$sqrtrhat<-c(out$sqrtrhat,potscalered.mcmc(em),potscalered.mcmc(ep))
   if(!quiet)
      cat("\tMax potential scale reduction (Gelman et al.'s sqrt(Rhat)) for all scalar estimands:",max(out$sqrtrhat[!is.nan(out$sqrtrhat)],na.rm=TRUE),"\n")
   }
   if(!quiet)
      cat("Preparing output.\n")
   #Whew, we're done with the MCMC.  Now, let's get that data together.
   out$net<-array(dim=c(reps*(slen-burntime),n,n))
   for(i in 1:reps)
      for(j in burntime:slen){
         out$net[(i-1)*(slen-burntime)+(j-burntime),,]<-a[i,j,,]
      }
   if(!quiet)
      cat("\tAggregated network variable draws\n")
   out$em<-em[1,(burntime+1):slen]
   out$ep<-ep[1,(burntime+1):slen]
   if(reps>=2)
      for(i in 2:reps){
         out$em<-c(out$em,em[i,(burntime+1):slen])
         out$ep<-c(out$ep,ep[i,(burntime+1):slen])
      }
   if(!quiet)
      cat("\tAggregated error parameters\n")
   #Mix up draws (keeping components together, of course!) to reduce dependence
   o<-sample(1:length(out$em))
   out$net<-out$net[o,,]
   out$em<-out$em[o]
   out$ep<-out$ep[o]
   if(!quiet)
      cat("Remixing draws\n")
   #Finish off the output and return it.
   out$anames<-anames
   out$onames<-onames
   out$nactors<-n
   out$nobservers<-m
   out$reps<-reps
   out$draws<-length(out$em)
   out$burntime<-burntime
   out$model<-"pooled"
   class(out)<-c("bbnam.pooled","bbnam")
   out
}


#bbnam.probtie - Probability of a given tie
bbnam.probtie<-function(dat,i,j,npriorij,em,ep){
   num<-npriorij
   denom<-1-npriorij
   num<-num*prod(dat[,i,j]*(1-em)+(1-dat[,i,j])*em,na.rm=TRUE)
   denom<-denom*prod(dat[,i,j]*ep+(1-dat[,i,j])*(1-ep),na.rm=TRUE)
   p<-num/(denom+num)
   p
}


#coef.lnam - Coefficient method for lnam
coef.lnam<-function(object, ...){
   coefs<-vector()
   cn<-vector()
   if(!is.null(object$rho1)){
      coefs<-c(coefs,object$rho1)
      cn<-c(cn,"rho1")
   }
   if(!is.null(object$rho2)){
      coefs<-c(coefs,object$rho2)
      cn<-c(cn,"rho2")
   }
   if(!is.null(object$beta)){
      coefs<-c(coefs,object$beta)
      cn<-c(cn,names(object$beta))
   }
   names(coefs)<-cn
   coefs
}


#consensus - Find a consensus structure, using one of several algorithms.  Note 
#that this is currently experimental, and that the routines are not guaranteed 
#to produce meaningful output
consensus<-function(dat,mode="digraph",diag=FALSE,method="central.graph",tol=0.01){
   n<-dim(dat)[2]
   m<-dim(dat)[1]
   #First, prepare the data
   if(m==1)
     dat<-array(dat,dim=c(1,n,n))
   if(mode=="graph")
      d<-upper.tri.remove(dat)
   else
      d<-dat
   if(!diag)
      d<-diag.remove(d)
   #Now proceed by method
   #First, use the central graph if called for
   if(method=="central.graph"){
      cong<-centralgraph(d)
   #Try the iterative reweighting algorithm....
   }else if(method=="iterative.reweight"){
      stop("Sorry, but iterative rewieghting is not currently supported.\n")
      oldrwv<-rep(1/m,m)
      gc<-gcor(d)
      gc[is.na(gc)]<-0
      diag(gc)<-1
      rwv<-apply(gc,1,sum)/sum(gc)
      while(sum(abs(rwv-oldrwv))>tol){
         cong<-apply(d*aperm(array(sapply(rwv,rep,n^2),dim=c(n,n,m)),c(3,2,1)),c(2,3),sum)
         oldrwv<-rwv
         rwv<-gcor(cong,d)
         rwv<-rwv/sum(rwv)
      }
   #Perform a single reweighting using mean correlation
   }else if(method=="single.reweight"){
      gc<-gcor(d)
      gc[is.na(gc)]<-0
      diag(gc)<-1
      rwv<-apply(gc,1,sum)
      rwv<-rwv/sum(rwv)
      cong<-apply(d*aperm(array(sapply(rwv,rep,n^2),dim=c(n,n,m)),c(3,2,1)),c(2,3),sum)
   #Perform a single reweighting using first component loadings
   }else if(method=="PCA.reweight"){
      gc<-gcor(d)
      gc[is.na(gc)]<-0
      diag(gc)<-1
      rwv<-eigen(gc)$vector[,1]
      cong<-apply(d*aperm(array(sapply(rwv,rep,n^2),dim=c(n,n,m)),c(3,2,1)),c(2,3),sum)
   #Use the Locally Aggregated Structure
   }else if(method=="LAS.intersection"){
      cong<-matrix(0,n,n)
      for(i in 1:n)
        for(j in 1:n)
          cong[i,j]<-as.numeric(d[i,i,j]&&d[j,i,j])
   }else if(method=="LAS.union"){
      cong<-matrix(0,n,n)
      for(i in 1:n)
        for(j in 1:n)
          cong[i,j]<-as.numeric(d[i,i,j]||d[j,i,j])
   }else if(method=="OR.row"){
      cong<-matrix(0,n,n)
      for(i in 1:n)
         cong[i,]<-d[i,i,]
   }else if(method=="OR.col"){
      cong<-matrix(0,n,n)
      for(i in 1:n)
         cong[,i]<-d[i,,i]
   }
   #Finish off and return the consensus graph
   if(mode=="graph")
      cong[upper.tri(cong)]<-t(cong)[upper.tri(cong)]
   if(!diag)
      diag(cong)<-0
   cong
}


#eval.edgeperturbation - Evaluate a function on a given graph with and without a
#given edge, returning the difference between the results in each case.
eval.edgeperturbation<-function(dat,i,j,FUN,...){
   #Get the function in question
   fun<-match.fun(FUN)
   #Set up the perturbation matrices
   present<-dat
   present[i,j]<-1
   absent<-dat
   absent[i,j]<-0
   #Evaluate the function across the perturbation and return the difference
   fun(present,...)-fun(absent,...)
}


#lnam - Fit a linear network autocorrelation model
#y = r1 * W1 %*% y + X %*% b + e, e = r2 * W2 %*% e + nu
#y =  (I-r1*W1)^-1%*%(X %*% b + e)
#y = (I-r1 W1)^-1 (X %*% b + (I-r2 W2)^-1 nu)
#e = (I-r2 W2)^-1 nu
#e = (I-r1 W1) y - X b
#nu = (I - r2 W2) [ (I-r1 W1) y - X b ]
#nu = (I-r2 W2) e
lnam<-function(y,x=NULL,W1=NULL,W2=NULL,theta.seed=NULL,null.model=c("meanstd","mean","std","none"),method="BFGS",control=list()){
   #Define the log-likelihood functions for each case
   lnLx<-function(theta,y,x,sigma.log=TRUE){ #theta=c(s,b)
      m<-length(theta)
      if(sigma.log)
        sig<-exp(theta[1])
      else
        sig<-theta[1]
      -2*sum(dnorm(y-x%*%(theta[2:m]),0,sig,log=TRUE))
   }
   lnL1<-function(theta,y,W1,sigma.log=TRUE){ #theta=c(r1,s)
      n<-length(y)
      if(sigma.log)
        sig<-exp(theta[2])
      else
        sig<-theta[2]
      -2*sum(dnorm((diag(n)-theta[1]*W1)%*%y,0,sig,log=TRUE))
   }
   lnL2<-function(theta,y,W2,sigma.log=TRUE){ #theta=c(r2,s)
      n<-length(y)
      if(sigma.log)
        sig<-exp(theta[2])
      else
        sig<-theta[2]
      -2*sum(dnorm((diag(n)-theta[1]*W2)%*%y,0,sig,log=TRUE))
   }
   lnLx1<-function(theta,y,x,W1,sigma.log=TRUE){ #theta=c(r1,s,b)
      n<-length(y)
      m<-length(theta)
      if(sigma.log)
        sig<-exp(theta[2])
      else
        sig<-theta[2]
      -2*sum(dnorm((diag(n)-theta[1]*W1)%*%y-x%*%theta[3:m],0,sig,log=TRUE))
   }
   lnLx2<-function(theta,y,x,W2,sigma.log=TRUE){ #theta=c(r2,s,b)
      n<-length(y)
      m<-length(theta)
      if(sigma.log)
        sig<-exp(theta[2])
      else
        sig<-theta[2]
      -2*sum(dnorm((diag(n)-theta[1]*W2)%*%(y-x%*%theta[3:m]),0,sig,log=TRUE))
   }
   lnL12<-function(theta,y,W1,W2,sigma.log=TRUE){ #theta=c(r1,r2,s)
      n<-length(y)
      if(sigma.log)
        sig<-exp(theta[3])
      else
        sig<-theta[3]
      -2*sum(dnorm((diag(n)-theta[2]*W2)%*%((diag(n)-theta[1]*W1)%*%y),0,sig,log=TRUE))
   }
   lnLx12<-function(theta,y,x,W1,W2,sigma.log=TRUE){ #theta=c(r1,r2,s,b)
      n<-length(y)
      m<-length(theta)
      if(sigma.log)
        sig<-exp(theta[3])
      else
        sig<-theta[3]
      -2*sum(dnorm((diag(n)-theta[2]*W2)%*%((diag(n)-theta[1]*W1)%*%y-x%*%theta[4:m]),0,sig,log=TRUE))
   }
   #How many data points are there?
   n<-length(y)
   #Fix x, if needed
   if(!is.null(x)&&is.vector(x))
     x<-as.matrix(x)
   #Determine the computation mode from the x,W1,W2 parameters
   comp.mode<-as.character(as.numeric(1*(!is.null(x))+10*(!is.null(W1))+100*(!is.null(W2))))
   if(comp.mode=="0")
      stop("At least one of x, W1, W2 must be specified.\n")
   #How many predictors?   
   m<-switch(comp.mode,
      "1"=dim(x)[2]+1,
      "10"=2,
      "100"=2,
      "11"=dim(x)[2]+2,
      "101"=dim(x)[2]+2,
      "110"=3,
      "111"=dim(x)[2]+3
   )
   #Initialize the parameter vector
   if(is.null(theta.seed)){
      theta<-rep(0,m)
   }else{
      theta<-theta.seed
      if(comp.mode=="1")           #Log the standard deviation parameter
         theta[1]<-log(theta[1])
      else if(comp.mode%in%c("10","100","11","101"))
         theta[2]<-log(theta[2])
      else
         theta[3]<-log(theta[3])
   }
   #Perform the MLE fit via a two-stage process
   fitted<-switch(comp.mode,
      "1"=optim(theta,lnLx,method=method,control=control,y=y,x=x),
      "10"=optim(theta,lnL1,method=method,control=control,y=y,W1=W1),
      "100"=optim(theta,lnL2,method=method,control=control,y=y,W2=W2),
      "11"=optim(theta,lnLx1,method=method,control=control,y=y,x=x,W1=W1),
      "101"=optim(theta,lnLx2,method=method,control=control,y=y,x=x,W2=W2),
      "110"=optim(theta,lnL12,method=method,control=control,y=y,W1=W1,W2=W2),
      "111"=optim(theta,lnLx12,method=method,control=control,y=y,x=x,W1=W1,W2=W2)
   )
   if(comp.mode=="1")           #De-log the standard deviation parameter
      fitted$par[1]<-exp(fitted$par[1])
   else if(comp.mode%in%c("10","100","11","101"))
      fitted$par[2]<-exp(fitted$par[2])
   else
      fitted$par[3]<-exp(fitted$par[3])
   theta<-fitted$par            #Prepare for the stage-2 fit
   fitted<-switch(comp.mode,
      "1"=optim(theta,lnLx,method=method,control=control,hessian=TRUE,y=y,x=x,sigma.log=FALSE),
      "10"=optim(theta,lnL1,method=method,control=control,hessian=TRUE,y=y,W1=W1,sigma.log=FALSE),
      "100"=optim(theta,lnL2,method=method,control=control,hessian=TRUE,y=y,W2=W2,sigma.log=FALSE),
      "11"=optim(theta,lnLx1,method=method,control=control,hessian=TRUE,y=y,x=x,W1=W1,sigma.log=FALSE),
      "101"=optim(theta,lnLx2,method=method,control=control,hessian=TRUE,y=y,x=x,W2=W2,sigma.log=FALSE),
      "110"=optim(theta,lnL12,method=method,control=control,hessian=TRUE,y=y,W1=W1,W2=W2,sigma.log=FALSE),
      "111"=optim(theta,lnLx12,method=method,control=control,hessian=TRUE,y=y,x=x,W1=W1,W2=W2,sigma.log=FALSE)
   )
   #Assemble and return the results
   o<-list()
   o$y<-y
   o$x<-x
   o$W1<-W1
   o$W2<-W2
   o$model<-comp.mode
   o$infomat<-fitted$hessian/2  
   o$acvm<-qr.solve(o$infomat)
   o$null.model<-match.arg(null.model)
   o$lnlik.null<-switch(match.arg(null.model),  #Fit a null model
      "meanstd"=sum(dnorm(y-mean(y),0,as.numeric(sqrt(var(y))),log=TRUE)),
      "mean"=sum(dnorm(y-mean(y),log=TRUE)),
      "std"=sum(dnorm(y,0,as.numeric(sqrt(var(y))),log=TRUE)),
      "none"=sum(dnorm(y,log=TRUE))
   )
   o$df.null.resid<-switch(match.arg(null.model),  #Find residual null df
      "meanstd"=n-2,
      "mean"=n-1,
      "std"=n-1,
      "none"=n
   )
   o$df.null<-switch(match.arg(null.model),  #Find null df
      "meanstd"=2,
      "mean"=1,
      "std"=1,
      "none"=0
   )
   o$null.param<-switch(match.arg(null.model),  #Find null params, if any
      "meanstd"=c(mean(y),sqrt(var(y))),
      "mean"=mean(y),
      "std"=sqrt(var(y)),
      "none"=NULL
   )
   o$lnlik.model<--fitted$value/2
   o$df.model<-m
   o$df.residual<-n-m
   o$df.total<-n
   o$rho1<-switch(comp.mode,   #Get the r1 parameter, if available
      "1"=NULL,
      "10"=fitted$par[1],
      "100"=NULL,
      "11"=fitted$par[1],
      "101"=NULL,
      "110"=fitted$par[1],
      "111"=fitted$par[1]
   )
   o$rho1.se<-switch(comp.mode,   #Get the r1 SE, if available
      "1"=NULL,
      "10"=sqrt(o$acvm[1,1]),
      "100"=NULL,
      "11"=sqrt(o$acvm[1,1]),
      "101"=NULL,
      "110"=sqrt(o$acvm[1,1]),
      "111"=sqrt(o$acvm[1,1])
   )
   o$rho2<-switch(comp.mode,   #Get the r2 parameter, if available
      "1"=NULL,
      "10"=NULL,
      "100"=fitted$par[1],
      "11"=NULL,
      "101"=fitted$par[1],
      "110"=fitted$par[2],
      "111"=fitted$par[2]
   )
   o$rho2.se<-switch(comp.mode,   #Get the r2 SE, if available
      "1"=NULL,
      "10"=NULL,
      "100"=sqrt(o$acvm[1,1]),
      "11"=NULL,
      "101"=sqrt(o$acvm[1,1]),
      "110"=sqrt(o$acvm[2,2]),
      "111"=sqrt(o$acvm[2,2])
   )
   o$sigma<-switch(comp.mode,   #Get the sigma parameter
      "1"=(fitted$par[1]),
      "10"=(fitted$par[2]),
      "100"=(fitted$par[2]),
      "11"=(fitted$par[2]),
      "101"=(fitted$par[2]),
      "110"=(fitted$par[3]),
      "111"=(fitted$par[3])
   )
   o$sigma.se<-switch(comp.mode,   #Get the sigma SE
      "1"=sqrt(o$acvm[1,1]),
      "10"=sqrt(o$acvm[2,2]),
      "100"=sqrt(o$acvm[2,2]),
      "11"=sqrt(o$acvm[2,2]),
      "101"=sqrt(o$acvm[2,2]),
      "110"=sqrt(o$acvm[3,3]),
      "111"=sqrt(o$acvm[3,3])
   )
   o$beta<-as.vector(switch(comp.mode,   #Get the beta parameters, if available
      "1"=fitted$par[2:m],
      "10"=NULL,
      "100"=NULL,
      "11"=fitted$par[3:m],
      "101"=fitted$par[3:m],
      "110"=NULL,
      "111"=fitted$par[4:m]
   ))
   o$beta.se<-as.vector(switch(comp.mode,   #Get the beta SE, if available
      "1"=sqrt(diag(o$acvm)[2:m]),
      "10"=NULL,
      "100"=NULL,
      "11"=sqrt(diag(o$acvm)[3:m]),
      "101"=sqrt(diag(o$acvm)[3:m]),
      "110"=NULL,
      "111"=sqrt(diag(o$acvm)[4:m])
   ))
   if(!is.null(colnames(x))){
      names(o$beta)<-colnames(x)
      names(o$beta.se)<-colnames(x)
   }else{
      names(o$beta)<-paste("X",1:dim(x)[2],sep="")
      names(o$beta.se)<-paste("X",1:dim(x)[2],sep="")
   }
   o$disturbances<-as.vector(switch(comp.mode,  #The estimated disturbances
      "1"=y-x%*%o$beta,
      "10"=(diag(n)-o$rho1*W1)%*%y,
      "100"=(diag(n)-o$rho2*W2)%*%y,
      "11"=(diag(n)-o$rho1*W1)%*%y-x%*%o$beta,
      "101"=(diag(n)-o$rho2*W2)%*%(y-x%*%o$beta),
      "110"=(diag(n)-o$rho2*W2)%*%((diag(n)-o$rho1*W1)%*%y),
      "111"=(diag(n)-o$rho2*W2)%*%((diag(n)-o$rho1*W1)%*%y-x%*%o$beta)
   ))
   o$fitted.values<-as.vector(switch(comp.mode,  #Compute the fitted values
      "1"=x%*%o$beta,
      "10"=rep(0,n),
      "100"=rep(0,n),
      "11"=qr.solve(diag(n)-o$rho1*W1,x%*%o$beta),
      "101"=x%*%o$beta,
      "110"=rep(0,n),
      "111"=qr.solve(diag(n)-o$rho1*W1,x%*%o$beta)
   ))
   o$residuals<-as.vector(y-o$fitted.values)
   o$call<-match.call()
   class(o)<-c("lnam")
   o
}


#netcancor - Canonical correlations for network variables.  NOTE: requires that 
#mva be loaded for R<2.0.0.
netcancor<-function(y,x,mode="digraph",diag=FALSE,nullhyp="cugtie",reps=1000){
   if(R.version$major<2)  #Only invoke mva if we're using an old R version
      require(mva)
   if(length(dim(y))>2){
      iy<-matrix(nrow=dim(y)[1],ncol=dim(y)[2]*dim(y)[3])
   }else{
      iy<-matrix(nrow=1,ncol=dim(y)[1]*dim(y)[2])
      temp<-y
      y<-array(dim=c(1,dim(temp)[1],dim(temp)[2]))
      y[1,,]<-temp
   }
   if(length(dim(x))>2){
      ix<-matrix(nrow=dim(x)[1],ncol=dim(x)[2]*dim(x)[3])
   }else{
      ix<-matrix(nrow=1,ncol=dim(x)[1]*dim(x)[2])
      temp<-x
      x<-array(dim=c(1,dim(temp)[1],dim(temp)[2]))
      x[1,,]<-temp
   }
   my<-dim(y)[1]
   mx<-dim(x)[1]
   n<-dim(y)[2]
   out<-list()
   out$xdist<-array(dim=c(reps,mx,mx))
   out$ydist<-array(dim=c(reps,my,my))
   #Convert the response first.
   for(i in 1:my){
      d<-y[i,,]
      #if(!diag){
      #   diag(d)<-NA
      #}
      #if(mode!="digraph")
      #   d[lower.tri(d)]<-NA
      iy[i,]<-as.vector(d)
   }
   #Now for the independent variables.
   for(i in 1:mx){
      d<-x[i,,]
      #if(!diag){
      #   diag(d)<-NA
      #}
      #if(mode!="digraph")
      #   d[lower.tri(d)]<-NA
      ix[i,]<-as.vector(d)
   }   
   #Run the initial model fit
   nc<-cancor(t(ix),t(iy))  #Had to take out na.action=na.omit, since it's not supported
   #Now, repeat the whole thing an ungodly number of times.
   out$cdist<-array(dim=c(reps,length(nc$cor)))
   for(i in 1:reps){
      #Clear out the internal structures
      iy<-matrix(nrow=dim(y)[1],ncol=dim(y)[2]*dim(y)[3])
      ix<-matrix(nrow=dim(x)[1],ncol=dim(x)[2]*dim(x)[3])
      #Convert (and mutate) the response first.
      for(j in 1:my){
         d<-switch(nullhyp,
            qap = rmperm(y[j,,]),
            cug = rgraph(n,1,mode=mode,diag=diag),
            cugden = rgraph(n,1,tprob=gden(y[j,,],mode=mode,diag=diag),mode=mode,diag=diag),
            cugtie = rgraph(n,1,mode=mode,diag=diag,tielist=y[j,,])
         )
         #if(!diag){
         #   diag(d)<-NA
         #}
         #if(mode!="digraph")
         #   d[lower.tri(d)]<-NA
         iy[j,]<-as.vector(d)
      }
      #Now for the independent variables.
      for(j in 1:mx){
         d<-switch(nullhyp,
            qap = rmperm(x[j,,]),
            cug = rgraph(n,1,mode=mode,diag=diag),
            cugden = rgraph(n,1,tprob=gden(x[j,,],mode=mode,diag=diag),mode=mode,diag=diag),
            cugtie = rgraph(n,1,mode=mode,diag=diag,tielist=x[j,,])
         )
         #if(!diag){
         #   diag(d)<-NA
         #}
         #if(mode!="digraph")
         #   d[lower.tri(d)]<-NA
         ix[j,]<-as.vector(d)
      }   
      #Finally, fit the test model
      tc<-cancor(t(ix),t(iy))         #Had to take out na.action=na.omit, since it's not supported
      #Gather the coefficients for use later...
      out$cdist[i,]<-tc$cor
      out$xdist[i,,]<-tc$xcoef
      out$ydist[i,,]<-tc$ycoef
   }
   #Find the p-values for our monte carlo null hypothesis tests
   out$cor<-nc$cor
   out$xcoef<-nc$xcoef
   out$ycoef<-nc$ycoef
   out$cpgreq<-vector(length=length(nc$cor))
   out$cpleeq<-vector(length=length(nc$cor))
   for(i in 1:length(nc$cor)){
      out$cpgreq[i]<-mean(out$cdist[,i]>=out$cor[i],na.rm=TRUE)
      out$cpleeq[i]<-mean(out$cdist[,i]<=out$cor[i],na.rm=TRUE)
   }
   out$xpgreq<-matrix(ncol=mx,nrow=mx)
   out$xpleeq<-matrix(ncol=mx,nrow=mx)
   for(i in 1:mx){
      for(j in 1:mx){ 
         out$xpgreq[i,j]<-mean(out$xdist[,i,j]>=out$xcoef[i,j],na.rm=TRUE)
         out$xpleeq[i,j]<-mean(out$xdist[,i,j]<=out$xcoef[i,j],na.rm=TRUE)
      }
   }
   out$ypgreq<-matrix(ncol=my,nrow=my)
   out$ypleeq<-matrix(ncol=my,nrow=my)
   for(i in 1:my){
      for(j in 1:my){ 
         out$ypgreq[i,j]<-mean(out$ydist[,i,j]>=out$ycoef[i,j],na.rm=TRUE)
         out$ypleeq[i,j]<-mean(out$ydist[,i,j]<=out$ycoef[i,j],na.rm=TRUE)
      }
   }
   #Having completed the model fit and MC tests, we gather useful information for
   #the end user.  This is a combination of cancor output and our own stuff.
   out$cnames<-as.vector(paste("cor",1:min(mx,my),sep=""))
   out$xnames<-as.vector(paste("x",1:mx,sep=""))
   out$ynames<-as.vector(paste("y",1:my,sep=""))
   out$xcenter<-nc$xcenter
   out$ycenter<-nc$ycenter
   out$nullhyp<-nullhyp
   class(out)<-c("netcancor")
   out
}


#netlm - OLS network regression routine using a QAP/CUG null hypotheses.  This routine is
#frighteningly slow, since it's essentially a front end to the builtin lm routine with a bunch of
#network hypothesis testing stuff thrown in for good measure.
#netlm2<-function(y,x,intercept=TRUE,mode="digraph",diag=FALSE,nullhyp="cugtie",reps=1000){
#  #Create the output list
#  out<-list()
#  out$r.squared.dist<-vector(length=reps)
#  out$adj.r.squared.dist<-vector(length=reps)
#  out$sigma.dist<-vector(length=reps)
#  #Perform the initial vectorization
#  vy<-as.vector(gvectorize(y,mode=mode,diag=diag))
#  vx<-gvectorize(x,mode=mode,diag=diag)
#  #Add an intercept, if needed
#  if(intercept)
#    vx<-cbind(rep(1,dim(vx)[1]),vx)
#  #Get some initial stats
#  n<-dim(y)[1]
#  p<-dim(vx)[2]
#  #Perform the initial model fit
#  xnam<-paste("vx[,",1:p,"]",sep="")
#  fmstr<-paste("vy ~ ",paste(xnam,collapse="+"))
#  fmla<-as.formula(fmstr)
#  nm<-lm(fmla,na.action=na.omit,singular.ok=TRUE)
#  #Now, repeat everything to test the appropriate null hypothesis
#  if(match.arg(nullhyp)=="qap"){ 
#  #QAP semi-partialling "plus"
#  #For Y ~ b0 + b1 X1 + b2 X2 + ... + bp Xp
#  #for(i in 1:p)
#  #  Fit Xi ~ b0* + b1* X1 + ... + bp* Xp (omit Xi)
#  #  Let ei = resid of above lm
#  #  for(j in 1:reps)
#  #    eij = rmperm (ei)
#  #    Fit Y ~ b0** + b1** X1 + ... + bi** eij + ... + bp** Xp
#  #Use resulting permutation distributions to test coefficients
#    #Identify rows without missing data (w/out permutations)
#    nona<-apply(!is.na(cbind(vy,vx)),1,all)
#    #Walk through the predictors
#    for(i in 1:p){
#      #Regress the appropriate X on its peers
#      xm<-lm.fit(vx[nona,-i,drop=FALSE],vx[nona,i])
#      #Convert the residuals of this regression back into matrix form
#      ex<-nona
#      ex[!nona]<-NA
#      ex[nona]<-xm$residuals
#      ex<-matrix(ex,n,n)
#      if(mode=="graph")   #If the matrix is symmetric, restore it
#        ex[upper.tri(ex)]<-t(ex)[upper.tri(ex)]
#      #Perform the QAP replications
#      for(j in 1:reps){
#        #Set up the new predictors
#	tx<-vx
#        tx[,i]<-gvectorize(rmperm(ex),mode=mode,diag=diag)
#        rnona<-apply(!is.na(cbind(vy,tx)),1,all)
#	#Fit the test model
#	tm<-lm.fit(vy[rnona],tx[rnona,])
#	#Gather the coefficient, for later use
#	out$dist[j,i]<-tm$coefficients[i]
#      }
#    }
#  }else{
#  }
#}
netlm<-function(y,x,mode="digraph",diag=FALSE,nullhyp="cugtie",reps=1000){
   out<-list()
   out$r.squared.dist<-vector(length=reps)
   out$adj.r.squared.dist<-vector(length=reps)
   out$sigma.dist<-vector(length=reps)
   iy<-vector()
   if(length(dim(x))>2){
      ix<-matrix(nrow=dim(x)[1],ncol=dim(x)[2]*dim(x)[3])
   }else{
      ix<-matrix(nrow=1,ncol=dim(x)[1]*dim(x)[2])
      temp<-x
      x<-array(dim=c(1,dim(temp)[1],dim(temp)[2]))
      x[1,,]<-temp
   }
   n<-dim(y)[1]
   m<-dim(x)[1]
   out$dist<-matrix(nrow=reps,ncol=m+1)
   #Convert the response first.
   d<-y
   if(!diag){
      diag(d)<-NA
   }
   if(mode!="digraph")
      d[lower.tri(d)]<-NA
   iy<-as.vector(d)
   #Now for the independent variables.
   for(i in 1:m){
      d<-x[i,,]
      if(!diag){
         diag(d)<-NA
      }
      if(mode!="digraph")
         d[lower.tri(d)]<-NA
      ix[i,]<-as.vector(d)
   }   
   #Run the initial model fit
   xnam <- paste("ix[", 1:m, ",]", sep="")
   fmla <- as.formula(paste("iy ~ ", paste(xnam, collapse= "+")))
   nm<-lm(fmla,na.action=na.omit,singular.ok=TRUE)
   #Now, repeat the whole thing an ungodly number of times.
   for(i in 1:reps){
      #Clear out the internal structures
      iy<-vector()
      ix<-matrix(nrow=dim(x)[1],ncol=dim(x)[2]*dim(x)[3])
      #Convert (and mutate) the response first.
      d<-switch(nullhyp,
         qap = rmperm(y),
         cug = rgraph(n,1,mode=mode,diag=diag),
         cugden = rgraph(n,1,tprob=gden(y,mode=mode,diag=diag),mode=mode,diag=diag),
         cugtie = rgraph(n,1,mode=mode,diag=diag,tielist=y)
      )
      if(!diag){
         diag(d)<-NA
      }
      if(mode!="digraph")
         d[lower.tri(d)]<-NA
      iy<-as.vector(d)
      #Now for the independent variables.
      for(j in 1:m){
         d<-switch(nullhyp,
            qap = rmperm(x[j,,]),
            cug = rgraph(n,1,mode=mode,diag=diag),
            cugden = rgraph(n,1,tprob=gden(x[j,,],mode=mode,diag=diag),mode=mode,diag=diag),
            cugtie = rgraph(n,1,mode=mode,diag=diag,tielist=x[j,,])
         )
         if(!diag){
            diag(d)<-NA
         }
         if(mode!="digraph")
            d[lower.tri(d)]<-NA
         ix[j,]<-as.vector(d)
      }   
      #Finally, fit the test model
      xnam <- paste("ix[", 1:m, ",]", sep="")
      fmla <- as.formula(paste("iy ~ ", paste(xnam, collapse= "+")))
      tm<-lm(fmla,na.action=na.omit,singular.ok=TRUE)
      #Gather the coefficients for use later...
      out$dist[i,]<-as.numeric(coef(tm))
      #Also grab R^2, sigma, and adjusted R^2s
      mss<-if(attr(tm$terms,"intercept"))
         sum((fitted(tm)-mean(fitted(tm)))^2)
      else
         sum(fitted(tm)^2)
      rss<-sum(resid(tm)^2)
      qn<-NROW(tm$qr$qr)
      df.int<-if(attr(tm$terms,"intercept")) 1
         else 0
      rdf<-qn-tm$rank
      out$r.squared.dist[i]<-mss/(mss+rss)
      out$adj.r.squared.dist[i]<-1-(1-out$r.squared.dist[i])*((qn-df.int)/rdf)
      out$sigma.dist[i]<-sqrt(rss/rdf)
   }
   #Find the p-values for our monte carlo null hypothesis tests
   out$coefficients<-nm$coefficients
   out$pgreq<-vector(length=m+1)
   out$pleeq<-vector(length=m+1)
   for(i in 1:(m+1)){
      out$pgreq[i]<-mean(out$dist[,i]>=out$coefficients[i],na.rm=TRUE)
      out$pleeq[i]<-mean(out$dist[,i]<=out$coefficients[i],na.rm=TRUE)
   }
   #Having completed the model fit and MC tests, we gather useful information for
   #the end user.  This is a combination of GLM output and our own stuff.
   out$names<-as.vector(c("(intercept)",paste("x",1:m,sep="")))
   out$nullhyp<-nullhyp
   out$residuals<-nm$residuals
   out$qr<-nm$qr
   out$fitted.values<-nm$fitted.values
   out$rank<-nm$rank
   out$terms<-nm$terms
   out$df.residual<-nm$df.residual
   class(out)<-c("netlm")
   out
}


#netlogit - God help me, it's a network regression routine using a 
#binomial/logit GLM.  It's also frighteningly slow, since it's essentially a 
#front end to the builtin GLM routine with a bunch of network hypothesis testing
#stuff thrown in for good measure.
netlogit<-function(y,x,mode="digraph",diag=FALSE,nullhyp="cugtie",reps=1000){
   out<-list()
   out$dist<-matrix(nrow=reps,ncol=dim(x)[1]+1)
   iy<-vector()
   if(length(dim(x))>2){
      ix<-matrix(nrow=dim(x)[1],ncol=dim(x)[2]*dim(x)[3])
   }else{
      ix<-matrix(nrow=1,ncol=dim(x)[1]*dim(x)[2])
      temp<-x
      x<-array(dim=c(1,dim(temp)[1],dim(temp)[2]))
      x[1,,]<-temp
   }
   n<-dim(y)[1]
   m<-dim(x)[1]
   out$dist<-matrix(nrow=reps,ncol=m+1)
   #Convert the response first.
   d<-y
   if(!diag){
      diag(d)<-NA
   }
   if(mode!="digraph")
      d[lower.tri(d)]<-NA
   iy<-as.vector(d)
   #Now for the independent variables.
   for(i in 1:m){
      d<-x[i,,]
      if(!diag){
         diag(d)<-NA
      }
      if(mode!="digraph")
         d[lower.tri(d)]<-NA
      ix[i,]<-as.vector(d)
   }   
   #Run the initial model fit
   xnam <- paste("ix[", 1:m, ",]", sep="")
   fmla <- as.formula(paste("iy ~ ", paste(xnam, collapse= "+")))
   nm<-glm(fmla,family=binomial,na.action=na.omit)
   out$ctable<-table(as.numeric(fitted.values(nm)>=0.5),iy[!is.na(iy)],dnn=c("Predicted","Actual"))  #Get the contingency table 
   #Now, repeat the whole thing an ungodly number of times.
   for(i in 1:reps){
      #Clear out the internal structures
      iy<-vector()
      ix<-matrix(nrow=dim(x)[1],ncol=dim(x)[2]*dim(x)[3])
      #Convert (and mutate) the response first.
      d<-switch(nullhyp,
         qap = rmperm(y),
         cug = rgraph(n,1,mode=mode,diag=diag),
         cugden = rgraph(n,1,tprob=gden(y,mode=mode,diag=diag),mode=mode,diag=diag),
         cugtie = rgraph(n,1,mode=mode,diag=diag,tielist=y)
      )
      if(!diag){
         diag(d)<-NA
      }
      if(mode!="digraph")
         d[lower.tri(d)]<-NA
      iy<-as.vector(d)
      #Now for the independent variables.
      for(j in 1:m){
         d<-switch(nullhyp,
            qap = rmperm(x[j,,]),
            cug = rgraph(n,1,mode=mode,diag=diag),
            cugden = rgraph(n,1,tprob=gden(x[j,,],mode=mode,diag=diag),mode=mode,diag=diag),
            cugtie = rgraph(n,1,mode=mode,diag=diag,tielist=x[j,,])
         )
         if(!diag){
            diag(d)<-NA
         }
         if(mode!="digraph")
            d[lower.tri(d)]<-NA
         ix[j,]<-as.vector(d)
      }   
      #Finally, fit the test model
      xnam <- paste("ix[", 1:m, ",]", sep="")
      fmla <- as.formula(paste("iy ~ ", paste(xnam, collapse= "+")))
      tm<-glm(fmla,family=binomial,na.action=na.omit)
      #Gather the coefficients for use later...
      out$dist[i,]<-as.numeric(coef(tm))
   }
   #Find the p-values for our monte carlo null hypothesis tests
   out$coefficients<-nm$coefficients
   out$pgreq<-vector(length=m+1)
   out$pleeq<-vector(length=m+1)
   for(i in 1:(m+1)){
      out$pgreq[i]<-mean(out$dist[,i]>=out$coefficients[i],na.rm=TRUE)
      out$pleeq[i]<-mean(out$dist[,i]<=out$coefficients[i],na.rm=TRUE)
   }
   #Having completed the model fit and MC tests, we gather useful information for
   #the end user.  This is a combination of GLM output and our own stuff.
   out$names<-as.vector(c("(intercept)",paste("x",1:m,sep="")))
   out$nullhyp<-nullhyp
   out$deviance<-nm$deviance
   out$df.residual<-nm$df.residual
   out$df.null<-nm$df.null
   out$aic<-nm$aic
   out$null.deviance<-nm$null.deviance
   class(out)<-c("netlogit","netglm","netlm")
   out
}


#npostpred - Take posterior predictive draws for functions of networks.
npostpred<-function(b,FUN,...){
   #Find the desired function
   fun<-match.fun(FUN)
   #Take the draws
   out<-apply(b$net,1,fun,...)
   out
}


#plot.bbnam - Plot method for bbnam
plot.bbnam<-function(x,mode="density",intlines=TRUE,...){
   UseMethod("plot",x)
}


#plot.bbnam.actor - Plot method for bbnam.actor
plot.bbnam.actor<-function(x,mode="density",intlines=TRUE,...){
   #Get the initial graphical settings, so we can restore them later
   oldpar<-par()
   #Change plotting params
   par(ask=TRUE)
   #Initial plot: global error distribution
   par(mfrow=c(2,1))
   if(mode=="density"){   #Approximate the pdf using kernel density estimation
      #Plot marginal population (i.e. across actors) density of p(false negative)
      plot(density(x$em),main=paste("Estimated Marginal Population Density of",expression(e^"-"),",",x$draws,"Draws"),xlab=expression({e^{"-"}}),xlim=c(0,1),...)
      #Plot interval lines if required.
      if(intlines)
         abline(v=quantile(x$em,c(0.05,0.5,0.95)),lty=c(3,2,3))
      #Plot marginal population (i.e. across actors) density of p(false positive)
      plot(density(x$ep),main=paste("Estimated Marginal Population Density of",expression(e^"+"),",",x$draws,"Draws"),xlab=expression({e^{"+"}}),xlim=c(0,1),...)
      #Plot interval lines if required.
      if(intlines)
         abline(v=quantile(x$ep,c(0.05,0.5,0.95)),lty=c(3,2,3))
   }else{     #Use histograms to plot the estimated density
      #Plot marginal population (i.e. across actors) density of p(false negative)
      hist(x$em,main=paste("Histogram of",expression(e^"-"),",",x$draws,"Draws"),xlab=expression({e^{"-"}}),xlim=c(0,1),...)
      #Plot interval lines if required.
      if(intlines)
         abline(v=quantile(x$em,c(0.05,0.5,0.95)),lty=c(3,2,3))
      #Plot marginal population (i.e. across actors) density of p(false positive)
      hist(x$ep,main=paste("Histogram of",expression(e^"+"),",",x$draws,"Draws"),xlab=expression({e^{"+"}}),xlim=c(0,1),...)
      #Plot interval lines if required.
      if(intlines)
         abline(v=quantile(x$ep,c(0.05,0.5,0.95)),lty=c(3,2,3))
   }
   #Plot e- next
   par(mfrow=c(floor(sqrt(x$nobservers)),ceiling(sqrt(x$nobservers))))
   for(i in 1:x$nobservers){
      if(mode=="density"){
         plot(density(x$em[,i]),main=paste("Estimated Density of",expression(e^"-"[i]),",",x$draws,"Draws"),xlab=expression({e^{"-"}}[i]),xlim=c(0,1),...)
         #Plot interval lines if required.
         if(intlines)
            abline(v=quantile(x$em[,i],c(0.05,0.5,0.95)),lty=c(3,2,3))
      }else{
         hist(x$em[,i],main=paste("Histogram of",expression(e^"-"[i]),",",x$draws,"Draws"),xlab=expression({e^{"-"}}[i]),xlim=c(0,1),...)
         #Plot interval lines if required.
         if(intlines)
            abline(v=quantile(x$em[,i],c(0.05,0.5,0.95)),lty=c(3,2,3))
      }
   }
   #Now plot e+
   par(mfrow=c(floor(sqrt(x$nobservers)),ceiling(sqrt(x$nobservers))))
   for(i in 1:x$nobservers){
      if(mode=="density"){
         plot(density(x$ep[,i]),main=paste("Estimated Density of",expression(e^"+"[i]),",",x$draws,"Draws"),xlab=expression({e^{"+"}}[i]),xlim=c(0,1),...)
         #Plot interval lines if required.
         if(intlines)
            abline(v=quantile(x$ep[,i],c(0.05,0.5,0.95)),lty=c(3,2,3))
      }else{
         hist(x$ep[,i],main=paste("Histogram of",expression(e^"+"[i]),",",x$draws,"Draws"),xlab=expression({e^{"+"}}[i]),xlim=c(0,1),...)
         #Plot interval lines if required.
         if(intlines)
            abline(v=quantile(x$ep[,i],c(0.05,0.5,0.95)),lty=c(3,2,3))
      }
   }
   #Finally, try to plot histograms of tie probabilities
   par(mfrow=c(1,1))
   plot.sociomatrix(apply(x$net,c(2,3),mean),labels=list(x$anames,x$anames),main="Marginal Posterior Tie Probability Distribution")
   #Clean up
   par(oldpar)
}


#plot.bbnam.fixed - Plot method for bbnam.fixed
plot.bbnam.fixed<-function(x,mode="density",intlines=TRUE,...){
   #Get the initial graphical settings, so we can restore them later
   oldpar<-par()
   #Perform matrix plot of tie probabilities
   par(mfrow=c(1,1))
   plot.sociomatrix(apply(x$net,c(2,3),mean),labels=list(x$anames,x$anames),main="Marginal Posterior Tie Probability Distribution")
   #Clean up
   par(oldpar)
}


#plot.bbnam.pooled - Plot method for bbnam.pooled
plot.bbnam.pooled<-function(x,mode="density",intlines=TRUE,...){
   #Get the initial graphical settings, so we can restore them later
   oldpar<-par()
   #Change plotting params
   par(ask=TRUE)
   #Initial plot: pooled error distribution
   par(mfrow=c(2,1))
   if(mode=="density"){   #Approximate the pdf using kernel density estimation
      #Plot marginal population (i.e. across actors) density of p(false negative)
      plot(density(x$em),main=paste("Estimated Marginal Posterior Density of",expression(e^"-"),",",x$draws,"Draws"),xlab=expression({e^{"-"}}),xlim=c(0,1),...)
      #Plot interval lines if required.
      if(intlines)
         abline(v=quantile(x$em,c(0.05,0.5,0.95)),lty=c(3,2,3))
      #Plot marginal population (i.e. across actors) density of p(false positive)
      plot(density(x$ep),main=paste("Estimated Marginal Posterior Density of",expression(e^"+"),",",x$draws,"Draws"),xlab=expression({e^{"+"}}),xlim=c(0,1),...)
      #Plot interval lines if required.
      if(intlines)
         abline(v=quantile(x$ep,c(0.05,0.5,0.95)),lty=c(3,2,3))
   }else{     #Use histograms to plot the estimated density
      #Plot marginal population (i.e. across actors) density of p(false negative)
      hist(x$em,main=paste("Histogram of",expression(e^"-"),",",x$draws,"Draws"),xlab=expression({e^{"-"}}),xlim=c(0,1),...)
      #Plot interval lines if required.
      if(intlines)
         abline(v=quantile(x$em,c(0.05,0.5,0.95)),lty=c(3,2,3))
      #Plot marginal population (i.e. across actors) density of p(false positive)
      hist(x$ep,main=paste("Histogram of",expression(e^"+"),",",x$draws,"Draws"),xlab=expression({e^{"+"}}),xlim=c(0,1),...)
      #Plot interval lines if required.
      if(intlines)
         abline(v=quantile(x$ep,c(0.05,0.5,0.95)),lty=c(3,2,3))
   }
   #Finally, try to plot histograms of tie probabilities
   par(mfrow=c(1,1))
   plot.sociomatrix(apply(x$net,c(2,3),mean),labels=list(x$anames,x$anames),main="Marginal Posterior Tie Probability Distribution")
   #Clean up
   par(oldpar)
}


#plot.lnam - Plot method for lnam
plot.lnam<-function(x,...){
   if(R.version$major<2)  #Only invoke mva if we're using an old R version
     require(mva)
   r<-residuals(x)
   f<-fitted(x)
   d<-x$disturbances
   sdr<-sd(r)
   ci<-c(-1.959964,1.959964)
   old.par <- par(no.readonly = TRUE)
   on.exit(par(old.par))
   par(mfrow=c(2,2))
   #Plot residual versus actual values
   plot(x$y,f,ylab=expression(hat(y)),xlab=expression(y),main="Fitted vs. Observed Values")
   abline(ci[1]*sdr,1,lty=3)
   abline(0,1,lty=2)
   abline(ci[2]*sdr,1,lty=3)
   #Plot disturbances versus fitted values
   plot(f,d,ylab=expression(hat(nu)),xlab=expression(hat(y)), ylim=c(min(ci[1]*x$sigma,d),max(ci[2]*x$sigma,d)),main="Fitted Values vs. Estimated Disturbances")
   abline(h=c(ci[1]*x$sigma,0,ci[2]*x$sigma),lty=c(3,2,3))
   #QQ-Plot the residuals
   qqnorm(r,main="Normal Q-Q Residual Plot")
   qqline(r)
   #Plot an influence diagram
   if(!(is.null(x$W1)&&is.null(x$W2))){
      inf<-matrix(0,nc=x$df.total,nr=x$df.total)
      if(!is.null(x$W1))
         inf<-inf+qr.solve(diag(x$df.total)-x$rho1*x$W1)
      if(!is.null(x$W2))
         inf<-inf+qr.solve(diag(x$df.total)-x$rho2*x$W2)
      syminf<-abs(inf)+abs(t(inf))
      diag(syminf)<-0
      infco<-cmdscale(as.dist(max(syminf)-syminf),k=2)
      diag(inf)<-NA
      stdinf<-inf-mean(inf,na.rm=TRUE)
      infsd<-sd(as.vector(stdinf),na.rm=TRUE)
      stdinf<-stdinf/infsd
      gplot(abs(stdinf),thresh=1.96,coord=infco,main="Net Influence Plot",edge.lty=1,edge.lwd=abs(stdinf)/2,edge.col=2+(inf>0)) 
   }
   #Restore plot settings
   invisible()
}


#potscalered.mcmc - Potential scale reduction (sqrt(Rhat)) for scalar estimands.
#Input must be a matrix whose columns correspond to replicate chains.  This, 
#clearly, doesn't belong here, but lacking a better place to put it I have 
#included it nonetheless.
potscalered.mcmc<-function(psi){
   #Use Gelman et al. notation, for convenience
   J<-dim(psi)[2]
   n<-dim(psi)[1]
   #Find between-group variance estimate
   mpsij<-apply(psi,2,mean)
   mpsitot<-mean(mpsij)
   B<-(n/(J-1))*sum((mpsij-mpsitot)^2)
   #Find within-group variance estimate
   s2j<-apply(psi,2,var)
   W<-mean(s2j)
   #Calculate the (estimated) marginal posterior variance of the estimand
   varppsi<-((n-1)/n)*W+(1/n)*B
   #Return the potential scale reduction estimate
   sqrt(varppsi/W)
}


#print.bayes.factor - A fairly generic routine for printing bayes factors, here used for the bbnam routine.
print.bayes.factor<-function(x,...){
   tab<-x$int.lik
   rownames(tab)<-x$model.names
   colnames(tab)<-x$model.names
   cat("Bayes Factors by Model:\n\n(Diagonals indicate raw integrated likelihood estimates.)\n\n")
   print(tab)
   cat("\n")   
}


#print.bbnam - Print method for bbnam
print.bbnam<-function(x,...){
   UseMethod("print",x)
}


#print.bbnam.actor - Print method for bbnam.actor
print.bbnam.actor<-function(x,...){
   cat("\nButts' Hierarchical Bayes Model for Network Estimation/Informant Accuracy\n\n")
   cat("Multiple Error Probability Model\n\n")
   #Dump marginal posterior network
   cat("Marginal Posterior Network Distribution:\n\n")
   d<-apply(x$net,c(2,3),mean)
   rownames(d)<-as.vector(x$anames)
   colnames(d)<-as.vector(x$anames)
   print.table(d,digits=2)
   cat("\n")
   #Dump summary of error probabilities
   cat("Marginal Posterior Global Error Distribution:\n\n")
   d<-matrix(ncol=2,nrow=6)
   d[1:3,1]<-quantile(x$em,c(0,0.25,0.5),names=FALSE,na.rm=TRUE)
   d[4,1]<-mean(x$em,na.rm=TRUE)
   d[5:6,1]<-quantile(x$em,c(0.75,1.0),names=FALSE,na.rm=TRUE)
   d[1:3,2]<-quantile(x$ep,c(0,0.25,0.5),names=FALSE,na.rm=TRUE)
   d[4,2]<-mean(x$ep,na.rm=TRUE)
   d[5:6,2]<-quantile(x$ep,c(0.75,1.0),names=FALSE,na.rm=TRUE)
   colnames(d)<-c("e^-","e^+")
   rownames(d)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
   print.table(d,digits=4)
   cat("\n")
}


#print.bbnam.fixed - Print method for bbnam.fixed
print.bbnam.fixed<-function(x,...){
   cat("\nButts' Hierarchical Bayes Model for Network Estimation/Informant Accuracy\n\n")
   cat("Fixed Error Probability Model\n\n")
   #Dump marginal posterior network
   cat("Marginal Posterior Network Distribution:\n\n")
   d<-apply(x$net,c(2,3),mean)
   rownames(d)<-as.vector(x$anames)
   colnames(d)<-as.vector(x$anames)
   print.table(d,digits=2)
   cat("\n")
}


#print.bbnam.pooled - Print method for bbnam.pooled
print.bbnam.pooled<-function(x,...){
   cat("\nButts' Hierarchical Bayes Model for Network Estimation/Informant Accuracy\n\n")
   cat("Pooled Error Probability Model\n\n")
   #Dump marginal posterior network
   cat("Marginal Posterior Network Distribution:\n\n")
   d<-apply(x$net,c(2,3),mean)
   rownames(d)<-as.vector(x$anames)
   colnames(d)<-as.vector(x$anames)
   print.table(d,digits=2)
   cat("\n")
   #Dump summary of error probabilities
   cat("Marginal Posterior Global Error Distribution:\n\n")
   d<-matrix(ncol=2,nrow=6)
   d[1:3,1]<-quantile(x$em,c(0,0.25,0.5),names=FALSE,na.rm=TRUE)
   d[4,1]<-mean(x$em,na.rm=TRUE)
   d[5:6,1]<-quantile(x$em,c(0.75,1.0),names=FALSE,na.rm=TRUE)
   d[1:3,2]<-quantile(x$ep,c(0,0.25,0.5),names=FALSE,na.rm=TRUE)
   d[4,2]<-mean(x$ep,na.rm=TRUE)
   d[5:6,2]<-quantile(x$ep,c(0.75,1.0),names=FALSE,na.rm=TRUE)
   colnames(d)<-c("e^-","e^+")
   rownames(d)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
   print.table(d,digits=4)
   cat("\n")
}


#print.lnam - Print method for lnam
print.lnam<-function(x, digits = max(3, getOption("digits") - 3), ...){
   cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
   cat("Coefficients:\n")
   print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
   cat("\n")
}


#print.netcancor - Print method for netcancor
print.netcancor<-function(x,...){
   cat("\nCanonical Network Correlation\n\n")

   cat("Canonical Correlations:\n\n")
   cmat<-matrix(data=x$cor,ncol=length(x$cor),nrow=1)
   rownames(cmat)<-""
   colnames(cmat)<-as.vector(x$cnames)
   print.table(cmat)
   cat("\n")
   cat("Pr(>=cor):\n\n")
   cmat <- matrix(data=format(x$cpgreq),ncol=length(x$cpgreq),nrow=1)
   colnames(cmat) <- as.vector(x$cnames)
   rownames(cmat)<- ""
   print.table(cmat)
   cat("\n")
   cat("Pr(<=cor):\n\n")
   cmat <- matrix(data=format(x$cpleeq),ncol=length(x$cpleeq),nrow=1)
   colnames(cmat) <- as.vector(x$cnames)
   rownames(cmat)<- ""
   print.table(cmat)
   cat("\n")

   cat("X Coefficients:\n\n")
   cmat <- format(x$xcoef)
   colnames(cmat) <- as.vector(x$xnames)
   rownames(cmat)<- as.vector(x$xnames)
   print.table(cmat)
   cat("\n")
   cat("Pr(>=xcoef):\n\n")
   cmat <- format(x$xpgreq)
   colnames(cmat) <- as.vector(x$xnames)
   rownames(cmat)<- as.vector(x$xnames)
   print.table(cmat)
   cat("\n")
   cat("Pr(<=xcoef):\n\n")
   cmat <- format(x$xpleeq)
   colnames(cmat) <- as.vector(x$xnames)
   rownames(cmat)<- as.vector(x$xnames)
   print.table(cmat)
   cat("\n")

   cat("Y Coefficients:\n\n")
   cmat <- format(x$ycoef)
   colnames(cmat) <- as.vector(x$ynames)
   rownames(cmat)<- as.vector(x$ynames)
   print.table(cmat)
   cat("\n")
   cat("Pr(>=ycoef):\n\n")
   cmat <- format(x$ypgreq)
   colnames(cmat) <- as.vector(x$ynames)
   rownames(cmat)<- as.vector(x$ynames)
   print.table(cmat)
   cat("\n")
   cat("Pr(<=ycoef):\n\n")
   cmat <- format(x$ypleeq)
   colnames(cmat) <- as.vector(x$ynames)
   rownames(cmat)<- as.vector(x$ynames)
   print.table(cmat)
   cat("\n")
}


#print.netlm - Print method for netlm
print.netlm<-function(x,...){
   cat("\nOLS Network Model\n\n")
   cat("Coefficients:\n\n")
   cmat <- as.vector(format(as.numeric(x$coefficients)))
   cmat <- cbind(cmat, as.vector(format(x$pgreq)))
   cmat <- cbind(cmat, as.vector(format(x$pleeq)))
   colnames(cmat) <- c("Estimate", "Pr(>=b)", "Pr(<=b)")
   rownames(cmat)<- as.vector(x$names)
   print.table(cmat)
   #Goodness of fit measures
   mss<-if(attr(x$terms,"intercept"))
      sum((fitted(x)-mean(fitted(x)))^2)
   else
      sum(fitted(x)^2)
   rss<-sum(resid(x)^2)
   qn<-NROW(x$qr$qr)
   df.int<-if(attr(x$terms,"intercept")) 1
      else 0
   rdf<-qn-x$rank
   resvar<-rss/rdf
   fstatistic<-c(value=(mss/(x$rank-df.int))/resvar,numdf=x$rank-df.int,dendf=rdf)
   r.squared<-mss/(mss+rss)
   adj.r.squared<-1-(1-r.squared)*((qn-df.int)/rdf)
   sigma<-sqrt(resvar)
   cat("\nResidual standard error:",format(sigma,digits=4),"on",rdf,"degrees of freedom\n")
   cat("F-statistic:",formatC(fstatistic[1],digits=4),"on",fstatistic[2],"and",fstatistic[3],"degrees of freedom, \tp-value:",formatC(1-pf(fstatistic[1],fstatistic[2],fstatistic[3]),dig=4),"\n")
   cat("Multiple R^2:",format(r.squared,digits=4),"\n")
   cat("Adjusted R^2:",format(adj.r.squared,digits=4),"\n")
   cat("\n")
}


#print.netlogit - Print method for netlogit
print.netlogit<-function(x,...){
   cat("\nNetwork Logit Model\n\n")
   cat("Coefficients:\n\n")
   cmat <- as.vector(format(as.numeric(x$coefficients)))
   cmat <- cbind(cmat, as.vector(format(x$pgreq)))
   cmat <- cbind(cmat, as.vector(format(x$pleeq)))
   colnames(cmat) <- c("Estimate", "Pr(>=b)", "Pr(<=b)")
   rownames(cmat)<- as.vector(x$names)
   print.table(cmat)
   cat("\nGoodness of Fit Statistics:\n")
   cat("\nNull deviance (-2*Ln(L)):",x$null.deviance,"on",x$df.null,"degrees of freedom\n")
   cat("Residual deviance (-2*Ln(L)):",x$deviance,"on",x$df.residual,"degrees of freedom\n")
   cat("Chi-Squared test of fit improvement:\n\t",x$null.deviance-x$deviance,"on",x$df.null-x$df.residual,"degrees of freedom, p-value",1-pchisq(x$null.deviance-x$deviance,df=x$df.null-x$df.residual),"\n") 
   cat("AIC:",x$aic,"\tBIC:",x$deviance+log(x$df.null+1)*(x$df.null-x$df.residual),"\nPseudo-R^2 Measures:\n\t(Dn-Dr)/(Dn-Dr+dfn):",(x$null.deviance-x$deviance)/(x$null.deviance-x$deviance+x$df.null),"\n\t(Dn-Dr)/Dn:",1-x$deviance/x$null.deviance,"\n")
   cat("\n")
}


#print.summary.bayes.factor - Printing for bayes factor summary objects
print.summary.bayes.factor<-function(x,...){
   cat("Bayes Factors by Model:\n\n(Diagonals indicate raw integrated likelihood estimates.)\n\n")
   print(x$int.lik)
   stdtab<-matrix(x$int.lik.std,nrow=1)
   colnames(stdtab)<-x$model.names
   cat("\n\nInverse Bayes Factors:\n\n(Diagonals indicate posterior probability of model under within-set choice constraints and uniform model priors.\n\n")
   print(x$inv.bf)
   cat("\n\nDiagnostics:\n\nReplications - ",x$reps,"\n\nStd deviations of integrated likelihood estimates:\n\n")
   print(x$int.lik.std)
   cat("\n\nVector of hyperprior parameters:\n\n")
   priortab<-matrix(x$prior.param,nrow=1,ncol=length(x$prior.param))
   colnames(priortab)<-x$prior.param.names
   print(priortab)
   cat("\n\n")   
}


#print.summary.bbnam - Print method for summary.bbnam
print.summary.bbnam<-function(x,...){
   UseMethod("print",x)
}


#print.summary.bbnam.actor - Print method for summary.bbnam.actor
print.summary.bbnam.actor<-function(x,...){
   cat("\nButts' Hierarchical Bayes Model for Network Estimation/Informant Accuracy\n\n")
   cat("Multiple Error Probability Model\n\n")
   #Dump marginal posterior network
   cat("Marginal Posterior Network Distribution:\n\n")
   d<-apply(x$net,c(2,3),mean)
   rownames(d)<-as.vector(x$anames)
   colnames(d)<-as.vector(x$anames)
   print.table(d,digits=2)
   cat("\n")
   #Dump summary of error probabilities
   cat("Marginal Posterior Global Error Distribution:\n\n")
   d<-matrix(ncol=2,nrow=6)
   d[1:3,1]<-quantile(x$em,c(0,0.25,0.5),names=FALSE,na.rm=TRUE)
   d[4,1]<-mean(x$em,na.rm=TRUE)
   d[5:6,1]<-quantile(x$em,c(0.75,1.0),names=FALSE,na.rm=TRUE)
   d[1:3,2]<-quantile(x$ep,c(0,0.25,0.5),names=FALSE,na.rm=TRUE)
   d[4,2]<-mean(x$ep,na.rm=TRUE)
   d[5:6,2]<-quantile(x$ep,c(0.75,1.0),names=FALSE,na.rm=TRUE)
   colnames(d)<-c("e^-","e^+")
   rownames(d)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
   print.table(d,digits=4)
   cat("\n")
   #Dump error probability estimates per observer
   cat("Marginal Posterior Error Distribution (by observer):\n\n")
   cat("Probability of False Negatives (e^-):\n\n")
   d<-matrix(ncol=6)
   for(i in 1:x$nobservers){
      dv<-matrix(c(quantile(x$em[,i],c(0,0.25,0.5),names=FALSE,na.rm=TRUE),mean(x$em[,i],na.rm=TRUE),quantile(x$em[,i],c(0.75,1.0),names=FALSE,na.rm=TRUE)),nrow=1,ncol=6)
      d<-rbind(d,dv)
   }
   d<-d[2:(x$nobservers+1),]
   rownames(d)<-as.vector(x$onames)
   colnames(d)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
   print.table(d,digits=4)
   cat("\n")
   cat("Probability of False Positives (e^+):\n\n")
   d<-matrix(ncol=6)
   for(i in 1:x$nobservers){
      dv<-matrix(c(quantile(x$ep[,i],c(0,0.25,0.5),names=FALSE,na.rm=TRUE),mean(x$ep[,i],na.rm=TRUE),quantile(x$ep[,i],c(0.75,1.0),names=FALSE,na.rm=TRUE)),nrow=1,ncol=6)
      d<-rbind(d,dv)
   }
   d<-d[2:(x$nobservers+1),]
   rownames(d)<-as.vector(x$onames)
   colnames(d)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
   print.table(d,digits=4)
   cat("\n")
   #Dump MCMC diagnostics
   cat("MCMC Diagnostics:\n\n")
   cat("\tReplicate Chains:",x$reps,"\n")
   cat("\tBurn Time:",x$burntime,"\n")
   cat("\tDraws per Chain:",x$draws/x$reps,"Total Draws:",x$draws,"\n")
   if("sqrtrhat" %in% names(x))
      cat("\tPotential Scale Reduction (G&R's sqrt(Rhat)):\n \t\tMax:",max(x$sqrtrhat[!is.nan(x$sqrtrhat)]),"\n\t\tMed:",median(x$sqrtrhat[!is.nan(x$sqrtrhat)]),"\n\t\tIQR:",IQR(x$sqrtrhat[!is.nan(x$sqrtrhat)]),"\n")
   cat("\n")
}


#print.summary.bbnam.fixed - Print method for summary.bbnam.fixed
print.summary.bbnam.fixed<-function(x,...){
   cat("\nButts' Hierarchical Bayes Model for Network Estimation/Informant Accuracy\n\n")
   cat("Fixed Error Probability Model\n\n")
   #Dump marginal posterior network
   cat("Marginal Posterior Network Distribution:\n\n")
   d<-apply(x$net,c(2,3),mean)
   rownames(d)<-as.vector(x$anames)
   colnames(d)<-as.vector(x$anames)
   print.table(d,digits=2)
   cat("\n")
   #Dump model diagnostics
   cat("Model Diagnostics:\n\n")
   cat("\tTotal Draws:",x$draws,"\n\t(Note: Draws taken directly from network posterior.)")
   cat("\n")
}


#print.summary.bbnam.pooled - Print method for summary.bbnam.pooled
print.summary.bbnam.pooled<-function(x,...){
   cat("\nButts' Hierarchical Bayes Model for Network Estimation/Informant Accuracy\n\n")
   cat("Pooled Error Probability Model\n\n")
   #Dump marginal posterior network
   cat("Marginal Posterior Network Distribution:\n\n")
   d<-apply(x$net,c(2,3),mean)
   rownames(d)<-as.vector(x$anames)
   colnames(d)<-as.vector(x$anames)
   print.table(d,digits=2)
   cat("\n")
   #Dump summary of error probabilities
   cat("Marginal Posterior Error Distribution:\n\n")
   d<-matrix(ncol=2,nrow=6)
   d[1:3,1]<-quantile(x$em,c(0,0.25,0.5),names=FALSE,na.rm=TRUE)
   d[4,1]<-mean(x$em,na.rm=TRUE)
   d[5:6,1]<-quantile(x$em,c(0.75,1.0),names=FALSE,na.rm=TRUE)
   d[1:3,2]<-quantile(x$ep,c(0,0.25,0.5),names=FALSE,na.rm=TRUE)
   d[4,2]<-mean(x$ep,na.rm=TRUE)
   d[5:6,2]<-quantile(x$ep,c(0.75,1.0),names=FALSE,na.rm=TRUE)
   colnames(d)<-c("e^-","e^+")
   rownames(d)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
   print.table(d,digits=4)
   cat("\n")
   #Dump MCMC diagnostics
   cat("MCMC Diagnostics:\n\n")
   cat("\tReplicate Chains:",x$reps,"\n")
   cat("\tBurn Time:",x$burntime,"\n")
   cat("\tDraws per Chain:",x$draws/x$reps,"Total Draws:",x$draws,"\n")
   if("sqrtrhat" %in% names(x))
      cat("\tPotential Scale Reduction (G&R's sqrt(Rhat)):\n \t\tMax:",max(x$sqrtrhat[!is.nan(x$sqrtrhat)]),"\n\t\tMed:",median(x$sqrtrhat[!is.nan(x$sqrtrhat)]),"\n\t\tIQR:",IQR(x$sqrtrhat[!is.nan(x$sqrtrhat)]),"\n")
   cat("\n")
}


#print.summary.lnam - Print method for summary.lnam
print.summary.lnam<-function(x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), ...){
   cat("\nCall:\n")
   cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
   cat("Residuals:\n")
   nam <- c("Min", "1Q", "Median", "3Q", "Max")
   resid<-x$residuals 
   rq <- if (length(dim(resid)) == 2) 
      structure(apply(t(resid), 1, quantile), dimnames = list(nam, dimnames(resid)[[2]]))
   else structure(quantile(resid), names = nam)
   print(rq, digits = digits, ...)
   cat("\nCoefficients:\n")
   cmat<-cbind(coef(x),se.lnam(x))
   cmat<-cbind(cmat,cmat[,1]/cmat[,2],(1-pnorm(abs(cmat[,1]),0,cmat[,2]))*2)
   colnames(cmat)<-c("Estimate","Std. Error","Z value","Pr(>|z|)")
   #print(format(cmat,digits=digits),quote=FALSE)
   printCoefmat(cmat,digits=digits,signif.stars=signif.stars,...)
   cat("\n")
   cmat<-cbind(x$sigma,x$sigma.se)
   colnames(cmat)<-c("Estimate","Std. Error")
   rownames(cmat)<-"Sigma"
   printCoefmat(cmat,digits=digits,signif.stars=signif.stars,...)
   cat("\nGoodness-of-Fit:\n")
   rss<-sum(x$residuals^2)
   mss<-sum((x$fitted-mean(x$fitted))^2)
   rdfns<-x$df.residual+1
   cat("\tResidual standard error: ",format(sqrt(rss/rdfns),digits=digits)," on ",rdfns," degrees of freedom (w/o Sigma)\n",sep="")
   cat("\tMultiple R-Squared: ",format(mss/(mss+rss),digits=digits),", Adjusted R-Squared: ",format(1-(1-mss/(mss+rss))*x$df.total/rdfns,digits=digits),"\n",sep="")
   cat("\tModel log likelihood:", format(x$lnlik.model,digits=digits), "on", x$df.resid, "degrees of freedom (w/Sigma)\n\tAIC:",format(-2*x$lnlik.model+2*x$df.model,digits=digits),"BIC:",format(-2*x$lnlik.model+log(x$df.total)*x$df.model,digits=digits),"\n")
   cat("\n\tNull model:",x$null.model,"\n")
   cat("\tNull log likelihood:", format(x$lnlik.null,digits=digits), "on", x$df.null.resid, "degrees of freedom\n\tAIC:",format(-2*x$lnlik.null+2*x$df.null,digits=digits),"BIC:",format(-2*x$lnlik.null+log(x$df.total)*x$df.null,digits=digits),"\n")
   cat("\tAIC difference (model versus null):",format(-2*x$lnlik.null+2*x$df.null+2*x$lnlik.model-2*x$df.model,digits=digits),"\n")
   cat("\tHeuristic Log Bayes Factor (model versus null): ",format(-2*x$lnlik.null+log(x$df.total)*x$df.null+2*x$lnlik.model-log(x$df.total)*x$df.model,digits=digits),"\n")
   cat("\n")
}


#print.summary.netcancor - Print method for summary.netcancor
print.summary.netcancor<-function(x,...){
   cat("\nCanonical Network Correlation\n\n")

   cat("Canonical Correlations:\n\n")
   cmat<-as.vector(x$cor)
   cmat<-rbind(cmat,as.vector((x$cor)^2))
   rownames(cmat)<-c("Correlation","Coef. of Det.")
   colnames(cmat)<-as.vector(x$cnames)
   print.table(cmat)
   cat("\n")
   cat("Pr(>=cor):\n\n")
   cmat <- matrix(data=format(x$cpgreq),ncol=length(x$cpgreq),nrow=1)
   colnames(cmat) <- as.vector(x$cnames)
   rownames(cmat)<- ""
   print.table(cmat)
   cat("\n")
   cat("Pr(<=cor):\n\n")
   cmat <- matrix(data=format(x$cpleeq),ncol=length(x$cpleeq),nrow=1)
   colnames(cmat) <- as.vector(x$cnames)
   rownames(cmat)<- ""
   print.table(cmat)
   cat("\n")

   cat("X Coefficients:\n\n")
   cmat <- format(x$xcoef)
   colnames(cmat) <- as.vector(x$xnames)
   rownames(cmat)<- as.vector(x$xnames)
   print.table(cmat)
   cat("\n")
   cat("Pr(>=xcoef):\n\n")
   cmat <- format(x$xpgreq)
   colnames(cmat) <- as.vector(x$xnames)
   rownames(cmat)<- as.vector(x$xnames)
   print.table(cmat)
   cat("\n")
   cat("Pr(<=xcoef):\n\n")
   cmat <- format(x$xpleeq)
   colnames(cmat) <- as.vector(x$xnames)
   rownames(cmat)<- as.vector(x$xnames)
   print.table(cmat)
   cat("\n")

   cat("Y Coefficients:\n\n")
   cmat <- format(x$ycoef)
   colnames(cmat) <- as.vector(x$ynames)
   rownames(cmat)<- as.vector(x$ynames)
   print.table(cmat)
   cat("\n")
   cat("Pr(>=ycoef):\n\n")
   cmat <- format(x$ypgreq)
   colnames(cmat) <- as.vector(x$ynames)
   rownames(cmat)<- as.vector(x$ynames)
   print.table(cmat)
   cat("\n")
   cat("Pr(<=ycoef):\n\n")
   cmat <- format(x$ypleeq)
   colnames(cmat) <- as.vector(x$ynames)
   rownames(cmat)<- as.vector(x$ynames)
   print.table(cmat)
   cat("\n")

   cat("Test Diagnostics:\n\n")
   cat("\tNull Hypothesis:")
   if(x$nullhyp=="qap")
      cat(" QAP\n")
   else
      cat(" CUG\n")
   cat("\tReplications:",dim(x$cdist)[1],"\n")
   cat("\tDistribution Summary for Correlations:\n\n")
   dmat<-apply(x$cdist,2,min,na.rm=TRUE)
   dmat<-rbind(dmat,apply(x$cdist,2,quantile,probs=0.25,names=FALSE,na.rm=TRUE))
   dmat<-rbind(dmat,apply(x$cdist,2,quantile,probs=0.5,names=FALSE,na.rm=TRUE))
   dmat<-rbind(dmat,apply(x$cdist,2,mean,na.rm=TRUE))
   dmat<-rbind(dmat,apply(x$cdist,2,quantile,probs=0.75,names=FALSE,na.rm=TRUE))
   dmat<-rbind(dmat,apply(x$cdist,2,max,na.rm=TRUE))
   colnames(dmat)<-as.vector(x$cnames)
   rownames(dmat)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
   print.table(dmat,digits=4)
   cat("\n")
}


#print.summary.netlm - Print method for summary.netlm
print.summary.netlm<-function(x,...){
   cat("\nOLS Network Model\n\n")
   cat("Coefficients:\n\n")
   cmat <- as.vector(format(as.numeric(x$coefficients)))
   cmat <- cbind(cmat, as.vector(format(x$pgreq)))
   cmat <- cbind(cmat, as.vector(format(x$pleeq)))
   colnames(cmat) <- c("Estimate", "Pr(>=b)", "Pr(<=b)")
   rownames(cmat)<- as.vector(x$names)
   print.table(cmat)
   #Goodness of fit measures
   mss<-if(attr(x$terms,"intercept"))
      sum((fitted(x)-mean(fitted(x)))^2)
   else
      sum(fitted(x)^2)
   rss<-sum(resid(x)^2)
   qn<-NROW(x$qr$qr)
   df.int<-if(attr(x$terms,"intercept")) 1
      else 0
   rdf<-qn-x$rank
   resvar<-rss/rdf
   fstatistic<-c(value=(mss/(x$rank-df.int))/resvar,numdf=x$rank-df.int,dendf=rdf)
   r.squared<-mss/(mss+rss)
   adj.r.squared<-1-(1-r.squared)*((qn-df.int)/rdf)
   sigma<-sqrt(resvar)
   cat("\nResidual standard error:",format(sigma,digits=4),"on",rdf,"degrees of freedom\n")
   cat("F-statistic:",formatC(fstatistic[1],digits=4),"on",fstatistic[2],"and",fstatistic[3],"degrees of freedom, \tp-value:",formatC(1-pf(fstatistic[1],fstatistic[2],fstatistic[3]),dig=4),"\n")
   cat("Multiple R^2:",format(r.squared,digits=4),"\n")
   cat("Adjusted R^2:",format(adj.r.squared,digits=4),"\n")
   #Test diagnostics
   cat("\n\nTest Diagnostics:\n\n")
   cat("\tNull Hypothesis:")
   if(x$nullhyp=="qap")
      cat(" QAP\n")
   else
      cat(" CUG\n")
   cat("\tReplications:",dim(x$dist)[1],"\n")
   cat("\tGoodness of Fit Distribution Summary:\n\n")
   gof<-cbind(x$sigma.dist,x$r.squared.dist,x$adj.r.squared.dist)
   dmat<-apply(gof,2,min,na.rm=TRUE)
   dmat<-rbind(dmat,apply(gof,2,quantile,probs=0.25,names=FALSE,na.rm=TRUE))
   dmat<-rbind(dmat,apply(gof,2,quantile,probs=0.5,names=FALSE,na.rm=TRUE))
   dmat<-rbind(dmat,apply(gof,2,mean,na.rm=TRUE))
   dmat<-rbind(dmat,apply(gof,2,quantile,probs=0.75,names=FALSE,na.rm=TRUE))
   dmat<-rbind(dmat,apply(gof,2,max,na.rm=TRUE))
   colnames(dmat)<-c("Sigma","R^2","Adj. R^2")
   rownames(dmat)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
   print.table(dmat,digits=4)
   cat("\n")
   cat("\tCoefficient Distribution Summary:\n\n")
   dmat<-apply(x$dist,2,min,na.rm=TRUE)
   dmat<-rbind(dmat,apply(x$dist,2,quantile,probs=0.25,names=FALSE,na.rm=TRUE))
   dmat<-rbind(dmat,apply(x$dist,2,quantile,probs=0.5,names=FALSE,na.rm=TRUE))
   dmat<-rbind(dmat,apply(x$dist,2,mean,na.rm=TRUE))
   dmat<-rbind(dmat,apply(x$dist,2,quantile,probs=0.75,names=FALSE,na.rm=TRUE))
   dmat<-rbind(dmat,apply(x$dist,2,max,na.rm=TRUE))
   colnames(dmat)<-as.vector(x$names)
   rownames(dmat)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
   print.table(dmat,digits=4)
   cat("\n")
}


#print.summary.netlogit - Print method for summary.netlogit
print.summary.netlogit<-function(x,...){
   cat("\nNetwork Logit Model\n\n")
   cat("Coefficients:\n\n")
   cmat <- as.vector(format(as.numeric(x$coefficients)))
   cmat <- cbind(cmat, as.vector(format(exp(as.numeric(x$coefficients)))))
   cmat <- cbind(cmat, as.vector(format(x$pgreq)))
   cmat <- cbind(cmat, as.vector(format(x$pleeq)))
   colnames(cmat) <- c("Estimate", "Exp(b)", "Pr(>=b)", "Pr(<=b)")
   rownames(cmat)<- as.vector(x$names)
   print.table(cmat)
   cat("\nGoodness of Fit Statistics:\n")
   cat("\nNull deviance (-2*Ln(L)):",x$null.deviance,"on",x$df.null,"degrees of freedom\n")
   cat("Residual deviance (-2*Ln(L)):",x$deviance,"on",x$df.residual,"degrees of freedom\n")
   cat("Chi-Squared test of fit improvement:\n\t",x$null.deviance-x$deviance,"on",x$df.null-x$df.residual,"degrees of freedom, p-value",1-pchisq(x$null.deviance-x$deviance,df=x$df.null-x$df.residual),"\n") 
   cat("AIC:",x$aic,"\tBIC:",x$deviance+log(x$df.null+1)*(x$df.null-x$df.residual),"\nPseudo-R^2 Measures:\n\t(Dn-Dr)/(Dn-Dr+dfn):",(x$null.deviance-x$deviance)/(x$null.deviance-x$deviance+x$df.null),"\n\t(Dn-Dr)/Dn:",1-x$deviance/x$null.deviance,"\n")
   cat("Contingency Table (predicted (rows) x actual (cols)):\n\n")
   print.table(x$ctable,print.gap=3)
   cat("\n\tTotal Fraction Correct:",(x$ctable[1,1]+x$ctable[2,2])/sum(x$ctable),"\n\tFraction 1s Correct:",x$ctable[2,2]/sum(x$ctable[2,]),"\n\tFraction 0s Correct:",x$ctable[1,1]/sum(x$ctable[1,]),"\n")
   cat("\nTest Diagnostics:\n\n")
   cat("\tNull Hypothesis:")
   if(x$nullhyp=="qap")
      cat(" QAP\n")
   else
      cat(" CUG\n")
   cat("\tReplications:",dim(x$dist)[1],"\n")
   cat("\tDistribution Summary:\n\n")
   dmat<-apply(x$dist,2,min,na.rm=TRUE)
   dmat<-rbind(dmat,apply(x$dist,2,quantile,probs=0.25,names=FALSE,na.rm=TRUE))
   dmat<-rbind(dmat,apply(x$dist,2,quantile,probs=0.5,names=FALSE,na.rm=TRUE))
   dmat<-rbind(dmat,apply(x$dist,2,mean,na.rm=TRUE))
   dmat<-rbind(dmat,apply(x$dist,2,quantile,probs=0.75,names=FALSE,na.rm=TRUE))
   dmat<-rbind(dmat,apply(x$dist,2,max,na.rm=TRUE))
   colnames(dmat)<-as.vector(x$names)
   rownames(dmat)<-c("Min","1stQ","Median","Mean","3rdQ","Max")
   print.table(dmat,digits=4)
   cat("\n")
}


#pstar - Perform an approximate p* analysis using the logistic regression 
#approximation.  Note that the result of this is returned as a GLM object, and 
#subsequent printing/summarizing/etc. should be treated accordingly.
pstar<-function(dat,effects=c("choice","mutuality","density","reciprocity","stransitivity","wtransitivity","stranstri","wtranstri","outdegree","indegree","betweenness","closeness","degcentralization","betcentralization","clocentralization","connectedness","hierarchy","lubness","efficiency"),attr=NULL,memb=NULL,diag=FALSE,mode="digraph"){
   #First, take care of various details
   n<-dim(dat)[1]
   m<-dim(dat)[2]
   o<-list()
   #Next, add NAs as needed
   d<-dat
   if(!diag)
      d<-diag.remove(d)
   if(mode=="graph")
      d<-upper.tri.remove(d)
   #Make sure that attr and memb are well-behaved
   if(!is.null(attr)){
      if(is.vector(attr))
         attr<-matrix(attr,ncol=1)
      if(is.null(colnames(attr)))
         colnames(attr)<-paste("Attribute",1:dim(attr)[2])
   }
   if(!is.null(memb)){
      if(is.vector(memb))
         memb<-matrix(memb,ncol=1)
      if(is.null(colnames(memb)))
         colnames(memb)<-paste("Membership",1:dim(memb)[2])
   }
   #Now, evaluate each specified effect given each possible perturbation
   tiedat<-vector()
   for(i in 1:n)
      for(j in 1:m)
         if(!is.na(d[i,j])){
            #Assess the effects
            td<-vector()
            if(!is.na(pmatch("choice",effects))){  #Compute a choice effect
               td<-c(td,1)  #Always constant
            }
            if(!is.na(pmatch("mutuality",effects))){  #Compute a mutuality effect
               td<-c(td,eval.edgeperturbation(d,i,j,"mutuality"))
            }
            if(!is.na(pmatch("density",effects))){  #Compute a density effect
               td<-c(td,eval.edgeperturbation(d,i,j,"gden",mode=mode,diag=diag))
            }
            if(!is.na(pmatch("reciprocity",effects))){  #Compute a reciprocity effect
               td<-c(td,eval.edgeperturbation(d,i,j,"grecip"))
            }
            if(!is.na(pmatch("stransitivity",effects))){  #Compute a strong transitivity effect
               td<-c(td,eval.edgeperturbation(d,i,j,"gtrans",mode=mode,diag=diag,measure="strong"))
            }
            if(!is.na(pmatch("wtransitivity",effects))){  #Compute a weak transitivity effect
               td<-c(td,eval.edgeperturbation(d,i,j,"gtrans",mode=mode,diag=diag,measure="weak"))
            }
            if(!is.na(pmatch("stranstri",effects))){  #Compute a strong trans census effect
               td<-c(td,eval.edgeperturbation(d,i,j,"gtrans",mode=mode,diag=diag,measure="strongcensus"))
            }
            if(!is.na(pmatch("wtranstri",effects))){  #Compute a weak trans census effect
               td<-c(td,eval.edgeperturbation(d,i,j,"gtrans",mode=mode,diag=diag,measure="weakcensus"))
            }
            if(!is.na(pmatch("outdegree",effects))){  #Compute outdegree effects
               td<-c(td,eval.edgeperturbation(d,i,j,"degree",cmode="outdegree",gmode=gmode,diag=diag))
            }
            if(!is.na(pmatch("indegree",effects))){  #Compute indegree effects
               td<-c(td,eval.edgeperturbation(d,i,j,"degree",cmode="indegree",gmode=gmode,diag=diag))
            }
            if(!is.na(pmatch("betweenness",effects))){  #Compute betweenness effects
               td<-c(td,eval.edgeperturbation(d,i,j,"betweenness",gmode=mode,diag=diag))
            }
            if(!is.na(pmatch("closeness",effects))){  #Compute closeness effects
               td<-c(td,eval.edgeperturbation(d,i,j,"closeness",gmode=mode,diag=diag))
            }
            if(!is.na(pmatch("degcentralization",effects))){  #Compute degree centralization effects
               td<-c(td,eval.edgeperturbation(d,i,j,"centralization","degree",mode=mode,diag=diag))
            }
            if(!is.na(pmatch("betcentralization",effects))){  #Compute betweenness centralization effects
               td<-c(td,eval.edgeperturbation(d,i,j,"centralization","betweenness",mode=mode,diag=diag))
            }
            if(!is.na(pmatch("clocentralization",effects))){  #Compute closeness centralization effects
               td<-c(td,eval.edgeperturbation(d,i,j,"centralization","closeness",mode=mode,diag=diag))
            }
            if(!is.na(pmatch("connectedness",effects))){  #Compute connectedness effects
               td<-c(td,eval.edgeperturbation(d,i,j,"connectedness"))
            }
            if(!is.na(pmatch("hierarchy",effects))){  #Compute hierarchy effects
               td<-c(td,eval.edgeperturbation(d,i,j,"hierarchy"))
            }
            if(!is.na(pmatch("lubness",effects))){  #Compute lubness effects
               td<-c(td,eval.edgeperturbation(d,i,j,"lubness"))
            }
            if(!is.na(pmatch("efficiency",effects))){  #Compute efficiency effects
               td<-c(td,eval.edgeperturbation(d,i,j,"efficiency",diag=diag))
            }
            #Add attribute differences, if needed
            if(!is.null(attr))
               td<-c(td,abs(attr[i,]-attr[j,]))
            #Add membership similarities, if needed
            if(!is.null(memb))
               td<-c(td,as.numeric(memb[i,]==memb[j,]))
            #Add this data to the aggregated tie data
            tiedat<-rbind(tiedat,c(d[i,j],td))
         }
   #Label the tie data matrix
   tiedat.lab<-"EdgeVal"
   if(!is.na(pmatch("choice",effects)))  #Label the choice effect
      tiedat.lab<-c(tiedat.lab,"Choice")
   if(!is.na(pmatch("mutuality",effects)))  #Label the mutuality effect
      tiedat.lab<-c(tiedat.lab,"Mutuality")
   if(!is.na(pmatch("density",effects)))  #Label the density effect
      tiedat.lab<-c(tiedat.lab,"Density")
   if(!is.na(pmatch("reciprocity",effects)))  #Label the reciprocity effect
      tiedat.lab<-c(tiedat.lab,"Reciprocity")
   if(!is.na(pmatch("stransitivity",effects)))  #Label the strans effect
      tiedat.lab<-c(tiedat.lab,"STransitivity")
   if(!is.na(pmatch("wtransitivity",effects)))  #Label the wtrans effect
      tiedat.lab<-c(tiedat.lab,"WTransitivity")
   if(!is.na(pmatch("stranstri",effects)))  #Label the stranstri effect
      tiedat.lab<-c(tiedat.lab,"STransTriads")
   if(!is.na(pmatch("wtranstri",effects)))  #Label the wtranstri effect
      tiedat.lab<-c(tiedat.lab,"WTransTriads")
   if(!is.na(pmatch("outdegree",effects)))  #Label the outdegree effect
      tiedat.lab<-c(tiedat.lab,paste("Outdegree",1:n,sep="."))
   if(!is.na(pmatch("indegree",effects)))  #Label the indegree effect
      tiedat.lab<-c(tiedat.lab,paste("Indegree",1:n,sep="."))
   if(!is.na(pmatch("betweenness",effects)))  #Label the betweenness effect
      tiedat.lab<-c(tiedat.lab,paste("Betweenness",1:n,sep="."))
   if(!is.na(pmatch("closeness",effects)))  #Label the closeness effect
      tiedat.lab<-c(tiedat.lab,paste("Closeness",1:n,sep="."))
   if(!is.na(pmatch("degcent",effects)))  #Label the degree centralization effect
      tiedat.lab<-c(tiedat.lab,"DegCentralization")
   if(!is.na(pmatch("betcent",effects)))  #Label the betweenness centralization effect
      tiedat.lab<-c(tiedat.lab,"BetCentralization")
   if(!is.na(pmatch("clocent",effects)))  #Label the closeness centralization effect
      tiedat.lab<-c(tiedat.lab,"CloCentralization")
   if(!is.na(pmatch("connectedness",effects)))  #Label the connectedness effect
      tiedat.lab<-c(tiedat.lab,"Connectedness")
   if(!is.na(pmatch("hierarchy",effects)))  #Label the hierarchy effect
      tiedat.lab<-c(tiedat.lab,"Hierarchy")
   if(!is.na(pmatch("lubness",effects)))  #Label the lubness effect
      tiedat.lab<-c(tiedat.lab,"LUBness")
   if(!is.na(pmatch("efficiency",effects)))  #Label the efficiency effect
      tiedat.lab<-c(tiedat.lab,"Efficiency")
   if(!is.null(attr))
      tiedat.lab<-c(tiedat.lab,colnames(attr))
   if(!is.null(memb))
      tiedat.lab<-c(tiedat.lab,colnames(memb))
   colnames(tiedat)<-tiedat.lab
   #Having had our fun, it's time to get serious.  Run a GLM on the resulting data.
   fmla<-as.formula(paste("EdgeVal ~ -1 + ",paste(colnames(tiedat)[2:dim(tiedat)[2]],collapse=" + ")))
   o<-glm(fmla,family="binomial",data=as.data.frame(tiedat))
   o$tiedata<-tiedat
   #Return the result
   o
}


#se.lnam - Standard error method for lnam
se.lnam<-function(object, ...){
   se<-vector()
   sen<-vector()
   if(!is.null(object$rho1.se)){
      se<-c(se,object$rho1.se)
      sen<-c(sen,"rho1")
   }
   if(!is.null(object$rho2.se)){
      se<-c(se,object$rho2.se)
      sen<-c(sen,"rho2")
   }
   if(!is.null(object$beta.se)){
      se<-c(se,object$beta.se)
      sen<-c(sen,names(object$beta.se))
   }
   names(se)<-sen
   se
}


#summary.bayes.factor - A fairly generic summary routine for bayes factors.  
#Clearly, this belongs in some other library than sna, but for the moment this 
#will have to do...
summary.bayes.factor<-function(object, ...){
   o<-object
   rownames(o$int.lik)<-o$model.names
   colnames(o$int.lik)<-o$model.names
   o$inv.bf<-1/o$int.lik
   for(i in 1:dim(o$int.lik)[1])
      o$inv.bf[i,i]<-o$int.lik[i,i]/sum(diag(o$int.lik))
   class(o)<-c("summary.bayes.factor","bayes.factor")
   o
}


#summary.bbnam - Summary method for bbnam
summary.bbnam<-function(object, ...){
   out<-object
   class(out)<-c("summary.bbnam",class(out))
   out
}


#summary.bbnam.actor - Summary method for bbnam.actor
summary.bbnam.actor<-function(object, ...){
   out<-object
   class(out)<-c("summary.bbnam.actor",class(out))
   out
}


#summary.bbnam.fixed - Summary method for bbnam.fixed
summary.bbnam.fixed<-function(object, ...){
   out<-object
   class(out)<-c("summary.bbnam.fixed",class(out))
   out
}


#summary.bbnam.pooled - Summary method for bbnam.pooled
summary.bbnam.pooled<-function(object, ...){
   out<-object
   class(out)<-c("summary.bbnam.pooled",class(out))
   out
}


#summary.lnam - Summary method lnam
summary.lnam<-function(object, ...){
   ans<-object
   class(ans)<-c("summary.lnam","lnam")
   ans
}


#summary.netcancor - Summary method for netcancor
summary.netcancor<-function(object, ...){
   out<-object
   class(out)<-c("summary.netcancor",class(out))
   out
}


#summary.netlm - Summary method for netlm
summary.netlm<-function(object, ...){
   out<-object
   class(out)<-c("summary.netlm",class(out))
   out
}


#summary.netlogit - Summary method for netlogit
summary.netlogit<-function(object, ...){
   out<-object
   class(out)<-c("summary.netlogit",class(out))
   out
}
