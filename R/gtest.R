######################################################################
#
# gtest.R
#
# copyright (c) 2004, Carter T. Butts <buttsc@uci.edu>
# Last Modified 11/25/04
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/sna package
#
# This file contains routines related to null hypothesis testing.
#
# Contents:
#  cugtest
#  plot.qaptest
#  print.cugtest
#  print.qaptest
#  print.summary.cugtest
#  print.summary.qaptest
#  qaptest
#  summary.cugtest
#  summary.qaptest
#
######################################################################


#cugtest - Generate, print, and plot CUG (conditional uniform graph) test 
#objects.
cugtest<-function(dat,FUN,reps=1000,gmode="digraph",cmode="density",diag=FALSE,g1=1,g2=2,...){
   out<-list()
   #First, find the test value for fun on dat
   fun<-match.fun(FUN)
   out$testval<-fun(dat,g1=g1,g2=g2,...)
   #Next, determine on what we are conditioning
   if(cmode=="density"){
      d<-c(gden(dat,g=g1,mode=gmode,diag=diag),gden(dat,g=g2,mode=gmode,diag=diag))
   }else{
      d<-c(0.5,0.5)
   }
   #Now, perform reps replications on random recreations of the data
   out$dist<-vector(mode="numeric",length=reps)
   for(i in 1:reps){
      if(cmode=="ties"){
         out$dist[i]<-fun(rgraph(dim(dat)[2],2,diag=diag,mode=gmode,tielist=dat),g1=1,g2=2,...)
      }else{
         out$dist[i]<-fun(rgraph(dim(dat)[2],2,tprob=d,diag=diag,mode=gmode),g1=1,g2=2,...)
      }
   }
   #Find p values
   out$pgreq<-mean(as.numeric(out$dist>=out$testval))
   out$pleeq<-mean(as.numeric(out$dist<=out$testval))
   class(out)<-c("cugtest","cug")
   out
}


#plot.cugtest - Plotting method for cugtest
plot.cugtest<-function(x,mode="density",...){
   if(mode=="density"){
      plot(density(x$dist),main="Estimated Density of CUG Replications",xlab="Test Statistic",...)
   }else{
      hist(x$dist,main="Histogram of CUG Replications",xlab="Test Statistic",...)
   }
   abline(v=x$testval,lty=2)
}


#plot.qaptest - Plotting method for qaptest
plot.qaptest<-function(x,mode="density",...){
   if(mode=="density"){
      plot(density(x$dist),main="Estimated Density of QAP Replications",xlab="Test Statistic",...)
   }else{
      hist(x$dist,main="Histogram of QAP Replications",xlab="Test Statistic",...)
   }
   abline(v=x$testval,lty=2)
}


#print.cugtest - Print method for cugtest
print.cugtest<-function(x,...){
   cat("\nCUG Test Results\n\n")
   cat("Estimated p-values:\n")
   cat("\tp(f(rnd) >= f(d)):",x$pgreq,"\n")
   cat("\tp(f(rnd) <= f(d)):",x$pleeq,"\n\n")      
}


#print.qaptest - Print method for qaptest
print.qaptest<-function(x,...){
   cat("\nQAP Test Results\n\n")
   cat("Estimated p-values:\n")
   cat("\tp(f(perm) >= f(d)):",x$pgreq,"\n")
   cat("\tp(f(perm) <= f(d)):",x$pleeq,"\n\n")      
}


#print.summary.cugtest - Print method for summary.cugtest
print.summary.cugtest<-function(x,...){
   cat("\nCUG Test Results\n\n")
   cat("Estimated p-values:\n")
   cat("\tp(f(rnd) >= f(d)):",x$pgreq,"\n")
   cat("\tp(f(rnd) <= f(d)):",x$pleeq,"\n")
   cat("\nTest Diagnostics:\n")
   cat("\tTest Value (f(d)):",x$testval,"\n")
   cat("\tReplications:",length(x$dist),"\n")
   cat("\tDistribution Summary:\n")
   cat("\t\tMin:\t",quantile(x$dist,probs=0,names=FALSE),"\n")
   cat("\t\t1stQ:\t",quantile(x$dist,probs=0.25,names=FALSE),"\n")
   cat("\t\tMed:\t",quantile(x$dist,probs=0.5,names=FALSE),"\n")
   cat("\t\tMean:\t",mean(x$dist),"\n")
   cat("\t\t3rdQ:\t",quantile(x$dist,probs=0.75,names=FALSE),"\n")
   cat("\t\tMax:\t",quantile(x$dist,probs=1,names=FALSE),"\n")
   cat("\n")
}


#print.summary.qaptest - Print method for summary.qaptest
print.summary.qaptest<-function(x,...){
   cat("\nQAP Test Results\n\n")
   cat("Estimated p-values:\n")
   cat("\tp(f(perm) >= f(d)):",x$pgreq,"\n")
   cat("\tp(f(perm) <= f(d)):",x$pleeq,"\n")
   cat("\nTest Diagnostics:\n")
   cat("\tTest Value (f(d)):",x$testval,"\n")
   cat("\tReplications:",length(x$dist),"\n")
   cat("\tDistribution Summary:\n")
   cat("\t\tMin:\t",quantile(x$dist,probs=0,names=FALSE),"\n")
   cat("\t\t1stQ:\t",quantile(x$dist,probs=0.25,names=FALSE),"\n")
   cat("\t\tMed:\t",quantile(x$dist,probs=0.5,names=FALSE),"\n")
   cat("\t\tMean:\t",mean(x$dist),"\n")
   cat("\t\t3rdQ:\t",quantile(x$dist,probs=0.75,names=FALSE),"\n")
   cat("\t\tMax:\t",quantile(x$dist,probs=1,names=FALSE),"\n")
   cat("\n")
}


#qaptest - Generate a QAP test object
qaptest<-function(dat,FUN,reps=1000,...){
   out<-list()
   #First, find the test value for fun on dat
   fun<-match.fun(FUN)
   out$testval<-fun(dat,...)
   #Now, perform reps replications on random permutations of the data
   out$dist<-vector(mode="numeric",length=reps)
   for(i in 1:reps){
      out$dist[i]<-fun(rmperm(dat),...)
   }
   #Find p values
   out$pgreq<-mean(as.numeric(out$dist>=out$testval))
   out$pleeq<-mean(as.numeric(out$dist<=out$testval))
   class(out)<-c("qaptest","qap")
   out
}


#summary.cugtest - Summary method for cugtest
summary.cugtest<-function(object, ...){
   out<-object
   class(out)<-c("summary.cugtest",class(out))
   out
}


#summary.qaptest - Summary method for qaptest
summary.qaptest<-function(object, ...){
   out<-object
   class(out)<-c("summary.qaptest",class(out))
   out
}
