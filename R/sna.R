#
# Various Useful Tools for Network Analysis in S
#
# By Carter Butts, buttsc@uci.edu
#
# Current Version: 0.44
#
# Last updated 5/24/04
#
#Contents:
#
#%c% - Composition of two adjacancy matrices
#addisolates - Add isolates to a graph set
#bbnam - Draw from Butts' (Hierarchical) Bayesian Network Accuracy Model
#bbnam.bf - Bayes factors for Butts' (Hierarchical) Bayesian Network Accuracy
#   Model
#bbnam.probtie - Internal routine for finding tie likelihoods under bbnam
#bbnam.jntlik - Internal routine for computing joint likelihoods under bbnam
#bbnam.jntlik.slice - Internal routine for joint likelihood computation under 
#   bbnam, by data slice
#betweenness - Find the betweenness centralities of network positions
#blockmodel - Generate blockmodels based on partitions of network positions
#blockmodel.expand - Generate a graph from a given blockmodel using particular 
#   expansion rules
#bonpow - Find the Bonacich power centrality scores of network positions
#centralgraph - Find the central graph of a graph stack
#centralization - Generic centralization routine (uses centrality routines)
#closeness - Find the closeness centralities of network positions
#component.dist - Find the distribution of (maximal) component sizes within a 
#   graph
#components - Find the number of (maximal) components within a given graph
#connectedness - Find the Krackhardt connectedness of a graph or graph stack
#consensus - Find a consensus structure, using one of several algorithms
#cugtest - Generic Conditional Uniform Graph (CUG) test for bigraphic statistics
#degree - Computes the degree centralities of network positions
#diag.remove - NAs the diagonals of graphs in a stack
#dyad.census - Compute the Holland and Leinhardt MAN dyad census for a graph or 
#   graph stack
#efficiency - Find the Krackhardt efficiency of a graph or graph stack
#equiv.clust - Find clusters of positions based on an equivalence relation
#eval.edgeperturbation - Evaluate a function on a given graph with and without a
#   given edge
#evcent - Find the eigenvector centralities of network positions
#event2dichot - Convert observed event matrices to dichotomous matrices
#gapply - Apply functions over vertex neighborhoods
#gclust.boxstats - Plot statistics associated with clusters
#gclust.centralgraph - Get central graphs associated with clusters
#gcor - Computes the correlation between graphs
#gcov - Computes the covariance between graphs
#gden - Computes the density of one or more graphs
#gdist.plotdiff - Plot differences in graph-level statistics against inter-graph
#   distances
#gdist.plotstats - Plot statistics associated with graphs against (projected) 
#   inter-graph distances
#geodist - Finds geodesic distances and shortest paths within a graph
#gliop - Return a binary operation on GLI values computed on two graphs (for 
#   test routines)
#gplot - Plots graphs using eigenvector projection, classical MDS, or custom 
#   coords
#graphcent - Find the graph centralities of network positions
#grecip - Computes the reciprocity of one or more graphs
#gtrans - Computes the transitivity of one or more graphs
#gscor - Computes the structural correlation between graphs
#gscov - Computes the structural covariance between graphs
#gvectorize - Vectorization of adjacency matrices
#hdist - Computes the Hamming distance between two labeled graphs
#hierarchy - Find the hierarchy score of a graph or graph stack
#infocent - Find the information centrality scores of network positions
#interval.graph - Construct one or more interval graphs (and exchangeability 
#   vectors) from a set of spells
#is.isolate - Determines whether a particular vertex is isolated
#isolates - Returns a list of isolates
#lab.optimize - Optimize a bivariate graph statistic over a set of labelings
#lnam - Fits a linear network autocorrelation model
#lower.tri.remove - NAs the lower triangles of graphs in graph stacks
#lubness - Find Krackhardt's Least Upper Boundedness of a graph or graph stack
#make.stochastic - Make a graph stack row, column, or row-column stochastic
#mutuality - Find the number of mutual (i.e., reciprocated) edges in a graph
#netcancor - Canonical correlations for network variables (requires mva)
#netlm - Fits a network OLS regression model (w/CUG or QAP tests)
#netlogit - Fits a network logit model (w/CUG or QAP tests)
#npostpred - Take posterior predictive draws for functions of graphs
#nties - Find the number of ties in a given graph
#numperm - Get the nth permutation vector, using a periodic numbering scheme
#plot.bbnam - Plotting for bbnam objects
#plot.blockmodel - Plotting for blockmodel objects
#plot.cugtest - Plotting for cugtest objects
#plot.equiv.clust - Plotting for equivalence clustering objects
#plot.lnam - Plotting for lnam objects
#plot.sociomatrix - Plotting of arbitrary matrices
#plot.qaptest - Plotting for qaptest objects
#potscalered.mcmc - Computes Gelman et al.'s potential scale reduction statistic
#   for scalar estimands
#prestige - Find actor prestige scores from one of several measures
#print.bayes.factor - Printing for Bayes factor objects
#print.bbnam - Printing for bbnam objects
#print.blockmodel - Printing for blockmodel objects
#print.cugtest - Printing for cugtest objects
#print.lnam - Printing for lnam objects
#print.netcancor - Printing for netcancor objects
#print.netlm - Printing for netlm objects
#print.netlogit - Printing for netlogit objects
#print.qaptest - Printing for qaptest objects
#print.summary.bayes.factor - Printing for Bayes factor summary objects 
#print.summary.bbnam - Printing for bbnam summary objects
#print.summary.blockmodel - Printing for blockmodel summary objects
#print.summary.cugtest - Printing for cugtest summary objects
#print.summary.lnam - Printing for lnam summary objects
#print.summary.netcancor - Printing for netcancor summary objects
#print.summary.netlm - Printing for netlm summary objects
#print.summary.netlogit - Printing for netlogit summary objects
#print.summary.qaptest - Printing for qaptest summary objects
#pstar - Fit a p* model using the logistic regression approximation.
#qaptest - Generic QAP test for bigraphic statistics
#reachability - Find the reachability matrix of a graph
#read.nos - Import data files in Neo-OrgStat format
#rgraph - Draws Bernoulli graphs
#rmperm - Randomly permutes the rows and columns of a graph stack
#rperm - Draw a random permutation vector with exchangability constraints
#sdmat - Estimate the matrix of structural distances among unlabeled graphs
#sedist - Find distances between positions based on structural equivalence
#sr2css - Convert a row-wise self-report matrix to a CSS matrix with missing 
#   observations
#stackcount -Find the number of matrices in a graph stack (matrix or array data 
#   acceptable)
#stresscent - Find the stress centralities of network positions
#structdist - Estimate the structural distance between unlabeled graphs
#summary.bayes.factor - Detailed printing for Bayes factor objects
#summary.bbnam - Detailed printing for bbnam objects
#summary.blockmodel - Detailed printing for blockmodel objects
#summary.cugtest - Detailed printing for qaptest objects
#summary.lnam - Detailed printing for lnam objects
#summary.netcancor - Detailed printing for netcancor objects
#summary.netlm - Detailed printing for netlm objects
#summary.netlogit - Detailed printing for netlogit objects
#summary.qaptest - Detailed printing for qaptest objects
#symmetrize - Symmetrize a graph or graph stack
#triad.census - Conduct a Davis and Leinhardt triad census for a graph or graph 
#   stack
#triad.classify - Return the Davis and Leinhardt classification of a given triad
#upper.tri.remove - NAs the upper triangles of graphs in graph stacks
#
#NOTES:
#
#License: This software is made available under the terms of the GNU Public 
#License (GPL), version 2.0 or later.  Please visit the Free Software Foundation
#(www.fsf.org) for more details.
#
#Lack of support: It must be noted that this code is provided AS IS, without 
#support or warranty (express or implied).  The author will generally try to fix
#bugs where possible, but the end user is generally on his/her own with respect 
#to this software; the author apologizes for this fact, but would like to point 
#out that he has numerous professional and other obligations which prevent his
#working full-time on this project.
#
#Implementation: This code has been written for the R implementation of the S 
#language, and is currently thought to work with version 1.2 and higher.  
#Systematic testing, however, has been limited: caveat emptor. 
#
#A general note on data frames: Graphs are assumed to be n x n matrices (for 
#single structures) or m x n x n arrays (for graph stacks).  By default, all 
#network tools will expect this sort of data; this winds up being important for 
#allowing things like general network hypothesis testing routines.
#
#
#CHANGELOG:
#
#v0.44 - New Functions, New Features, Changes, and Bug Fixes
#   New Functions:
#      %c% - Composition of two adjacency matrices
#      gapply - Apply functions over vertex neighborhoods.
#      gplot.layout.* - Layout functions for gplot
#      lnam - Fit a linear network autocorrelation model
#      plot.lnam - Plotting for lnam objects
#      print.lnam - Printing for lnam objects
#      print.summary.lnam - Printing for lnam summary objects
#      summary.lnam - Detailed printing for lnam objects
#   New Features:
#      consensus now supports union/intersection LAS and column/row raw report
#        methods (useful for replicating "classic" CSS work)
#      gplot sports a wide range of new features, including:
#        An interactive mode, which allows for manual repositioning
#          of vertices.  (It's not pretty, but it works.)  
#        Silent return of the x,y coordinates for all vertices.  This is 
#          useful for adding features to an existing plot, or for saving
#          a given layout for later reuse.
#        Adjacency matrices can now be used as parameters to set edge widths,
#          line types, and colors [submitted by Alex Montgomery]
#        Overplotting support [submitted by Alex Montgomery]
#        Support for loops (i.e., self-ties) [submitted by Alex Montgomery]
#        Curved edges, with adjustable curvature [submitted by Alex Montgomery]
#        Support for user-supplied layout methods
#        Expansion of existing layout methods, including new arguments
#   Changes:
#      triad.census now automatically removes the diagonal before calculating
#        the triad census; a note on valued/unvalued data has also been added
#        to the man page [suggested by Skye Bender-deMoll and Dan McFarland]
#      Vertex layouts for gplot are now generated externally, via the
#        gplot.layout.* functions.  Arbitrary parameters may be passed to
#        the layout functions via an argument list; further, since layout
#        functions are identified with match.fun(), user-added functions
#        can also be employed.  The default layout method is now segeo (spring
#        embedder results were sometimes very good, but too uneven).
#   Bug Fixes:
#      degree was reporting incorrect tmaxdev values in some cases
#      triad.classify switched 111D and 111U [submitted by Dan McFarland and
#        Skye Bender-deMoll]
#
#v0.43 - Minor Changes, Updates, and New Features
#   Changes:
#      Contact URL has been updated
#   Updates:
#      In keeping with the new rigorousness regarding data.frame structures in
#        1.8.0, many data.frames have been changed to lists.  This should be
#        transparent to the end user, but will avoid the generation of errors
#        under the new system.
#      Removed references to (deprecated) plot.hclust
#   New Features:
#      gplot now supports spring embedding.  Providing unsupported layout modes
#        now produces an error, rather than undefined behavior (as before).  
#        Options have also been added for suppressing the printing of axes,
#        and for placing opaque boxes behind vertex labels.
#v0.42 - Minor Changes
#   Changes:
#      Author contact information has been updated
#      plot.matrix is now plot.sociomatrix, in order to ensure compatibility
#         with the new method standards; all code should be updated to reflect
#	 this fact
#v0.41 - Updates, New Features, and Bug Fixes
#   Updates:
#      Deprecated function .Alias removed (was used in netlm, netlogit)
#      Changed keyword "network" to "math" in all help pages [as requested
#         by the Keepers of R]
#      Various internal changes to plot/print/summary methods, in order
#         to maintain consistency with their generic equivalents; these
#         will (hopefully) have no visible effect
#   New Features:
#      component.dist now supports weak, unilateral, strong, and recursive
#         component definitions
#   Bug Fixes:
#      component.dist was calculating recursively connected components for
#         connected="strong", instead of strongly connected components
#      pstar was dumping (internal) edge perturbation data to the screen,
#         which was harmless but very annoying; names for pstar coefficients
#         were not being recognized by glm
#
#v0.4 - New Features, Changes, and Fixes
#   New Functions:
#      connectedness - Find the Krackhardt connectedness of a graph or graph 
#         stack
#      dyad.census - Compute the Holland and Leinhardt MAN dyad census for a 
#         graph or graph stack
#      efficiency - Find the Krackhardt efficiency of a graph or graph stack
#      hierarchy - Find the hierarchy score of a graph or graph stack.
#      infocent - Find the information centrality scores of network positions
#         [submitted by David Barron]
#      lubness - Find Krackhardt's Least Upper Boundedness of a graph or graph 
#         stack
#      reachability - Find the reachability matrix of a graph.
#      triad.census - Conduct a Davis and Leinhardt triad census for a graph or 
#         graph stack
#      triad.classify - Return the Davis and Leinhardt classification of a given
#         triad
#   New Features:
#      gplot now adjusts line width for valued graphs, via the edge.lwd 
#         parameter, and allows users to set the vertex point type (using
#         vertex.pch.  gmode=="graph" now sets usearrows<-FALSE, and new
#         gmode "twomode" automagically plots two-mode data.  New display modes
#         have also been added: geodist (MDS of proximities); adj (MDS with 
#         adjacency as similarity); and seham (MDS of SE dist using Hamming 
#         method).
#      grecip now supports a "measure" option, which allows the user to choose
#         between "dyadic" reciprocity (aRb iff bRa) and "edgewise" reciprocity
#         (aRb if bRa).  The old version (and default) is the "dyadic" option,
#         which (as the above implies) takes null dyads into account; the 
#         "edgewise" definition omits null dyads from the calculation.
#      gtrans now supports a "measure" option, which allows the user to choose
#         between "weak" transitivity (aRc if aRb & bRc) and "strong" 
#         transitivity (aRc iff aRb & bRc).  The old version was strong-only, 
#         but much of the field prefers the weak version.  Each of these options
#         has a respective census variant ("weakcensus", "strongcensus") which 
#         returns a count of legal triads rather than a rate.
#      pstar now supports separate options for strong/weak transitivity scores,
#         and for strong/weak transitive triad counts.
#   Bug Fixes:
#      Labeling optimizers were not pre-sorting to guard against
#         mixed-up exchange lists (these are now checked, too).
#      Various centrality measures were not returning the correct
#         absolute maximum deviation with tmaxdev==TRUE.
#      gtrans was ignoring settings and counting diagonal entries [submitted by
#         David Barron]
#      pstar behaved badly when external predictors were unnamed [submitted by
#         Gindo Tampubolon]
#   Changes:
#      Comparable labelings are now _enforced_ where applicable.
#      lab.optimize.gumbel now just tries to scare you off, rather
#         than refusing you altogether.
#      Path-based indices (betweenness, closeness, stresscent,
#         graphcent, etc.) now automatically assume 
#         cmode=="undirected" whenever gmode="graph".
#      The default mode for gtrans is now "weak", to match the usage of W&F
#
#v0.3 - New Features, Changes, and Fixes
#   General:
#      All standard functions are now documented
#      R package format is now supported
#   New Functions:
#      component.dist - Find the distribution of (maximal) component sizes 
#         within a graph
#      components - Find the number of (maximal) components within a given graph
#      eval.edgeperturbation - Evaluate a function on a given graph with and 
#         without a given edge, returning the difference between the results in 
#         each case.
#      interval.graph - Construct one or more interval graphs (and 
#         exchangeability vectors) from a set of spells
#      mutuality - Find the number of mutual (i.e., reciprocated) edges in a 
#         graph
#      gscor - Computes the structural correlation between graphs
#      gscov - Computes the structural covariance between graphs
#      gtrans - Compute the transitivity of an input graph or graph stack.
#      lab.optimize - Optimize a bivariate graph statistic over a set of 
#         labelings
#      pstar - Fits a p* model using the logistic regression approximation 
#      read.nos - Reads input files in Neo-OrgStat format
#      rperm - Draw a random permutation vector with exchangability constraints
#   Features and Modifications:
#      diag.remove, upper.tri.remove, and lower.tri.remove now allow replacement
#         with any value
#      gplot now provides a slew of new parameters to change color, size, etc. 
#         of vertices, edges, and labels
#      gscor and gscov now delegate to lab.optimize
#      gscor, gscov, structdist now support exchange lists (via lab.optimize)
#
#v0.2 - New Features and Some Bug Fixes
#   New Functions:
#      blockmodel - Generate blockmodels based on partitions of network 
#         positions
#      blockmodel.expand - Generate a graph from a given blockmodel using 
#         particular expansion rules
#      bonpow - Find the Bonacich power centrality scores of network positions
#      equiv.clust - Find clusters of positions based on an equivalence relation
#      evcent - Find the eigenvector centralities of network positions
#      gdist.plotdiff - Plot differences in graph-level statistics against 
#         inter-graph distances
#      gdist.plotstats - Plot statistics associated with graphs against 
#         (projected) inter-graph distances
#      make.stochastic - Make a graph stack row, column, or row-column 
#         stochastic
#      plot.blockmodel - Plotting for blockmodel objects
#      plot.equiv.clust - Plotting for equivalence clustering objects
#      prestige - Find actor prestige scores from one of several measures
#      print.blockmodel - Printing for blockmodel objects
#      print.summary.blockmodel - Printing for blockmodel summary objects
#      sedist - Find distances between positions based on structural equivalence
#      stackcount -Find the number of matrices in a graph stack (matrix or array
#         data acceptable)
#      summary.blockmodel - Detailed printing for blockmodel objects
#   Features and Modifications:
#      All centrality routines can now be rescaled (division by maximum realized
#         value); default is always FALSE   
#   Bug Fixes:
#      Various centrality routines once again return values for the selected 
#         graph and vertices
#      
#

#%c% - The composition operator for graph adjacency matrices
"%c%"<-function(x,y)
  round((x%*%y)>0)


#make.stochastic - Make a graph stack row, column, or row-column stochastic

make.stochastic<-function(dat,mode="rowcol",tol=0.005,maxiter=prod(dim(dat))*100,anneal.decay=0.01,errpow=1){
   #Organize the data
   m<-stackcount(dat)
   if(m==1){
      n<-dim(dat)[1]
      o<-dim(dat)[2]
      d<-array(dim=c(m,n,o))
      d[1,,]<-dat
   }else{
      n<-dim(dat)[2]
      o<-dim(dat)[3]
      d<-dat
   }
   #Stochasticize
   if(mode=="row"){
      for(i in 1:m)
         d[i,,]<-d[i,,]/t(sapply(apply(d[i,,],1,sum),rep,o))
   }else if(mode=="col"){
      for(i in 1:m)
         d[i,,]<-d[i,,]/sapply(apply(d[i,,],2,sum),rep,n)
   }else if(mode=="rowcol"){
      for(i in 1:m){
         f<-d[i,,]/t(sapply(apply(d[i,,],1,sum),rep,o))   #Seed with the row-stochastic form
         f<-f/sapply(apply(f,2,sum),rep,n)   #Col-stochasticize for good measure (sometimes this works)
         edgelist<-cbind(rep(1:n,o),rep(1:o,rep(n,o)))
         edgelist<-edgelist[d[i,,][edgelist]>0,]   #Skip edges which are forced to be zero-valued
         err<-sum(abs(apply(f,2,sum)-rep(1,o))^errpow,abs(apply(f,1,sum)-rep(1,n))^errpow)
         iter<-0
         while((err>(n+o)*tol)&(iter<maxiter)){  #Right now, use an annealer to find an approximate solution
            edge<-sample(1:dim(edgelist)[1],1)
            x<-edgelist[edge,1]
            y<-edgelist[edge,2]
            draw<-max(0,min(rnorm(1,f[x,y],d[i,x,y]/10),d[i,x,y]))
            nerr<-err-abs(sum(f[x,])-1)^errpow-abs(sum(f[,y])-1)^errpow+abs(sum(f[x,][-y])+draw-1)^errpow+abs(sum(f[,y][-x])+draw-1)^errpow
            if((nerr<err)|(runif(1,0,1)<exp(-anneal.decay*iter))){
               f[x,y]<-draw
               err<-nerr
            }
            iter<-iter+1
         }
         d[i,,]<-f
         if(err>(n+o)*tol)
            warning(paste("Annealer unable to reduce total error below apx",round(err,digits=7),"in matrix",i,". Hope that's OK....\n"))
      }
   }else if(mode=="total"){
         for(i in 1:m)
            d[i,,]<-d[i,,]/sum(d[i,,])
   }
   #Patch NaN values
   d[is.nan(d)]<-0
   #Return the output
   if(m==1)
      d[1,,]
   else
      d
}


#stackcount -How many matrices in a given stack?

stackcount<-function(d){
   if(length(dim(d))>2)
      dim(d)[1]
   else
      1
}


#gvectorize - Vectorization of adjacency matrices

gvectorize<-function(mats,mode="digraph",diag=FALSE,censor.as.na=TRUE){
   #Build the input data structures
   if(length(dim(mats))>2){
      m<-dim(mats)[1]
      n<-dim(mats)[2]
      n<-dim(mats)[3]
      d<-mats
   }else{
      m<-1
      n<-dim(mats)[1]
      o<-dim(mats)[2]
      d<-array(dim=c(1,n,o))
      d[1,,]<-mats
   }
   #If using NA censoring, turn unused parts of the matrices to NAs and vectorize
   if(censor.as.na){
      if(mode=="graph")
         d<-upper.tri.remove(d)
      if(!diag)
         d<-diag.remove(d)
      out<-apply(d,1,as.vector)
   }else{   #Otherwise, vectorize only the useful parts
      mask<-apply(d,1,lower.tri,diag=diag)
      out<-apply(d,1,as.vector)
      if(m==1)
         out<-out[mask]
      else
         out<-matrix(out[mask],ncol=m)      
   }
   out
}


#gclust.centralgraph - Get central graphs associated with clusters

gclust.centralgraph<-function(h,k,mats,...){
   #h must be an hclust object, k the number of groups, and mats the array of graphs
   out<-array(dim=c(k,dim(mats)[2],dim(mats)[2]))
   gmat<-matrix(nrow=dim(mats)[1],ncol=2)
   gmat[,1]<-c(1:dim(mats)[1])
   gmat[,2]<-cutree(h,k=k)
   for(i in 1:k)
      out[i,,]<-centralgraph(mats[gmat[gmat[,2]==i,1],,],...)
   out
}


#gclust.boxstats - Plot statistics associated with clusters

gclust.boxstats<-function(h,k,meas,...){
   #h must be an hclust object, k the number of groups, and meas the group measurement vector
   out<-matrix(nrow=length(meas),ncol=k)
   gmat<-matrix(nrow=length(meas),ncol=2)
   gmat[,1]<-c(1:length(meas))
   gmat[,2]<-cutree(h,k=k)
   for(i in 1:k){
      out[1:length(meas[gmat[gmat[,2]==i,1]]),i]<-meas[gmat[gmat[,2]==i,1]]
   }
   boxplot(data.frame(out),...)
}


#gdist.plotdiff - Plot differences in graph-level statistics against inter-graph distances

gdist.plotdiff<-function(d,meas,method="manhattan",jitter=TRUE,xlab="Inter-Graph Distance",ylab="Measure Distance",lm.line=FALSE,...){
   #Note that d must be a matrix of distances, and meas must be a matrix of measures
   #Create a matrix of differences in graph-level statistics
   require(mva)
   md<-dist(meas,method=method)
   #Vectorize
   dv<-as.vector(as.dist(d))
   mdv<-as.vector(md)
   if(jitter){
      dv<-jitter(dv)
      mdv<-jitter(mdv)
   }
   #Plot the results
   plot(dv,mdv,xlab=xlab,ylab=ylab,...)   
   #If needed, add a line fit
   abline(lm(mdv~dv),col="red")
}


#gdist.plotstats - Plot statistics associated with graphs against (projected) inter-graph distances

gdist.plotstats<-function(d,meas,siz.lim=c(0,0.15),rescale="quantile",display.scale="radius",display.type="circleray",cex=0.5,pch=1,labels=NULL,pos=1,labels.cex=1,legend=NULL,legend.xy=NULL,legend.cex=1,...){
   #d must be a matrix of distances (e.g., from structdist or hdist) and meas a matrix or vector of measures
   #Perform an MDS on the distances
   require(mva)
   xy<-cmdscale(as.dist(d))
   n<-dim(xy)[1]
   #Adjust and rescale the measure(s) prior to display
   if(is.null(dim(meas)))
      m<-matrix(meas,ncol=1)
   else
      m<-meas
   nm<-dim(m)[2]
   if(rescale=="quantile"){   #Rescale by (empirical) quantiles (ordinal rescaling)
      m<-apply(m,2,order)
      m<-sweep(m,2,apply(m,2,min))
      m<-sweep(m,2,apply(m,2,max),"/")
   }else if(rescale=="affine"){   #Rescale the measure(s) by affine transformation to the [0,1] interval
      m<-sweep(m,2,apply(m,2,min))
      m<-sweep(m,2,apply(m,2,max),"/")
   }else if(rescale=="normalize"){   #Rescale the measure(s) by their maximum values
      m<-sweep(m,2,apply(m,2,max),"/")
   }
   #Determine how large our drawn symbols are to be
   if(display.scale=="area")  #If we're using area scaling, we need to take square roots
      m<-sqrt(m)
   msize<-m*siz.lim[2]+siz.lim[1]   #Now, express the scaled measures as fractions of the plotting range
   pwid<-max(xy)-min(xy)  #Grab width of plotting range
   msize<-msize*pwid          #Adjust the msize matrix.
   #Plot the data
   plot(xy,xlab=expression(lambda[1]),ylab=expression(lambda[2]),cex=cex,pch=pch,xlim=c(min(xy),max(xy)),ylim=c(min(xy),max(xy)),...)   #Plot the graphs' MDS positions
   if(display.type=="poly"){  #Plot measures using polygons
      for(i in 1:nm){
         for(j in 1:n){
            x<-xy[j,1]+sin(2*pi*((0:nm)/nm))*msize[j,i]
            y<-xy[j,2]+cos(2*pi*((0:nm)/nm))*msize[j,i]
            lines(x,y,col=i)
         }
      }
   }else if(display.type=="circle"){  #Plot measures using circles
      for(i in 1:nm){
         for(j in 1:n){
            x<-xy[j,1]+sin(2*pi*((0:500)/500))*msize[j,i]
            y<-xy[j,2]+cos(2*pi*((0:500)/500))*msize[j,i]
            lines(x,y,col=i)
         }
      }
   }else if(display.type=="ray"){  #Plot measures using rays
      for(i in 1:nm){
         for(j in 1:n){
            lines(c(xy[j,1],xy[j,1]+sin(2*pi*((i-1)/nm))*msize[j,i]),c(xy[j,2],xy[j,2]+cos(2*pi*((i-1)/nm))*msize[j,i]),col=i)
         }
      }      
   }else if(display.type=="polyray"){  #Plot measures using polys and rays
      for(i in 1:nm){
         for(j in 1:n){
            x<-xy[j,1]+sin(2*pi*((0:nm)/nm))*msize[j,i]
            y<-xy[j,2]+cos(2*pi*((0:nm)/nm))*msize[j,i]
            lines(x,y,col=i)
            lines(c(xy[j,1],xy[j,1]+sin(2*pi*((i-1)/nm))*msize[j,i]),c(xy[j,2],xy[j,2]+cos(2*pi*((i-1)/nm))*msize[j,i]),col=i)
         }
      }      
   }else if(display.type=="circleray"){  #Plot measures using circles and rays
      for(i in 1:nm){
         for(j in 1:n){
            x<-xy[j,1]+sin(2*pi*((0:500)/500))*msize[j,i]
            y<-xy[j,2]+cos(2*pi*((0:500)/500))*msize[j,i]
            lines(x,y,col=i)
            lines(c(xy[j,1],xy[j,1]+sin(2*pi*((i-1)/nm))*msize[j,i]),c(xy[j,2],xy[j,2]+cos(2*pi*((i-1)/nm))*msize[j,i]),col=i)
         }
      }      
   }
   #Add labels, if needed
   if(!is.null(labels))
      text(xy[,1],xy[,2],labels,pos=pos,cex=labels.cex)
   #Add legend?
   if(!is.null(legend)){
      if(is.null(legend.xy))
         legend.xy<-c(min(xy),max(xy))
      legend(legend.xy[1],legend.xy[2],legend=legend,fill=1:nm,cex=legend.cex)
   }
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


#addisolates - Add isolates to one or more graphs

addisolates<-function(dat,n){
   if(length(dim(dat))>2){
      d<-array(dim=c(dim(dat)[1],dim(dat)[2]+n,dim(dat)[3]+n))
      d[,,]<-0
      for(i in 1:dim(dat)[1])
         d[i,1:dim(dat)[2],1:dim(dat)[2]]<-dat[i,,]
   }
   else{
      d<-matrix(nrow=dim(dat)[1]+n,ncol=dim(dat)[2]+n)
      d[,]<-0
      d[1:dim(dat)[2],1:dim(dat)[2]]<-dat
   }   
   d
}


#diag.remove - NA the diagonals of adjacency matrices in a graph stack

diag.remove<-function(dat,remove.val=NA){
   if(length(dim(dat))>2){
      d<-dat
      for(i in 1:dim(dat)[1])
         diag(d[i,,])<-remove.val
   }
   else{
      d<-dat
      diag(d)<-remove.val
   }   
   d
}


#upper.tri.remove - NA the upper triangles of adjacency matrices in a graph stack

upper.tri.remove<-function(dat,remove.val=NA){
   if(length(dim(dat))>2){
      d<-dat
      for(i in 1:dim(dat)[1]){
         temp<-d[i,,]
         temp[upper.tri(temp,diag=FALSE)]<-remove.val
         d[i,,]<-temp
      }
   }
   else{
      d<-dat
      d[upper.tri(d,diag=FALSE)]<-remove.val
   }   
   d
}


#lower.tri.remove - NA the lower triangles of adjacency matrices in a graph stack

lower.tri.remove<-function(dat,remove.val=NA){
   if(length(dim(dat))>2){
      d<-dat
      for(i in 1:dim(dat)[1]){
         temp<-d[i,,]
         temp[lower.tri(temp,diag=FALSE)]<-remove.val
         d[i,,]<-temp
      }
   }
   else{
      d<-dat
      d[lower.tri(d,diag=FALSE)]<-remove.val
   }   
   d
}


#nties - Find the number of ties in a given graph or stack

nties<- function(dat,mode="digraph",diag=FALSE){
   #Did someone send us a stack?
   if(length(dim(dat))>2)
      shiftit<-1
   else
      shiftit<-0
   #Get size params
   n<-dim(dat)[1+shiftit]
   m<-dim(dat)[2+shiftit]
   #Sanity check for hypergraphs
   if(mode=="hgraph")
      diag<-TRUE
   #Get the initial count
   count<-switch(mode,
      digraph = n*n,
      graph = (n*n-n)/2+n,
      hgraph = n*m
   )
   #Modify for diag, if needed
   if(!diag)
      count<-count-n
   #Return the needed info
   if(shiftit)
      rep(count,dim(dat)[1])
   else
      count                
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


#rmperm - Randomly permutes the rows and columns of an input matrix.

rmperm<-function(m){
   if(length(dim(m))==2){
      #Only a single matrix is included
      o<-sample(1:dim(m)[1])
      p<-matrix(data=m[o,o],nrow=dim(m)[1],ncol=dim(m)[2])
   }else{
      #Here, we assume a stack of matrices
      p<-array(dim=c(dim(m)[1],dim(m)[2],dim(m)[3]))
      for(i in 1:dim(m)[1]){
         o<-sample(1:dim(m)[2])
         p[i,,]<-array(m[i,o,o])
      }
   }
   p
}


#qaptest - Generate, print, and plot QAP test objects.

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

summary.qaptest<-function(object, ...){
   out<-object
   class(out)<-c("summary.qaptest",class(out))
   out
}

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

print.qaptest<-function(x,...){
   cat("\nQAP Test Results\n\n")
   cat("Estimated p-values:\n")
   cat("\tp(f(perm) >= f(d)):",x$pgreq,"\n")
   cat("\tp(f(perm) <= f(d)):",x$pleeq,"\n\n")      
}

plot.qaptest<-function(x,mode="density",...){
   if(mode=="density"){
      plot(density(x$dist),main="Estimated Density of QAP Replications",xlab="Test Statistic",...)
   }else{
      hist(x$dist,main="Histogram of QAP Replications",xlab="Test Statistic",...)
   }
   abline(v=x$testval,lty=2)
}


#lab.optimize - Optimize a function over the accessible permutation groups of two or more graphs.  This routine is a front end for various method-specific functions, and is in turn intended to be called from structural distance/covariance routines and the like.  The methods supported at this time include "exhaustive" (exhaustive search - I hope these are _small_ graphs!), "mc" (simple monte carlo search), " 

lab.optimize<-function(d1,d2,FUN,exchange.list=0,seek="min",opt.method=c("anneal","exhaustive","mc","hillclimb","gumbel"),...){
   meth<-match.arg(opt.method)
   if(meth=="anneal")
      lab.optimize.anneal(d1,d2,FUN,exchange.list,seek,...)
   else if(meth=="exhaustive")
      lab.optimize.exhaustive(d1,d2,FUN,exchange.list,seek,...)
   else if(meth=="mc")
      lab.optimize.mc(d1,d2,FUN,exchange.list,seek,...)
   else if(meth=="hillclimb")
      lab.optimize.hillclimb(d1,d2,FUN,exchange.list,seek,...)
   else if(meth=="gumbel"){
      warning("Warning, gumbel method not yet supported. Try at your own risk.\n")
      lab.optimize.gumbel(d1,d2,FUN,exchange.list,seek,...)
   }
}

lab.optimize.exhaustive<-function(d1,d2,FUN,exchange.list=0,seek="min",...){
   #Find the data set size
   n<-dim(d1)[2]
   #If exchange list is a single number or vector, expand it via replication in a reasonable manner
   if(is.null(dim(exchange.list))){       #Exchange list was given as a single number or vector
      if(length(exchange.list)==1){                 #Single number case
         el<-matrix(rep(exchange.list,2*n),nrow=2,ncol=n)
      }else{                                                    #Vector case
         el<-sapply(exchange.list,rep,2)
      }  
   }else                         #Exchange list was given as a matrix; keep it.
      el<-exchange.list
   #Initialize various things
   fun<-match.fun(FUN)       #Find the function to be optimized
   d1<-d1[order(el[1,]),order(el[1,])]  #Reorder d1
   d2<-d2[order(el[2,]),order(el[2,])]  #Reorder d2
   el[1,]<-el[1,order(el[1,])]  #Reorder the exchange lists to match
   el[2,]<-el[2,order(el[2,])]
   if(any(el[1,]!=el[2,]))  #Make sure the exlist is legal
      stop("Illegal exchange list; lists must be comparable!\n")
   best<-fun(d1,d2,...)  #Take the seed value (this has to be legal)
   #Search exhaustively - I hope you're not in a hurry!
   if(any(duplicated(el[1,])))  #If we're dealing with the labeled case, don't bother.
      for(k in 0:(gamma(n+1)-1)){
         o<-numperm(n,k)
         if(all(el[1,]==el[2,o])){
            if(seek=="min")
               best<-min(best,fun(d1,d2[o,o],...))
            else
               best<-max(best,fun(d1,d2[o,o],...))
         }
      }
   #Report the results
   best
}


lab.optimize.mc<-function(d1,d2,FUN,exchange.list=0,seek="min",draws=1000,...){
   #Find the data set size
   n<-dim(d1)[2]
   #If exchange list is a single number or vector, expand it via replication in a reasonable manner
   if(is.null(dim(exchange.list))){       #Exchange list was given as a single number or vector
      if(length(exchange.list)==1){                 #Single number case
         el<-matrix(rep(exchange.list,2*n),nrow=2,ncol=n)
      }else{                                                    #Vector case
         el<-sapply(exchange.list,rep,2)
      }  
   }else                         #Exchange list was given as a matrix; keep it.
      el<-exchange.list
   #Initialize various things
   fun<-match.fun(FUN)       #Find the function to be optimized
   d1<-d1[order(el[1,]),order(el[1,])]  #Reorder d1
   d2<-d2[order(el[2,]),order(el[2,])]  #Reorder d2
   el[1,]<-el[1,order(el[1,])]  #Reorder the exchange lists to match
   el[2,]<-el[2,order(el[2,])]
   if(any(el[1,]!=el[2,]))  #Make sure the exlist is legal
      stop("Illegal exchange list; lists must be comparable!\n")
   best<-fun(d1,d2,...)  #Take the seed value (this has to be legal)
   #Search via blind monte carlo - slow, yet ineffectual
   if(any(duplicated(el[1,])))  #If we're dealing with the labeled case, don't bother.
      for(i in 1:draws){
         o<-rperm(el[2,])
         if(seek=="min")
            best<-min(best,fun(d1,d2[o,o],...))
         else
            best<-max(best,fun(d1,d2[o,o],...))
      }
   #Report the results
   best
}


lab.optimize.gumbel<-function(d1,d2,FUN,exchange.list=0,seek="min",draws=500,tol=1e-5,estimator="median",...){
   #Find the data set size
   n<-dim(d1)[2]
   #If exchange list is a single number or vector, expand it via replication in a reasonable manner
   if(is.null(dim(exchange.list))){       #Exchange list was given as a single number or vector
      if(length(exchange.list)==1){                 #Single number case
         el<-matrix(rep(exchange.list,2*n),nrow=2,ncol=n)
      }else{                                                    #Vector case
         el<-sapply(exchange.list,rep,2)
      }  
   }else                         #Exchange list was given as a matrix; keep it.
      el<-exchange.list
   #Initialize various things
   fun<-match.fun(FUN)       #Find the function to be optimized
   fg<-vector()                    #Set up the function
   d1<-d1[order(el[1,]),order(el[1,])]  #Reorder d1
   d2<-d2[order(el[2,]),order(el[2,])]  #Reorder d2
   el[1,]<-el[1,order(el[1,])]  #Reorder the exchange lists to match
   el[2,]<-el[2,order(el[2,])]
   if(any(el[1,]!=el[2,]))  #Make sure the exlist is legal
      stop("Illegal exchange list; lists must be comparable!\n")
   #Approximate the distribution using Monte Carlo
   for(i in 1:draws){
      o<-rperm(el[2,])
      fg[i]<-fun(d1,d2[o,o],...)
   }
   #Use the approximated distribution to fit a Gumbel model for the extreme values;
   #this is only approximate, since the extreme value model assumes an unbounded, continuous underlying
   #distribution.  Also, these results are "unproven," in the sense that no actual permutation has been
   #found by the algorithm which results in the predicted value (unlike the other methods); OTOH, in 
   #a world of approximations, this one may not be any worse than the others....
   b<-1
   b.old<-1
   bdiff<-Inf
   mfg<-mean(fg)
   print(quantile(fg))
   while(bdiff>tol){         #Solve iteratively for bhat
      cat("bold=",b.old,"b=",b,"bdiff=",bdiff,"\n")
      b.old<-b
      b<-mfg-sum(fg*exp(-fg/b))/sum(exp(-fg/b))
      bdiff<-abs(b.old-b)
   }
   a<--b*log(sum(exp(-fg/b))/draws)     #Given this, ahat is a function of bhat and the data
   #Report the results
   cat("a=",a,"b=",b,"\n")
   switch(estimator,
      mean=a-b*digamma(1),
      mode=a,
      median=a-b*log(log(2))
   )
}


lab.optimize.anneal<-function(d1,d2,FUN,exchange.list=0,seek="min",prob.init=1,prob.decay=0.99,freeze.time=1000,full.neighborhood=TRUE,...){
   #Find the data set size
   n<-dim(d1)[2]
   #If exchange list is a single number or vector, expand it via replication in a reasonable manner
   if(is.null(dim(exchange.list))){       #Exchange list was given as a single number or vector
      if(length(exchange.list)==1){                 #Single number case
         el<-matrix(rep(exchange.list,2*n),nrow=2,ncol=n)
      }else{                                                    #Vector case
         el<-sapply(exchange.list,rep,2)
      }  
   }else                         #Exchange list was given as a matrix; keep it.
      el<-exchange.list
   #Initialize various things
   fun<-match.fun(FUN)       #Find the function to be optimized
   d1<-d1[order(el[1,]),order(el[1,])]  #Reorder d1
   d2<-d2[order(el[2,]),order(el[2,])]  #Reorder d2
   el[1,]<-el[1,order(el[1,])]  #Reorder the exchange lists to match
   el[2,]<-el[2,order(el[2,])]
   if(any(el[1,]!=el[2,]))  #Make sure the exlist is legal
      stop("Illegal exchange list; lists must be comparable!\n")
   best<-fun(d1,d2,...)  #Take the seed value (this has to be legal)
   o<-1:n                              #Set the initial ordering
   global.best<-best              #Set global best values
   global.o<-o
   prob<-prob.init                  #Set acceptance prob
   ftime<-freeze.time            #Set time until freezing occurs
   nc<-choose(n,2)               #How many candidate steps?
   candp<-sapply(o,rep,choose(n,2))   #Build the candidate permutation matrix
   ccount<-1
   for(i in 1:n)
      for(j in i:n)
         if(i!=j){                                          #Perform binary exchanges 
            temp<-candp[ccount,i]
            candp[ccount,i]<-candp[ccount,j]
            candp[ccount,j]<-temp
            ccount<-ccount+1
         }         
   #Run the annealer
   flag<-FALSE
   if(any(duplicated(el[2,])))                #If we're dealing with the labeled case, don't bother.
      while((!flag)|(ftime>0)){           #Until we both freeze _and_ reach an optimum...
         #cat("Best: ",o," Perf: ",best," Global best: ",global.o," Global perf: ",global.best," Temp: ",prob,"\n")
         #cat("Perf: ",best," Global perf: ",global.best," Temp: ",prob,"\n")
         flag<-TRUE
         if(full.neighborhood){ #Full neighborhood search method - much slower, but more likely to find the optimum
            candperf<-vector()
            for(i in 1:nc)                           #Use candidate permutation matrix to produce new candidates
               if(all(el[2,]==el[2,o[candp[i,]]]))   #Is this legal?
                  candperf[i]<-fun(d1,d2[o[candp[i,]],o[candp[i,]]],...)
               else
                  candperf[i]<-NA   #If not, put the results in as missing data
            if(seek=="min"){
               bestcand<-(1:nc)[candperf==min(candperf,na.rm=TRUE)]   #Find the best candidate
               bestcand<-bestcand[!is.na(bestcand)]            
               if(length(bestcand)>1)
                  bestcand<-sample(bestcand,1)  #If we have multiple best candidates, choose one at random
               #cat(min(candperf,na.rm=TRUE),bestcand,candperf[bestcand],"\n")
               if(candperf[bestcand]<best){  #If this is better, move on and keep looking...        
                  o<-o[candp[bestcand,]]
                  best<-candperf[bestcand]
                  flag<-FALSE
                  if(best<global.best){   #Check to see if this is better than the global best
                     global.best<-best
                     global.o<-o
                  }
               }else if((ftime>0)&(runif(1,0,1)<prob)){   #...but if not frozen and no better option, take a chance.
                  #print(candperf)
                  bestcand<-sample(1:nc,1)               #Choose randomly from the available options
                  while(!all(el[2,]==el[2,o[candp[bestcand,]]]))  #Make sure we have a legal one...
                     bestcand<-sample(1:nc,1)
                  #cat("Wildcard - perm ",bestcand," perf ",candperf[bestcand],"\n")
                  o<-o[candp[bestcand,]]                #Accept the new candidate
                  best<-candperf[bestcand]
               }
            }else{
               bestcand<-(1:nc)[candperf==max(candperf,na.rm=TRUE)]   #Find the best candidate
               bestcand<-bestcand[!is.na(bestcand)]            
               if(length(bestcand)>1)
                  bestcand<-sample(bestcand,1)          #If we have multiple best candidates, choose one at random
               if((candperf[bestcand]>best)|(runif(1,0,1)<prob)){  #If this is better, move on and keep looking...        
                  o<-o[candp[bestcand,]]
                  best<-candperf[bestcand]
                  flag<-FALSE
                  if(best>global.best){   #Check to see if this is better than the global best
                     global.best<-best
                     global.o<-o
                  }
               }else if((ftime>0)&(runif(1,0,1)<prob)){   #...but if not frozen and no better option, take a chance.
                  bestcand<-sample(1:nc,1)               #Choose randomly from the available options
                  while(!all(el[2,]==el[2,o[candp[bestcand,]]]))  #Make sure we have a legal one...
                     bestcand<-sample(1:nc,1)
                  o<-o[candp[bestcand,]]                #Accept the new candidate
                  best<-candperf[bestcand]
               }
            }
         }else{      #Single candidate method.  Much faster, but less likely to find the optimum.
            #Use candidate permutation matrix to produce new candidates
            i<-sample(1:nc,1)       
            while(!all(el[2,]==el[2,o[candp[i,]]]))   #Is this legal?
               i<-sample(1:nc,1)       #Keep trying till we get it right.
            #Assess candidate performance
            candperf<-fun(d1,d2[o[candp[i,]],o[candp[i,]]],...)
            #Make a decision
            if(seek=="min"){
               if(candperf<best){  #If this is better, move on and keep looking...        
                  o<-o[candp[i,]]
                  best<-candperf
                  flag<-FALSE
                  if(best<global.best){   #Check to see if this is better than the global best
                     global.best<-best
                     global.o<-o
                  }
               }else if((ftime>0)&(runif(1,0,1)<prob)){   #...but if not frozen and no better option, take a chance.
                  i<-sample(1:nc,1)               #Choose randomly from the available options
                  while(!all(el[2,]==el[2,o[candp[i,]]]))  #Make sure we have a legal one...
                     i<-sample(1:nc,1)
                  o<-o[candp[i,]]                #Accept the new candidate
                  best<-candperf
                  if(best<global.best){   #Check to see if this is better than the global best
                     global.best<-best
                     global.o<-o
                  }
               }
            }else{
               if(candperf>best){  #If this is better, move on and keep looking...        
                  o<-o[candp[i,]]
                  best<-candperf
                  flag<-FALSE
                  if(best>global.best){   #Check to see if this is better than the global best
                     global.best<-best
                     global.o<-o
                  }
               }else if((ftime>0)&(runif(1,0,1)<prob)){   #...but if not frozen and no better option, take a chance.
                  i<-sample(1:nc,1)               #Choose randomly from the available options
                  while(!all(el[2,]==el[2,o[candp[i,]]]))  #Make sure we have a legal one...
                     i<-sample(1:nc,1)
                  o<-o[candp[i,]]                #Accept the new candidate
                  best<-candperf
                  if(best>global.best){   #Check to see if this is better than the global best
                     global.best<-best
                     global.o<-o
                  }
               }
            }
         }
         #Set things up for the next iteration (if there is one)
         ftime<-ftime-1               #Continue the countdown to the freezing point
         prob<-prob*prob.decay   #Cool things off a bit
      }
   #Report the results
   global.best
}

lab.optimize.hillclimb<-function(d1,d2,FUN,exchange.list=0,seek="min",...){
   #Find the data set size
   n<-dim(d1)[2]
   #If exchange list is a single number or vector, expand it via replication in a reasonable manner
   if(is.null(dim(exchange.list))){       #Exchange list was given as a single number or vector
      if(length(exchange.list)==1){                 #Single number case
         el<-matrix(rep(exchange.list,2*n),nrow=2,ncol=n)
      }else{                                                    #Vector case
         el<-sapply(exchange.list,rep,2)
      }  
   }else                         #Exchange list was given as a matrix; keep it.
      el<-exchange.list
   #Initialize various things
   fun<-match.fun(FUN)       #Find the function to be optimized
   d1<-d1[order(el[1,]),order(el[1,])]  #Reorder d1
   d2<-d2[order(el[2,]),order(el[2,])]  #Reorder d2
   el[1,]<-el[1,order(el[1,])]  #Reorder the exchange lists to match
   el[2,]<-el[2,order(el[2,])]
   if(any(el[1,]!=el[2,]))  #Make sure the exlist is legal
      stop("Illegal exchange list; lists must be comparable!\n")
   best<-fun(d1,d2,...)  #Take the seed value (this has to be legal)
   o<-1:n                              #Set the initial ordering
   nc<-choose(n,2)               #How many candidate steps?
   candp<-sapply(o,rep,choose(n,2))   #Build the candidate permutation matrix
   ccount<-1
   for(i in 1:n)
      for(j in i:n)
         if(i!=j){                                          #Perform binary exchanges 
            temp<-candp[ccount,i]
            candp[ccount,i]<-candp[ccount,j]
            candp[ccount,j]<-temp
            ccount<-ccount+1
         }         
   #Run the hill climber
   flag<-FALSE
   while(!flag){           #Until we reach an optimum...
      #cat("Best: ",o," Perf: ",best,"\n")
      flag<-TRUE
      candperf<-vector()
      for(i in 1:nc)                           #Use candidate permutation matrix to produce new candidates
         if(all(el[2,]==el[2,o[candp[i,]]]))   #Is this legal?
            candperf[i]<-fun(d1,d2[o[candp[i,]],o[candp[i,]]],...)
         else
            candperf[i]<-NA   #If not, put the results in as missing data
      if(seek=="min"){
         bestcand<-(1:nc)[candperf==min(candperf,na.rm=TRUE)]   #Find the best candidate
         if(length(bestcand)>1)
            bestcand<-sample(bestcand,1)          #If we have multiple best candidates, choose one at random
         if(candperf[bestcand]<best){          #If this is better, move on and keep looking...        
            o<-o[candp[bestcand,]]
            best<-candperf[bestcand]
            flag<-FALSE
         }
      }else{
         bestcand<-(1:nc)[candperf==max(candperf,na.rm=TRUE)]   #Find the best candidate
         if(length(bestcand)>1)
            bestcand<-sample(bestcand,1)          #If we have multiple best candidates, choose one at random
         if(candperf[bestcand]>best){          #If this is better, move on and keep looking...        
            o<-o[candp[bestcand,]]
            best<-candperf[bestcand]
            flag<-FALSE
         }
      }
   }
   #Report the results
   best
}


#hdist - Find the Hamming distances between two or more graphs.

hdist<-function(dat,dat2=NULL,g1=c(1:dim(dat)[1]),g2=c(1:dim(dat)[1]),normalize=FALSE,diag=FALSE,mode="digraph"){
   #Collect data and various parameters
   if(!is.null(dat2)){
      if(length(dim(dat))>2)
         temp1<-dat
      else{
         temp1<-array(dim=c(1,dim(dat)[2],dim(dat)[2]))
         temp1[1,,]<-dat
      }
      if(length(dim(dat2))>2)
         temp2<-dat2
      else{
         temp2<-array(dim=c(1,dim(dat2)[2],dim(dat2)[2]))
         temp2[1,,]<-dat2
      }
      if(dim(temp1)[2]>dim(temp2)[2])
         temp2<-addisolates(temp2,dim(temp1)[2]-dim(temp2)[2])
      if(dim(temp2)[2]>dim(temp1)[2])
         temp1<-addisolates(temp1,dim(temp2)[2]-dim(temp1)[2])
      n<-dim(temp1)[2]
      gn<-dim(temp1)[1]+dim(temp2)[1]
      gn1<-dim(temp1)[1]
      gn2<-dim(temp2)[1]
      d<-array(dim=c(gn,n,n))
      d[1:gn1,,]<-temp1
      d[(gn1+1):(gn2+gn1),,]<-temp2
      g1<-1:gn1
      g2<-(gn1+1):(gn1+gn2)
   }else{
      d<-dat
      n<-dim(dat)[2]
      gn<-dim(dat)[1]
      gn1<-length(g1)
      gn2<-length(g2)
   }
   #Scrap the diagonals, if required
   if(!diag)
      d<-diag.remove(d)
   #Now, get rid of the upper triangle if these are simple graphs
   if(mode=="graph")
      d<-upper.tri.remove(d)
   #Compute the raw distance matrix
   hd<-matrix(nrow=gn1,ncol=gn2)
   rownames(hd)<-g1
   colnames(hd)<-g2
   for(i in 1:gn1)
      for(j in 1:gn2)
         hd[i,j]<-sum(abs(d[g1[i],,]-d[g2[j],,]),na.rm=TRUE)
   #Normalize if need be
   if(normalize)
      hd<-hd/nties(dat[1,,],mode=mode,diag=diag)
   #If only one comparison requested, return as an element
   if((gn1==1)&(gn2==1))
      hd[1,1]
   else
      hd
}


#rperm - Draw a random permutation vector with exchangability constraints

rperm<-function(exchange.list){
   #Note that exchange.list should be a vector whose entries correspond to the class identity
   #of the respective element.  It doesn't matter what the values are, so long as elements have
   #the same value iff they are exchangeable.
   n<-length(exchange.list)   #Get the length of the output vector
   grp<-unique(exchange.list) #Get the groups
   o<-1:n  #Create the initial ordering
   #Randomly scramble orders within groups
   for(i in grp){
      v<-(1:n)[exchange.list==i]
      if(length(v)>1)      #Need this test, because sample is too smart for its own good...
         o[v]<-sample(v)
   }
   #Return the permutation
   o
}


#numperm - Get the nth permutation vector by periodic placement

numperm<-function(olength,permnum){
   if((permnum>gamma(olength+1)-1)|(permnum<0)){
      cat("permnum must be an integer in [0,olength!-1]\n")
   }
   o<-vector(length=olength)
   o[]<--1
   pnum<-permnum
   for(i in 1:olength){
      relpos<-pnum%%(olength-i+1)
      flag<-FALSE
      p<-1
      while(!flag)
         if(o[p]==-1){
            if(relpos==0){
               o[p]<-i
               flag<-TRUE
            }else{
               p<-p+1
               relpos<-relpos-1
            }
         }else
            p<-p+1      
      pnum<-pnum%/%(olength-i+1)
   }
   o
}


#structdist - Estimate the structural distance between two or more unlabeled graphs.

structdist<-function(dat,g1=c(1:dim(dat)[1]),g2=c(1:dim(dat)[1]),normalize=FALSE,diag=FALSE,mode="digraph",method="anneal",reps=1000,prob.init=0.9,prob.decay=0.85,freeze.time=25,full.neighborhood=TRUE,mut=0.05,pop=20,trials=5,exchange.list=rep(0,dim(dat)[2])){
   #Collect data and various parameters
   d<-dat
   n<-dim(dat)[2]
   gn<-dim(dat)[1]
   gn1<-length(g1)
   gn2<-length(g2)
   #Scrap the diagonals, if required
   if(!diag)
      d<-diag.remove(d)
   #Now, get rid of the upper triangle if these are simple graphs
   if(mode=="graph")
      d<-upper.tri.remove(d)
   #If exchange list is a single number or vector, expand it via replication in a reasonable manner
   if(is.null(dim(exchange.list))){       #Exchange list was given as a single number or vector
      if(length(exchange.list)==1){             #Single number case
         el<-matrix(rep(exchange.list,gn*n),nrow=gn,ncol=n)
      }else{                                    #Vector case
         el<-sapply(exchange.list,rep,gn)
      }  
   }else                #Exchange list was given as a matrix; keep it.
      el<-exchange.list
   #Compute the distance matrix
   hd<-matrix(nrow=gn1,ncol=gn2)
   rownames(hd)<-g1
   colnames(hd)<-g2
   if(method=="exhaustive"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            hd[i,j]<-lab.optimize.exhaustive(d[g1[i],,],d[g2[j],,],function(m1,m2){sum(abs(m1-m2),na.rm=TRUE)},exchange.list=el[c(g1[i],g2[j]),],seek="min")
   }else if(method=="anneal"){
      for(i in 1:gn1)
         for(j in 1:gn2){
            hd[i,j]<-lab.optimize.anneal(d[g1[i],,],d[g2[j],,],function(m1,m2){sum(abs(m1-m2),na.rm=TRUE)},exchange.list=el[c(g1[i],g2[j]),],seek="min",prob.init=prob.init,prob.decay=prob.decay,freeze.time=freeze.time,full.neighborhood=full.neighborhood)
         }
   }else if(method=="hillclimb"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            hd[i,j]<-lab.optimize.hillclimb(d[g1[i],,],d[g2[j],,],function(m1,m2){sum(abs(m1-m2),na.rm=TRUE)},exchange.list=el[c(g1[i],g2[j]),],seek="min")
   }else if(method=="mc"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            hd[i,j]<-lab.optimize.mc(d[g1[i],,],d[g2[j],,],function(m1,m2){sum(abs(m1-m2),na.rm=TRUE)},exchange.list=el[c(g1[i],g2[j]),],seek="min",draws=reps)
   }else if(method=="ga"){
      #This is broken right now - exit with a warning
      stop("Sorry, GA mode is not currently supported.\n")
   }else if(method=="none"){
      for(i in 1:gn1)
         for(j in 1:gn2){
            d1<-d[g1[i],order(el[1,]),order(el[1,])]  #Reorder d1
            d2<-d[g2[j],order(el[2,]),order(el[2,])]  #Reorder d2
            if(any(el[1,]!=el[2,]))  #Make sure the exlist is legal
               stop("Illegal exchange list; lists must be comparable!\n")
            hd[i,j]<-sum(abs(d1-d2),na.rm=TRUE)
         }
   }else{
      cat("Method",method,"not implemented yet.\n")
   }
   #Normalize if need be
   if(normalize)
      hd<-hd/nties(dat[1,,],mode=mode,diag=diag)
   #If only one comparison requested, return as an element
   if((gn1==1)&(gn2==1))
      hd[1,1]
   else
      hd
}

#sdmat - Estimate the matrix of structural distances among a set of unlabeled graphs.
#NOTE: This is redundant, as the included functionality is now included within hdist and
#structdist, but the function is left for reasons of compatibility as well as speed (currently,
#the distance functions don't check for duplicate calculations when building distance matrices).

sdmat<-function(dat,normalize=FALSE,diag=FALSE,mode="digraph",output="matrix",method="mc",exchange.list=rep(0,dim(dat)[2]),...){
   m<-dim(dat)[1]
   sdm<-matrix(nrow=m,ncol=m)
   for(i in 1:m)
      if(i<m)
         for(j in i:m)
         sdm[i,j]<-structdist(dat,g1=i,g2=j,normalize=normalize,diag=diag,mode=mode,method=method,exchange.list=exchange.list,...)
   diag(sdm)<-0
   for(i in 1:m)
      sdm[i,c(1:i)[-i]]<-sdm[c(1:i)[-i],i]
   if(output=="dist"){
      #This is for compatibility with the mva library
      require(mva)
      sdm<-as.dist(sdm)
   }
   sdm
}


#gcor - Correlation between two or more graphs.

gcor<-function(dat,dat2=NULL,g1=c(1:dim(dat)[1]),g2=c(1:dim(dat)[1]),diag=FALSE,mode="digraph"){
   #Collect data and various parameters
   if(!is.null(dat2)){
      if(length(dim(dat))>2)
         temp1<-dat
      else{
         temp1<-array(dim=c(1,dim(dat)[2],dim(dat)[2]))
         temp1[1,,]<-dat
      }
      if(length(dim(dat2))>2)
         temp2<-dat2
      else{
         temp2<-array(dim=c(1,dim(dat2)[2],dim(dat2)[2]))
         temp2[1,,]<-dat2
      }
      if(dim(temp1)[2]>dim(temp2)[2])
         temp2<-addisolates(temp2,dim(temp1)[2]-dim(temp2)[2])
      if(dim(temp2)[2]>dim(temp1)[2])
         temp1<-addisolates(temp1,dim(temp2)[2]-dim(temp1)[2])
      n<-dim(temp1)[2]
      gn<-dim(temp1)[1]+dim(temp2)[1]
      gn1<-dim(temp1)[1]
      gn2<-dim(temp2)[1]
      d<-array(dim=c(gn,n,n))
      d[1:gn1,,]<-temp1
      d[(gn1+1):(gn2+gn1),,]<-temp2
      g1<-1:gn1
      g2<-(gn1+1):(gn1+gn2)
   }else{
      d<-dat
      n<-dim(dat)[2]
      gn<-dim(dat)[1]
      gn1<-length(g1)
      gn2<-length(g2)
   }
   #Scrap the diagonals, if required
   if(!diag)
      d<-diag.remove(d)
   #Now, get rid of the upper triangle if these are simple graphs
   if(mode=="graph")
      d<-upper.tri.remove(d)
   #Compute the graph correlation matrix
   gd<-matrix(nrow=gn1,ncol=gn2)
   rownames(gd)<-g1
   colnames(gd)<-g2
   for(i in 1:gn1)
      for(j in 1:gn2)
         gd[i,j]<-cor(as.vector(d[g1[i],,]),as.vector(d[g2[j],,]),use="complete.obs")
   #If only one comparison requested, return as an element
   if((gn1==1)&(gn2==1))
      gd[1,1]
   else
      gd
}


#gcov - Covariance between two or more graphs.

gcov<-function(dat,dat2=NULL,g1=c(1:dim(dat)[1]),g2=c(1:dim(dat)[1]),diag=FALSE,mode="digraph"){
   #Collect data and various parameters
   if(!is.null(dat2)){
      if(length(dim(dat))>2)
         temp1<-dat
      else{
         temp1<-array(dim=c(1,dim(dat)[2],dim(dat)[2]))
         temp1[1,,]<-dat
      }
      if(length(dim(dat2))>2)
         temp2<-dat2
      else{
         temp2<-array(dim=c(1,dim(dat2)[2],dim(dat2)[2]))
         temp2[1,,]<-dat2
      }
      if(dim(temp1)[2]>dim(temp2)[2])
         temp2<-addisolates(temp2,dim(temp1)[2]-dim(temp2)[2])
      if(dim(temp2)[2]>dim(temp1)[2])
         temp1<-addisolates(temp1,dim(temp2)[2]-dim(temp1)[2])
      n<-dim(temp1)[2]
      gn<-dim(temp1)[1]+dim(temp2)[1]
      gn1<-dim(temp1)[1]
      gn2<-dim(temp2)[1]
      d<-array(dim=c(gn,n,n))
      d[1:gn1,,]<-temp1
      d[(gn1+1):(gn2+gn1),,]<-temp2
      g1<-1:gn1
      g2<-(gn1+1):(gn1+gn2)
   }else{
      d<-dat
      n<-dim(dat)[2]
      gn<-dim(dat)[1]
      gn1<-length(g1)
      gn2<-length(g2)
   }
   #Scrap the diagonals, if required
   if(!diag)
      d<-diag.remove(d)
   #Now, get rid of the upper triangle if these are simple graphs
   if(mode=="graph")
      d<-upper.tri.remove(d)
   #Compute the graph covariance matrix
   gd<-matrix(nrow=gn1,ncol=gn2)
   rownames(gd)<-g1
   colnames(gd)<-g2
   for(i in 1:gn1)
      for(j in 1:gn2)
         gd[i,j]<-cov(as.vector(d[g1[i],,]),as.vector(d[g2[j],,]),use="complete.obs")
   #If only one comparison requested, return as an element
   if((gn1==1)&(gn2==1))
      gd[1,1]
   else
      gd
}


#gscor - Structural correlation between two or more graphs.

gscor<-function(dat,dat2=NULL,g1=c(1:dim(dat)[1]),g2=c(1:dim(dat)[1]),diag=FALSE,mode="digraph",method="anneal",reps=1000,prob.init=0.9,prob.decay=0.85,freeze.time=25,full.neighborhood=TRUE,exchange.list=rep(0,dim(dat)[2])){
   #Collect data and various parameters
   if(!is.null(dat2)){
      if(length(dim(dat))>2)
         temp1<-dat
      else{
         temp1<-array(dim=c(1,dim(dat)[2],dim(dat)[2]))
         temp1[1,,]<-dat
      }
      if(length(dim(dat2))>2)
         temp2<-dat2
      else{
         temp2<-array(dim=c(1,dim(dat2)[2],dim(dat2)[2]))
         temp2[1,,]<-dat2
      }
      if(dim(temp1)[2]>dim(temp2)[2])
         temp2<-addisolates(temp2,dim(temp1)[2]-dim(temp2)[2])
      if(dim(temp2)[2]>dim(temp1)[2])
         temp1<-addisolates(temp1,dim(temp2)[2]-dim(temp1)[2])
      n<-dim(temp1)[2]
      gn<-dim(temp1)[1]+dim(temp2)[1]
      gn1<-dim(temp1)[1]
      gn2<-dim(temp2)[1]
      d<-array(dim=c(gn,n,n))
      d[1:gn1,,]<-temp1
      d[(gn1+1):(gn2+gn1),,]<-temp2
      g1<-1:gn1
      g2<-(gn1+1):(gn1+gn2)
   }else{
      d<-dat
      n<-dim(dat)[2]
      gn<-dim(dat)[1]
      gn1<-length(g1)
      gn2<-length(g2)
   }
   #Scrap the diagonals, if required
   if(!diag)
      d<-diag.remove(d)
   #Now, get rid of the upper triangle if these are simple graphs
   if(mode=="graph")
      d<-upper.tri.remove(d)
   #If exchange list is a single number or vector, expand it via replication in a reasonable manner
   if(is.null(dim(exchange.list))){       #Exchange list was given as a single number or vector
      if(length(exchange.list)==1){                 #Single number case
         el<-matrix(rep(exchange.list,gn*n),nrow=gn,ncol=n)
      }else{                                                    #Vector case
         el<-sapply(exchange.list,rep,gn)
      }  
   }else                         #Exchange list was given as a matrix; keep it.
      el<-exchange.list
   #Compute the structural correlation matrix
   gd<-matrix(nrow=gn1,ncol=gn2)
   rownames(gd)<-g1
   colnames(gd)<-g2
   if(method=="none"){
      for(i in 1:gn1)
         for(j in 1:gn2){
            d1<-d[g1[i],order(el[1,]),order(el[1,])]  #Reorder d1
            d2<-d[g2[j],order(el[2,]),order(el[2,])]  #Reorder d2
            if(any(el[1,]!=el[2,]))  #Make sure the exlist is legal
               stop("Illegal exchange list; lists must be comparable!\n")
            gd[i,j]<-cor(as.vector(d1),as.vector(d2),use="complete.obs")
         }
   }else if(method=="exhaustive"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            gd[i,j]<-lab.optimize.exhaustive(d[g1[i],,],d[g2[j],,],function(m1,m2){cor(as.vector(m1),as.vector(m2),use="complete.obs")},exchange.list=el[c(g1[i],g2[j]),],seek="max")
   }else if(method=="anneal"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            gd[i,j]<-lab.optimize.anneal(d[g1[i],,],d[g2[j],,],function(m1,m2){cor(as.vector(m1),as.vector(m2),use="complete.obs")},exchange.list=el[c(g1[i],g2[j]),],seek="max",prob.init=prob.init,prob.decay=prob.decay,freeze.time=freeze.time,full.neighborhood=full.neighborhood)
   }else if(method=="hillclimb"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            gd[i,j]<-lab.optimize.hillclimb(d[g1[i],,],d[g2[j],,],function(m1,m2){cor(as.vector(m1),as.vector(m2),use="complete.obs")},exchange.list=el[c(g1[i],g2[j]),],seek="max")
   }else if(method=="mc"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            gd[i,j]<-lab.optimize.mc(d[g1[i],,],d[g2[j],,],function(m1,m2){cor(as.vector(m1),as.vector(m2),use="complete.obs")},exchange.list=el[c(g1[i],g2[j]),],seek="max",draws=reps)
   }
   #If only one comparison requested, return as an element
   if((gn1==1)&(gn2==1))
      gd[1,1]
   else
      gd
}


#gscov - Structural covariance between two or more graphs.

gscov<-function(dat,dat2=NULL,g1=c(1:dim(dat)[1]),g2=c(1:dim(dat)[1]),diag=FALSE,mode="digraph",method="anneal",reps=1000,prob.init=0.9,prob.decay=0.85,freeze.time=25,full.neighborhood=TRUE,exchange.list=rep(0,dim(dat)[2])){
   #Collect data and various parameters
   if(!is.null(dat2)){
      if(length(dim(dat))>2)
         temp1<-dat
      else{
         temp1<-array(dim=c(1,dim(dat)[2],dim(dat)[2]))
         temp1[1,,]<-dat
      }
      if(length(dim(dat2))>2)
         temp2<-dat2
      else{
         temp2<-array(dim=c(1,dim(dat2)[2],dim(dat2)[2]))
         temp2[1,,]<-dat2
      }
      if(dim(temp1)[2]>dim(temp2)[2])
         temp2<-addisolates(temp2,dim(temp1)[2]-dim(temp2)[2])
      if(dim(temp2)[2]>dim(temp1)[2])
         temp1<-addisolates(temp1,dim(temp2)[2]-dim(temp1)[2])
      n<-dim(temp1)[2]
      gn<-dim(temp1)[1]+dim(temp2)[1]
      gn1<-dim(temp1)[1]
      gn2<-dim(temp2)[1]
      d<-array(dim=c(gn,n,n))
      d[1:gn1,,]<-temp1
      d[(gn1+1):(gn2+gn1),,]<-temp2
      g1<-1:gn1
      g2<-(gn1+1):(gn1+gn2)
   }else{
      d<-dat
      n<-dim(dat)[2]
      gn<-dim(dat)[1]
      gn1<-length(g1)
      gn2<-length(g2)
   }
   #Scrap the diagonals, if required
   if(!diag)
      d<-diag.remove(d)
   #Now, get rid of the upper triangle if these are simple graphs
   if(mode=="graph")
      d<-upper.tri.remove(d)
   #If exchange list is a single number or vector, expand it via replication in a reasonable manner
   if(is.null(dim(exchange.list))){       #Exchange list was given as a single number or vector
      if(length(exchange.list)==1){                 #Single number case
         el<-matrix(rep(exchange.list,gn*n),nrow=gn,ncol=n)
      }else{                                                    #Vector case
         el<-sapply(exchange.list,rep,gn)
      }  
   }else                         #Exchange list was given as a matrix; keep it.
      el<-exchange.list
   #Compute the structural covariance matrix
   gd<-matrix(nrow=gn1,ncol=gn2)
   rownames(gd)<-g1
   colnames(gd)<-g2
   if(method=="none"){
      for(i in 1:gn1)
         for(j in 1:gn2){
            d1<-d[g1[i],order(el[1,]),order(el[1,])]  #Reorder d1
            d2<-d[g2[j],order(el[2,]),order(el[2,])]  #Reorder d2
            if(any(el[1,]!=el[2,]))  #Make sure the exlist is legal
               stop("Illegal exchange list; lists must be comparable!\n")
            gd[i,j]<-cov(as.vector(d1),as.vector(d2),use="complete.obs")
         }
   }else if(method=="exhaustive"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            gd[i,j]<-lab.optimize.exhaustive(d[g1[i],,],d[g2[j],,],function(m1,m2){cov(as.vector(m1),as.vector(m2),use="complete.obs")},exchange.list=el[c(g1[i],g2[j]),],seek="max")
   }else if(method=="anneal"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            gd[i,j]<-lab.optimize.anneal(d[g1[i],,],d[g2[j],,],function(m1,m2){cov(as.vector(m1),as.vector(m2),use="complete.obs")},exchange.list=el[c(g1[i],g2[j]),],seek="max",prob.init=prob.init,prob.decay=prob.decay,freeze.time=freeze.time,full.neighborhood=full.neighborhood)
   }else if(method=="hillclimb"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            gd[i,j]<-lab.optimize.hillclimb(d[g1[i],,],d[g2[j],,],function(m1,m2){cov(as.vector(m1),as.vector(m2),use="complete.obs")},exchange.list=el[c(g1[i],g2[j]),],seek="max")
   }else if(method=="mc"){
      for(i in 1:gn1)
         for(j in 1:gn2)
            gd[i,j]<-lab.optimize.mc(d[g1[i],,],d[g2[j],,],function(m1,m2){cov(as.vector(m1),as.vector(m2),use="complete.obs")},exchange.list=el[c(g1[i],g2[j]),],seek="max",draws=reps)
   }
   #If only one comparison requested, return as an element
   if((gn1==1)&(gn2==1))
      gd[1,1]
   else
      gd
}


#gliop - Return a binary operation on GLI values computed on two graphs (for test routines).

gliop<-function(dat,GFUN,OP="-",g1=1,g2=2,...){
   fun<-match.fun(GFUN)
   op<-match.fun(OP)
   op(fun(dat[g1,,],...),fun(dat[g2,,],...))
}


#cugtest - Generate, print, and plot CUG (conditional uniform graph) test objects.

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

summary.cugtest<-function(object, ...){
   out<-object
   class(out)<-c("summary.cugtest",class(out))
   out
}

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

print.cugtest<-function(x,...){
   cat("\nCUG Test Results\n\n")
   cat("Estimated p-values:\n")
   cat("\tp(f(rnd) >= f(d)):",x$pgreq,"\n")
   cat("\tp(f(rnd) <= f(d)):",x$pleeq,"\n\n")      
}

plot.cugtest<-function(x,mode="density",...){
   if(mode=="density"){
      plot(density(x$dist),main="Estimated Density of CUG Replications",xlab="Test Statistic",...)
   }else{
      hist(x$dist,main="Histogram of CUG Replications",xlab="Test Statistic",...)
   }
   abline(v=x$testval,lty=2)
}

#geodist - Find the numbers and lengths of geodesics among nodes in a graph using a BFS, a la
#Brandes (2000).  (Thanks, Ulrik!)

geodist<-function(dat,inf.replace=dim(dat)[2]){
   n<-dim(dat)[2]
   #Initialize the matrices
   sigma<-matrix(0,nrow=n,ncol=n)
   gd<-matrix(inf.replace,nrow=n,ncol=n)  #Use the infinite replace value
   #Cycle through each node, performing a BFS
   for(v in 1:n){
      visited<-rep(0,n)   #No nodes have been visited
      i<-v                      #Start with the source node
      visited[i]<-1
      sigma[v,v]<-1
      gd[v,v]<-0
      #cat("\nBeginning trace for node",v,"\n")
      while(any(visited==1)){
         while(any(visited==1)){
            i<-match(1,visited)   #Increment to the next visitable node 
            #cat(i,"->")
            visited[i]<-3     #Mark this node as visited
            for(j in (1:n)[((visited==0)|(visited==2))&dat[i,]]){
               #Walk through the unvisited neighborhood
               #cat(j)
               if(visited[j]==0)                #If we've never seen this node, we'll visit it next time
                  visited[j]<-2
               gd[v,j]<-gd[v,i]+1          #Geodesic distance is this node's+1
               sigma[v,j]<-sigma[v,j]+sigma[v,i]     #Add the accumulated paths to the total path count
            }
            #cat("\n")
            if(any(visited==1))
               i<-match(1,visited)   #Increment to the next visitable node 
         }   #Continue until we run out of nodes at this level
         visited[visited==2]<-1     #Mark the to-be-visited nodes as visitable
      }   #Keep going until there's no one left at all
   }
   #Return the results
   o<-list()
   o$counts<-sigma
   o$gdist<-gd
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


#component.dist - Returns a data frame containing a vector of length n such that the ith element
#contains the number of components of G having size i, and a vector of length n giving component
#membership.  Component strength is determined by the rule which is used to symmetrize the matrix;
#this controlled by the eponymous parameter given to the symmetrize command.

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
   #Cycle through each node, performing a BFS
   for(v in 1:n)
      if(membership[v]==0){   #Ignore nodes whose membership is already certain
         comp<-max(membership)+1  #Which component are we on?
         visited<-rep(0,n)   #No nodes have been visited
         i<-v                      #Start with the source node
         visited[i]<-1
         membership[i]<-comp   #Set initial membership
         #cat("\nBeginning trace for node",v,"\n")
         while(any(visited==1)){
            while(any(visited==1)){
               i<-match(1,visited)   #Increment to the next visitable node 
               #cat(i,"->")
               visited[i]<-3     #Mark this node as visited
               membership[i]<-comp  #Set membership to current component
               visited[(visited==0)&dat[i,]]<-2  #We'll visit these next time
            }  #Continue until we run out of nodes at this level
            #cat("\n")
            visited[visited==2]<-1  #Mark the to-be-visited nodes as visitable
         }  #Keep going until there's no one left at all   
      }   
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


#betweenness - Find the betweenness centralities of network positions

betweenness<-function(dat,g=1,nodes=c(1:dim(dat)[2]),gmode="digraph",diag=FALSE,tmaxdev=FALSE,cmode="directed",geodist.precomp=NULL,rescale=FALSE){
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
      for(i in 1:n){
         for(j in 1:n)
            for(k in 1:n)
               if((j!=i)&(k!=i)&(j!=k)&(gd$gdist[j,k]>=gd$gdist[j,i]+gd$gdist[i,k]))
                  bet[i]<-bet[i]+gd$counts[j,i]*gd$counts[i,k]/gd$counts[j,k]
      }
      if(cmode=="undirected")
         bet<-bet/2
      #Return the results
      if(rescale)
         bet<-bet/sum(bet)
      bet<-bet[nodes]
   }
   bet
}


#stresscent - Find the stress centralities of network positions

stresscent<-function(dat,g=1,nodes=c(1:dim(dat)[2]),gmode="digraph",diag=FALSE,tmaxdev=FALSE,cmode="directed",geodist.precomp=NULL,rescale=FALSE){
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
      for(i in 1:n){
         for(j in 1:n)
            for(k in 1:n)
               if((j!=i)&(k!=i)&(j!=k)&(gd$gdist[j,k]>=gd$gdist[j,i]+gd$gdist[i,k]))
                  str[i]<-str[i]+gd$counts[j,i]*gd$counts[i,k]
      }
      if(cmode=="undirected")
         str<-str/2
      #Return the results
      if(rescale)
         str<-str/sum(str)
      str<-str[nodes]
   }
   str
}


#closeness - Find the closeness centralities of network positions

closeness<-function(dat,g=1,nodes=c(1:dim(dat)[2]),gmode="digraph",diag=FALSE,tmaxdev=FALSE,cmode="directed",geodist.precomp=NULL,rescale=FALSE){
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


#graphcent - Find the graph centralities of network positions

graphcent<-function(dat,g=1,nodes=c(1:dim(dat)[2]),gmode="digraph",diag=FALSE,tmaxdev=FALSE,cmode="directed",geodist.precomp=NULL,rescale=FALSE){
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


#degree - Find the degree centralities of network positions

degree<-function(dat,g=1,nodes=c(1:dim(dat)[2]),gmode="digraph",diag=FALSE,tmaxdev=FALSE,cmode="freeman",rescale=FALSE){
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

evcent<-function(dat,g=1,nodes=c(1:dim(dat)[2]),gmode="digraph",diag=FALSE,tmaxdev=FALSE,rescale=FALSE){
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


#bonpow - Find the Bonacich power centrality scores of network positions

bonpow<-function(dat,g=1,nodes=c(1:dim(dat)[2]),gmode="digraph",diag=FALSE,tmaxdev=FALSE,exponent=1,rescale=FALSE,tol=1e-7){
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


#prestige - Find actor prestige scores from one of several measures

prestige<-function(dat,g=1,nodes=c(1:dim(dat)[2]),gmode="digraph",diag=FALSE,cmode="indegree",tmaxdev=FALSE,rescale=FALSE,tol=1e-7){
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


# infocent - Find actor information centrality scores
# Wasserman & Faust pp. 192-197; based on code generously submitted by David
# Barron (thanks!) and tweaked by myself to enable compatibility with the
# centralization() routine.
infocent <- function(dat,g=1,nodes=c(1:dim(dat)[2]),gmode="digraph",diag=FALSE,cmode="weak",tmaxdev=FALSE,rescale=FALSE,tol=1e-20){
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


#centralization - Find the centralization of a graph (for some arbitrary centrality measure)

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


#efficiency - Find the Krackhardt efficiency of a graph or graph stack

efficiency<-function(dat,g=1:stackcount(dat),diag=FALSE){
   #Define an internal function, for convenience
   inteff<-function(g,diag){
      comsz<-component.dist(g,connected="weak")$csize
      reqedge<-sum(comsz-1)     #Get required edges
      if(!diag)
         g<-diag.remove(g)
      edgec<-sum(g,na.rm=TRUE)  #Get count of actual edges
      1-(edgec-reqedge)/(prod(dim(g))-(!diag)*dim(g)[1]-reqedge)
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
         if(length(vi)>2){  #Don't bother unless we have at least three vertices
            #Accumulate violations of LUBness
            for(j in 1:length(vi))     #Check each dyad
               for(k in j:length(vi))
                  if(k>j){
                     ub<-vi[r[vi,vi[j]]*r[vi,vi[k]]>0]  #Get upper bounds
                     if((length(ub)==0)||(!any(apply(r[ub,ub,drop=FALSE],2,prod))))  #Any least upper bounds?
                        nolub<-nolub+1
                  }
            #Also accumulate maximum violations
            maxnolub<-maxnolub+(length(vi)-1)*(length(vi)-2)/2 
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


#gplot - Graph visualization
gplot<-function(dat,g=1,gmode="digraph",diag=FALSE,label=c(1:dim(dat)[2]),coord=NULL,jitter=TRUE,thresh=0,usearrows=TRUE,mode="segeo",displayisolates=TRUE,interactive=FALSE,boxed.labels=TRUE,xlab=NULL,ylab=NULL,pad=0.2,vertex.pch=19,label.cex=1,vertex.cex=1,label.col=1,edge.col=1,vertex.col=1,arrowhead.length=0.2,edge.lty=1,edge.lwd=0,edge.len=0.5,edge.curve=0.1,edge.steps=50,diag.size=0.025,uselen=FALSE,usecurve=FALSE,suppress.axes=TRUE,new=TRUE,layout.par=NULL,...){
   #Turn the annoying locator bell off
   bellstate<-options()$locatorBell
   on.exit(options(locatorBell=bellstate))
   options(locatorBell=FALSE)
   #Create a useful interval inclusion operator
   "%iin%"<-function(x,int) (x>=int[1])&(x<=int[2])
   #Extract the graph to be displayed
   if(length(dim(dat))>2)
      d<-dat[g,,]
   else
      d<-dat
   #Make adjustments for gmode, if required
   if(gmode=="graph"){
      usearrows<-FALSE
      n<-dim(d)[1]
   }else if(gmode=="twomode"){
      n<-sum(dim(d))
      temp<-matrix(0,nrow=n,ncol=n)
      temp[1:dim(d)[1],(dim(d)[1]+1):n]<-d
      d<-temp
      if(all(label==1:dim(dat)[2]))
         label<-1:n
   }else 
      n<-dim(d)[1]
   #Replace NAs with 0s
   d[is.na(d)]<-0
   #Save a copy of d, in case values are needed
   d.raw<-d
   #Dichotomize d
   d<-matrix(as.numeric(d>thresh),n,n)
   #Determine coordinate placement
   if(!is.null(coord)){      #If the user has specified coords, override all other considerations
      x<-coord[,1]
      y<-coord[,2]
   }else{   #Otherwise, use the specified layout function
     layout.fun<-try(match.fun(paste("gplot.layout.",mode,sep="")),silent=TRUE)
     if(class(layout.fun)=="try-error")
       stop("Error in gplot: no layout function for mode",mode)
     temp<-layout.fun(d,layout.par)
     x<-temp[,1]
     y<-temp[,2]
   }
   #Jitter the coordinates if need be
   if(jitter){
      x<-jitter(x)
      y<-jitter(y)
   }
   #Which nodes should we use?
   use<-displayisolates|(!is.isolate(d,ego=1:dim(d)[1]))   
   #Deal with axis labels
   if(is.null(xlab))
     xlab=""
   if(is.null(ylab))
     ylab=""
   #Create the base plot, if needed
   xlim<-c(min(x[use])-pad,max(x[use])+pad)  #Save x, y limits
   ylim<-c(min(y[use])-pad,max(y[use])+pad)
   if(new){  #If new==FALSE, we add to the existing plot; else create a new one
     if((length(x)>0)&(!all(use==FALSE)))
        plot(x[use],y[use],xlim=xlim,ylim=ylim,type="p",pch=vertex.pch,xlab=xlab,ylab=ylab,col=vertex.col,cex=vertex.cex,axes=!suppress.axes,...)
     else
        plot(0,0,type="n",pch=vertex.pch,xlab=expression(lambda[1]),ylab=expression(lambda[2]),col=vertex.col,cex=vertex.cex,axes=!suppress.axes,...)
   }
   #Generate the edges and their attributes
   px0<-vector()   #Create position vectors (tail, head)
   py0<-vector()
   px1<-vector()
   py1<-vector()
   e.lwd<-vector() #Create edge attribute vectors
   e.curv<-vector()
   e.type<-vector()
   e.col<-vector()
   e.diag<-vector() #Indicator for self-ties
   if(!is.array(edge.col))   #Coerce edge.col/edge.lty to array form
     edge.col<-array(edge.col,dim=dim(d))
   if(!is.array(edge.lty))
     edge.lty<-array(edge.lty,dim=dim(d))
   dist<-as.matrix(dist(cbind(x,y)))  #Get the inter-point distances for curves
   tl<-d.raw*dist   #Get rescaled edge lengths
   tl.max<-max(tl)  #Get maximum edge length   
   for(i in (1:n)[use])    #Plot edges for displayed vertices
     for(j in (1:n)[use])
       if(d[i,j]){       #Perform for actually existing edges
         px0<-c(px0,as.real(x[i]))  #Store endpoint coordinates
         py0<-c(py0,as.real(y[i]))
         px1<-c(px1,as.real(x[j]))
         py1<-c(py1,as.real(y[j]))
         e.col<-c(e.col,edge.col[i,j])    #Store other edge attributes
         e.type<-c(e.type,edge.lty[i,j])
         if(!is.array(edge.lwd)){
           if(edge.lwd>0)
             e.lwd<-c(e.lwd,edge.lwd*d.raw[i,j])
           else
             e.lwd<-c(e.lwd,1)
         }else
           e.lwd<-c(e.lwd,edge.lwd[i,j])
         e.diag<-c(e.diag,i==j)  #Is this a loop?
         if(uselen){   #Should we base curvature on interpoint distances?
           if(tl[i,j]>0){ 
             e.len<-dist[i,j]*tl.max/tl[i,j]
             e.curv<-c(e.curv,edge.len*sqrt((e.len/2)^2-(dist[i,j]/2)^2))
           }else{      
             e.curv<-c(e.curv,0)   
           }
         }else{        #Otherwise, use prespecified edge.curve
           if(!is.array(edge.curve)){
             if(!is.null(edge.curve))  #If it's a scalar, multiply by edge str
               e.curv<-c(e.curv,edge.curve*d.raw[i,j])
             else
               e.curv<-c(e.curv,0)
           }else{
            e.curv<-c(e.curv,edge.curve[i,j])
           }
         }
       }
   #Plot loops for the diagonals, if diag==TRUE
   if(diag&&(length(px0)>0)&&sum(e.diag>0)){  #Are there any loops present?
     for(i in (1:length(px0))[e.diag]){  #Only walk the loops
       xctr<-(xlim[2]+xlim[1])/2  #Find the center of the plot
       yctr<-(ylim[2]+ylim[1])/2
       xoff<-(px0[i]-xctr)/(xlim[2]-xctr)  #Create offsets
       yoff<-(py0[i]-yctr)/(ylim[2]-yctr)
       hoff<-sqrt(xoff^2+yoff^2)
       dx<-diag.size*(xlim[2]-xlim[1])
       dy<-diag.size*(ylim[2]-ylim[1])
       xoff<-dx*xoff/hoff
       yoff<-dy*yoff/hoff
       t2<-edge.steps
       for(t in 0:(t2-1)){                #Walk through the curve
         segments(xoff+px0[i]+dx*sin(2*pi*t/t2), yoff+py0[i]+dy*cos(2*pi*t/t2), xoff+px1[i]+dx*sin(2*pi*(t+1)/t2), yoff+py1[i]+dy*cos(2*pi*(t+1)/t2), col=e.col[i], lty=e.type[i], lwd=e.lwd[i])
       }
       if(usearrows){                     #Add arrowheads, if needed
         arrows(xoff+px0[i]+dx*sin(2*pi*t/t2), yoff+py0[i]+dy*cos(2*pi*t/t2), xoff+px1[i]+dx*sin(2*pi*(t+1)/t2), yoff+py1[i]+dy*cos(2*pi*(t+1)/t2), col=e.col[i], lty=e.type[i], lwd=e.lwd[i], length=arrowhead.length,angle=15)
       }
     }
   }
   #Plot standard (i.e., non-loop) edges
   if(length(px0)>0){  #If edges are present, remove loops from consideration
     px0<-px0[!e.diag] 
     py0<-py0[!e.diag]
     px1<-px1[!e.diag]
     py1<-py1[!e.diag]
     e.curv<-e.curv[!e.diag]
     e.lwd<-e.lwd[!e.diag]
     e.type<-e.type[!e.diag]
     e.col<-e.col[!e.diag]
   }
   if(!usecurve&!uselen){   #Straight-line edge case
     if(usearrows&(length(px0)>0))
        arrows(as.vector(px0),as.vector(py0),as.vector(px1),as.vector(py1),length=arrowhead.length,angle=15,col=e.col,lty=e.type,lwd=e.lwd)
     else if(length(px0)>0)
        segments(as.vector(px0),as.vector(py0),as.vector(px1),as.vector(py1),col=e.col,lty=e.type,lwd=e.lwd)
   }else{   #Curved edge case
     if(length(px0)>0){
       for(i in 1:length(px0)){
         h<-e.curv[i]
         t2<-edge.steps
         dx<-(px1[i]-px0[i])/t2                 #Compute incremental changes
         dy<-(py1[i]-py0[i])/t2
         dzx<--sign(dy)*abs(dy/sqrt(dx^2+dy^2)) #Compute normals
         dzy<-dx/sqrt(dx^2+dy^2)
         dz1<-0
         for(t in 0:(t2-1)){                    #Walk through the curve
           dz2<-(h-(h/(t2/2)^2)*(t+1-(t2/2))^2)
           segments(px0[i]+dx*t+dz1*dzx, py0[i]+dy*t+dz1*dzy, px0[i]+dx*(t+1)+dz2*dzx, py0[i]+dy*(t+1)+dz2*dzy, col=e.col[i], lty=e.type[i], lwd=e.lwd[i])
           if(t<(t2-1))
             dz1<-dz2
         }
         if(usearrows)                         #Add arrowheads, if needed
           arrows(px0[i]+dx*t+dz1*dzx, py0[i]+dy*t+dz1*dzy, px0[i]+dx*(t+1)+dz2*dzx, py0[i]+dy*(t+1)+dz2*dzy, col=e.col[i], lty=e.type[i], lwd=e.lwd[i], length=arrowhead.length, angle=15)
       }
       
     }
   }
   #Plot vertex labels, if needed
   if((!all(label==""))&(!all(use==FALSE))){
      if(boxed.labels){
        lw<-strwidth(label[use],cex=label.cex)/2
        lh<-strheight(label[use],cex=label.cex)/2
	os<-c(0.2,0.4)*par()$cxy*label.cex
        rect(x[use]-lw-os[1],y[use]-3.5*lh-os[2],x[use]+lw+os[1],y[use]-2.5*lh+os[2],col="white")
      }	
      text(x[use],y[use]-3*lh,label[use],cex=label.cex,col=label.col)
   }
   #If interactive, allow the user to mess with things
   if(interactive&&((length(x)>0)&&(!all(use==FALSE)))){
     #Set up the text offset increment
     os<-c(0.2,0.4)*par()$cxy
     #Get the location for text messages, and write to the screen
     textloc<-c(min(x[use])-pad,max(y[use])+pad)
     tm<-"Select a vertex to move, or click \"Finished\" to end."
     tmh<-strheight(tm)
     tmw<-strwidth(tm)
     text(textloc[1],textloc[2],tm,adj=c(0,0.5)) #Print the initial instruction
     fm<-"Finished"
     finx<-c(textloc[1],textloc[1]+strwidth(fm))
     finy<-c(textloc[2]-3*tmh-strheight(fm)/2,textloc[2]-3*tmh+strheight(fm)/2)
     finbx<-finx+c(-os[1],os[1])
     finby<-finy+c(-os[2],os[2])
     rect(finbx[1],finby[1],finbx[2],finby[2],col="white")
     text(finx[1],mean(finy),fm,adj=c(0,0.5))
     #Get the click location
     clickpos<-unlist(locator(1))
     #If the click is in the "finished" box, end our little game.  Otherwise,
     #relocate a vertex and redraw.
     if((clickpos[1]%iin%finbx)&&(clickpos[2]%iin%finby)){
       cl<-match.call()                #Get the args of the current function
       cl$interactive<-FALSE           #Turn off interactivity
       cl$coord<-cbind(x,y)            #Set the coordinates
       cl$dat<-dat                     #"Fix" the data array
       return(eval(cl))     #Execute the function and return
     }else{
       #Figure out which vertex was selected
       clickdis<-sqrt((clickpos[1]-x[use])^2+(clickpos[2]-y[use])^2)
       selvert<-match(min(clickdis),clickdis)
       #Create usable labels, if the current ones aren't
       if(all(label==""))
         label<-1:n
       #Clear out the old message, and write a new one
       rect(textloc[1],textloc[2]-tmh/2,textloc[1]+tmw,textloc[2]+tmh/2,border="white",col="white")
       tm<-"Where should I move this vertex?"
       tmh<-strheight(tm)
       tmw<-strwidth(tm)
       text(textloc[1],textloc[2],tm,adj=c(0,0.5))
       fm<-paste("Vertex",label[use][selvert],"selected")
       finx<-c(textloc[1],textloc[1]+strwidth(fm))
       finy<-c(textloc[2]-3*tmh-strheight(fm)/2,textloc[2]-3*tmh+strheight(fm)/2)
       finbx<-finx+c(-os[1],os[1])
       finby<-finy+c(-os[2],os[2])
       rect(finbx[1],finby[1],finbx[2],finby[2],col="white")
       text(finx[1],mean(finy),fm,adj=c(0,0.5))
       #Get the destination for the new vertex
       clickpos<-unlist(locator(1))
       #Set the coordinates accordingly
       x[use][selvert]<-clickpos[1]
       y[use][selvert]<-clickpos[2]
       #Iterate (leaving interactivity on)
       cl<-match.call()                #Get the args of the current function
       cl$coord<-cbind(x,y)            #Set the coordinates
       cl$dat<-dat                     #"Fix" the data array
       return(eval(cl))     #Execute the function and return
     }
   }
   #Return the vertex positions, should they be needed
   invisible(cbind(x,y))
}

#gplot.layout.* - Layout functions for gplot
#gplot.layout.princoord - Place using the eigenstructure of the correlation 
#matrix among concatenated rows/columns (principal coordinates by position
#similarity)
gplot.layout.princoord<-function(d,layout.par){     
  #Determine the vectors to be related
  if(is.null(layout.par$var))
    vm<-rbind(d,t(d))
  else
    vm<-switch(layout.par$var,
      rowcol=rbind(d,t(d)),
      col=d,
      row=t(d),
      rcsum=d+t(d),
      rcdiff=d-t(d),
      user=layout.par$vm
    )
  #Find the correlation/covariance matrix
  if(is.null(layout.par$cor)||layout.par$cor)
    cd<-cor(vm,use="pairwise.complete.obs")
  else    
    cd<-cov(vm,use="pairwise.complete.obs")
  cd<-replace(cd,is.na(cd),0)
  #Obtain the eigensolution
  e<-eigen(cd,symmetric=TRUE)
  x<-Re(e$vectors[,1])
  y<-Re(e$vectors[,2])
  cbind(x,y)
}
#gplot.layout.eigen - Place vertices based on the first two eigenvectors of
#an adjacency matrix
gplot.layout.eigen<-function(d,layout.par){     
  #Determine the matrix to be used
  if(is.null(layout.par$var))
    vm<-d
  else
    vm<-switch(layout.par$var,
      symupper=symmetrize(d,rule="uppper"),
      symlower=symmetrize(d,rule="lower"),
      symstrong=symmetrize(d,rule="strong"),
      symweak=symmetrize(d,rule="weak"),
      user=layout.par$mat,
      raw=d
    )
  #Pull the eigenstructure
  e<-eigen(vm)
  if(is.null(layout.par$evsel))
    coord<-Re(e$vectors[,1:2])
  else
    coord<-switch(layout.par$evsel,
      first=Re(e$vectors[,1:2]),
      size=Re(e$vectors[,rev(order(abs(e$values)))[1:2]])
    )
  #Return the result
  coord
}
#gplot.layout.mds - Place vertices based on metric multidimensional scaling
#of a distance matrix
gplot.layout.mds<-function(d,layout.par){     
  #Determine the raw inputs for the scaling
  if(is.null(layout.par$var))
    vm<-cbind(d,t(d))
  else
    vm<-switch(layout.par$var,
      rowcol=cbind(d,t(d)),
      col=t(d),
      row=d,
      rcsum=d+t(d),
      rcdiff=t(d)-d,
      invadj=max(d)-d,
      geodist=geodist(d)$gdist,
      user=layout.par$vm
    )
  #If needed, construct the distance matrix
  if(is.null(layout.par$dist))
    dm<-as.matrix(dist(vm))
  else
    dm<-switch(layout.par$dist,
      euclidean=as.matrix(dist(vm)),
      maximum=as.matrix(dist(vm,method="maximum")),
      manhattan=as.matrix(dist(vm,method="manhattan")),
      canberra=as.matrix(dist(vm,method="canberra")),
      none=vm
    )
  #Transform the distance matrix, if desired
  if(is.null(layout.par$exp))
    dm<-dm^2
  else
    dm<-dm^layout.par$exp
  #Perform the scaling and return
  cmdscale(dm,2)
}
gplot.layout.rmds<-function(d,layout.par){
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$var="row"
  layout.par$dist="euclidean"
  layout.par$exp=1
  gplot.layout.mds(d,layout.par)
}
gplot.layout.geodist<-function(d,layout.par){
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$var="geodist"
  layout.par$dist="none"
  layout.par$exp=1
  gplot.layout.mds(d,layout.par)
}
gplot.layout.adj<-function(d,layout.par){
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$var="invadj"
  layout.par$dist="none"
  layout.par$exp=1
  gplot.layout.mds(d,layout.par)
}
gplot.layout.seham<-function(d,layout.par){
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$var="rowcol"
  layout.par$dist="manhattan"
  layout.par$exp=1
  gplot.layout.mds(d,layout.par)
}
gplot.layout.segeo<-function(d,layout.par){
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$var="geodist"
  layout.par$dist="euclidean"
  gplot.layout.mds(d,layout.par)
}
#gplot.layout.random
gplot.layout.random<-function(d,layout.par){     
  n<-dim(d)[1]
  #Determine the distribution
  if(is.null(layout.par$dist))
    temp<-matrix(runif(2*n,-1,1),n,2)
  else if (layout.par$dist=="unif")
    temp<-matrix(runif(2*n,-1,1),n,2)
  else if (layout.par$dist=="uniang"){
    tempd<-rnorm(n,1,0.25)
    tempa<-runif(n,0,2*pi)
    temp<-cbind(tempd*sin(tempa),tempd*cos(tempa))
  }else if (layout.par$dist=="normal")
    temp<-matrix(rnorm(2*n),n,2)
  #Return the result
  temp
}
gplot.layout.circrand<-function(d,layout.par){ 
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$dist="uniang"
  gplot.layout.random(d,layout.par)
}
#gplot.layout.circle - Place vertices in a circular layout
gplot.layout.circle<-function(d,layout.par){ 
  n<-dim(d)[1]
  cbind(sin(2*pi*((0:(n-1))/n)),cos(2*pi*((0:(n-1))/n)))
}
#gplot.layout.spring - Place vertices using a spring embedder
gplot.layout.spring<-function(d,layout.par){
  #Set up the embedder params
  ep<-vector()
  if(is.null(layout.par$mass))  #Mass is in "quasi-kilograms"
    ep[1]<-0.1
  else
    ep[1]<-layout.par$mass
  if(is.null(layout.par$equil)) #Equilibrium extension is in "quasi-meters"
    ep[2]<-1
  else
    ep[2]<-layout.par$equil
  if(is.null(layout.par$k)) #Spring coefficient is in "quasi-Newtons/qm"
    ep[3]<-0.001
  else
    ep[3]<-layout.par$k
  if(is.null(layout.par$repeqdis)) #Repulsion equilibrium is in qm
    ep[4]<-0.1
  else
    ep[4]<-layout.par$repeqdis
  if(is.null(layout.par$kfr)) #Base coef of kinetic friction is in qn-qkg
    ep[5]<-0.01
  else
    ep[5]<-layout.par$kfr
  if(is.null(layout.par$repulse))
    repulse<-FALSE
  else
    repulse<-layout.par$repulse
  #Create initial condidions
  n<-dim(d)[1]
  f.x<-rep(0,n)       #Set initial x/y forces to zero
  f.y<-rep(0,n)
  v.x<-rep(0,n)       #Set initial x/y velocities to zero
  v.y<-rep(0,n)
  tempa<-sample((0:(n-1))/n) #Set initial positions randomly on the circle
  x<-n/(2*pi)*sin(2*pi*tempa)
  y<-n/(2*pi)*cos(2*pi*tempa)
  ds<-symmetrize(d,"weak")            #Symmetrize/dichotomize the graph
  kfr<-ep[5]                          #Set initial friction level
  niter<-1                            #Set the iteration counter
  #Simulate, with increasing friction, until motion stops    
  repeat{
    niter<-niter+1                    #Update the iteration counter
    dis<-as.matrix(dist(cbind(x,y)))  #Get inter-point distances
    #Get angles relative to the positive x direction
    theta<-acos(t(outer(x,x,"-"))/dis)*sign(t(outer(y,y,"-"))) 
    #Compute spring forces; note that we assume a base spring coefficient
    #of ep[3] units ("pseudo-Newtons/quasi-meter"?), with an equilibrium
    #extension of ep[2] units for all springs
    f.x<-apply(ds*cos(theta)*ep[3]*(dis-ep[2]),1,sum,na.rm=TRUE)
    f.y<-apply(ds*sin(theta)*ep[3]*(dis-ep[2]),1,sum,na.rm=TRUE)
    #If node repulsion is active, add a force component for this
    #as well.  We employ an inverse cube law which is equal in power
    #to the attractive spring force at distance ep[4]
    if(repulse){
      f.x<-f.x-apply(cos(theta)*ep[3]/(dis/ep[4])^3,1,sum,na.rm=TRUE)
      f.y<-f.y-apply(sin(theta)*ep[3]/(dis/ep[4])^3,1,sum,na.rm=TRUE)
    }
    #Adjust the velocities (assume a mass of ep[1] units); note that the
    #motion is roughly modeled on the sliding of flat objects across
    #a uniform surface (e.g., spring-connected cylinders across a table).
    #We assume that the coefficients of static and kinetic friction are
    #the same, which should only trouble you if you are under the 
    #delusion that this is a simulation rather than a graph drawing
    #exercise (in which case you should be upset that I'm not using
    #Runge-Kutta or the like!).
    v.x<-v.x+f.x/ep[1]         #Add accumulated spring/repulsion forces
    v.y<-v.y+f.y/ep[1]
    spd<-sqrt(v.x^2+v.y^2)     #Determine frictional forces
    fmag<-pmin(spd,kfr)  #We can't let friction _create_ motion!
    theta<-acos(v.x/spd)*sign(v.y)  #Calculate direction of motion
    f.x<-fmag*cos(theta)        #Decompose frictional forces
    f.y<-fmag*sin(theta)
    f.x[is.nan(f.x)]<-0         #Correct for any 0/0 problems
    f.y[is.nan(f.y)]<-0
    v.x<-v.x-f.x                #Apply frictional forces (opposing motion -
    v.y<-v.y-f.y                #note that mass falls out of equation)
    #Adjust the positions (yep, it's primitive linear updating time!)
    x<-x+v.x
    y<-y+v.y
    #Check for cessation of motion, and increase friction
    mdist<-mean(dis)
    if(all(v.x<mdist*1e-5)&&all(v.y<mdist*1e-5))
      break
    else
      kfr<-ep[5]*exp(0.1*niter)
  }
  #Return the result
  cbind(x,y)
}
gplot.layout.springrepulse<-function(d,layout.par){
  if(is.null(layout.par))
    layout.par<-list()
  layout.par$repulse<-TRUE
  gplot.layout.spring(d,layout.par)
}

#netlogit - God help me, it's a network regression routine using a binomial/logit GLM.  It's also
#frighteningly slow, since it's essentially a front end to the builtin GLM routine with a bunch of
#network hypothesis testing stuff thrown in for good measure.

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

summary.netlogit<-function(object, ...){
   out<-object
   class(out)<-c("summary.netlogit",class(out))
   out
}

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

summary.netlm<-function(object, ...){
   out<-object
   class(out)<-c("summary.netlm",class(out))
   out
}

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


#netcancor - Canonical correlations for network variables.  NOTE: requires that mva be loaded.

netcancor<-function(y,x,mode="digraph",diag=FALSE,nullhyp="cugtie",reps=1000){
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

summary.netcancor<-function(object, ...){
   out<-object
   class(out)<-c("summary.netcancor",class(out))
   out
}

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


#centralgraph - Find the central graph of a graph stack

centralgraph<-function(dat,normalize=FALSE){
   #Check to see if someone foolishly called this with one graph
   if(length(dim(dat))==2)
      out<-dat
   else{
      if(normalize)
         out<-apply(dat,c(2,3),mean,na.rm=TRUE)
      else
         out<-matrix(data=as.numeric(apply(dat,c(2,3),mean,na.rm=TRUE)>=0.5),nrow=dim(dat)[2],ncol=dim(dat)[2])
   }
   out
}


#plot.sociomatrix - An odd sort of plotting routine; plots a matrix (e.g., a Bernoulli graph density, or a set of
#adjacencies) as an image.  Very handy for visualizing large valued matrices...

plot.sociomatrix<-function(x,labels=list(seq(1:dim(x)[1]),seq(1:dim(x)[2])),drawlab=TRUE,diaglab=TRUE,...){       
   n<-dim(x)[1]
   o<-dim(x)[2]
   d<-1-(x-min(x))/(max(x)-min(x))
   plot(1,1,xlim=c(0,o+1),ylim=c(n+1,0),type="n",axes=FALSE,...)
   for(i in 1:n)
      for(j in 1:o)
         rect(j-0.5,i+0.5,j+0.5,i-0.5,col=gray(d[i,j]),xpd=TRUE)
   if(drawlab){
      text(rep(0,n),1:n,labels[[1]])
      text(1:o,rep(0,o),labels[[2]])
   }
   if((n==o)&(drawlab)&(diaglab))
      if(labels[[1]]==labels[[2]])
         text(1:o,1:n,labels[[1]])
}


#bbnam.bf - Estimate Bayes Factors for the Butts Bayesian Network Accuracy Model.  This implementation
#relies on monte carlo integration to estimate the BFs, and tests the fixed probability, pooled, and pooled by
#actor models.

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

bbnam.jntlik<-function(dat,log=FALSE,...){
   p<-sum(sapply(1:dim(dat)[1],bbnam.jntlik.slice,dat=dat,log=TRUE,...))
   if(!log)
      exp(p)
   else
      p
}

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


#print.bayes.factor - A fairly generic routine for printing bayes factors, here used for the bbnam routine.

print.bayes.factor<-function(x,...){
   tab<-x$int.lik
   rownames(tab)<-x$model.names
   colnames(tab)<-x$model.names
   cat("Bayes Factors by Model:\n\n(Diagonals indicate raw integrated likelihood estimates.)\n\n")
   print(tab)
   cat("\n")   
}


#summary.bayes.factor - A fairly generic summary routine for bayes factors.  Clearly, this belongs in
#some other library than sna, but for the moment this will have to do...

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


#potscalered.mcmc - Potential scale reduction (sqrt(Rhat)) for scalar estimands.  Input must be a 
#matrix whose columns correspond to replicate chains.  This, clearly, doesn't belong here, but lacking
#a better place to put it I have included it nonetheless.

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


#bbnam - Draw from Butts' Bayesian Network Accuracy Model.  This version uses a Gibbs' Sampler,
#and assumes error rates to be drawn from conditionally independent betas for each actor.  Note
#that dat MUST be an n x n x n array, and that the data in question MUST be dichotomous.  Priors
#are also assumed to be in the right form (n x 2 matrices of alpha, beta pairs for em and ep, and
#an n x n probability matrix for the network itself), and are not checked; default behavior if no
#priors are provided is the uninformative case.

#Probability of a given tie
bbnam.probtie<-function(dat,i,j,npriorij,em,ep){
   num<-npriorij
   denom<-1-npriorij
   num<-num*prod(dat[,i,j]*(1-em)+(1-dat[,i,j])*em,na.rm=TRUE)
   denom<-denom*prod(dat[,i,j]*ep+(1-dat[,i,j])*(1-ep),na.rm=TRUE)
   p<-num/(denom+num)
   p
}

#Wrapper function for the various bbnam models
bbnam<-function(dat,model="actor",...){
   if(model=="actor")
      bbnam.actor(dat,...)
   else if(model=="pooled")
      bbnam.pooled(dat,...)
   else if(model=="fixed")
      bbnam.fixed(dat,...)
}

#Draw from the fixed probability error model
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

#Draw from the pooled error model
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


#Draw from the error-prob-by-actor model
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

print.bbnam<-function(x,...){
   UseMethod("print",x)
}

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

summary.bbnam<-function(object, ...){
   out<-object
   class(out)<-c("summary.bbnam",class(out))
   out
}

summary.bbnam.fixed<-function(object, ...){
   out<-object
   class(out)<-c("summary.bbnam.fixed",class(out))
   out
}

summary.bbnam.pooled<-function(object, ...){
   out<-object
   class(out)<-c("summary.bbnam.pooled",class(out))
   out
}

summary.bbnam.actor<-function(object, ...){
   out<-object
   class(out)<-c("summary.bbnam.actor",class(out))
   out
}

print.summary.bbnam<-function(x,...){
   UseMethod("print",x)
}

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

plot.bbnam<-function(x,mode="density",intlines=TRUE,...){
   UseMethod("plot",x)
}

plot.bbnam.fixed<-function(x,mode="density",intlines=TRUE,...){
   #Get the initial graphical settings, so we can restore them later
   oldpar<-par()
   #Perform matrix plot of tie probabilities
   par(mfrow=c(1,1))
   plot.sociomatrix(apply(x$net,c(2,3),mean),labels=list(x$anames,x$anames),main="Marginal Posterior Tie Probability Distribution")
   #Clean up
   par(oldpar)
}

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


#npostpred - Take posterior predictive draws for functions of networks.

npostpred<-function(b,FUN,...){
   #Find the desired function
   fun<-match.fun(FUN)
   #Take the draws
   out<-apply(b$net,1,fun,...)
   out
}


#sr2css - Convert a row-wise self-report matrix to a CSS matrix with missing observations.

sr2css<-function(net){
   n<-dim(net)[1]
   css<-array(dim=c(n,n,n))
   for(i in 1:n){
      css[i,,]<-NA
      css[i,i,]<-net[i,]
   }
   css
}


#event2dichot - Convert an observed event matrix to a dichotomous matrix.  Methods are quantile,
#mean, rquantile, rmean, cquantile, cmean, absolute, rank, rrank, and crank.  Thresh specifies the
#cutoff, in terms of whatever method is used (if applicable).

event2dichot<-function(m,method="quantile",thresh=0.5,leq=FALSE){
   if(method=="quantile"){
      q<-quantile(m,thresh,na.rm=TRUE, names=FALSE)
      out<-as.numeric(m>q)
   } else if(method=="rquantile"){
      q<-quantile(m[1,],thresh,na.rm=TRUE, names=FALSE)
      out<-as.numeric(m[1,]>q)
      for(i in 2:dim(m)[1]){      
         q<-quantile(m[i,],thresh,na.rm=TRUE, names=FALSE)
         out<-rbind(out,as.numeric(m[i,]>q))
      }
   } else if(method=="cquantile"){
      q<-quantile(m[,1],thresh,na.rm=TRUE, names=FALSE)
      out<-as.numeric(m[,1]>q)
      for(i in 2:dim(m)[2]){      
         q<-quantile(m[,i],thresh,na.rm=TRUE, names=FALSE)
         out<-cbind(out,as.numeric(m[,i]>q))
      }
   } else if(method=="mean"){
      q<-mean(m)
      out<-as.numeric(m>q)
   } else if(method=="rmean"){
      q<-mean(m[1,])
      out<-as.numeric(m[1,]>q)
      for(i in 2:dim(m)[1]){      
         q<-mean(m[i,])
         out<-rbind(out,as.numeric(m[i,]>q))
      }
   } else if(method=="cmean"){
      q<-mean(m[,1])
      out<-as.numeric(m[,1]>q)
      for(i in 2:dim(m)[2]){      
         q<-mean(m[,i])
         out<-rbind(out,as.numeric(m[,i]>q))
      }
   } else if(method=="absolute"){
      out<-as.numeric(m>thresh)
   } else if(method=="rank"){
      o<-order(m)
      out<-as.numeric((max(o)-o+1)<thresh)
   } else if(method=="rrank"){
      o<-order(m[1,])
      out<-as.numeric((max(o)-o+1)<thresh)
      for(i in 2:dim(m)[1]){      
         o<-order(m[i,])
         out<-rbind(out,as.numeric((max(o)-o+1)<thresh))
      }
   } else if(method=="crank"){
      o<-order(m[,1])
      out<-as.numeric((max(o)-o+1)<thresh)
      for(i in 2:dim(m)[2]){      
         o<-order(m[,i])
         out<-cbind(out,as.numeric((max(o)-o+1)<thresh))
      }
   }
   if(leq==TRUE)
      out<-1-out
   if(is.null(dim(out))!=is.null(dim(m)))
      out<-array(out,dim=dim(m))
   else
      if(dim(out)!=dim(m))
         out<-array(out,dim=dim(m))
   out
}


#symmetrize - Convert a graph or graph stack to a symmetric form.  Current rules for symmetrizing
#include "upper" and "lower" diagonals, "weak" connectedness rule, and a "strong" connectedness rule.

symmetrize<-function(mats,rule="weak"){
   #Build the input data structures
   if(length(dim(mats))>2){
      m<-dim(mats)[1]
      n<-dim(mats)[2]
      o<-dim(mats)[3]
      d<-mats
   }else{
      m<-1
      n<-dim(mats)[1]
      o<-dim(mats)[2]
      d<-array(dim=c(1,n,o))
      d[1,,]<-mats
   }
   #Apply the symmetry rule
   for(i in 1:m){
      if(rule=="upper"){
         temp<-d[i,,]
         for(j in 1:n)
            temp[j:n,j]<-temp[j,j:n]
         d[i,,]<-temp
      }else if(rule=="lower"){
         temp<-d[i,,]
         for(j in 1:n)
            temp[j,j:n]<-temp[j:n,j]
         d[i,,]<-temp
      }else if(rule=="weak"){
         d[i,,]<-matrix(as.numeric(d[i,,]|t(d[i,,])),nrow=n,ncol=o)
      }else if(rule=="strong"){
         d[i,,]<-matrix(as.numeric(d[i,,]&t(d[i,,])),nrow=n,ncol=o)
      }
   }
   #Return the symmetrized matrix
   if(m==1)
      out<-d[1,,]
   else
      out<-d
   out
}


#consensus - Find a consensus structure, using one of several algorithms.  Note that this is currently
#experimental, and that the routines are not guaranteed to produce meaningful output

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


#sedist - Find a matrix of distances between positions based on structural equivalence

sedist<-function(dat,g=c(1:dim(dat)[1]),method="hamming",joint.analysis=FALSE,mode="digraph",diag=FALSE,code.diss=FALSE){
   #First, prepare the data
   if(length(dim(dat))>2){
      n<-dim(dat)[2]
      m<-length(g)
      d<-dat[g,,]
   }else{
      n<-dim(dat)[2]
      m<-1
      d<-array(dim=c(m,n,n))
      d[1,,]<-dat
   }
   if(mode=="graph")
      d<-upper.tri.remove(d)
   if(!diag)
      d<-diag.remove(d)
   #Are we conducting a joint analysis?
   if(joint.analysis){
      o<-array(dim=c(1,n,n))
      #Build the data matrix
      v<-vector()
      for(i in 1:n)
         v<-cbind(v,c(as.vector(d[,i,]),as.vector(d[,,i])))
      #Proceed by method
      if(method=="correlation"){
         o[1,,]<-cor(v,use="pairwise")
         #Reverse code?
         if(code.diss)
            o<--o
      }else if(method=="euclidean"){
         for(i in 1:n)
            for(j in 1:n)
               o[1,i,j]<-sqrt(sum((v[,i]-v[,j])^2,na.rm=TRUE))         
      }else if(method=="hamming"){
         for(i in 1:n)
            for(j in 1:n)
               o[1,i,j]<-sum(abs(v[,i]-v[,j]),na.rm=TRUE)
      }else if(method=="gamma"){
         for(i in 1:n)
            for(j in 1:n){
               concord<-sum(as.numeric(v[,i]==v[,j]),na.rm=TRUE)
               discord<-sum(as.numeric(v[,i]!=v[,j]),na.rm=TRUE)
               o[1,i,j]<-(concord-discord)/(concord+discord)
            }                  
         #Reverse code?
         if(code.diss)
            o<--o
      }else if(method=="exact"){
         for(i in 1:n)
            for(j in 1:n)
               o[1,i,j]<-as.numeric(any(v[!(is.na(v[,i])|is.na(v[,j])),i]!=v[!(is.na(v[,i])|is.na(v[,j])),j]))
      }
   }else{  #Analyze each graph seperately
      o<-array(dim=c(m,n,n))
      for(k in 1:m){
         #Build the data matrix
         v<-vector()
         for(i in 1:n)
            v<-cbind(v,c(as.vector(d[k,i,]),as.vector(d[k,,i])))
         #Proceed by method
         if(method=="correlation"){
            o[k,,]<-cor(v,use="pairwise")
            o[k,,][is.na(o[k,,])]<-0
            #Reverse code?
            if(code.diss)
               o[k,,]<--o[k,,]
         }else if(method=="euclidean"){
            for(i in 1:n)
               for(j in 1:n)
                  o[k,i,j]<-sqrt(sum((v[,i]-v[,j])^2,na.rm=TRUE))         
         }else if(method=="hamming"){
            for(i in 1:n)
               for(j in 1:n)
                  o[k,i,j]<-sum(abs(v[,i]-v[,j]),na.rm=TRUE)
         }else if(method=="gamma"){
            for(i in 1:n)
               for(j in 1:n){
                  concord<-sum(as.numeric(v[,i]==v[,j]),na.rm=TRUE)
                  discord<-sum(as.numeric(v[,i]!=v[,j]),na.rm=TRUE)
                  o[k,i,j]<-(concord-discord)/(concord+discord)
               }                  
            #Reverse code?
            if(code.diss)
               o[k,,]<--o[k,,]
         }else if(method=="exact"){
            for(i in 1:n)
               for(j in 1:n)
                  o[k,i,j]<-as.numeric(any(v[!(is.na(v[,i])|is.na(v[,j])),i]!=v[!(is.na(v[,i])|is.na(v[,j])),j]))
         }
      }
   }
   #Dump the output
   if(dim(o)[1]==1)
      as.matrix(o[1,,])
   else
      o
}


#equiv.clust - Find clusters of positions based on an equivalence relation

equiv.clust<-function(dat,g=c(1:dim(dat)[1]),equiv.fun="sedist",method="hamming",mode="digraph",diag=FALSE,cluster.method="complete",glabels=dimnames(dat)[[1]][g],plabels=dimnames(dat)[[2]],...){
   #First, find the equivalence distances using the appropriate function and method
   equiv.dist.fun<-match.fun(equiv.fun)
   equiv.dist<-equiv.dist.fun(dat,g=g,method=method,joint.analysis=TRUE,mode=mode,diag=diag,code.diss=TRUE,...)
   #Load the mva package, if it's not already loaded
   require(mva)
   #Generate the output object
   o<-list()
   #Produce the hierarchical clustering
   o$cluster<-hclust(as.dist(equiv.dist),method=cluster.method)
   #Set the output class and take care of other details
   o$metric<-method
   o$equiv.fun<-equiv.fun
   o$cluster.method<-cluster.method
   if((length(dim(dat))==1)&(length(glabels)>1))
      o$glabels<-glabels[1]
   else
      o$glabels<-glabels
   o$plabels<-plabels
   class(o)<-"equiv.clust"
   #Produce the output
   o
}


#plot.equiv.clust - Plotting for equivalence clustering objects

plot.equiv.clust<-function(x,labels=x$plabels,...){
   require(mva)
   class(x)<-"hclust"
   if(is.null(labels))
      plot(x$cluster,...)
   else
      plot(x$cluster,labels=labels,...)
}


#blockmodel - Generate blockmodels based on partitions of network positions

blockmodel<-function(dat,ec,k=NULL,h=NULL,block.content="density",plabels=ec$plabels,glabels=ec$glabels,rlabels=NULL,mode="digraph",diag=FALSE){
   require(mva)
   #First, extract the blocks
   b<-cutree(ec$cluster,k,h)
   #Prepare the data
   n<-dim(dat)[2]
   if(length(dim(dat))>2)
      d<-dat
   else{
      d<-array(dim=c(1,n,n))
      d[1,,]<-dat
   }
   if(!diag)
      d<-diag.remove(d)
   #Now, construct a model
   rn<-max(b)
   rm<-dim(d)[1]
   if(is.null(rlabels))      #Add labels for roles if needed
      rlabels<-paste("Block",1:rn)
   bm<-array(dim=c(rm,rn,rn))
   for(i in 1:rm)
      for(j in 1:rn)
         for(k in 1:rn){
            if(block.content=="density")
               bm[i,j,k]<-mean(d[i,b==j,b==k,drop=FALSE],na.rm=TRUE)
            else if(block.content=="meanrowsum"){
               bm[i,j,k]<-mean(apply(d[i,b==j,b==k,drop=FALSE],2,sum,na.rm=TRUE))
            }else if(block.content=="meancolsum"){
               bm[i,j,k]<-mean(apply(d[i,b==j,b==k,drop=FALSE],3,sum,na.rm=TRUE))
            }else if(block.content=="sum"){
               bm[i,j,k]<-sum(d[i,b==j,b==k,drop=FALSE],na.rm=TRUE)
            }else if(block.content=="median"){
               bm[i,j,k]<-median(d[i,b==j,b==k,drop=FALSE],na.rm=TRUE)
            }else if(block.content=="min"){
               bm[i,j,k]<-min(d[i,b==j,b==k,drop=FALSE],na.rm=TRUE)
            }else if(block.content=="max"){
               bm[i,j,k]<-max(d[i,b==j,b==k,drop=FALSE],na.rm=TRUE)
            }else if(block.content=="types"){
               temp<-mean(d[i,b==j,b==k,drop=FALSE],na.rm=TRUE)
               if(is.nan(temp))    #Is this a nan block (due to having only one actor)?
                  bm[i,j,k]<-"NA"
               else if(temp==0)    #Is this a null block?
                  bm[i,j,k]<-"null"
               else if(temp==1)   #How about a complete block?
                  bm[i,j,k]<-"complete"
               else if(all(apply(d[i,b==j,b==k,drop=FALSE],2,sum,na.rm=TRUE)>0,apply(d[i,b==j,b==k,drop=FALSE],3,sum,na.rm=TRUE)>0))
                  bm[i,j,k]<-"1 covered"   #1 covered block
               else if(all(apply(d[i,b==j,b==k,drop=FALSE],2,sum,na.rm=TRUE)>0))
                  bm[i,j,k]<-"1 row-covered"   #1 row-covered block
               else if(all(apply(d[i,b==j,b==k,drop=FALSE],3,sum,na.rm=TRUE)>0))
                  bm[i,j,k]<-"1 col-covered"   #1 col-covered block
               else
                  bm[i,j,k]<-"other"   #other block
            }
         }
   #Prepare the output object
   o<-list()
   o$block.membership<-b[ec$cluster$order]
   o$order.vector<-ec$cluster$order
   o$block.content<-block.content
   if(length(dim(dat))>2){
      o$blocked.data<-dat[,ec$cluster$order,ec$cluster$order]
      dimnames(o$blocked.data)<-list(glabels,plabels[ec$cluster$order],plabels[ec$cluster$order])
   }else{
      o$blocked.data<-dat[ec$cluster$order,ec$cluster$order]      
      dimnames(o$blocked.data)<-list(plabels[ec$cluster$order],plabels[ec$cluster$order])
   }
   if(dim(bm)[1]==1){
      o$block.model<-bm[1,,]
      rownames(o$block.model)<-rlabels
      colnames(o$block.model)<-rlabels
   }else{
      o$block.model<-bm
      dimnames(o$block.model)<-list(glabels,rlabels,rlabels)
   }
   o$plabels<-plabels[ec$cluster$order]
   o$glabels<-glabels
   o$rlabels<-rlabels
   o$cluster.method<-ec$cluster.method
   o$equiv.fun<-ec$equiv.fun
   o$equiv.metric<-ec$metric
   class(o)<-"blockmodel"
   o   
}


#print.blockmodel - Printing for blockmodel objects

print.blockmodel<-function(x,...){
   cat("\nBlockmodel by Equivalence Clustering:\n\n")
   cat("Block membership:\n\n")
   if(is.null(x$plabels))                    #Get position labels
      plab<-(1:length(x$block.membership))[x$order.vector]
   else
      plab<-x$plabels
   temp<-matrix(x$block.membership,nrow=1)
   dimnames(temp)<-list("",plab)
   print(temp[1,order(x$order.vector)])  #Print in original order
   cat("\nReduced form blockmodel:\n\n")
   if(length(dim(x$block.model))>2){
      for(i in 1:dim(x$block.model)[1]){
         temp<-x$block.model[i,,]
         dimnames(temp)<-list(x$rlabels,x$rlabels)
         cat("\t",x$glabels[i],"\n") 
         print(temp)
         cat("\n")
      }
   }else{
         temp<-x$block.model
         dimnames(temp)<-list(x$rlabels,x$rlabels)
         cat("\t",x$glabels[i],"\n") 
         print(temp)
   }
}


#summary.blockmodel - Detailed printing for blockmodel objects

summary.blockmodel<-function(object, ...){
   o<-object
   class(o)<-"summary.blockmodel"
   o
}


#print.summary.blockmodel - Printing for blockmodel summary objects

print.summary.blockmodel<-function(x,...){
   cat("\nBlockmodel by Equivalence Clustering:\n\n")

   cat("\nGeneral information:\n\n")
   cat("\tEquivalence function: ",x$equiv.fun,"\n")
   cat("\tEquivalence metric: ",x$equiv.metric,"\n")
   cat("\tClustering method: ",x$cluster.method,"\n")
   cat("\tBlockmodel content: ",x$block.content,"\n")

   cat("\n\nBlock membership by actor:\n\n")
   if(is.null(x$plabels))                    #Get position labels
      plab<-(1:length(x$block.membership))[x$order.vector]
   else
      plab<-x$plabels
   temp<-matrix(x$block.membership,nrow=1)
   dimnames(temp)<-list("",plab)
   print(temp[1,order(x$order.vector)])  #Print in original order

   cat("\n\nBlock membership by block:\n\n")
   for(i in 1:max(x$block.membership))
      cat("\t",x$rlabels[i],":",plab[x$block.membership==i],"\n")
   
   cat("\n\nReduced form blockmodel:\n\n")
   if(length(dim(x$block.model))>2){
      for(i in 1:dim(x$block.model)[1]){
         temp<-x$block.model[i,,]
         dimnames(temp)<-list(x$rlabels,x$rlabels)
         cat("\t",x$glabels[i],"\n") 
         print(temp)
         cat("\n")
      }
   }else{
         temp<-x$block.model
         dimnames(temp)<-list(x$rlabels,x$rlabels)
         cat("\t",x$glabels[i],"\n") 
         print(temp)
   }

   cat("\n\nBlocked data:\n\n")
   if(length(dim(x$block.model))>2){
      for(i in 1:dim(x$block.model)[1]){
         temp<-x$blocked.data[i,,]
         dimnames(temp)<-list(plab,plab)
         cat("\t",x$glabels[i],"\n") 
         print(temp)
         cat("\n")
      }
   }else{
         temp<-x$blocked.data
         dimnames(temp)<-list(plab,plab)
         cat("\t",x$glabels[i],"\n") 
         print(temp)
   }

}


#plot.blockmodel - Plotting for blockmodel objects

plot.blockmodel<-function(x,...){
   #Save old settings
   oldpar<-par()
   #Get new settings from data
   n<-dim(x$blocked.data)[2]
   m<-stackcount(x$blocked.data)
   if(!is.null(x$plabels))
      plab<-x$plabels
   else
      plab<-(1:n)[x$order.vector]
   if(!is.null(x$glabels))
      glab<-x$glabels
   else
      glab<-1:m
   print(glab)
   #Now, plot the blocked data
   par(mfrow=c(floor(sqrt(m)),ceiling(m/floor(sqrt(m)))))
   if(m>1)
      for(i in 1:m){
         plot.sociomatrix(x$blocked.data[i,,],labels=list(plab,plab),main=paste("Relation - ",glab[i]))
         for(j in 2:n)
            if(x$block.membership[j]!=x$block.membership[j-1])
               abline(v=j-0.5,h=j-0.5,lty=3)
      }
   else{
      plot.sociomatrix(x$blocked.data,labels=list(plab,plab),main=paste("Relation - ",glab[1]))
      for(j in 2:n)
         if(x$block.membership[j]!=x$block.membership[j-1])
            abline(v=j-0.5,h=j-0.5,lty=3)
   }
   #Fix the display settings
   par(oldpar)
}


#blockmodel.expand - Generate a graph (or stack) from a given blockmodel using particular expansion rules

blockmodel.expand<-function(b,ev,mode="digraph",diag=FALSE){
   #First, get some useful parameters and such
   en<-sum(ev)
   el<-length(ev)
   bn<-max(b$block.membership)
   bm<-stackcount(b$block.model)
   if(bm>1)
      block.model<-b$block.model
   else{
      block.model<-array(dim=c(1,bn,bn))
      block.model[1,,]<-b$block.model
   }
   #Now, perform the expansion)
   expanded<-array(dim=c(bm,en,en))
   for(i in 1:bm){
      if(b$block.content=="density"){
         tp<-matrix(nrow=en,ncol=en)
         for(j in 1:el)
            for(k in 1:el)
               tp[(cumsum(ev)[j]-ev[j]+1):(cumsum(ev)[j]),(cumsum(ev)[k]-ev[k]+1):(cumsum(ev)[k])]<-block.model[i,j,k]
         tp[is.na(tp)|is.nan(tp)]<-0   #Fill in any NA or NaN blocks with zero
         expanded[i,,]<-rgraph(en,1,tprob=tp,mode=mode,diag=diag)
      }else
         stop(paste("\nContent type",b$block.content,"not supported yet.\n"))
   }
   #Return the output data
   if(dim(expanded)[1]>1)
      expanded
   else
      expanded[1,,]
}


#interval.graph - Construct one or more interval graphs (and exchangeability vectors) from a set of spells

interval.graph<-function(slist,type="simple",diag=FALSE){
   #Note that each slice of slist must have one spell per row, with col 1 containing the spell type,
   #col 2 containing the spell onset, and col 3 containing the spell termination.  If there are multiple
   #slices present, they must be indexed by the first dimension of the array.
   #First, the preliminaries
   o<-list()
   m<-stackcount(slist)          #Get the number of stacks
   if(m==1){
      d<-array(dim=c(m,dim(slist)[1],dim(slist)[2]))
      d[1,,]<-slist
   }else
      d<-slist
   ns<-dim(d)[2]                     #Get the number of spells
   o$exchange.list<-d[,,1]   #Exchange list is just the vector of spell types
   #Now, for the graph itself...
   o$graph<-array(dim=c(m,ns,ns))
   for(i in 1:ns)
      for(j in 1:ns)
         o$graph[,i,j]<-switch(type,
            simple=as.numeric((d[,i,2]<=d[,j,3])&(d[,i,3]>=d[,j,2])),  #"Start before the end, end after the beginning"
            overlap=pmax(pmin(d[,i,3],d[,j,3])-pmax(d[,i,2],d[,j,2]),0),
            fracxy=pmax(pmin(d[,i,3],d[,j,3])-pmax(d[,i,2],d[,j,2]),0)/(d[,i,3]-d[,i,2]),
            fracyx=pmax(pmin(d[,i,3],d[,j,3])-pmax(d[,i,2],d[,j,2]),0)/(d[,j,3]-d[,j,2]),
            jntfrac=2*pmax(pmin(d[,i,3],d[,j,3])-pmax(d[,i,2],d[,j,2]),0)/(d[,i,3]-d[,i,2]+d[,j,3]-d[,j,2])
         )
   #Patch up those loose ends.
   if(m==1)
      o$graph<-o$graph[1,,]
   if(!diag)
      o$graph<-diag.remove(o$graph,remove.val=0)
   #Return the data structure
   o
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


#eval.edgeperturbation - Evaluate a function on a given graph with and without a given edge, returning
#the difference between the results in each case.

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


#pstar - Perform an approximate p* analysis using the logistic regression approximation.  Note that the
#result of this is returned as a GLM object, and subsequent printing/summarizing/etc. should be treated
#accordingly.

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


#read.nos - Read an input file in Neo-OrgStat format.  At this time, only the graph stack is read; any 
#coloring information is ignored.

read.nos<-function(file){
   #Get the formatting information
   f<-sapply(readLines(file,n=2),strsplit," ")
   #Parse the formatting information
   m<-as.numeric((f[[1]])[1])
   n<-as.numeric((f[[2]])[1])
   o<-as.numeric((f[[2]])[2])
   #Read the input data
   dat<-scan(file,skip=3)
   #Unpack the data in the proper order
   gstack<-array(dim=c(m,n,o))
   for(i in 1:m)
      for(j in 1:n)
         for(k in 1:o)
            gstack[i,j,k]<-dat[(i-1)*n*o+(j-1)*o+k]   
   #Return the stack
   gstack
}


#dyad.census - Return the Holland and Leinhardt MAN dyad census for a given graph 
#or graph stack
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


#triad.classify - Return the Davis and Leinhardt classification of a given triad
triad.classify<-function(dat,g=1,tri=c(1,2,3)){
   #Zeroth step: extract the triad
   if(length(dim(dat))==2)
      d<-dat[tri,tri]
   else
      d<-dat[g,tri,tri]
   #First, classify as NA if any entries are missing
   if(any(is.na(d[upper.tri(d)|lower.tri(d)])))
      man<-NA
   else{
      man<-dyad.census(d)  #Start with the dyad census
      #Refine the classification using configural properties
      if(all(man==c(0,2,1))){   #The two asym/one null triad
         ind<-apply(d,2,sum)  
         outd<-apply(d,1,sum)
         if(any(ind==2))
            man<-c(man,"U")   #"Up" variant
         else if(any(outd==2))
            man<-c(man,"D")   #"Down" variant
         else
            man<-c(man,"C")   #"Cyclic" variant
      }else if(all(man==c(1,1,1))){   #The one mut/one asym/one null triad
         ind<-apply(d,2,sum)  
         if(any(ind==2))
            man<-c(man,"D")   #"Down" variant
         else
            man<-c(man,"U")   #"Up" variant
      }else if(all(man==c(0,3,0))){   #The three asym triad
         ind<-apply(d,2,sum)  
         if(any(ind==2))
            man<-c(man,"T")   #"Transitive" variant
         else
            man<-c(man,"C")   #"Cyclic" variant
      }else if(all(man==c(1,2,0))){   #The one mut/two asym triad
         ind<-apply(d,2,sum)  
         outd<-apply(d,1,sum)
         if(any(ind==0))
            man<-c(man,"D")   #"Down" variant
         else if(any(outd==0))
            man<-c(man,"U")   #"Up" variant
         else
            man<-c(man,"C")   #"Cyclic" variant
      }         
   }
   #Return the classification
   paste(man,collapse="")
}


#triad.census - Conduct a Davis and Leinhardt triad census for a graph or graph stack
triad.census<-function(dat,g=1:stackcount(dat)){
   #First, define the triad class vector
   tc<-c("003","012","102","021D","021U","021C","111D","111U","030T","030C","201","120D","120U","120C","210","300")
   #Define an internal census function
   intcalc<-function(m,tcv){
      n<-dim(m)[1]
      tricent<-rep(0,length(tcv))
      for(i in 1:n)
         for(j in i:n)
            for(k in j:n)
               if((i!=j)&&(j!=k)&&(i!=k)){
                  tric<-triad.classify(m,tri=c(i,j,k))
                  tricent[tcv==tric]<-tricent[tcv==tric]+1
               }
      tricent
   }
   #Organize the data
   if(length(dim(dat))>2)
      d<-dat[g,,,drop=FALSE]
   else{
      d<-array(dim=c(1,dim(dat)[1],dim(dat)[2]))
      d[1,,]<-dat
   }
   d<-diag.remove(d,remove.val=0)  #Remove any diagonals
   #Perform the census
   census<-t(apply(d,1,intcalc,tc))
   colnames(census)<-tc
   #Return the result
   census
}

#gapply - Apply a function to vertex neighborhoods within a graph
gapply<-function(X,MARGIN,STATS,FUN,...,mode="digraph",diag=FALSE,distance=1,thresh=0,simplify=TRUE){
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

print.lnam<-function(x, digits = max(3, getOption("digits") - 3), ...){
   cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
   cat("Coefficients:\n")
   print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
   cat("\n")
}

summary.lnam<-function(object, ...){
   ans<-object
   class(ans)<-c("summary.lnam","lnam")
   ans
}

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

plot.lnam<-function(x,...){
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
