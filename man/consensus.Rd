\name{consensus}
\alias{consensus}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Estimate a Consensus Structure from Multiple Observations }
\description{
  \code{consensus} estimates a central or consensus structure given multiple observations, using one of several algorithms.
}
\usage{
consensus(dat, mode="digraph", diag=FALSE, method="central.graph", 
    tol=0.01)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ An m x n x n graph stack}
  \item{mode}{ ``digraph'' for directed data, else ``graph'' }
  \item{diag}{ A boolean indicating whether the diagonals (loops) should be treated as data }
  \item{method}{ One of ``central.graph'', ``single.reweight'', ``PCA.reweight''  }
  \item{tol}{ Tolerance for the iterative reweighting algorithm (not currently supported)}
}
\details{

The term ``consensus structure'' is used by a number of authors to reflect a notion of shared or common perceptions of social structure among a set of observers.  As there are many interpretations of what is meant by ``consensus'' (and as to how best to estimate it), several algorithms are employed here:

\enumerate{
\item \code{central.graph}: Estimate the consensus structure using the central graph.  This correponds to a ``median response'' notion of consensus.
\item \code{single.reweight}: Estimate the consensus structure using subject responses, reweighted by mean graph correlation.  This corresponds to an ``expertise-weighted vote'' notion of consensus.
\item \code{PCA.reweight}: Estimate the consensus using the (scores on the) first component of a network PCA.  This corresponds to a ``shared theme'' or ``common element'' notion of consensus.
}

Note that the reweighted algorithms are not dichotomized by default; since these return valued graphs, dichotomization may be desirable prior to use.

It should be noted that a model for estimating an underlying criterion structure from multiple informant reports is provided in \code{\link{bbnam}}; if your goal is to reconstruct an ``objective'' network from informant reports, this may prove more useful.
}
\value{
An adjacency matrix representing the consensus structure
}
\references{ 
Banks, D.L., and Carley, K.M.  (1994).  ``Metric Inference for Social Networks.''  \emph{Journal of Classification.'' 

Butts, C.T., and Carley, K.M.  (2001).  ``Multivariate Methods for Inter-Structural Analysis.''  CASOS Working Paper, Carnegie Mellon University.

Krackhardt, D.  (1987).  ``Cognitive Social Structures.'' \emph{Social Networks,} 9, 109-134. 
}
\author{ Carter T. Butts \email{ctb@andrew.cmu.edu} }

\note{ Eventually, this routine will also support the (excellent) consensus methods of Romney and Batchelder; since these are similar in many respects to the \code{\link{bbnam}} model, users may wish to try this alternative for now. 
}

\seealso{ \code{\link{bbnam}}, \code{\link{centralgraph}} }

\examples{

#Generate some test data
g<-rgraph(5)
g.pobs<-g*0.9+(1-g)*0.5
g.obs<-rgraph(5,5,tprob=g.pobs)

#Find some consensus structures
consensus(g.obs)                           #Central graph
consensus(g.obs,method="single.reweight")  #Single reweighting
consensus(g.obs,method="PCA.reweight")     #1st component in network PCA
}
\keyword{ network}%-- one or more ...
