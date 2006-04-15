\name{consensus}
\alias{consensus}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Estimate a Consensus Structure from Multiple Observations }
\description{
  \code{consensus} estimates a central or consensus structure given multiple observations, using one of several algorithms.
}
\usage{
consensus(dat, mode="digraph", diag=FALSE, method="central.graph", 
    tol=1e-06, maxiter=1e3, verbose=TRUE, no.bias=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ a set of input graphs (must have same order). }
  \item{mode}{ \code{"digraph"} for directed data, else \code{"graph"}. }
  \item{diag}{ logical; should diagonals (loops) be treated as data? }
  \item{method}{ one of \code{"central.graph"}, \code{"single.reweight"}, \code{"iterative.reweight"}, \code{"romney.batchelder"}, \code{"PCA.reweight"}, \code{"LAS.intersection"}, \code{"LAS.union"}, \code{"OR.row"}, or \code{"OR.col"}.  }
  \item{tol}{ convergence tolerance for the iterative reweighting and B-R algorithms.}
  \item{maxiter}{ maximum number of iterations to take (regardless of convergence) for the iterative reweighting and B-R algorithms.}
  \item{verbose}{ logical; should bias and competency parameters be reported (where computed)?}
  \item{no.bias}{ logical; should responses be assumed to be unbiased? } 
}
\details{

The term ``consensus structure'' is used by a number of authors to reflect a notion of shared or common perceptions of social structure among a set of observers.  As there are many interpretations of what is meant by ``consensus'' (and as to how best to estimate it), several algorithms are employed here:

\enumerate{
\item \code{central.graph}: Estimate the consensus structure using the central graph.  This corresponds to a ``median response'' notion of consensus.
\item \code{single.reweight}: Estimate the consensus structure using subject responses, reweighted by mean graph correlation.  This corresponds to an ``expertise-weighted vote'' notion of consensus.
\item \code{iterative.reweight}: Similar to \code{single.reweight}, but the consensus structure and accuracy parameters are estimated via an iterated proportional fitting scheme.  The implementation employed here uses both bias and competency parameters.
\item \code{romney.batchelder}: Fits a Romney-Batchelder informant accuracy model using IPF.  This is very similar to \code{iterative.reweight}, but can be interpreted as the result of a process in which each informant report is correct with a probability equal to the informant's competency score, and otherwise equal to a Bernoulli trial with parameter equal to the informant's bias score.
\item \code{PCA.reweight}: Estimate the consensus using the (scores on the) first component of a network PCA.  This corresponds to a ``shared theme'' or ``common element'' notion of consensus.
\item \code{LAS.intersection}: Estimate the consensus structure using the locally aggregated structure (intersection rule).  In this model, an i->j edge exists iff i \emph{and} j agree that it exists.
\item \code{LAS.union}: Estimate the consensus structure using the locally aggregated structure (union rule).  In this model, an i->j edge exists iff i \emph{or} j agree that it exists.
\item \code{OR.row}: Estimate the consensus structure using own report.  Here, we take each informant's outgoing tie reports to be correct.
\item \code{OR.col}: Estimate the consensus structure using own report.  Here, we take each informant's incoming tie reports to be correct.
}

Note that the results returned by the single weighting algorithms are not dichotomized by default; since some algorithms thus return valued graphs, dichotomization may be desirable prior to use.

It should be noted that a model for estimating an underlying criterion structure from multiple informant reports is provided in \code{\link{bbnam}}; if your goal is to reconstruct an ``objective'' network from informant reports, this (or the R-B model) may prove more useful than the ad-hoc solutions.
}
\value{
An adjacency matrix representing the consensus structure
}
\references{ 
Banks, D.L., and Carley, K.M.  (1994).  ``Metric Inference for Social Networks.''  \emph{Journal of Classification,}  11(1), 121-49.

Butts, C.T., and Carley, K.M.  (2001).  ``Multivariate Methods for Inter-Structural Analysis.''  CASOS Working Paper, Carnegie Mellon University.

Krackhardt, D.  (1987).  ``Cognitive Social Structures.'' \emph{Social Networks,} 9, 109-134. 

Romney, A.K.; Weller, S.C.; and Batchelder, W.H.  (1986).  ``Culture as Consensus: A Theory of Culture and Informant Accuracy.''  \emph{American Anthropologist,} 88(2), 313-38.
}
\author{ Carter T. Butts \email{buttsc@uci.edu} }

%\note{ Eventually, this routine will also support the (excellent) consensus methods of Romney and Batchelder; since these are similar in many respects to the \code{\link{bbnam}} model, users may wish to try this alternative for now. 
%}

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
\keyword{ graphs }%-- one or more ...
