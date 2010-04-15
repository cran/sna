\name{hdist}
\alias{hdist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Find the Hamming Distances Between Two or More Graphs }
\description{
  \code{hdist} returns the Hamming distance between the labeled graphs \code{g1} and \code{g2} in set \code{dat} for dichotomous data, or else the absolute (manhattan) distance.  If \code{normalize} is true, this distance is divided by its dichotomous theoretical maximum (conditional on |V(G)|).
}
\usage{
hdist(dat, dat2=NULL, g1=NULL, g2=NULL, normalize=FALSE, 
    diag=FALSE, mode="digraph")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ a stack of input graphs. }
  \item{dat2}{ a second graph stack (optional). }
  \item{g1}{ a vector indicating which graphs to compare (by default, all elements of \code{dat}). }
  \item{g2}{ a vector indicating against which the graphs of \code{g1} should be compared (by default, all graphs). }
  \item{normalize}{ divide by the number of available dyads? }
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{mode}{ string indicating the type of graph being evaluated.  "digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  \code{mode} is set to "digraph" by default. }
}
\details{
  The Hamming distance between two labeled graphs \eqn{G_1} and \eqn{G_2} is equal to \eqn{|\{e : (e \in E(G_1), e \not\in E(G_2)) \wedge (e \not\in E(G_1), e \in E(G_2))\}|}{|\{e : (e in E(G_1) and e not in E(G_2)) or (e not in E(G_1) and e in E(G_2))\}|}.  In more prosaic terms, this may be thought of as the number of addition/deletion operations required to turn the edge set of \eqn{G_1} into that of \eqn{G_2}.  The Hamming distance is a highly general measure of structural similarity, and forms a metric on the space of graphs (simple or directed).  Users should be reminded, however, that the Hamming distance is extremely sensitive to nodal labeling, and should not be employed directly when nodes are interchangeable.  The structural distance (Butts and Carley (2001)), implemented in \code{\link{structdist}}, provides a natural generalization of the Hamming distance to the more general case of unlabeled graphs.

Null hypothesis testing for Hamming distances is available via \code{\link{cugtest}}, and \code{\link{qaptest}}; graphs which minimize the Hamming distances to all members of a graph set can be found by \code{\link{centralgraph}}.  For an alternative means of comparing the similarity of graphs, consider \code{\link{gcor}}.
}
\value{
A matrix of Hamming distances
}
\references{ 
Banks, D., and Carley, K.M.  (1994).  \dQuote{Metric Inference for Social Networks.}  \emph{Journal of Classification}, 11(1), 121-49.

Butts, C.T. and Carley, K.M.  (2005).  \dQuote{Some Simple Algorithms for Structural Comparison.}  \emph{Computational and Mathematical Organization Theory,} 11(4), 291-305.

Butts, C.T., and Carley, K.M.  (2001).  \dQuote{Multivariate Methods for Interstructural Analysis.}  CASOS Working Paper, Carnegie Mellon University. 

Hamming, R.W. (1950). \dQuote{Error Detecting and Error Correcting Codes.} \emph{Bell System Technical Journal,} 29, 147-160.
}

\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ For non-dichotomous data, the distance which is returned is simply the sum of the absolute edge-wise differences. }


 %~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{sdmat}}, \code{\link{structdist}} }

\examples{
#Get some random graphs
g<-rgraph(5,5,tprob=runif(5,0,1))

#Find the Hamming distances
hdist(g)
}
\keyword{ univar }
\keyword{ multivariate }%-- one or more ...
\keyword{ graphs }
