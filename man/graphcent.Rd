\name{graphcent}
\alias{graphcent}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the (Harary) Graph Centrality Scores of Network Positions }
\description{
  \code{graphcent} takes a graph stack (\code{dat}) and returns the Harary graph centralities of positions within one graph (indicated by \code{nodes} and \code{g}, respectively).  Depending on the specified mode, graph centrality on directed or undirected geodesics will be returned; this function is compatible with \code{\link{centralization}}, and will return the theoretical maximum absolute deviation (from maximum) conditional on size (which is used by \code{\link{centralization}} to normalize the observed centralization score).
}
\usage{
graphcent(dat, g=1, nodes=c(1:dim(dat)[2]), gmode="digraph", 
    diag=FALSE, tmaxdev=FALSE, cmode="directed", 
    geodist.precomp=NULL, rescale=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ Data array to be analyzed.  By assumption, the first dimension of the array indexes the graph, with the next two indexing the actors. Alternately, this can be an n x n matrix (if only one graph is involved). }
  \item{g}{ Integer indicating the index of the graph for which centralities are to be calculated.  By default, \code{g==1}. }
  \item{nodes}{ List indicating which nodes are to be included in the calculation.  By default, all nodes are included. }
  \item{gmode}{ String indicating the type of graph being evaluated.  "digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  \code{gmode} is set to "digraph" by default. }
  \item{diag}{ Boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{tmaxdev}{ Boolean indicating whether or not the theoretical maximum absolute deviation from the maximum nodal centrality should be returned.  By default, \code{tmaxdev==FALSE}. }
  \item{cmode}{ String indicating the type of graph centrality being computed (directed or undirected geodesics). }
  \item{geodist.precomp}{ A \code{\link{geodist}} object precomputed for the graph to be analyzed (optional) }
  \item{rescale}{ If true, centrality scores are rescaled such that they sum to 1. }
}
\details{
  The Harary graph centrality of a vertex v is equal to \eqn{\frac{1}{\max_u d(v,u)}}{1/(max_u d(v,u))}, where \eqn{d(v,u)}{d(v,u)} is the geodesic distance from v to u.  Vertices with low graph centrality scores are likely to be near the ``edge'' of a graph, while those with high scores are likely to be near the ``middle.''  Compare this with \code{\link{closeness}}, which is based on the reciprocal of the sum of distances to all other vertices (rather than simply the maximum).
}
\value{
  A vector containing the centrality scores
}
\references{  }
\author{ Carter T. Butts \email{ctb@andrew.cmu.edu} }
\note{ Judicious use of \code{geodist.precomp} can save a great deal of time when computing multiple path-based indices on the same network. }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{centralization}} }

\examples{
g<-rgraph(10)     #Draw a random graph with 10 members
graphcent(g)    #Compute centrality scores
}
\keyword{ univar }%-- one or more ...
\keyword{ math }
