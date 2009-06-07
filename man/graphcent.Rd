\name{graphcent}
\alias{graphcent}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the (Harary) Graph Centrality Scores of Network Positions }
\description{
  \code{graphcent} takes one or more graphs (\code{dat}) and returns the Harary graph centralities of positions (selected by \code{nodes}) within the graphs indicated by \code{g}.  Depending on the specified mode, graph centrality on directed or undirected geodesics will be returned; this function is compatible with \code{\link{centralization}}, and will return the theoretical maximum absolute deviation (from maximum) conditional on size (which is used by \code{\link{centralization}} to normalize the observed centralization score).
}
\usage{
graphcent(dat, g=1, nodes=NULL, gmode="digraph", diag=FALSE, 
    tmaxdev=FALSE, cmode="directed", geodist.precomp=NULL, 
    rescale=FALSE, ignore.eval)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{g}{ integer indicating the index of the graph for which centralities are to be calculated (or a vector thereof).  By default, \code{g==1}. }
  \item{nodes}{ list indicating which nodes are to be included in the calculation.  By default, all nodes are included. }
  \item{gmode}{ string indicating the type of graph being evaluated.  "digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  \code{gmode} is set to "digraph" by default. }
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{tmaxdev}{ boolean indicating whether or not the theoretical maximum absolute deviation from the maximum nodal centrality should be returned.  By default, \code{tmaxdev==FALSE}. }
  \item{cmode}{ string indicating the type of graph centrality being computed (directed or undirected geodesics). }
  \item{geodist.precomp}{ a \code{\link{geodist}} object precomputed for the graph to be analyzed (optional) }
  \item{rescale}{ if true, centrality scores are rescaled such that they sum to 1. }
  \item{ignore.eval}{ logical; should edge values be ignored when calculating geodesics?}
}
\details{
  The Harary graph centrality of a vertex v is equal to \eqn{\frac{1}{\max_u d(v,u)}}{1/(max_u d(v,u))}, where \eqn{d(v,u)}{d(v,u)} is the geodesic distance from v to u.  Vertices with low graph centrality scores are likely to be near the ``edge'' of a graph, while those with high scores are likely to be near the ``middle.''  Compare this with \code{\link{closeness}}, which is based on the reciprocal of the sum of distances to all other vertices (rather than simply the maximum).
}
\value{
  A vector, matrix, or list containing the centrality scores (depending on the number and size of the input graphs).
}
\references{ Hage, P. and Harary, F.  (1995).  ``Eccentricity and Centrality in Networks.''  \emph{Social Networks}, 17:57-63. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ Judicious use of \code{geodist.precomp} can save a great deal of time when computing multiple path-based indices on the same network. }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{centralization}} }

\examples{
g<-rgraph(10)     #Draw a random graph with 10 members
graphcent(g)    #Compute centrality scores
}
\keyword{ univar }%-- one or more ...
\keyword{ math }
