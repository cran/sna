\name{closeness}
\alias{closeness}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the Closeness Centrality Scores of Network Positions }
\description{
 \code{closeness} takes a graph stack (\code{dat}) and returns the closeness centralities of positions within one graph (indicated by \code{nodes} and \code{g}, respectively).  Depending on the specified mode, closeness on directed or undirected geodesics will be returned; this function is compatible with \code{\link{centralization}}, and will return the theoretical maximum absolute deviation (from maximum) conditional on size (which is used by \code{\link{centralization}} to normalize the observed centralization score).
}
\usage{
closeness(dat, g=1, nodes=c(1:dim(dat)[2]), gmode="digraph", 
    diag=FALSE, tmaxdev=FALSE, cmode="directed", 
    geodist.precomp=NULL, rescale=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ Data array to be analyzed.  By assumption, the first dimension of the array indexes the graph, with the next two indexing the actors. Alternately, this can be an n x n matrix (if only one graph is involved). }
  \item{g}{ Integer indicating the index of the graph for which centralities are to be calculated.  By default, \code{g}=1. }
  \item{nodes}{ List indicating which nodes are to be included in the calculation.  By default, all nodes are included. }
  \item{gmode}{ String indicating the type of graph being evaluated.  "digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  \code{gmode} is set to "digraph" by default. }
  \item{diag}{ Boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{tmaxdev}{ Boolean indicating whether or not the theoretical maximum absolute deviation from the maximum nodal centrality should be returned.  By default, \code{tmaxdev==FALSE}. }
  \item{cmode}{ String indicating the type of closeness centrality being computed (distances on directed or undirected geodesics). }
  \item{geodist.precomp}{ A \code{\link{geodist}} object precomputed for the graph to be analyzed (optional) }
  \item{rescale}{ If true, centrality scores are rescaled such that they sum to 1. }
}
\details{
  The closeness of a vertex v is defined as

\deqn{C_C(v) = \frac{\left|V\left(G\right)\right|-1}{\sum_{i : i \neq v} d(v,i)}}{%
C_C(v) = (|V(G)|-1)/sum( d(v,i), i in V(G), i!=v )}

where \eqn{d(i,j)}{d(i,j)} is the geodesic distance between i and j (where defined).  Closeness is ill-defined on disconnected graphs; in such cases, this routine substitutes \code{Inf}.  It should be understood that this modification is not canonical (though it is common), but can be avoided by not attempting to measure closeness on disconnected graphs in the first place!  Intuitively, closeness provides an index of the extent to which a given vertex has short paths to all other vertices in the graph; this is one reasonable measure of the extent to which a vertex is in the ``middle'' of a given structure.
}
\value{
A vector containing the closeness scores.
}
\references{Freeman, L.C.  (1979).  ``Centrality in Social Networks I: Conceptual Clarification.'' \emph{Social Networks}, 1, 215-239. }
\author{ Carter T. Butts, \email{buttsc@uci.edu} }
\note{ Judicious use of \code{geodist.precomp} can save a great deal of time when computing multiple path-based indices on the same network.  Before \code{sna} version 0.50, \code{closeness} substituted one plus the maximum possible path length for disconnected vertex pairs.}



\seealso{ \code{\link{centralization}} }

\examples{
g<-rgraph(10)     #Draw a random graph with 10 members
closeness(g)      #Compute closeness scores

}
\keyword{univar }%-- one or more ...
\keyword{ math }
\keyword{graphs}
