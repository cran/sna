\name{closeness}
\alias{closeness}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the Closeness Centrality Scores of Network Positions }
\description{
 \code{closeness} takes one or more graphs (\code{dat}) and returns the closeness centralities of positions (selected by \code{nodes}) within the graphs indicated by \code{g}.  Depending on the specified mode, closeness on directed or undirected geodesics will be returned; this function is compatible with \code{\link{centralization}}, and will return the theoretical maximum absolute deviation (from maximum) conditional on size (which is used by \code{\link{centralization}} to normalize the observed centralization score).
}
\usage{
closeness(dat, g=1, nodes=NULL, gmode="digraph", diag=FALSE, 
    tmaxdev=FALSE, cmode="directed", geodist.precomp=NULL, 
    rescale=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{g}{ integer indicating the index of the graph for which centralities are to be calculated (or a vector thereof).  By default, \code{g}=1. }
  \item{nodes}{ list indicating which nodes are to be included in the calculation.  By default, all nodes are included. }
  \item{gmode}{ string indicating the type of graph being evaluated.  "digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  \code{gmode} is set to "digraph" by default. }
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{tmaxdev}{ boolean indicating whether or not the theoretical maximum absolute deviation from the maximum nodal centrality should be returned.  By default, \code{tmaxdev==FALSE}. }
  \item{cmode}{ string indicating the type of closeness centrality being computed (distances on directed or undirected geodesics). }
  \item{geodist.precomp}{ a \code{\link{geodist}} object precomputed for the graph to be analyzed (optional) }
  \item{rescale}{ if true, centrality scores are rescaled such that they sum to 1. }
}
\details{
  The closeness of a vertex v is defined as

\deqn{C_C(v) = \frac{\left|V\left(G\right)\right|-1}{\sum_{i : i \neq v} d(v,i)}}{%
C_C(v) = (|V(G)|-1)/sum( d(v,i), i in V(G), i!=v )}

where \eqn{d(i,j)}{d(i,j)} is the geodesic distance between i and j (where defined).  Closeness is ill-defined on disconnected graphs; in such cases, this routine substitutes \code{Inf}.  It should be understood that this modification is not canonical (though it is common), but can be avoided by not attempting to measure closeness on disconnected graphs in the first place!  Intuitively, closeness provides an index of the extent to which a given vertex has short paths to all other vertices in the graph; this is one reasonable measure of the extent to which a vertex is in the ``middle'' of a given structure.
}
\value{
A vector, matrix, or list containing the closeness scores (depending on the number and size of the input graphs).
}
\references{Freeman, L.C.  (1979).  ``Centrality in Social Networks I: Conceptual Clarification.'' \emph{Social Networks}, 1, 215-239. }
\author{ Carter T. Butts, \email{buttsc@uci.edu} }
\note{ Judicious use of \code{geodist.precomp} can save a great deal of time when computing multiple path-based indices on the same network. }



\seealso{ \code{\link{centralization}} }

\examples{
g<-rgraph(10)     #Draw a random graph with 10 members
closeness(g)      #Compute closeness scores

}
\keyword{univar }%-- one or more ...
\keyword{ math }
\keyword{graphs}
