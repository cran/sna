\name{stresscent}
\alias{stresscent}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the Stress Centrality Scores of Network Positions }
\description{
  \code{stresscent} takes a graph stack (\code{dat}) and returns the stress centralities of positions within one graph (indicated by \code{nodes} and \code{g}, respectively).  Depending on the specified mode, stress on directed or undirected geodesics will be returned; this function is compatible with \code{\link{centralization}}, and will return the theoretical maximum absolute deviation (from maximum) conditional on size (which is used by \code{\link{centralization}} to normalize the observed centralization score).
}
\usage{
stresscent(dat, g=1, nodes=c(1:dim(dat)[2]), gmode="digraph", 
    diag=FALSE, tmaxdev=FALSE, cmode="directed", 
    geodist.precomp=NULL, rescale=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ Data array to be analyzed.  By assumption, the first dimension of the array indexes the graph, with the next two indexing the actors. Alternately, this can be an n x n matrix (if only one graph is involved). }
  \item{g}{ Integer indicating the index of the graph for which centralities are to be calculated.  By default, \code{g==1}. }
  \item{nodes}{ List indicating which nodes are to be included in the calculation.  By default, all nodes are included. }
  \item{gmode}{ String indicating the type of graph being evaluated.  \code{"digraph"} indicates that edges should be interpreted as directed; \code{"graph"} indicates that edges are undirected.  \code{gmode} is set to \code{"digraph"} by default.}
  \item{diag}{ Boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{tmaxdev}{ Boolean indicating whether or not the theoretical maximum absolute deviation from the maximum nodal centrality should be returned.  By default, \code{tmaxdev==FALSE}. }
  \item{cmode}{ String indicating the type of betweenness centrality being computed (directed or undirected geodesics). }
  \item{geodist.precomp}{A \code{\link{geodist}} object precomputed for the graph to be analyzed (optional)  }
  \item{rescale}{ If true, centrality scores are rescaled such that they sum to 1. }
}
\details{
 The stress of a vertex, v, is given by

\deqn{C_S(v) = \sum_{i,j : i \neq j,i \neq v,j \neq v} g_{ivj}}{%
C_S(v) = sum( g_ivj, i,j: i!=j,i!=v,j!=v)}

where \eqn{g_{ijk}}{g_ijk} is the number of geodesics from i to k through j.  Conceptually, high-stress vertices lie on a large number of shortest paths between other vertices; they can thus be thought of as ``bridges'' or ``boundary spanners.''  Compare this with \code{\link{betweenness}}, which weights shortest paths by the inverse of their redundancy.
}
\value{
A vector of centrality scores
}
%\references{ }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ Judicious use of \code{geodist.precomp} can save a great deal of time when computing multiple path-based indices on the same network. }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{centralization}}}

\examples{
g<-rgraph(10)     #Draw a random graph with 10 members
stresscent(g)     #Compute stress scores
}
\keyword{univar}
\keyword{ math }
