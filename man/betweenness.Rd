\name{betweenness}
\alias{betweenness}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the Betweenness Centrality Scores of Network Positions }
\description{
   \code{betweenness} takes a graph stack (\code{dat}) and returns the betweenness centralities of positions within one graph (indicated by \code{nodes} and \code{g}, respectively).  Depending on the specified mode, betweenness on directed or undirected geodesics will be returned; this function is compatible with \code{\link{centralization}}, and will return the theoretical maximum absolute deviation (from maximum) conditional on size (which is used by \code{\link{centralization}} to normalize the observed centralization score).
}
\usage{
betweenness(dat, g=1, nodes=c(1:dim(dat)[2]), gmode="digraph", 
    diag=FALSE, tmaxdev=FALSE, cmode="directed", 
    geodist.precomp=NULL, rescale=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ Data array to be analyzed.  By assumption, the first dimension of the array indexes the graph, with the next two indexing the actors. Alternately, this can be an n x n matrix (if only one graph is involved). }
  \item{g}{ Integer indicating the index of the graph for which centralities are to be calculated.  By default, \code{g}=1. }
  \item{nodes}{ List indicating which nodes are to be included in the calculation.  By default, all nodes are included. }
  \item{gmode}{ String indicating the type of graph being evaluated.  "digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  \code{dmode} is set to "digraph" by default. }
  \item{diag}{ Boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{tmaxdev}{ Boolean indicating whether or not the theoretical maximum absolute deviation from the maximum nodal centrality should be returned.  By default, \code{tmaxdev}==\code{FALSE}. }
  \item{cmode}{ String indicating the type of betweenness centrality being computed (directed or undirected geodesics). }
  \item{geodist.precomp}{ A \code{\link{geodist}} object precomputed for the graph to be analyzed (optional) }
  \item{rescale}{ If true, centrality scores are rescaled such that they sum to 1. }
}
\details{
The betweenness of a vertex, v, is given by

\deqn{
C_B(v) = \sum_{i,j : i \neq j, i \neq v, j \neq v} \frac{g_{ivj}}{g_{ij}}}{%
C_B(v) = sum( g_ivj / g_ij, i,j: i!=j,i!=v,j!=v )}

where \eqn{g_{ijk}}{g_ijk} is the number of geodesics from i to k through j.  Conceptually, high-betweenness vertices lie on a large number of non-redundant shortest paths between other vertices; they can thus be thought of as ``bridges'' or ``boundary spanners.'' 
}
\value{
A vector containing the betweenness scores.
}
\references{ Freeman, L.C.  (1979).  ``Centrality in Social Networks I: Conceptual Clarification.'' \emph{Social Networks}, 1, 215-239.}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ Judicious use of \code{geodist.precomp} can save a great deal of time when computing multiple path-based indices on the same network. }

\section{Warning }{Rescale may cause unexpected results if all actors have zero betweenness.}

\seealso{ \code{\link{centralization}} }

\examples{
g<-rgraph(10)     #Draw a random graph with 10 members
betweenness(g)    #Compute betweenness scores
}

\keyword{univar}
\keyword{ graphs }
