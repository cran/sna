\name{reachability}
\alias{reachability}
\alias{reachability_R}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Find the Reachability Matrix of a Graph }
\description{
  \code{reachability} takes one or more (possibly directed) graphs as input, producing the associated reachability matrices.
}
\usage{
reachability(dat, geodist.precomp=NULL, return.as.edgelist=FALSE, na.omit=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more graphs (directed or otherwise). }
  \item{geodist.precomp}{ optionally, a precomputed \code{\link{geodist}} object. }
  \item{return.as.edgelist}{ logical; return the result as an sna edgelist? }
  \item{na.omit}{ logical; omit missing edges when computing reach? }
}
\details{
For a digraph \eqn{G=(V,E)} with vertices \eqn{i} and \eqn{j}, let \eqn{P_{ij}}{P_ij} represent a directed \eqn{ij} path.  Then the (di)graph

\deqn{
R = \left(V\left(G\right),\left\{\left(i,j\right):i,j \in V\left(G\right), P_{ij} \in G\right\}\right)}{%
R = ( V(G), \{ (i,j): i,j in V(G), P_ij in G \} )}

is said to be the \emph{reachability graph} of \eqn{G}, and the adjacency matrix of \eqn{R} is said to be \eqn{G}'s \emph{reachability matrix}.  (Note that when \eqn{G} is undirected, we simply take each undirected edge to be bidirectional.)  Vertices which are adjacent in the reachability graph are connected by one or more directed paths in the original graph; thus, structural equivalence classes in the reachability graph are synonymous with strongly connected components in the original structure.

Bear in mind that -- as with all matters involving connectedness -- reachability is strongly related to size and density.  Since, for any given density, almost all structures of sufficiently large size are connected, reachability graphs associated with large structures will generally be complete.  Measures based on the reachability graph, then, will tend to become degenerate in the large \eqn{|V(G)|} limit (assuming constant positive density).

By default, \code{reachability} will try to build the reachability graph using an internal sparse graph approximation; this is no help on fully connected graphs (but not a lot worse than using an adjacency matrix), but will result in considerable savings for large graphs that are heavily fragmented.  (The intended design tradeoff is thus that one pays a small cost on the usually cheap cases, in exchange for much greater efficiency on the cases that would otherwise be prohibitively expensive.)  If \code{geodist.precomp} is given, however, the \eqn{O(N^2)} cost of an adjacency matrix representation has already been paid, and we simply employ what we are given -- so, if you want to force the internal use of adjacency matrices, just pass a \code{\link{geodist}} object.  Because the internal representation used is otherwise list based, using \code{return.as.edgelist=TRUE} will save resources; if you are using \code{reachability} as part of a more complex series of calls, it is thus recommended that you both pass and return sna edgelists unless you have a good reason not to do so.

When set, \code{na.omit} results in missing edges (i.e., edges with \code{NA} values) being removed prior to computation.  Since paths are not recomputed when \code{geodist.precomp} is passed, this option is only active when \code{geodist.precomp==NULL}; if this behavior is desired and precomputed distances are being used, such edges should be removed prior to the \code{\link{geodist}} call.
}
\value{
A reachability matrix, or the equivalent edgelist representation
}
\references{ Wasserman, S., and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }

\seealso{ \code{\link{geodist}} }

\examples{
#Find the reachability matrix for a sparse random graph
g<-rgraph(10,tprob=0.15)
rg<-reachability(g)
g  #Compare the two structures
rg

#Compare to the output of geodist
all(rg==(geodist(g)$counts>0))
}
\keyword{ algebra }% __ONLY ONE__ keyword per line
\keyword{ graphs }
