\name{reachability}
\alias{reachability}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Find the Reachability Matrix of a Graph }
\description{
  \code{reachability} takes a (possibly directed) graph as input, producing the associated reachability matrix.
}
\usage{
reachability(dat, geodist.precomp=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ A graph adjacency matrix (directed or otherwise) }
  \item{geodist.precomp}{ Optionally, a precomputed \code{\link{geodist}} object }
}
\details{
For a digraph \eqn{G=(V,E)} with vertices \eqn{i} and \eqn{j}, let \eqn{P_{ij}}{P_ij} represent a directed \eqn{ij} path.  Then the graph

\deqn{
R = \left(V\left(G\right),\left\{\left(i,j\right):i,j \in V\left(G\right), P_{ij} \in G\right\}\right)}{%
R = ( V(G), \{ (i,j): i,j in V(G), P_ij in G \} )}

is said to be the \emph{reachability graph} of \eqn{G}, and the adjacency matrix of \eqn{R} is said to be \eqn{G}'s \emph{reachability matrix}.  (Note that when \eqn{G} is undirected, we simply take each undirected edge to be bidirectional.)  Vertices which are adjacent in the reachability graph are connected by one or more directed paths in the original graph; thus, structural equivalence classes in the reachability graph are synonymous with strongly connected components in the original structure.

Bear in mind that -- as with all matters involving connectedness -- reachability is strongly related to size and density.  Since, for any given density, almost all structures of sufficiently large size are connected, reachability graphs associated with large structures will generally be complete.  Measures based on the reachability graph, then, will tend to become degenerate in the large \eqn{|V(G)|} limit (assuming constant positive density).
}
\value{
A reachability matrix
}
\references{ Wasserman, S., and Faust, K.  (1994).  ``Social Network Analysis: Methods and Applications.''  Cambridge: Cambridge University Press. }
\author{ Carter T. Butts \email{ctb@andrew.cmu.edu} }

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
