\name{equiv.clust}
\alias{equiv.clust}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Find Clusters of Positions Based on an Equivalence Relation }
\description{
  \code{equiv.clust} uses a definition of approximate equivalence (\code{equiv.fun}) to form a hierarchical clustering of network positions.  Where \code{dat} consists of multiple relations, all specified relations are considered jointly in forming the equivalence clustering.
}
\usage{
equiv.clust(dat, g=c(1:dim(dat)[1]), equiv.fun="sedist", 
   method="hamming", mode="digraph", diag=FALSE, 
   cluster.method="complete", glabels=dimnames(dat)[[1]][g], 
   plabels=dimnames(dat)[[2]], ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ A graph adjacency matrix, or stack thereof }
  \item{g}{ The elements of \code{dat} to use in clustering the vertices; by default, all structures are used }
  \item{equiv.fun}{ The distance function to use in clustering vertices (defaults to \code{\link{sedist}}) }
  \item{method}{ \code{method} parameter to be passed to \code{equiv.fun} }
  \item{mode}{ ``graph'' or ``digraph,'' as appropriate }
  \item{diag}{ A boolean indicating whether or not matrix diagonals (loops) should be interpreted as useful data }
  \item{cluster.method}{ The hierarchical clustering method to use (see \code{\link{hclust}}) }
  \item{glabels}{ Labels for the various graphs in \code{dat} }
  \item{plabels}{ Labels for the vertices of \code{dat} }
  \item{\dots}{ Additional arguments to \code{equiv.dist} }
}
\details{
This routine is essentially a joint front-end to \code{\link{hclust}} and various positional distance functions, though it defaults to structural equivalence in particular.  Taking the specified graphs as input, \code{equiv.clust} computes the distances between all pairs of positions using \code{equiv.dist}, and then performs a cluster analysis of the result.  The return value is an object of class \code{equiv.clust}, for which various secondary analysis methods exist.
}
\value{
An object of class \code{equiv.clust}
}
\references{ 
Breiger, R.L.; Boorman, S.A.; and Arabie, P.  (1975).  ``An Algorithm for Clustering Relational Data with Applications to Social Network Analysis and Comparison with Multidimensional Scaling.''  \emph{Journal of Mathematical Psychology}, 12, 328-383.

Burt, R.S.  (1976).  ``Positions in Networks.''  \emph{Social Forces}, 55, 93-122.

Wasserman, S., and Faust, K.  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ See \code{\link{sedist}} for an example of a distance function compatible with \code{equiv.clust}. }

%\section{Requires}{\code{mva}}

\seealso{ \code{\link{sedist}}, \code{\link{blockmodel}}}

\examples{

#Create a random graph with _some_ edge structure
g.p<-sapply(runif(20,0,1),rep,20)  #Create a matrix of edge 
                                   #probabilities
g<-rgraph(20,tprob=g.p)            #Draw from a Bernoulli graph 
                                   #distribution

#Cluster based on structural equivalence
eq<-equiv.clust(g)
plot(eq)
}
\keyword{ cluster }%-- one or more ...
\keyword{ math }
