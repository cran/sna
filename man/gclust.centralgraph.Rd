\name{gclust.centralgraph}
\alias{gclust.centralgraph}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Get Central Graphs Associated with Graph Clusters }
\description{
  Calculates central graphs associated with particular graph clusters (as indicated by the \code{k} partition of \code{h}).
}
\usage{
gclust.centralgraph(h, k, dat, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{h}{ an \code{\link{hclust}} object, based on a graph stack in \code{dat}. }
  \item{k}{ the number of groups to evaluate. }
  \item{dat}{ one or more graphs (on which the clustering was performed). }
  \item{\dots}{ additional arguments to \code{\link{centralgraph}}. }
}
\details{
\code{gclust.centralgraph} uses \code{\link{cutree}} to cut the hierarchical clustering in \code{h} into \code{k} groups.  \code{\link{centralgraph}} is then called on each cluster, and the results are returned as a graph stack.  This is a useful tool for interpreting clusters of (labeled) graphs, with the resulting central graphs being subsequently analyzed using standard SNA methods.
}
\value{
An array containing the stack of central graph adjacency matrices
}
\references{ Butts, C.T., and Carley, K.M.  (2001).  ``Multivariate Methods for Interstructural Analysis.''  CASOS working paper, Carnegie Mellon University. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

 %~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{hclust}}, \code{\link{centralgraph}}, \code{\link{gclust.boxstats}}, \code{\link{gdist.plotdiff}}, \code{\link{gdist.plotstats}} }

\examples{
#Create some random graphs
g<-rgraph(10,20,tprob=c(rbeta(10,15,2),rbeta(10,2,15)))

#Find the Hamming distances between them
g.h<-hdist(g)

#Cluster the graphs via their Hamming distances
g.c<-hclust(as.dist(g.h))

#Now find central graphs by cluster for a two cluster solution
g.cg<-gclust.centralgraph(g.c,2,g)

#Plot the central graphs
gplot(g.cg[1,,])
gplot(g.cg[2,,])
}
\keyword{ cluster }%-- one or more ...
\keyword{ math }
\keyword{ graphs }
