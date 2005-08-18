\name{gclust.boxstats}
\alias{gclust.boxstats}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Plot Statistics Associated with Graph Clusters }
\description{
  \code{gclust.boxstats} creates side-by-side boxplots of graph statistics based on a hierarchical clustering of networks (cut into \code{k} sets).
}
\usage{
gclust.boxstats(h, k, meas, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{h}{ an \code{\link{hclust}} object, presumably formed by clustering a set of structural distances. }
  \item{k}{ the number of groups to evaluate. }
  \item{meas}{ a vector of length equal to the number of graphs in \code{h}, containing a GLI to be evaluated. }
  \item{\dots}{ additional parameters to \code{\link{boxplot}}. }
}
\details{
  \code{gclust.boxstats} simply takes the \code{\link{hclust}} object in \code{h}, applies \code{\link{cutree}} to form \code{k} groups, and then uses \code{\link{boxplot}} on the distribution of \code{meas} by group.  This can be quite handy for assessing graph clusters.
}
\value{
None
}
\references{ Butts, C.T., and Carley, K.M.  (2001).  ``Multivariate Methods for Interstructural Analysis.''  CASOS working paper, Carnegie Mellon University. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ Actually, this function will work with any \code{\link{hclust}} object and measure matrix; the data need not originate with social networks.  For this reason, the clever may also employ this function in conjunction with \code{\link{sedist}} or \code{\link{equiv.clust}} to plot NLIs against clusters of positions within a graph.}

 %~Make other sections like WARNING with \section{WARNING }{....} ~

%\section{Requires}{\code{mva}}

\seealso{ \code{\link{gclust.centralgraph}}, \code{\link{gdist.plotdiff}}, \code{\link{gdist.plotstats}}  }

\examples{
#Create some random graphs
g<-rgraph(10,20,tprob=c(rbeta(10,15,2),rbeta(10,2,15)))

#Find the Hamming distances between them
g.h<-hdist(g)

#Cluster the graphs via their Hamming distances
g.c<-hclust(as.dist(g.h))

#Now display boxplots of density by cluster for a two cluster solution
gclust.boxstats(g.c,2,gden(g))
}
\keyword{ hplot }%-- one or more ...
