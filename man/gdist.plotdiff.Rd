\name{gdist.plotdiff}
\alias{gdist.plotdiff}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Plot Differences in Graph-level Statistics Against Inter-graph Distances}
\description{
  For a given graph set, \code{gdist.plotdiff} plots the distances between graphs against their distances (or differences) on a set of graph-level measures.
}
\usage{
gdist.plotdiff(d, meas, method="manhattan", jitter=TRUE, 
    xlab="Inter-Graph Distance", ylab="Measure Distance", 
    lm.line=FALSE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{d}{ A matrix containing the inter-graph distances }
  \item{meas}{ An n x m matrix containing the graph-level indices; rows of this matrix must correspond to graphs, and columns to indices }
  \item{method}{ The distance method used by \code{\link{dist}} to establish differences/distances between graph GLI values.  By default, absolute ("manhattan") differences are used. }
  \item{jitter}{ Should values be jittered prior to display? }
  \item{xlab}{ A label for the X axis }
  \item{ylab}{ A label for the Y axis }
  \item{lm.line}{ Include a least-squares line? }
  \item{\dots}{ Additional arguments to \code{\link{plot}} }
}
\details{
  \code{gdist.plotdiff} works by taking the distances between all graphs on \code{meas} and then plotting these distances against \code{d} for all pairs of graphs (with, optionally, an added least-squares line for reference value).  This can be a useful exploratory tool for relating inter-graph distances (e.g., Hamming distances) to differences on other attributes.  
}
\value{
None
}
\references{ Butts, C.T., and Carley, K.M.  (2001).  ``Multivariate Methods for Interstructural Analysis.''  CASOS working paper, Carnegie Mellon University. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ This function is actually quite generic, and can be used with node-level -- or even non-network -- data as well. }

 %~Make other sections like WARNING with \section{WARNING }{....} ~

%\section{Requires}{\code{mva}}

\seealso{ \code{\link{gdist.plotstats}}, \code{\link{gclust.boxstats}}, \code{\link{gclust.centralgraph}} }

\examples{
#Generate some random graphs with varying densities
g<-rgraph(10,20,tprob=runif(20,0,1))

#Find the Hamming distances between graphs
g.h<-hdist(g)

#Plot the relationship between distance and differences in density
gdist.plotdiff(g.h,gden(g),lm.line=TRUE)
}
\keyword{ hplot }%-- one or more ...
