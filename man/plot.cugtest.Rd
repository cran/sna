\name{plot.cugtest}
\alias{plot.cugtest}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Plotting for cugtest Objects }
\description{
  Plots the distribution of a CUG test statistic.
}
\usage{
\method{plot}{cugtest}(x, mode="density", ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ A \code{\link{cugtest}} object }
  \item{mode}{ ``density'' for kernel density estimation, ``hist'' for histogram }
  \item{\dots}{ Additional arguments to \code{\link{plot}} }
}
\details{
  In addition to the quantiles associated with a CUG test, it is often useful to examine the form of the distribution of the test statistic.  \code{plot.cugtest} facilitates this.
}
\value{
None
}
\references{ Anderson, B.S.; Butts, C.T.; and Carley, K.M. (1999). ``The Interaction of Size and Density with Graph-Level Indices.'' \emph{Social Networks}, 21(3), 239-267. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

 %~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{cugtest}} }

\examples{
#Draw two random graphs, with different tie probabilities
dat<-rgraph(20,2,tprob=c(0.2,0.8))

#Is their correlation higher than would be expected, conditioning 
#only on size?
cug<-cugtest(dat,gcor,cmode="order")
summary(cug)
plot(cug)

#Now, let's try conditioning on density as well.
cug<-cugtest(dat,gcor)
plot(cug)
}
\keyword{ hplot }%-- one or more ...
