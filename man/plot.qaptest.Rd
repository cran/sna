\name{plot.qaptest}
\alias{plot.qaptest}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Plotting for qaptest Objects }
\description{
  Plots the Distribution of a QAP Test Statistic.
}
\usage{
\method{plot}{qaptest}(x, mode="density", ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ A \code{\link{qaptest}} object }
  \item{mode}{ ``density'' for kernel density estimation, ``hist'' for histogram }
  \item{\dots}{ Additional arguments to \code{\link{plot}} }
}
\details{
  In addition to the quantiles associated with a QAP test, it is often useful to examine the form of the distribution of the test statistic.  \code{plot.qaptest} facilitates this.
}
\value{
  None
}
\references{
Hubert, L.J., and Arabie, P.  (1989).  ``Combinatorial Data Analysis: Confirmatory Comparisons Between Sets of Matrices.''  \emph{Applied Stochastic Models and Data Analysis}, 5, 273-325.

Krackhardt, D.  (1987).  ``QAP Partialling as a Test of Spuriousness.'' \emph{Social Networks}, 9 171-186.

Krackhardt, D.  (1988).  ``Predicting With Networks: Nonparametric Multiple Regression Analyses of Dyadic Data.''  \emph{Social Networks}, 10, 359-382.}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{qaptest}} }

\examples{

#Generate three graphs
g<-array(dim=c(3,10,10))
g[1,,]<-rgraph(10)
g[2,,]<-rgraph(10,tprob=g[1,,]*0.8)
g[3,,]<-1; g[3,1,2]<-0              #This is nearly a clique

#Perform qap tests of graph correlation
q.12<-qaptest(g,gcor,g1=1,g2=2)
q.13<-qaptest(g,gcor,g1=1,g2=3)

#Examine the results
summary(q.12)
plot(q.12)
summary(q.13)
plot(q.13)
}
\keyword{ hplot }%-- one or more ...
