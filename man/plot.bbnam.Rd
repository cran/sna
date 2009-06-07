\name{plot.bbnam}
\alias{plot.bbnam}
\alias{plot.bbnam.fixed}
\alias{plot.bbnam.pooled}
\alias{plot.bbnam.actor}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Plotting for bbnam Objects}
\description{
  Generates various plots of posterior draws from the \code{\link{bbnam}} model.
}
\usage{
\method{plot}{bbnam}(x, mode="density", intlines=TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ A \code{bbnam} object }
  \item{mode}{ ``density'' for kernel density estimators of posterior marginals; otherwise, histograms are used }
  \item{intlines}{ Plot lines for the 0.9 central posterior probability intervals? }
  \item{\dots}{ Additional arguments to \code{\link{plot}} }
}
\details{
  \code{plot.bbnam} provides plots of the estimated posterior marginals for the criterion graph and error parameters (as appropriate).  Plotting may run into difficulties when dealing with large graphs, due to the problem of getting all of the various plots on the page; the routine handles these issues reasonably intelligently, but there is doubtless room for improvement.
}
\value{
None
}
\references{Butts, C.T. (1999). ``Informant (In)Accuracy and Network Estimation: A Bayesian Approach.'' CASOS Working Paper, Carnegie Mellon University. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

 %~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{bbnam}}}

\examples{
#Create some random data
g<-rgraph(5)
g.p<-0.8*g+0.2*(1-g)
dat<-rgraph(5,5,tprob=g.p)

#Define a network prior
pnet<-matrix(ncol=5,nrow=5)
pnet[,]<-0.5
#Define em and ep priors
pem<-matrix(nrow=5,ncol=2)
pem[,1]<-3
pem[,2]<-5
pep<-matrix(nrow=5,ncol=2)
pep[,1]<-3
pep[,2]<-5

#Draw from the posterior
b<-bbnam(dat,model="actor",nprior=pnet,emprior=pem,epprior=pep,
    burntime=100,draws=100)
#Print a summary of the posterior draws
summary(b)
#Plot the result
plot(b)
}
\keyword{ hplot}%-- one or more ...
