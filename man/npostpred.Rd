\name{npostpred}
\alias{npostpred}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Take Posterior Predictive Draws for Functions of Networks }
\description{
  \code{npostpred} takes a list or data frame, \code{b}, and applies the function \code{FUN} to each element of \code{b}'s \code{net} member.
}
\usage{
npostpred(b, FUN, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{b}{ A list or data frame containing posterior network draws; these draws must take the form of a graph stack, and must be the member of \code{b} referenced by "\code{net}" }
  \item{FUN}{ Function for which posterior predictive is to be estimated }
  \item{\dots}{ Additional arguments to \code{FUN} }
}
\details{
Although created to work with \code{\link{bbnam}}, \code{npostpred} is quite generic.  The form of the posterior draws will vary with the output of \code{FUN}; since invocation is handled by \code{\link{apply}}, check there if unsure.
}
\value{
A series of posterior predictive draws
}
\references{ Gelman, A.; Carlin, J.B.; Stern, H.S.; and Rubin, D.B.  (1995).  \emph{Bayesian Data Analysis.}  London: Chapman and Hall. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{bbnam}} }

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
#Plot a summary of the posterior predictive of reciprocity
hist(npostpred(b,grecip))
}
\keyword{ models }%-- one or more ...
\keyword{ math }
