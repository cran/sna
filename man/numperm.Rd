\name{numperm}
\alias{numperm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Get the nth Permutation Vector by Periodic Placement }
\description{
  \code{numperm} implicitly numbers all permutations of length \code{olength}, returning the \code{permnum}th of these.
}
\usage{
numperm(olength, permnum)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{olength}{ The number of items to permute }
  \item{permnum}{ The number of the permutation to use (in \code{1:olength!}) }
}
\details{
  The n! permutations on n items can be deterministically ordered via a factorization process in which there are n slots for the first element, n-1 for the second, and n-i for the ith.  This fact is quite handy if you want to visit each permutation in turn, or if you wish to sample without replacement from the set of permutations on some number of elements: one just enumerates or samples from the integers on [1,n!], and then find the associated permutation.  \code{numperm} performs exactly this last operation, returning the \code{permnum}th permutation on \code{olength} items.
}
\value{
A permutation vector
}
%\references{  }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ Permutation search is central to the estimation of structural distances, correlations, and covariances on partially labeled graphs.  \code{numperm} is hence used by \code{\link{structdist}}, \code{\link{gscor}}, \code{\link{gscov}}, etc. }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{rperm}}, \code{\link{rmperm}} }

\examples{
#Draw a graph
g<-rgraph(5)

#Permute the rows and columns
p.1<-numperm(5,1)
p.2<-numperm(5,2)
p.3<-numperm(5,3)
g[p.1,p.1]
g[p.2,p.2]
g[p.3,p.3]
}
\keyword{ array }%-- one or more ...
\keyword{ math }
