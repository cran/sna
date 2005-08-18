\name{rmperm}
\alias{rmperm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Randomly Permute the Rows and Columns of an Input Matrix }
\description{
Given an input matrix (or stack thereof), \code{rmperm} performs a (random) simultaneous row/column permutation of the input data.
}
\usage{
rmperm(m)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{m}{ a matrix, or stack thereof (or a graph set, for that matter). }
}
\details{
Random matrix permutations are the essence of the QAP test; see \code{\link{qaptest}} for details.
}
\value{
The permuted matrix (or matrices)
}
%\references{  }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

 %~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{rperm}} }

\examples{
#Generate an input matrix
g<-rgraph(5)
g             #Examine it

#Examine a random permutation
rmperm(g)
}
\keyword{ array }
\keyword{ distribution }%-- one or more ...
