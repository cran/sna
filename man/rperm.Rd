\name{rperm}
\alias{rperm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Draw a Random Permutation Vector with Exchangeability Constraints }
\description{
  Draws a random permutation on \code{1:length(exchange.list)} such that no two elements whose corresponding \code{exchange.list} values are different are interchanged.
}
\usage{
rperm(exchange.list)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{exchange.list}{ A vector such that the permutation vector may exchange the ith and jth positions iff \code{exchange.list[i]==exchange.list[j]} }
}
\details{
\code{rperm} draws random permutation vectors given the constraints of exchangeability described above.  Thus, \code{rperm(c(0,0,0,0))} returns a random permutation of four elements in which all exchanges are allowed, while \code{rperm(c(1,1,"a","a")} (or similar) returns a random permutation of four elements in which only the first/second and third/fourth elements may be exchanged.  This turns out to be quite useful for searching permutation spaces with exchangeability constraints (e.g., for structural distance estimation).
}
\value{
A random permutation vector satisfying the given constraints
}
%\references{ }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

 %~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{rmperm}} }

\examples{
rperm(c(0,0,0,0))  #All elements may be exchanged
rperm(c(0,0,0,1))  #Fix the fourth element
rperm(c(0,0,1,1))  #Allow \{1,2\} and \{3,4\} to be swapped
rperm(c("a",4,"x",2))  #Fix all elements (the identity permutation)
}
\keyword{ distribution }
\keyword{ array }%-- one or more ...
