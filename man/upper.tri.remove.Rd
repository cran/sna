\name{upper.tri.remove}
\alias{upper.tri.remove}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Remove the Upper Triangles of Adjacency Matrices in a Graph Stack }
\description{
  Returns the input graph stack, with the upper triangle entries removed/replaced as indicated.
}
\usage{
upper.tri.remove(dat, remove.val=NA)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ a graph or graph stack. }
  \item{remove.val}{ the value with which to replace the existing upper triangles. }
}
\details{
  \code{upper.tri.remove} is simply a convenient way to apply \code{g[upper.tri(g)]<-remove.val} to an entire stack of adjacency matrices at once.
}
\value{
  The updated graph stack.
}
%\references{ }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{}

 %~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{upper.tri}}, \code{\link{lower.tri.remove}}, \code{\link{diag.remove}} }

\examples{
#Generate a random graph stack
g<-rgraph(3,5)
#Remove the upper triangles
g<-upper.tri.remove(g)
}
\keyword{ manip }
\keyword{ array }%-- one or more ...
\keyword{ math }
\keyword{ graphs }
