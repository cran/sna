\name{lower.tri.remove}
\alias{lower.tri.remove}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Remove the Lower Triangles of Adjacency Matrices in a Graph Stack }
\description{
  Returns the input graph set, with the lower triangle entries removed/replaced as indicated.
}
\usage{
lower.tri.remove(dat, remove.val=NA)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{remove.val}{ the value with which to replace the existing lower triangles. }
}
\details{
  \code{lower.tri.remove} is simply a convenient way to apply \code{g[lower.tri(g)]<-remove.val} to an entire stack of adjacency matrices at once.
}
\value{
  The updated graph set.
}
%\references{  }
\author{ Carter T. Butts \email{buttsc@uci.edu}}
%\note{  }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{lower.tri}}, \code{\link{upper.tri.remove}}, \code{\link{diag.remove}} }

\examples{
#Generate a random graph stack
g<-rgraph(3,5)
#Remove the lower triangles
g<-lower.tri.remove(g)

}
\keyword{ manip }%-- one or more ...
\keyword{ math }
\keyword{ graphs }
