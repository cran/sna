\name{diag.remove}
\alias{diag.remove}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Remove the Diagonals of Adjacency Matrices in a Graph Stack }
\description{
  Returns the input graphs, with the diagonal entries removed/replaced as indicated.
}
\usage{
diag.remove(dat, remove.val=NA)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more graphs. }
  \item{remove.val}{ the value with which to replace the existing diagonals }
}
\details{
  \code{diag.remove} is simply a convenient way to apply \code{\link{diag}} to an entire collection of adjacency matrices/\code{network} objects at once.
}
\value{
  The updated graphs.
}

\author{ Carter T. Butts \email{buttsc@uci.edu} }


\seealso{ \code{\link{diag}}, \code{\link{upper.tri.remove}}, \code{\link{lower.tri.remove}} }

\examples{
#Generate a random graph stack
g<-rgraph(3,5)
#Remove the diagonals
g<-diag.remove(g)
}
\keyword{ manip }%-- one or more ...
\keyword{ math }
\keyword{ graphs }
