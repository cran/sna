\name{plot.matrix}
\alias{plot.matrix}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Plot Matrices Using a Color/Intensity Grid }
\description{
  Plots a matrix, \code{m}, associating the magnitude of the i,jth cell of \code{m} with the color of the i,jth cell of an \code{nrow(m)} by \code{ncol(m)} grid.
}
\usage{
plot.matrix(m, labels=list(seq(1:dim(m)[1]), seq(1:dim(m)[2])), 
    drawlab=TRUE, diaglab=TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{m}{ An input matrix }
  \item{labels}{ A list containing the vectors of row and column labels (respectively) }
  \item{drawlab}{ Add labels to the plot? }
  \item{diaglab}{ Label the diagonal? }
  \item{\dots}{ Additional arguments to \code{\link{plot}} }
}
\details{
  \code{plot.matrix} is particularly valuable for examining large adjacency matrices, whose structure can be non-obvious otherwise.
}
\value{
None
}
%\references{ }
\author{ Carter T. Butts \email{ctb@andrew.cmu.edu} }
%\note{ }

 %~Make other sections like WARNING with \section{WARNING }{....} ~

%\seealso{  }

\examples{
#Plot a small adjacency matrix
plot.matrix(rgraph(5))

#Plot a much larger one
plot.matrix(rgraph(100),drawlab=F,diaglab=F)
}
\keyword{ hplot }%-- one or more ...
\keyword{ network }