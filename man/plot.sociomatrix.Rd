\name{plot.sociomatrix}
\alias{plot.sociomatrix}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Plot Matrices Using a Color/Intensity Grid }
\description{
  Plots a matrix, \code{m}, associating the magnitude of the i,jth cell of \code{m} with the color of the i,jth cell of an \code{nrow(m)} by \code{ncol(m)} grid.
}
\usage{
plot.sociomatrix(x, labels=list(seq(1:dim(x)[1]), seq(1:dim(x)[2])), 
    drawlab=TRUE, diaglab=TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ An input matrix }
  \item{labels}{ A list containing the vectors of row and column labels (respectively) }
  \item{drawlab}{ Add labels to the plot? }
  \item{diaglab}{ Label the diagonal? }
  \item{\dots}{ Additional arguments to \code{\link{plot}} }
}
\details{
  \code{plot.sociomatrix} is particularly valuable for examining large adjacency matrices, whose structure can be non-obvious otherwise.
}
\value{
None
}
%\references{ }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ }

 %~Make other sections like WARNING with \section{WARNING }{....} ~

%\seealso{  }

\examples{
#Plot a small adjacency matrix
plot.sociomatrix(rgraph(5))

#Plot a much larger one
plot.sociomatrix(rgraph(100),drawlab=FALSE,diaglab=FALSE)
}
\keyword{ hplot }%-- one or more ...
