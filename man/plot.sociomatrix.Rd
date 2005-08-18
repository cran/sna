\name{plot.sociomatrix}
\alias{plot.sociomatrix}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Plot Matrices Using a Color/Intensity Grid }
\description{
  Plots a matrix, \code{m}, associating the magnitude of the i,jth cell of \code{m} with the color of the i,jth cell of an \code{nrow(m)} by \code{ncol(m)} grid.
}
\usage{
plot.sociomatrix(x, labels=NULL, drawlab=TRUE, diaglab=TRUE, 
    drawlines=TRUE, xlab=NULL, ylab=NULL, cex.lab=1, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ an input graph. }
  \item{labels}{ a list containing the vectors of row and column labels (respectively); defaults to numerical labels. }
  \item{drawlab}{ logical; add row/column labels to the plot? }
  \item{diaglab}{ logical; label the diagonal? }
  \item{drawlines}{ logical; draw lines to mark cell boundaries? }
  \item{xlab}{ x axis label. }
  \item{ylab}{ y axis label. }
  \item{cex.lab}{ optional expansion factor for labels. }
  \item{\dots}{ additional arguments to \code{\link{plot}}. }
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

\seealso{ \code{\link{plot.blockmodel}} }

\examples{
#Plot a small adjacency matrix
plot.sociomatrix(rgraph(5))

#Plot a much larger one
plot.sociomatrix(rgraph(100),drawlab=FALSE,diaglab=FALSE)
}
\keyword{ hplot }%-- one or more ...
\keyword{graphs}
