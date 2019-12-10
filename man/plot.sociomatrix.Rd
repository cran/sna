\name{plot.sociomatrix}
\alias{plot.sociomatrix}
\alias{sociomatrixplot}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Plot Matrices Using a Color/Intensity Grid }
\description{
  Plots a matrix, \code{m}, associating the magnitude of the i,jth cell of \code{m} with the color of the i,jth cell of an \code{nrow(m)} by \code{ncol(m)} grid.
}
\usage{
\method{plot}{sociomatrix}(x, labels=NULL, drawlab=TRUE, diaglab=TRUE, 
    drawlines=TRUE, xlab=NULL, ylab=NULL, cex.lab=1, font.lab=1, col.lab=1,
    scale.values=TRUE, cell.col=gray, ...)

sociomatrixplot(x, labels=NULL, drawlab=TRUE, diaglab=TRUE, 
    drawlines=TRUE, xlab=NULL, ylab=NULL, cex.lab=1, font.lab=1, col.lab=1,
    scale.values=TRUE, cell.col=gray, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ an input graph. }
  \item{labels}{ a list containing the vectors of row and column labels (respectively); defaults to the row/column labels of \code{x} (if specified), or otherwise sequential numerical labels. }
  \item{drawlab}{ logical; add row/column labels to the plot? }
  \item{diaglab}{ logical; label the diagonal? }
  \item{drawlines}{ logical; draw lines to mark cell boundaries? }
  \item{xlab}{ x axis label. }
  \item{ylab}{ y axis label. }
  \item{cex.lab}{ optional expansion factor for labels. }
  \item{font.lab}{ optional font specification for labels. }
  \item{col.lab}{ optional color specification for labels. }
  \item{scale.values}{ logical; should cell values be affinely scaled to the [0,1] interval?  (Defaults to \code{TRUE}.) }
  \item{cell.col}{ function taking a vector of cell values as an argument and returning a corresponding vector of colors; defaults to \code{\link{gray}}. }
  \item{\dots}{ additional arguments to \code{\link{plot}}. }
}
\details{
  \code{plot.sociomatrix} is particularly valuable for examining large adjacency matrices, whose structure can be non-obvious otherwise.  \code{sociomatrixplot} is an alias to \code{plot.sociomatrix}, and may eventually supersede it.

  The \code{cell.col} argument can be any function that takes input cell values and returns legal colors; while \code{\link{gray}} will produce an error for cell values outside the [0,1] interval, user-specified functions can be employed to get other effects (see examples below).  Note that, by default, all input cell values are affinely scaled to the [0,1] interval before colors are computed, so \code{scale.values} must be set to \code{FALSE} to allow access to the raw inputs.
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
plot.sociomatrix(rgraph(100), drawlab=FALSE, diaglab=FALSE)

#Example involving a signed, valued graph and custom colors
mycolfun <- function(z){   #Custom color function
    ifelse(z<0, rgb(1,0,0,alpha=1-1/(1-z)), ifelse(z>0, 
        rgb(0,0,1,alpha=1-1/(1+z)), rgb(0,0,0,alpha=0)))
}
sg <- rgraph(25) * matrix(rnorm(25^2),25,25)
plot.sociomatrix(sg, scale.values=FALSE, cell.col=mycolfun)  #Blue pos/red neg
}
\keyword{ hplot }%-- one or more ...
\keyword{graphs}
