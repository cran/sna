\name{gplot.loop}
\alias{gplot.loop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Add Loops to a Plot }
\description{
  \code{gplot.loop} draws a "loop" at a specified location; this is used to designate self-ties in \code{\link{gplot}}.
}
\usage{
gplot.loop(x0, y0, length = 0.1, angle = 10, width = 0.01, col = 1,
    border = 1, lty = 1, offset = 0, sides = 10, arrowhead = TRUE,
    ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x0}{ A vector of x coordinates for points of origin }
  \item{y0}{ A vector of y coordinates for points of origin }
  \item{length}{ Arrowhead length, in current plotting units }
  \item{angle}{ Arrowhead angle (in degrees) }
  \item{width}{ Width for loop body, in current plotting units (can be a vector) }
  \item{col}{ Loop body color (can be a vector) }
  \item{border}{ Loop border color (can be a vector) }
  \item{lty}{ Loop border line type (can be a vector) }
  \item{offset}{ Offset for origin point (can be a vector) }
  \item{sides}{ Number of steps to use in approximating curves (can be a vector) }
  \item{arrowhead}{ Boolean; should arrowheads be used?  (Can be a vector) }
  \item{\dots}{ Additional arguments to \code{\link{polygon}} }
}
\details{
  \code{gplot.loop} is the companion to \code{\link{gplot.arrow}}; like the latter, plot elements produced by \code{gplot.loop} are drawn using \code{\link{polygon}}, and as such are scaled based on the current plotting device.  By default, loops are drawn so as to encircle the top left quadrant of a circular region of radius \code{offset} with origin at \code{x0,y0}.  This is useful for functions like \code{\link{gplot}}, which need to draw loops incident to vertices of varying radii.
}
\value{
  None.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{gplot.arrow}}, \code{\link{gplot}}, \code{\link{polygon}} }
\examples{
#Plot a few polygons with loops
plot(0,0,type="n",xlim=c(-2,2),ylim=c(-2,2),asp=1)
gplot.loop(c(0,0),c(1,-1),col=c(3,2),width=0.05,length=0.4,
  offset=sqrt(2)/4,angle=20,arrowhead=TRUE)
polygon(c(0.25,-0.25,-0.25,0.25,NA,0.25,-0.25,-0.25,0.25), 
    c(1.25,1.25,0.75,0.75,NA,-1.25,-1.25,-0.75,-0.75),col=c(2,3))

}
\keyword{ aplot }% at least one, from doc/KEYWORDS
\keyword{ graphs }% __ONLY ONE__ keyword per line
