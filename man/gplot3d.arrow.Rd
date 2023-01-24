\name{gplot3d.arrow}
\alias{gplot3d.arrow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Add Arrows a Three-Dimensional Plot }
\description{
  \code{gplot3d.arrow} draws an arrow between two pairs of points.
}
\usage{
gplot3d.arrow(a, b, radius, color = "white", alpha = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{ a vector or three-column matrix containing origin X,Y,Z coordinates. }
  \item{b}{ a vector or three-column matrix containing origin X,Y,Z coordinates. }
  \item{radius}{ the arrow radius, in current plotting units. May be a vector, if multiple arrows are to be drawn. }
  \item{color }{ the arrow color.  May be a vector, if multiple arrows are being drawn. }
  \item{alpha}{ alpha (transparency) value(s) for arrows.  (May be a vector.) }
}
\details{
  \code{gplot3d.arrow} draws one or more three-dimensional \dQuote{arrows} from the points given in \code{a} to those given in \code{b}.  Note that the \dQuote{arrows} are really cones, narrowing in the direction of the destination point.  
}
\value{
  None.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{gplot3d}}, \code{\link{gplot3d.loop}} }
%\examples{
%}
\keyword{ aplot }% at least one, from doc/KEYWORDS
\keyword{ graphs }% __ONLY ONE__ keyword per line
