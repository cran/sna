\name{gplot3d.loop}
\alias{gplot3d.loop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Add Loops to a Three-Dimensional Plot }
\description{
  \code{gplot3d.loop} draws a "loop" at a specified location; this is used to designate self-ties in \code{\link{gplot3d}}.
}
\usage{
gplot3d.loop(a, radius, color = "white", alpha = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{ a vector or three-column matrix containing origin X,Y,Z coordinates. }
  \item{radius}{ the loop radius, in current plotting units. May be a vector, if multiple loops are to be drawn. }
  \item{color }{ the loop color.  May be a vector, if multiple loops are being drawn. }
  \item{alpha}{ alpha (transparency) value(s) for loops.  (May be a vector.) }
}
\details{
  \code{gplot3d.loop} is the companion to \code{\link{gplot3d.arrow}}.  The "loops" produced by this routine currently look less like loops than like "hats" -- they are noticable as spike-like structures which protrude from vertices.  Eventually, something more attractice will be produced by this routine.
}
\value{
  None.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{gplot3d.arrow}}, \code{\link{gplot3d}}, \code{\link[rgl]{rgl-package}} }
%\examples{
%}
\keyword{ aplot }% at least one, from doc/KEYWORDS
\keyword{ graphs }% __ONLY ONE__ keyword per line
