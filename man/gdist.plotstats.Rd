\name{gdist.plotstats}
\alias{gdist.plotstats}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Plot Various Graph Statistics Over a Network MDS }
\description{
  Plots a two-dimensional metric MDS of \code{d}, with the corresponding values of \code{meas} indicated at each point.  Various options are available for controlling how \code{meas} is to be displayed.
}
\usage{
gdist.plotstats(d, meas, siz.lim=c(0, 0.15), rescale="quantile", 
    display.scale="radius", display.type="circleray", cex=0.5, pch=1,
    labels=NULL, pos=1, labels.cex=1, legend=NULL, legend.xy=NULL, 
    legend.cex=1, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{d}{ A matrix containing the inter-graph distances }
  \item{meas}{ An nxm matrix containing the graph-level measures; each row must correspond to a graph, and each column must correspond to an index }
  \item{siz.lim}{ The minimum and maximum sizes (respectively) of the plotted symbols, given as fractions of the total plotting range }
  \item{rescale}{ One of ``quantile'' for ordinal scaling, ``affine'' for max-min scaling, and ``normalize'' for rescaling by maximum value; these determine the scaling rule to be used in sizing the plotting symbols }
  \item{display.scale}{ One of ``area'' or ``radius''; this controls the attribute of the plotting symbol which is rescaled by the value of \code{meas}}
  \item{display.type}{ One of ``circle'', ``ray'', ``circleray'', ``poly'', or ``polyray''; this determines the type of plotting symbol used (circles, rays, polygons, or come combination of these) }
  \item{cex}{ Character expansion coefficient }
  \item{pch}{ Point types for the base plotting symbol (not the expanded symbols which are used to indicate \code{meas} values) }
  \item{labels}{ Point labels, if desired }
  \item{pos}{ Relative position of labels (see \code{\link{par}}) }
  \item{labels.cex}{ Character expansion factor for labels }
  \item{legend}{ Add a legend? }
  \item{legend.xy}{ x,y coordinates for legend }
  \item{legend.cex}{ Character expansion factor for legend }
  \item{\dots}{ Additional arguments to \code{\link{plot}} }
}
\details{
\code{gdist.plotstats} works by performing an MDS (using \code{\link{cmdscale}}) on \code{d}, and then using the values in \code{meas} to determine the shape of the points at each MDS coordinate.  Typically, these shapes involve rays of varying color and length indicating \code{meas} magnitude, with circles and polygons of the appropriate radius and/or error being options as well.  Various options are available (described above) to govern the details of the data display; some tinkering may be needed in order to produce an aesthetically pleasing visualization.

The primary use of \code{gdist.plotstats} is to explore broad relationships between graph properties and inter-graph distances.  This routine complements others in the \code{gdist} and \code{gclust} family of interstructural visualization tools.
}
\value{
None
}
\references{ Butts, C.T., and Carley, K.M.  (2001).  ``Multivariate Methods for Interstructural Analysis.'' CASOS working paper, Carnegie Mellon University. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ This routine does not actually depend on the data's being graphic in origin, and can be used with any distance matrix/measure matrix combination. }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

%\section{Requires}{\code{mva}}

\seealso{ \code{\link{gdist.plotdiff}}, \code{\link{gclust.boxstats}}, \code{\link{gclust.centralgraph}} }

\examples{
#Generate random graphs with varying density
g<-rgraph(10,20,tprob=runif(20,0,1))

#Get Hamming distances between graphs
g.h<-hdist(g)

#Plot the association of distance, density, and reciprocity
gdist.plotstats(g.h,cbind(gden(g),grecip(g)))
}
\keyword{ hplot }%-- one or more ...
