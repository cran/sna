\name{gplot}
\alias{gplot}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Two-Dimensional Visualization of Graphs }
\description{
  \code{gplot} produces a simple two-dimensional plot of graph \code{g} in stack \code{dat}.  A variety of options are available to control vertex placement, display details, color, etc.
}
\usage{
gplot(dat, g=1, gmode="digraph", diag=FALSE, label=c(1:dim(dat)[2]), 
    coord=NULL, jitter=TRUE, thresh=0, usearrows=TRUE, mode="mds", 
    displayisolates=TRUE, pad=0, vertex.pch = 19, label.cex=1, 
    vertex.cex=1, label.col=1, edge.col=1, vertex.col=1, 
    arrowhead.length=0.2, edge.type=1, edge.lwd = 0, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ A graph or graph stack.  This data may be valued. }
  \item{g}{ Integer indicating the index of the graph which is to be plotted.  By default, \code{g==1}. }
  \item{gmode}{ String indicating the type of graph being evaluated.  \code{"digraph"} indicates that edges should be interpreted as directed; \code{"graph"} indicates that edges are undirected; \code{"twomode"} indicates that data should be interpreted as bimodal (i.e., rows and columns are distinct vertex sets).  \code{gmode} is set to \code{"digraph"} by default.}
  \item{diag}{ Boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{label}{ A vector of vertex labels, if desired; defaults to the vertex index number. }
  \item{coord}{ User-specified vertex coordinates, in an nx2 matrix.  Where this is specified, it will override the \code{mode} setting. }
  \item{jitter}{ Should the output be jittered? }
  \item{thresh}{ Real number indicating the lower threshold for tie values.  Only ties of value >\code{thresh} are displayed.  By default, \code{thresh}=0. }
  \item{usearrows}{ Should arrows (rather than line segments) be used to indicate edges? }
  \item{mode}{ The vertex placement algorithm.  This must be one of \code{"princoord"}, \code{"eigen"}, \code{"mds"}, \code{"random"}, \code{"circle"}, \code{"circrand"}, \code{"rmds"}, \code{"geodist"}, \code{"adj"}, or \code{"seham"}. }
  \item{displayisolates}{ Should we display isolates? }
  \item{pad}{ Amount to pad the plotting range; useful if labels are being clipped }
  \item{vertex.pch}{ Character type for vertices; may be given as a vector, if vertices are to be of different types }
  \item{label.cex}{ Character expansion factor for label text }
  \item{vertex.cex}{ Character expansion factor for vertices }
  \item{label.col}{ A vector of colors for labels }
  \item{edge.col}{ A vector of colors for edges }
  \item{vertex.col}{ A vector of colors for vertices }
  \item{arrowhead.length}{ If arrows are used for edges, how long should the arrowhead be? }
  \item{edge.type}{ Line type for edges }
  \item{edge.lwd}{ Line width scale for edges; if set greater than 0, edge widths are given by \code{edge.lwd*dat} }
  \item{\dots}{ Additional arguments to \code{\link{plot}} }
}
\details{
\code{gplot} is the standard network visualization tool within the \code{sna} library.  By means of clever selection of display parameters, a fair amount of display flexibility can be obtained.  
}
\value{
None.
}
\references{ Wasserman, S., and Faust, K.  (1994).  ``Social Network Analysis: Methods and Applications.''  Cambridge: Cambridge University Press. }
\author{ Carter T. Butts \email{ctb@andrew.cmu.edu} }
%\note{  }

 %~Make other sections like WARNING with \section{WARNING }{....} ~
\section{Requires}{\code{mva}}

\seealso{\code{\link{plot}} }

\examples{
gplot(rgraph(5))  #Plot a random graph
}
\keyword{ hplot }%-- one or more ...
