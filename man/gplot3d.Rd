\name{gplot3d}
\alias{gplot3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Three-Dimensional Visualization of Graphs }
\description{
  \code{gplot3d} produces a three-dimensional plot of graph \code{g} in set \code{dat}.  A variety of options are available to control vertex placement, display details, color, etc.
}
\usage{
gplot3d(dat, g = 1, gmode = "digraph", diag = FALSE, 
    label = NULL, coord = NULL, jitter = TRUE, thresh = 0,
    mode = "fruchtermanreingold", displayisolates = TRUE, 
    displaylabels = !missing(label), xlab = NULL, ylab = NULL, 
    zlab = NULL, vertex.radius = NULL, absolute.radius = FALSE, 
    label.col = "gray50", edge.col = "black", vertex.col = NULL, 
    edge.alpha = 1, vertex.alpha = 1, edge.lwd = NULL, suppress.axes = TRUE,
    new = TRUE, bg.col = "white", layout.par = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{ a graph or set thereof.  This data may be valued. }
  \item{g}{ integer indicating the index of the graph (from \code{dat}) which is to be displayed. }
  \item{gmode}{ string indicating the type of graph being evaluated. \code{"digraph"} indicates that edges should be interpreted as directed; \code{"graph"} indicates that edges are undirected;\code{"twomode"} indicates that data should be interpreted as two-mode (i.e., rows and columns are distinct vertex sets). }
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  }
  \item{label}{ a vector of vertex labels; setting this to a zero-length string (e.g., \code{""}) omits }
  \item{coord}{ user-specified vertex coordinates, in an \code{NCOL(dat)}x3 matrix.  Where this is specified, it will override the \code{mode} setting. }
  \item{jitter}{ boolean; should vertex positions be jittered? }
  \item{thresh}{ real number indicating the lower threshold for tie values.  Only ties of value >\code{thresh} are displayed.  }
  \item{mode}{ the vertex placement algorithm; this must correspond to a \code{gplot3d.layout} function. }
  \item{displayisolates}{ boolean; should isolates be displayed? }
  \item{displaylabels}{ boolean; should vertex labels be displayed? }
  \item{xlab}{ X axis label. }
  \item{ylab}{ Y axis label. }
  \item{zlab}{ Z axis label. }
  \item{vertex.radius}{ vertex radius, relative to the baseline (which is set based on layout features); may be given as a vector, if radii vary across vertices. }
  \item{absolute.radius}{ vertex radius, specified in absolute terms; this may be given as a vector. }
  \item{label.col}{ color for vertex labels; may be given as a vector, if labels are to be of different colors. }
  \item{edge.col}{ color for edges; may be given as a vector or adjacency matrix, if edges are to be of different colors. }
  \item{vertex.col}{ color for vertices; may be given as a vector, if vertices are to be of different colors.   By default, red is used (or red and blue, if \code{gmode=="twomode"}).}
  \item{edge.alpha}{ alpha (transparency) values for edges; may be given as a vector or adjacency matrix, if edge transparency is to vary. }
  \item{vertex.alpha}{ alpha (transparency) values for vertices; may be given as a vector, if vertex transparency is to vary. }
  \item{edge.lwd}{ line width scale for edges; if set greater than 0, edge widths are rescaled by \code{edge.lwd*dat}.  May be given as a vector or adjacency matrix, if edges are to have different line widths. }
  \item{suppress.axes}{ boolean; suppress plotting of axes? }
  \item{new}{ boolean; create a new plot?  If \code{new==FALSE}, the RGL device will not be cleared prior to adding vertices and edges. }
  \item{bg.col}{ background color for display. }
  \item{layout.par}{ list of parameters to the \code{\link{gplot.layout}} function specified in \code{mode}. }
}
\details{
\code{gplot3d} is the three-dimensional companion to \code{gplot}.  As with the latter, clever manipulation of parameters can allow for a great deal of flexibility in the resulting display.  (Displays produced by \code{gplot3d} are also interactive, to the extent supported by \code{\link[rgl:rgl-package]{rgl}}.)  If vertex positions are not specified directly using \code{coord}, vertex layout is determined via one of the various available algorithms.  These should be specified via the \code{mode} argument; see \code{\link{gplot3d.layout}} for a full list.  User-supplied layout functions are also possible - see the aforementioned man page for details.

Note that where \code{gmode=="twomode"}, the supplied two-mode graph is converted to bipartite form prior to computing coordinates (assuming it is not in this form already).  It may be desirable to use parameters such as \code{vertex.col} to differentiate row and column vertices; by default, row vertices are colored red, and column vertices blue.
}
\value{
  A three-column matrix containing vertex coordinates
}
\references{ 
Wasserman, S. and Faust, K.  (1994)  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press.
}
\author{Carter T. Butts \email{buttsc@uci.edu}}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Requires }{\code{\link[rgl:rgl-package]{rgl}}}

\seealso{ \code{\link{gplot}}, \code{\link{gplot3d.layout}}, \code{\link[rgl:rgl-package]{rgl}} }
\examples{

\dontrun{
#A three-dimensional grid...
gplot3d(rgws(1,5,3,1,0))

#...rewired...
gplot3d(rgws(1,5,3,1,0.05))

#...some more!
gplot3d(rgws(1,5,3,1,0.2))
}
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ graphs }% __ONLY ONE__ keyword per line
