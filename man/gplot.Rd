\name{gplot}
\alias{gplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Two-Dimensional Visualization of Graphs }
\description{
  \code{gplot} produces a simple two-dimensional plot of graph \code{g} in stack \code{dat}.  A variety of options are available to control vertex placement, display details, color, etc.
}
\usage{
gplot(dat, g = 1, gmode = "digraph", diag = FALSE, label = c(1:dim(dat)[2]),
    coord = NULL, jitter = TRUE, thresh = 0, usearrows = TRUE, mode = "segeo",
    displayisolates = TRUE, interactive = FALSE, boxed.labels = TRUE, 
    xlab = NULL, ylab = NULL, pad = 0.2, vertex.pch = 19, label.cex = 1,
    vertex.cex = 1, label.col = 1, edge.col = 1, vertex.col = 1,
    arrowhead.length = 0.2, edge.lty = 1, edge.lwd = 0, edge.len = 0.5,
    edge.curve = 0.1, edge.steps = 50, diag.size = 0.025, uselen = FALSE,
    usecurve = FALSE, suppress.axes = TRUE, new = TRUE, layout.par = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{ a graph or graph stack.  This data may be valued. }
  \item{g}{ integer indicating the index of the graph which is to be plotted.  By default, \code{g==1}. }
  \item{gmode}{ String indicating the type of graph being evaluated.  \code{"digraph"} indicates that edges should be interpreted as directed; \code{"graph"} indicates that edges are undirected; \code{"twomode"} indicates that data should be interpreted as bimodal (i.e., rows and columns are distinct vertex sets).  \code{gmode} is set to \code{"digraph"} by default. }
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{label}{ a vector of vertex labels, if desired; defaults to the vertex index number. }
  \item{coord}{ user-specified vertex coordinates, in an NCOL(dat)x2 matrix.  Where this is specified, it will override the \code{mode} setting. }
  \item{jitter}{ boolean; should the output be jittered? }
  \item{thresh}{ real number indicating the lower threshold for tie values.  Only ties of value >\code{thresh} are displayed.  By default, \code{thresh}=0.}
  \item{usearrows}{ boolean; should arrows (rather than line segments) be used to indicate edges? }
  \item{mode}{ the vertex placement algorithm; this must correspond to a \code{\link{gplot.layout}} function. }
  \item{displayisolates}{ boolean; should isolates be displayed? }
  \item{interactive}{ boolean; should interactive adjustment of vertex placement be attempted? }
  \item{boxed.labels}{ boolean; place vertex labels within boxes? }
  \item{xlab}{ X axis label. }
  \item{ylab}{ Y axis label. }
  \item{pad}{ amount to pad the plotting range; useful if labels are being clipped. }
  \item{vertex.pch}{ character type for vertices; may be given as a vector, if vertices are to be of different types. }
  \item{label.cex}{ character expansion factor for label text. }
  \item{vertex.cex}{ character expansion factor for vertices; may be given as a vector, if vertices are to be of different types. }
  \item{label.col}{ color for vertex labels; may be given as a vector, if labels are to be of different colors. }
  \item{edge.col}{ color for edges; may be given as a vector or adjacency matrix, if edges are to be of different colors. }
  \item{vertex.col}{ color for vertices; may be given as a vector, if vertices are to be of different colors. }
  \item{arrowhead.length}{ arrowhead length to use with \code{arrow}, if \code{usearrows==TRUE}. }
  \item{edge.lty}{ line type for edges; may be given as a vector or adjacency matrix, if edges are to have different line types. }
  \item{edge.lwd}{ line width scale for edges; if set greater than 0, edge widths are given by \code{edge.lwd*dat}.  May be given as a vector or adjacency matrix, if edges are to have different line widths. }
  \item{edge.len}{ if \code{uselen==TRUE}, curved edge lengths are scaled by \code{edge.len}. }
  \item{edge.curve}{ if \code{usecurve==TRUE}, the extent of edge curvature is controlled by \code{edge.curv}.  May be given as a fixed value, vector, or adjacency matrix, if edges are to have different levels of curvature.  }
  \item{edge.steps}{ for curved edges (including loops), the number of line segments to use for the curve approximation. }
  \item{diag.size}{ loop size, if \code{diag==TRUE}. }
  \item{uselen}{ boolean; should we use \code{edge.len} to rescale edge lengths? }
  \item{usecurve}{ boolean; should we use \code{edge.curve}? }
  \item{suppress.axes}{ boolean; suppress plotting of axes? }
  \item{new}{ boolean; create a new plot?  If \code{new==FALSE}, vertices and edges will be added to the existing plot. }
  \item{layout.par}{ parameters to the \code{\link{gplot.layout}} function specified in \code{mode}. }
  \item{\dots}{ additional arguments to \code{\link{plot}}. }
}
\details{
\code{gplot} is the standard network visualization tool within the \code{sna} library.  By means of clever selection of display parameters, a fair amount of display flexibility can be obtained.  Graph layout -- if not specified directly using \code{coord} -- is determined via one of the various available algorithms.  These should be specified via the \code{mode} argument; see \code{\link{gplot.layout}} for a full list.  User-supplied layout functions are also possible -- see the aforementioned man page for details.

Note that where \code{gmode=="twomode"}, the supplied two-mode matrix is converted to bipartite adjacency form prior to computing coordinates.  If \code{interactive==TRUE}, then the user may modify the initial graph layout by selecting an individual vertex and then clicking on the location to which this vertex is to be moved; this process may be repeated until the layout is satisfactory.
}
\value{
A two-column matrix containing the vertex positions as x,y coordinates.  
}
\references{ Wasserman, S. and Faust, K.  (1994)  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press. }
\author{ 
Carter T. Butts \email{buttsc@uci.edu} 

Alex Montgomery \email{ahm@stanford.edu}
}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{plot}}, \code{\link{gplot.layout}} }
\examples{
gplot(rgraph(5))               #Plot a random graph
gplot(rgraph(5),usecurv=TRUE)  #This time, use curved edges
gplot(rgraph(5),mode="mds")    #Try an alternative layout scheme
}
\keyword{ graphs }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
