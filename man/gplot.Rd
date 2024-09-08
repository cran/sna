\name{gplot}
\alias{gplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Two-Dimensional Visualization of Graphs }
\description{
  \code{gplot} produces a two-dimensional plot of graph \code{g} in collection \code{dat}.  A variety of options are available to control vertex placement, display details, color, etc.
}
\usage{
gplot(dat, g = 1, gmode = "digraph", diag = FALSE, 
    label = NULL, coord = NULL, jitter = TRUE, thresh = 0, 
    thresh.absval=TRUE, usearrows = TRUE, mode = "fruchtermanreingold", 
    displayisolates = TRUE, interactive = FALSE, interact.bycomp = FALSE,
    xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, pad = 0.2, 
    label.pad = 0.5, displaylabels = !is.null(label), boxed.labels = FALSE, 
    label.pos = 0, label.bg = "white", vertex.enclose = FALSE, 
    vertex.sides = NULL, vertex.rot = 0, arrowhead.cex = 1, label.cex = 1, 
    loop.cex = 1, vertex.cex = 1, edge.col = 1, label.col = 1, 
    vertex.col = NULL, label.border = 1, vertex.border = 1, edge.lty = NULL,
    edge.lty.neg=2, label.lty = NULL, vertex.lty = 1, edge.lwd = 0, 
    label.lwd = par("lwd"), edge.len = 0.5, edge.curve = 0.1, 
    edge.steps = 50, loop.steps = 20, object.scale = 0.01, uselen = FALSE, 
    usecurve = FALSE, suppress.axes = TRUE, vertices.last = TRUE, 
    new = TRUE, layout.par = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{ a graph or set thereof.  This data may be valued. }
  \item{g}{ integer indicating the index of the graph which is to be plotted.  By default, \code{g==1}. }
  \item{gmode}{ String indicating the type of graph being evaluated.  \code{"digraph"} indicates that edges should be interpreted as directed; \code{"graph"} indicates that edges are undirected; \code{"twomode"} indicates that data should be interpreted as two-mode (i.e., rows and columns are distinct vertex sets).  \code{gmode} is set to \code{"digraph"} by default. }
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{label}{ a vector of vertex labels, if desired; defaults to the vertex index number. }
  \item{coord}{ user-specified vertex coordinates, in an NCOL(dat)x2 matrix.  Where this is specified, it will override the \code{mode} setting. }
  \item{jitter}{ boolean; should the output be jittered? }
  \item{thresh}{ real number indicating the lower threshold for tie values.  Only ties of value >\code{thresh} (by default in absolute value - see \code{thresh.absval})are displayed.  By default, \code{thresh}=0.}
  \item{thresh.absval}{ boolean; should the absolute value of edge weights be used when thresholding?  (Defaults to TRUE; setting to FALSE leads to thresholding by signed weights.)}
  \item{usearrows}{ boolean; should arrows (rather than line segments) be used to indicate edges? }
  \item{mode}{ the vertex placement algorithm; this must correspond to a \code{\link{gplot.layout}} function. }
  \item{displayisolates}{ boolean; should isolates be displayed? }
  \item{interactive}{ boolean; should interactive adjustment of vertex placement be attempted? }
  \item{interact.bycomp}{ boolean; if \code{interactive==TRUE}, should all vertices in the component be moved? }
  \item{xlab}{ x axis label. }
  \item{ylab}{ y axis label. }
  \item{xlim}{ the x limits (min, max) of the plot. }
  \item{ylim}{ the y limits of the plot. } 
  \item{pad}{ amount to pad the plotting range; useful if labels are being clipped. }
  \item{label.pad}{ amount to pad label boxes (if \code{boxed.labels==TRUE}), in character size units. }
  \item{displaylabels}{ boolean; should vertex labels be displayed? }
  \item{boxed.labels}{ boolean; place vertex labels within boxes? }
  \item{label.pos}{ position at which labels should be placed, relative to vertices.  \code{0} results in labels which are placed away from the center of the plotting region; \code{1}, \code{2}, \code{3}, and \code{4} result in labels being placed below, to the left of, above, and to the right of vertices (respectively); and \code{label.pos=5} results in labels which are plotted with no offset (i.e., at the vertex positions). }
  \item{label.bg}{ background color for label boxes (if \code{boxed.labels==TRUE}); may be a vector, if boxes are to be of different colors.}
  \item{vertex.enclose}{ boolean; should vertices be enclosed within circles?  (Can increase legibility for polygonal vertices.)}
  \item{vertex.sides}{ number of polygon sides for vertices; may be given as a vector, if vertices are to be of different types.  By default, 50 sides are used (or 50 and 4, for two-mode data).}
  \item{vertex.rot}{ angle of rotation for vertices (in degrees); may be given as a vector, if vertices are to be rotated differently. }
  \item{arrowhead.cex}{ expansion factor for edge arrowheads.}
  \item{label.cex}{ character expansion factor for label text. }
  \item{loop.cex}{ expansion factor for loops; may be given as a vector, if loops are to be of different sizes. }
  \item{vertex.cex}{ expansion factor for vertices; may be given as a vector, if vertices are to be of different sizes. }
  \item{edge.col}{ color for edges; may be given as a vector or adjacency matrix, if edges are to be of different colors. }
  \item{label.col}{ color for vertex labels; may be given as a vector, if labels are to be of different colors. }
  \item{vertex.col}{ color for vertices; may be given as a vector, if vertices are to be of different colors.  By default, red is used (or red and blue, for two-mode data).}
  \item{label.border}{ label border colors (if \code{boxed.labels==TRUE}); may be given as a vector, if label boxes are to have different colors. }
  \item{vertex.border}{ border color for vertices; may be given as a vector, if vertex borders are to be of different colors. }
  \item{edge.lty}{ line type for (positive weight) edges; may be given as a vector or adjacency matrix, if edges are to have different line types. }
  \item{edge.lty.neg}{ line type for negative weight edges, if any; may be given as per \code{edge.lty}.} 
  \item{label.lty}{ line type for label boxes (if \code{boxed.labels==TRUE}); may be given as a vector, if label boxes are to have different line types. }
  \item{vertex.lty}{ line type for vertex borders; may be given as a vector or adjacency matrix, if vertex borders are to have different line types. }
  \item{edge.lwd}{ line width scale for edges; if set greater than 0, edge widths are scaled by \code{edge.lwd*dat}.  May be given as a vector or adjacency matrix, if edges are to have different line widths. }
  \item{label.lwd}{ line width for label boxes (if \code{boxed.labels==TRUE}); may be given as a vector, if label boxes are to have different line widths. }
  \item{edge.len}{ if \code{uselen==TRUE}, curved edge lengths are scaled by \code{edge.len}. }
  \item{edge.curve}{ if \code{usecurve==TRUE}, the extent of edge curvature is controlled by \code{edge.curv}.  May be given as a fixed value, vector, or adjacency matrix, if edges are to have different levels of curvature.  }
  \item{edge.steps}{ for curved edges (excluding loops), the number of line segments to use for the curve approximation. }
  \item{loop.steps}{ for loops, the number of line segments to use for the curve approximation. }
  \item{object.scale}{ base length for plotting objects, as a fraction of the linear scale of the plotting region. Defaults to 0.01.}
  \item{uselen}{ boolean; should we use \code{edge.len} to rescale edge lengths? }
  \item{usecurve}{ boolean; should we use \code{edge.curve}? }
  \item{suppress.axes}{ boolean; suppress plotting of axes? }
  \item{vertices.last}{ boolean; plot vertices after plotting edges? }
  \item{new}{ boolean; create a new plot?  If \code{new==FALSE}, vertices and edges will be added to the existing plot. }
  \item{layout.par}{ parameters to the \code{\link{gplot.layout}} function specified in \code{mode}. }
  \item{\dots}{ additional arguments to \code{\link{plot}}. }
}
\details{
\code{gplot} is the standard network visualization tool within the \code{sna} library.  By means of clever selection of display parameters, a fair amount of display flexibility can be obtained.  Graph layout -- if not specified directly using \code{coord} -- is determined via one of the various available algorithms.  These should be specified via the \code{mode} argument; see \code{\link{gplot.layout}} for a full list.  User-supplied layout functions are also possible -- see the aforementioned man page for details.

Note that where \code{gmode=="twomode"}, the supplied two-mode network is converted to bipartite form prior to computing coordinates (if not in that form already).  \code{vertex.col} or other settings may be used to differentiate row and column vertices -- by default, row vertices are drawn as red circles, and column vertices are rendered as blue squares.  If \code{interactive==TRUE}, then the user may modify the initial graph layout by selecting an individual vertex and then clicking on the location to which this vertex is to be moved; this process may be repeated until the layout is satisfactory.  If \code{interact.bycomp==TRUE} as well, the vertex and all other vertices in the same component as that vertex are moved together.
}
\value{
A two-column matrix containing the vertex positions as x,y coordinates.  
}
\references{ Wasserman, S. and Faust, K.  (1994)  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press. }
\author{ 
Carter T. Butts \email{buttsc@uci.edu} 

Alex Montgomery \email{ahm@reed.edu}
}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{plot}}, \code{\link{gplot.layout}} }
\examples{
gplot(rgraph(5))               #Plot a random graph
gplot(rgraph(5),usecurv=TRUE)  #This time, use curved edges
gplot(rgraph(5),mode="mds")    #Try an alternative layout scheme

#A colorful demonstration...
gplot(rgraph(5,diag=TRUE),diag=TRUE,vertex.cex=1:5,vertex.sides=3:8,
    vertex.col=1:5,vertex.border=2:6,vertex.rot=(0:4)*72,
    displaylabels=TRUE,label.bg="gray90")
    
}
\keyword{ graphs }% at least one, from doc/KEYWORDS
\keyword{ hplot }% __ONLY ONE__ keyword per line
