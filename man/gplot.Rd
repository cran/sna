\name{gplot}
\alias{gplot}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Two-Dimensional Visualization of Graphs }
\description{
  \code{gplot} produces a simple two-dimensional plot of graph \code{g} in stack \code{dat}.  A variety of options are available to control vertex placement, display details, color, etc.
}
\usage{
gplot(dat, g=1, gmode="digraph", diag=FALSE, label=c(1:dim(dat)[2]), 
    coord=NULL, jitter=TRUE, thresh=0, usearrows=TRUE, 
    mode="springrepulse", displayisolates=TRUE, boxed.labels=TRUE, 
    xlab=NULL, ylab=NULL, pad=0.1, vertex.pch = 19, label.cex=1, 
    vertex.cex=1, label.col=1, edge.col=1, vertex.col=1, 
    arrowhead.length=0.2, edge.type=1, edge.lwd = 0, 
    suppress.axes=TRUE, embedder.params=c(0.001,1,0.01,0.2,0.001), ...)
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
  \item{mode}{ The vertex placement algorithm.  This must be one of \code{"princoord"}, \code{"eigen"}, \code{"mds"}, \code{"random"}, \code{"circle"}, \code{"circrand"}, \code{"rmds"}, \code{"geodist"}, \code{"adj"}, \code{"seham"}, \code{"spring"}, or \code{"springrepulse"}. }
  \item{displayisolates}{ Should we display isolates? }
  \item{boxed.labels}{ Place vertex labels within boxes? }
  \item{xlab}{ X axis label }
  \item{ylab}{ Y axis label }
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
  \item{suppress.axes}{ Suppress plotting of axes? }
  \item{embedder.params}{ A vector of parameters for the \code{spring} and \code{springrepulse} modes; see below }
  \item{\dots}{ Additional arguments to \code{\link{plot}} }
}
\details{
\code{gplot} is the standard network visualization tool within the \code{sna} library.  By means of clever selection of display parameters, a fair amount of display flexibility can be obtained.  Graph layout -- if not specified directly using \code{coord} -- is determined via one of the various available algorithms.  These are (briefly) as follows:
\enumerate{

\item \code{random}: Vertices are placed (uniformly) randomly within a square region about the origin.

\item \code{circle}: Vertices are placed evenly about the unit circle.

\item \code{circrand}: Vertices are placed in a ``gaussian donut,'' with distance from the origin following a normal distribution and angle relative to the X axis chosen (uniformly) randomly.

\item \code{eigen}, \code{princoord}: Vertices are placed via (the real components of) the first two eigenvectors of: 
  \enumerate{
  \item \code{eigen}: the matrix of correlations among (concatenated) rows/columns of the adjacency matrix

  \item \code{princoord}: the raw adjacency matrix.
  }

\item \code{mds}, \code{rmds}, \code{geodist}, \code{adj}, \code{seham}: Vertices are placed by a metric MDS.  The distance matrix used is given by:
  \enumerate{
  \item \code{mds}: absolute row/column differences within the adjacency matrix

  \item \code{rmds}: euclidean distances between rows of the adjacency matrix
  
  \item \code{geodist}: geodesic distances between vertices within the graph
  
  \item \code{adj}: \eqn{(\max A)-A}{(max A)-A}, where \eqn{A}{A} is the raw adjacency matrix
  
  \item \code{seham}: structural (dis)equivalence distances (i.e., as per \code{\link{sedist}}) based on the Hamming metric
  }
  
\item \code{spring}, \code{springrepulse}: Vertices are placed using a simple spring embedder.  Parameters for the embedding model are given by \code{embedder.params}, in the following order: vertex mass; equilibrium extension; spring coefficient; repulsion equilibrium distance; and base coefficient of friction.  Initial vertex positions are in random order around a circle, and simulation proceeds -- increasing the coefficient of friction by the specified base value per unit time -- until ``motion'' within the system ceases.  If \code{springrepulse} is specified, then an inverse-cube repulsion force between vertices is also simulated; this force is calibrated so as to be exactly equal to the force of a unit spring extension at a distance specified by the repulsion equilibrium distance.
}

Note that where \code{gmode=="twomode"}, the supplied two-mode matrix is converted to bipartite adjacency form prior to computing coordinates.
}
\value{
None.
}
\references{ Wasserman, S., and Faust, K.  (1994).  ``Social Network Analysis: Methods and Applications.''  Cambridge: Cambridge University Press. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

 %~Make other sections like WARNING with \section{WARNING }{....} ~
\section{Requires}{\code{mva}}

\seealso{\code{\link{plot}} }

\examples{
gplot(rgraph(5))  #Plot a random graph
}
\keyword{ hplot }%-- one or more ...
