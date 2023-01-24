\name{gplot.layout}
\alias{gplot.layout}
\alias{gplot.layout.adj}
\alias{gplot.layout.circle}
\alias{gplot.layout.circrand}
\alias{gplot.layout.eigen}
\alias{gplot.layout.fruchtermanreingold}
\alias{gplot.layout.geodist}
\alias{gplot.layout.hall}
\alias{gplot.layout.kamadakawai}
\alias{gplot.layout.mds}
\alias{gplot.layout.princoord}
\alias{gplot.layout.random}
\alias{gplot.layout.rmds}
\alias{gplot.layout.segeo}
\alias{gplot.layout.seham}
\alias{gplot.layout.spring}
\alias{gplot.layout.springrepulse}
\alias{gplot.layout.target}
\alias{gplot_layout_fruchtermanreingold_R}
\alias{gplot_layout_fruchtermanreingold_old_R}
\alias{gplot_layout_kamadakawai_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Vertex Layout Functions for gplot }
\description{
  Various functions which generate vertex layouts for the \code{\link{gplot}} visualization routine.
}
\usage{
gplot.layout.adj(d, layout.par)
gplot.layout.circle(d, layout.par)
gplot.layout.circrand(d, layout.par)
gplot.layout.eigen(d, layout.par)
gplot.layout.fruchtermanreingold(d, layout.par)
gplot.layout.geodist(d, layout.par)
gplot.layout.hall(d, layout.par)
gplot.layout.kamadakawai(d, layout.par)
gplot.layout.mds(d, layout.par)
gplot.layout.princoord(d, layout.par)
gplot.layout.random(d, layout.par)
gplot.layout.rmds(d, layout.par)
gplot.layout.segeo(d, layout.par)
gplot.layout.seham(d, layout.par)
gplot.layout.spring(d, layout.par)
gplot.layout.springrepulse(d, layout.par)
gplot.layout.target(d, layout.par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{ an adjacency matrix, as passed by \code{\link{gplot}}. }
  \item{layout.par}{ a list of parameters. }
}
\details{
  Vertex layouts for network visualization pose a difficult problem -- there is no single, ``good'' layout algorithm, and many different approaches may be valuable under different circumstances.  With this in mind, \code{\link{gplot}} allows for the use of arbitrary vertex layout algorithms via the \code{gplot.layout.*} family of routines.  When called, \code{\link{gplot}} searches for a \code{gplot.layout} function whose third name matches its \code{mode} argument (see \code{\link{gplot}} help for more information); this function is then used to generate the layout for the resulting plot.  In addition to the routines documented here, users may add their own layout functions as needed.  The requirements for a \code{gplot.layout} function are as follows:
  \enumerate{
  \item the first argument, \code{d}, must be the (dichotomous) graph adjacency matrix;
  \item the second argument, \code{layout.par}, must be a list of parameters (or \code{NULL}, if no parameters are specified); and
  \item the return value must be a real matrix of dimension \code{c(2,NROW(d))}, whose rows contain the vertex coordinates.
  }
  Other than this, anything goes.  (In particular, note that \code{layout.par} could be used to pass additional matrices, if needed.)  
  
  The \code{graph.layout} functions currently supplied by default are as follows:
  \describe{
  \item{circle}{ This function places vertices uniformly in a circle; it takes no arguments.}
  \item{eigen}{ This function places vertices based on the eigenstructure of the adjacency matrix.  It takes the following arguments:
    \describe{
    \item{\code{layout.par$var}}{ This argument controls the matrix to be used for the eigenanalysis.  \code{"symupper"}, \code{"symlower"}, \code{"symstrong"}, \code{"symweak"} invoke \code{\link{symmetrize}} on \code{d} with the respective symmetrizing rule.  \code{"user"} indicates a user-supplied matrix (see below), while \code{"raw"} indicates that \code{d} should be used as-is.  (Defaults to \code{"raw"}.)}
    \item{\code{layout.par$evsel}}{ If \code{"first"}, the first two eigenvectors are used; if \code{"size"}, the two eigenvectors whose eigenvalues have the largest magnitude are used instead. Note that only the real portion of the associated eigenvectors is used.  (Defaults to \code{"first"}.)}
    \item{\code{layout.par$mat}}{  If \code{layout.par$var=="user"}, this matrix is used for the eigenanalysis. (No default.)}
    }
  }
  \item{fruchtermanreingold}{ This function generates a layout using a variant of Fruchterman and Reingold's force-directed placement algorithm.  It takes the following arguments:
    \describe{
    \item{\code{layout.par$niter}}{ This argument controls the number of iterations to be employed.  Larger values take longer, but will provide a more refined layout.  (Defaults to 500.) }
    \item{\code{layout.par$max.delta}}{ Sets the maximum change in position for any given iteration.   (Defaults to \code{n}.)}
    \item{\code{layout.par$area}}{  Sets the \dQuote{area} parameter for the F-R algorithm. (Defaults to \code{n^2}.)}
    \item{\code{layout.par$cool.exp}}{  Sets the cooling exponent for the annealer. (Defaults to 3.)}
    \item{\code{layout.par$repulse.rad}}{  Determines the radius at which vertex-vertex repulsion cancels out attraction of adjacent vertices.  (Defaults to \code{area*log(n)}.)}
    \item{\code{layout.par$ncell}}{  To speed calculations on large graphs, the plot region is divided at each iteration into \code{ncell} by \code{ncell} \dQuote{cells}, which are used to define neighborhoods for force calculation.  Moderate numbers of cells result in fastest performance; too few cells (down to 1, which produces \dQuote{pure} F-R results) can yield odd layouts, while too many will result in long layout times.  (Defaults to \code{n^0.5}.)}
    \item{\code{layout.par$cell.jitter}}{  Jitter factor (in units of cell width) used in assigning vertices to cells.  Small values may generate \dQuote{grid-like} anomalies for graphs with many isolates.  (Defaults to \code{0.5}.)}
    \item{\code{layout.par$cell.pointpointrad}}{  Squared \dQuote{radius} (in units of cells) such that exact point interaction calculations are used for all vertices belonging to any two cells less than or equal to this distance apart.  Higher values approximate the true F-R solution, but increase computational cost.  (Defaults to \code{0}.)}
    \item{\code{layout.par$cell.pointcellrad}}{  Squared \dQuote{radius} (in units of cells) such that approximate point/cell interaction calculations are used for all vertices belonging to any two cells less than or equal to this distance apart (and not within the point/point radius).  Higher values provide somewhat better approximations to the true F-R solution at slightly increased computational cost.  (Defaults to \code{18}.)}
    \item{\code{layout.par$cell.cellcellrad}}{  Squared \dQuote{radius} (in units of cells) such that approximate cell/cell interaction calculations are used for all vertices belonging to any two cells less than or equal to this distance apart (and not within the point/point or point/cell radii).  Higher values provide somewhat better approximations to the true F-R solution at slightly increased computational cost.  Note that cells beyond this radius (if any) do not interact, save through edge attraction. (Defaults to \code{ncell^2}.)}
    \item{\code{layout.par$seed.coord}}{ A two-column matrix of initial vertex coordinates.  (Defaults to a random circular layout.) }  
    }
  }
  \item{hall}{ This function places vertices based on the last two eigenvectors of the Laplacian of the input matrix (Hall's algorithm).  It takes no arguments.}
  \item{kamadakawai}{ This function generates a vertex layout using a version of the Kamada-Kawai force-directed placement algorithm.  It takes the following arguments:
    \describe{
    \item{\code{layout.par$niter}}{ This argument controls the number of iterations to be employed.  (Defaults to 1000.) }
    \item{\code{layout.par$sigma}}{ Sets the base standard deviation of position change proposals.   (Defaults to \code{NROW(d)/4}.)}
    \item{\code{layout.par$initemp}}{  Sets the initial "temperature" for the annealing algorithm. (Defaults to 10.)}
    \item{\code{layout.par$cool.exp}}{  Sets the cooling exponent for the annealer. (Defaults to 0.99.)}
    \item{\code{layout.par$kkconst}}{  Sets the Kamada-Kawai vertex attraction constant.  (Defaults to \code{NROW(d)^2}.)}
    \item{\code{layout.par$elen}}{  Provides the matrix of interpoint distances to be approximated.  (Defaults to the geodesic distances of \code{d} after symmetrizing, capped at \code{sqrt(NROW(d))}.)}
    \item{\code{layout.par$seed.coord}}{ A two-column matrix of initial vertex coordinates.  (Defaults to a gaussian layout.) }  
    }
  }
  \item{mds}{ This function places vertices based on a metric multidimensional scaling of a specified distance matrix.  It takes the following arguments:
    \describe{
    \item{\code{layout.par$var}}{ This argument controls the raw variable matrix to be used for the subsequent distance calculation and scaling.  \code{"rowcol"}, \code{"row"}, and \code{"col"} indicate that the rows and columns (concatenated), rows, or columns (respectively) of \code{d} should be used.  \code{"rcsum"} and \code{"rcdiff"} result in the sum or difference of \code{d} and its transpose being employed.  \code{"invadj"} indicates that \code{max{d}-d} should be used, while \code{"geodist"} uses \code{\link{geodist}} to generate a matrix of geodesic distances from \code{d}. Alternately, an arbitrary matrix can be provided using \code{"user"}.  (Defaults to \code{"rowcol"}.)}
    \item{\code{layout.par$dist}}{ The distance function to be calculated on the rows of the variable matrix.  This must be one of the \code{method} parameters to \code{\link{dist}} (\code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, or \code{"canberra"}), or else \code{"none"}.  In the latter case, no distance function is calculated, and the matrix in question must be square (with dimension \code{dim(d)}) for the routine to work properly.  (Defaults to \code{"euclidean"}.)}
    \item{\code{layout.par$exp}}{ The power to which distances should be raised prior to scaling.  (Defaults to 2.)}
    \item{\code{layout.par$vm}}{ If \code{layout.par$var=="user"}, this matrix is used for the distance calculation. (No default.)}
    }
  
  Note: the following layout functions are based on \code{mds}:
    \describe{
    \item{adj}{ scaling of the raw adjacency matrix, treated as similarities (using \code{"invadj"}).}
    \item{geodist}{ scaling of the matrix of geodesic distances.}
    \item{rmds}{ euclidean scaling of the rows of \code{d}.}
    \item{segeo}{ scaling of the squared euclidean distances between row-wise geodesic distances (i.e., approximate structural equivalence).}
    \item{seham}{ scaling of the Hamming distance between rows/columns of \code{d} (i.e., another approximate structural equivalence scaling).}
    }
  }
  \item{princoord}{ This function places vertices based on the eigenstructure of a given correlation/covariance matrix. It takes the following arguments:
    \describe{
    \item{\code{layout.par$var}}{ The matrix of variables to be used for the correlation/covariance calculation.  \code{"rowcol"}, \code{"col"}, and \code{"row"} indicate that the rows/cols, columns, or rows (respectively) of \code{d} should be employed.  \code{"rcsum"} \code{"rcdiff"} result in the sum or difference of \code{d} and \code{t(d)} being used.  \code{"user"} allows for an arbitrary variable matrix to be supplied.  (Defaults to \code{"rowcol"}.)}
    \item{\code{layout.par$cor}}{ Should the correlation matrix (rather than the covariance matrix) be used?  (Defaults to \code{TRUE}.)}
    \item{\code{layout.par$vm}}{ If \code{layout.par$var=="user"}, this matrix is used for the correlation/covariance calculation. (No default.)}
    }
  }
  \item{random}{ This function places vertices randomly.  It takes the following argument:
    \describe{
    \item{\code{layout.par$dist}}{ The distribution to be used for vertex placement.  Currently, the options are \code{"unif"} (for uniform distribution on the square), \code{"uniang"} (for a ``gaussian donut'' configuration), and \code{"normal"} (for a straight Gaussian distribution).  (Defaults to \code{"unif"}.)}
    }
  Note: \code{circrand}, which is a frontend to the \code{"uniang"} option, is based on this function.
  }
  \item{spring}{ This function places vertices using a spring embedder.  It takes the following arguments:
    \describe{
    \item{\code{layout.par$mass}}{ The vertex mass (in ``quasi-kilograms'').  (Defaults to \code{0.1}.)}
    \item{\code{layout.par$equil}}{ The equilibrium spring extension (in ``quasi-meters''). (Defaults to \code{1}.)}
    \item{\code{layout.par$k}}{ The spring coefficient (in ``quasi-Newtons per quasi-meter''). (Defaults to \code{0.001}.)}
    \item{\code{layout.par$repeqdis}}{ The point at which repulsion (if employed) balances out the spring extension force (in ``quasi-meters''). (Defaults to \code{0.1}.)}
    \item{\code{layout.par$kfr}}{ The base coefficient of kinetic friction (in ``quasi-Newton quasi-kilograms''). (Defaults to \code{0.01}.)}
    \item{\code{layout.par$repulse}}{ Should repulsion be used?  (Defaults to \code{FALSE}.)}
    }
  Note: \code{springrepulse} is a frontend to \code{spring}, with repulsion turned on.
  }
  \item{target}{ This function produces a "target diagram" or "bullseye" layout, using a Brandes et al.'s force-directed placement algorithm.  (See also \code{\link{gplot.target}}.)  It takes the following arguments:
    \describe{
    \item{\code{layout.par$niter}}{ This argument controls the number of iterations to be employed.  (Defaults to 1000.) }
    \item{\code{layout.par$radii}}{ This argument should be a vector of length \code{NROW(d)} containing vertex radii.  Ideally, these should lie in the [0,1] interval (and odd behavior may otherwise result).  (Defaults to the affine-transformed Freeman \code{\link{degree}} centrality scores of \code{d}.) }
    \item{\code{layout.par$minlen}}{ Sets the minimum edge length, below which edge lengths are to be adjusted upwards.  (Defaults to 0.05.) }
    \item{\code{layout.par$area}}{  Sets the initial "temperature" for the annealing algorithm. (Defaults to 10.)}
    \item{\code{layout.par$cool.exp}}{  Sets the cooling exponent for the annealer. (Defaults to 0.99.)}
    \item{\code{layout.par$maxdelta}}{  Sets the maximum angular distance for vertex moves.  (Defaults to \code{pi}.)}
    \item{\code{layout.par$periph.outside}}{  Boolean; should "peripheral" vertices (in the Brandes et al. sense) be placed together outside the main target area?  (Defaults to \code{FALSE}.)}
    \item{\code{layout.par$periph.outside.offset}}{  Radius at which to place "peripheral" vertices if \code{periph.outside==TRUE}.  (Defaults to 1.2.)}
    \item{\code{layout.par$disconst}}{  Multiplier for the Kamada-Kawai-style distance potential.  (Defaults to 1.)}
    \item{\code{layout.par$crossconst}}{  Multiplier for the edge crossing potential.  (Defaults to 1.)}
    \item{\code{layout.par$repconst}}{  Multiplier for the vertex-edge repulsion potential.  (Defaults to 1.)}
    \item{\code{layout.par$minpdis}}{  Sets the "minimum distance" parameter for vertex repulsion.   (Defaults to 0.05.)}
    }
  }

  }
}
\value{
A matrix whose rows contain the x,y coordinates of the vertices of \code{d}.
}
\references{
Brandes, U.; Kenis, P.; and Wagner, D.  (2003).  \dQuote{Communicating Centrality in Policy Network Drawings.} \emph{IEEE Transactions on Visualization and Computer Graphics,} 9(2):241-253.

Fruchterman, T.M.J. and Reingold, E.M.  (1991).  \dQuote{Graph Drawing by Force-directed Placement.} \emph{Software - Practice and Experience,} 21(11):1129-1164.

Kamada, T. and Kawai, S.  (1989).  \dQuote{An Algorithm for Drawing General Undirected Graphs.} \emph{Information Processing Letters,} 31(1):7-15.
  }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{gplot}}, \code{\link{gplot.target}}, \code{\link{gplot3d.layout}}, \code{\link{cmdscale}}, \code{\link{eigen}} }
%\examples{
%}
\keyword{ graphs }% at least one, from doc/KEYWORDS
\keyword{ dplot }% __ONLY ONE__ keyword per line
