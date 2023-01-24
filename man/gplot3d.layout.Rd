\name{gplot3d.layout}
\alias{gplot3d.layout}
\alias{gplot3d.layout.adj}
\alias{gplot3d.layout.eigen}
\alias{gplot3d.layout.fruchtermanreingold}
\alias{gplot3d.layout.geodist}
\alias{gplot3d.layout.hall}
\alias{gplot3d.layout.kamadakawai}
\alias{gplot3d.layout.mds}
\alias{gplot3d.layout.princoord}
\alias{gplot3d.layout.random}
\alias{gplot3d.layout.rmds}
\alias{gplot3d.layout.segeo}
\alias{gplot3d.layout.seham}
\alias{gplot3d_layout_fruchtermanreingold_R}
\alias{gplot3d_layout_kamadakawai_R}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Vertex Layout Functions for gplot3d }
\description{
  Various functions which generate vertex layouts for the \code{\link{gplot3d}} visualization routine.
}
\usage{
gplot3d.layout.adj(d, layout.par)
gplot3d.layout.eigen(d, layout.par)
gplot3d.layout.fruchtermanreingold(d, layout.par)
gplot3d.layout.geodist(d, layout.par)
gplot3d.layout.hall(d, layout.par)
gplot3d.layout.kamadakawai(d, layout.par)
gplot3d.layout.mds(d, layout.par)
gplot3d.layout.princoord(d, layout.par)
gplot3d.layout.random(d, layout.par)
gplot3d.layout.rmds(d, layout.par)
gplot3d.layout.segeo(d, layout.par)
gplot3d.layout.seham(d, layout.par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{ an adjacency matrix, as passed by \code{\link{gplot3d}}. }
  \item{layout.par}{ a list of parameters. }
}
\details{
  Like \code{\link{gplot}}, \code{\link{gplot3d}} allows for the use of arbitrary vertex layout algorithms via the \code{gplot3d.layout.*} family of routines.  When called, \code{\link{gplot3d}} searches for a \code{gplot3d.layout} function whose third name matches its \code{mode} argument (see \code{\link{gplot3d}} help for more information); this function is then used to generate the layout for the resulting plot.  In addition to the routines documented here, users may add their own layout functions as needed.  The requirements for a \code{gplot3d.layout} function are as follows:
  \enumerate{
  \item the first argument, \code{d}, must be the (dichotomous) graph adjacency matrix;
  \item the second argument, \code{layout.par}, must be a list of parameters (or \code{NULL}, if no parameters are specified); and
  \item the return value must be a real matrix of dimension \code{c(3,NROW(d))}, whose rows contain the vertex coordinates.
  }
  Other than this, anything goes.  (In particular, note that \code{layout.par} could be used to pass additional matrices, if needed.)  
  
  The \code{gplot3d.layout} functions currently supplied by default are as follows:
  \describe{
  \item{eigen}{ This function places vertices based on the eigenstructure of the adjacency matrix.  It takes the following arguments:
    \describe{
    \item{\code{layout.par$var}}{ This argument controls the matrix to be used for the eigenanalysis.  \code{"symupper"}, \code{"symlower"}, \code{"symstrong"}, \code{"symweak"} invoke \code{\link{symmetrize}} on \code{d} with the respective symmetrizing rule.  \code{"user"} indicates a user-supplied matrix (see below), while \code{"raw"} indicates that \code{d} should be used as-is.  (Defaults to \code{"raw"}.)}
    \item{\code{layout.par$evsel}}{ If \code{"first"}, the first three eigenvectors are used; if \code{"size"}, the three eigenvectors whose eigenvalues have the largest magnitude are used instead. Note that only the real portion of the associated eigenvectors is used.  (Defaults to \code{"first"}.)}
    \item{\code{layout.par$mat}}{  If \code{layout.par$var=="user"}, this matrix is used for the eigenanalysis. (No default.)}
    }
  }
  \item{fruchtermanreingold}{ This function generates a layout using a variant of Fruchterman and Reingold's force-directed placement algorithm.  It takes the following arguments:
    \describe{
    \item{\code{layout.par$niter}}{ This argument controls the number of iterations to be employed.  (Defaults to 300.) }
    \item{\code{layout.par$max.delta}}{ Sets the maximum change in position for any given iteration.   (Defaults to \code{NROW(d)}.)}
    \item{\code{layout.par$volume}}{  Sets the "volume" parameter for the F-R algorithm. (Defaults to \code{NROW(d)^3}.)}
    \item{\code{layout.par$cool.exp}}{  Sets the cooling exponent for the annealer. (Defaults to 3.)}
    \item{\code{layout.par$repulse.rad}}{  Determines the radius at which vertex-vertex repulsion cancels out attraction of adjacent vertices.  (Defaults to \code{volume*NROW(d)}.)}
    \item{\code{layout.par$seed.coord}}{ A three-column matrix of initial vertex coordinates.  (Defaults to a random spherical layout.) }  
    }
  }
  \item{hall}{ This function places vertices based on the last three eigenvectors of the Laplacian of the input matrix (Hall's algorithm).  It takes no arguments.}
  \item{kamadakawai}{ This function generates a vertex layout using a version of the Kamada-Kawai force-directed placement algorithm.  It takes the following arguments:
    \describe{
    \item{\code{layout.par$niter}}{ This argument controls the number of iterations to be employed.  (Defaults to 1000.) }
    \item{\code{layout.par$sigma}}{ Sets the base standard deviation of position change proposals.   (Defaults to \code{NROW(d)/4}.)}
    \item{\code{layout.par$initemp}}{  Sets the initial "temperature" for the annealing algorithm. (Defaults to 10.)}
    \item{\code{layout.par$cool.exp}}{  Sets the cooling exponent for the annealer. (Defaults to 0.99.)}
    \item{\code{layout.par$kkconst}}{  Sets the Kamada-Kawai vertex attraction constant.  (Defaults to \code{NROW(d)^3}.)}
    \item{\code{layout.par$elen}}{  Provides the matrix of interpoint distances to be approximated.  (Defaults to the geodesic distances of \code{d} after symmetrizing, capped at \code{sqrt(NROW(d))}.)}
    \item{\code{layout.par$seed.coord}}{ A three-column matrix of initial vertex coordinates.  (Defaults to a gaussian layout.) }  
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
    \item{\code{layout.par$dist}}{ The distribution to be used for vertex placement.  Currently, the options are \code{"unif"} (for uniform distribution on the unit cube), \code{"uniang"} (for a ``gaussian sphere'' configuration), and \code{"normal"} (for a straight Gaussian distribution).  (Defaults to \code{"unif"}.)}
    }
  }
}
}
\value{
A matrix whose rows contain the x,y,z coordinates of the vertices of \code{d}.
}
\references{

Fruchterman, T.M.J. and Reingold, E.M.  (1991).  \dQuote{Graph Drawing by Force-directed Placement.} \emph{Software - Practice and Experience,} 21(11):1129-1164.

Kamada, T. and Kawai, S.  (1989).  \dQuote{An Algorithm for Drawing General Undirected Graphs.} \emph{Information Processing Letters,} 31(1):7-15.
  }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{gplot3d}}, \code{\link{gplot}}, \code{\link{gplot.layout}}, \code{\link{cmdscale}}, \code{\link{eigen}} }
%\examples{
%}
\keyword{ graphs }% at least one, from doc/KEYWORDS
\keyword{ dplot }% __ONLY ONE__ keyword per line
