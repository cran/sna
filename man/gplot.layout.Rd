\name{gplot.layout}
\alias{gplot.layout}
\alias{gplot.layout.princoord}
\alias{gplot.layout.mds}
\alias{gplot.layout.eigen}
\alias{gplot.layout.rmds}
\alias{gplot.layout.geodist}
\alias{gplot.layout.adj}
\alias{gplot.layout.seham}
\alias{gplot.layout.segeo}
\alias{gplot.layout.random}
\alias{gplot.layout.circrand}
\alias{gplot.layout.circle}
\alias{gplot.layout.spring}
\alias{gplot.layout.springrepulse}
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
gplot.layout.geodist(d, layout.par)
gplot.layout.mds(d, layout.par)
gplot.layout.princoord(d, layout.par)
gplot.layout.random(d, layout.par)
gplot.layout.rmds(d, layout.par)
gplot.layout.segeo(d, layout.par)
gplot.layout.seham(d, layout.par)
gplot.layout.spring(d, layout.par)
gplot.layout.springrepulse(d, layout.par)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{ an adjacency matrix, as passed by \code{\link{gplot}}. }
  \item{layout.par}{ a list of parameters. }
}
\details{
  Vertex layouts for network visualization pose a difficult problem -- there is no single, ``good'' layout algorithm, and many different approaches may be valuable under different circumstances.  With this in mind, \code{\link{gplot}} allows for the use of arbitrary vertex layout algorithms via the \code{gplot.layout.*} family of routines.  When called, \code{\link{gplot}} searches for a \code{graph.layout} function whose third name matches its \code{mode} argument (see \code{\link{gplot}} help for more information); this function is then used to generate the layout for the resulting plot.  In addition to the routines documented here, users may add their own layout functions as needed.  The requirements for a \code{graph.layout} function are as follows:
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
    \item{layout.par\$var}{ This argument controls the matrix to be used for the eigenanalysis.  \code{"symupper"}, \code{"symlower"}, \code{"symstrong"}, \code{"symweak"} invoke \code{\link{symmetrize}} on \code{d} with the respective symmetrizing rule.  \code{"user"} indicates a user-supplied matrix (see below), while \code{"raw"} indicates that \code{d} should be used as-is.  (Defaults to \code{"raw"}.)}
    \item{layout.par\$evsel}{ If \code{"first"}, the first two eigenvectors are used; if \code{"size"}, the two eigenvectors whose eigenvalues have the largest magnitude are used instead. Note that only the real portion of the associated eigenvectors is used.  (Defaults to \code{"first"}.)}
    \item{layout.par\$mat}{  If \code{layout.par\$var=="user"}, this matrix is used for the eigenanalysis. (No default.)}
    }
  }
  \item{mds}{ This function places vertices based on a metric multidimensional scaling of a specified distance matrix.  It takes the following arguments:
    \describe{
    \item{layout.par\$var}{ This argument controls the raw variable matrix to be used for the subsequent distance calculation and scaling.  \code{"rowcol"}, \code{"row"}, and \code{"col"} indicate that the rows and columns (concatenated), rows, or columns (respectively) of \code{d} should be used.  \code{"rcsum"} and \code{"rcdiff"} result in the sum or difference of \code{d} and its transpose being employed.  \code{"invadj"} indicates that \code{max{d}-d} should be used, while \code{"geodist"} uses \code{\link{geodist}} to generate a matrix of geodesic distances from \code{d}. Alternately, an arbitrary matrix can be provided using \code{"user"}.  (Defaults to \code{"rowcol"}.)}
    \item{layout.par\$dist}{ The distance function to be calculated on the rows of the variable matrix.  This must be one of the \code{method} parameters to \code{\list{dist}} (\code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, or \code{"canberra"}), or else \code{"none"}.  In the latter case, no distance function is calculated, and the matrix in question must be square (with dimension \code{dim(d)}) for the routine to work properly.  (Defaults to \code{"euclidean"}.)}
    \item{layout.par\$exp}{ The power to which distances should be raised prior to scaling.  (Defaults to 2.)}
    \item{layout.par\$vm}{ If \code{layout.par\$var=="user"}, this matrix is used for the distance calculation. (No default.)}
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
    \item{layout.par\$var}{ The matrix of variables to be used for the correlation/covariance calculation.  \code{"rowcol"}, \code{"col"}, and \code{"row"} indicate that the rows/cols, columns, or rows (respectively) of \code{d} should be employed.  \code{"rcsum"} \code{"rcdiff"} result in the sum or difference of \code{d} and \code{t(d)} being used.  \code{"user"} allows for an arbitrary variable matrix to be supplied.  (Defaults to \code{"rowcol"}.)}
    \item{layout.par\$cor}{ Should the correlation matrix (rather than the covariance matrix) be used?  (Defaults to \code{TRUE}.)}
    \item{layout.par\$vm}{ If \code{layout.par\$var=="user"}, this matrix is used for the correlation/covariance calculation. (No default.)}
    }
  }
  \item{random}{ This function places vertices randomly.  It takes the following argument:
    \describe{
    \item{layout.par\$dist}{ The distribution to be used for vertex placement.  Currently, the options are \code{"unif"} (for uniform distribution on the square), \code{"uniang"} (for a ``gaussian donut'' configuration), and \code{"normal"} (for a straight Gaussian distribution).  (Defaults to \code{"unif"}.)}
    }
  Note: \code{circrand}, which is a frontend to the \code{"uniang"} option, is based on this function.
  }
  \item{spring}{ This function places vertices using a spring embedder.  It takes the following arguments:
    \describe{
    \item{layout.par\$mass}{ The vertex mass (in ``quasi-kilograms'').  (Defaults to \code{0.1}.)}
    \item{layout.par\$equil}{ The equilibrium spring extension (in ``quasi-meters''). (Defaults to \code{1}.)}
    \item{layout.par\$k}{ The spring coefficient (in ``quasi-Newtons per quasi-meter''). (Defaults to \code{0.001}.)}
    \item{layout.par\$repeqdis}{ The point at which repulsion (if employed) balances out the spring extension force (in ``quasi-meters''). (Defaults to \code{0.1}.)}
    \item{layout.par\$kfr}{ The base coefficient of kinetic friction (in ``quasi-Newton quasi-kilograms''). (Defaults to \code{0.01}.)}
    \item{layout.par\$repulse}{ Should repulsion be used?  (Defaults to \code{FALSE}.)}
    }
  Note: \code{springrepulse} is a frontend to \code{spring}, with repulsion turned on.
  }
  }
}
\value{
A matrix whose rows contain the x,y coordinates of the vertices of \code{d}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{gplot}}, \code{\link{cmdscale}}, \code{\link{eigen}} }
\examples{
}
\keyword{ graphs }% at least one, from doc/KEYWORDS
\keyword{ dplot }% __ONLY ONE__ keyword per line
