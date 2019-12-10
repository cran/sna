\name{geodist}
\alias{geodist}
\alias{geodist_R}
\alias{geodist_adj_R}
\alias{geodist_val_R}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Fund the Numbers and Lengths of Geodesics Among Nodes in a Graph }
\description{
  \code{geodist} uses a BFS to find the number and lengths of geodesics between all nodes of \code{dat}.  Where geodesics do not exist, the value in \code{inf.replace} is substituted for the distance in question.
}
\usage{
geodist(dat, inf.replace=Inf, count.paths=TRUE, predecessors=FALSE,
    ignore.eval=TRUE, na.omit=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{inf.replace}{ the value to use for geodesic distances between disconnected nodes; by default, this is equal \code{Inf}. }
  \item{count.paths}{ logical; should a count of geodesics be included in the returned object?}
  \item{predecessors}{ logical; should a predecessor list be included in the returned object?}
  \item{ignore.eval}{ logical; should edge values be ignored when computing geodesics?}
  \item{na.omit}{ logical; should \code{NA}-valued edges be removed?}
}
\details{
  This routine is used by a variety of other functions; many of these will allow the user to provide manually precomputed \code{geodist} output so as to prevent expensive recomputation.  Note that the choice of infinite path length for disconnected vertex pairs is non-canonical (albeit common), and some may prefer to simply treat these as missing values.  \code{geodist} (without loss of generality) treats all paths as directed, a fact which should be kept in mind when interpreting \code{geodist} output.
  
  By default, \code{geodist} ignores edge values (except for \code{NA}ed edges, which are dropped when \code{na.omit==TRUE}).  Setting \code{ignore.eval=FALSE} will change this behavior, with edge values being interpreted as distances; where edge values reflect proximity or tie strength, transformation may be necessary.  Edge values should also be non-negative.  Because the valued-case algorithm is significantly slower than the unvalued-case algorithm, \code{ignore.eval} should be set to \code{TRUE} wherever possible.
}
\value{
A list containing:

  \item{counts}{If \code{count.paths==TRUE}, a matrix containing the number of geodesics between each pair of vertices}
  \item{gdist}{A matrix containing the geodesic distances between each pair of vertices}
  \item{predecessors}{If \code{predecessors}, a list whose ith element is a list of vectors, the jth of which contains the intervening vertices on all shortest paths from i to j}
}
\references{ 
Brandes, U.  (2000).  ``Faster Evaluation of Shortest-Path Based Centrality Indices.'' \emph{Konstanzer Schriften in Mathematik und Informatik}, 120.

West, D.B.  (1996).  \emph{Introduction to Graph Theory.}  Upper Saddle River, N.J.: Prentice Hall. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }

 %~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{component.dist}}, \code{\link{components}} }

\examples{
#Find geodesics on a random graph
gd<-geodist(rgraph(15))

#Examine the number of geodesics
gd$counts

#Examine the geodesic distances
gd$gdist
}
\keyword{ graphs }%-- one or more ...
\keyword{ math}
