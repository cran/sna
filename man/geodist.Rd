\name{geodist}
\alias{geodist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Fund the Numbers and Lengths of Geodesics Among Nodes in a Graph }
\description{
  \code{geodist} uses a BFS to find the number and lengths of geodesics between all nodes of \code{dat}.  Where geodesics do not exist, the value in \code{inf.replace} is substituted for the distance in question.
}
\usage{
geodist(dat, inf.replace=dim(dat)[2])
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ An adjacency matrix }
  \item{inf.replace}{ The value to use for geodesic distances between disconnected nodes; by default, this is equal to one greater than the maximum path length }
}
\details{
  This routine is used by a variety of other functions; many of these will allow the user to provide manually precomputed \code{geodist} output so as to prevent expensive recomputation.  Note that the choice of one greater than the maximum path length for disconnected vertex pairs is non-canonical, and some may prefer to simply treat these as missing values.  \code{geodist} (without loss of generality) treats all paths as directed, a fact which should be kept in mind when interpreting \code{geodist} output.
}
\value{
A data frame containing:

  \item{counts}{A matrix containing the number of geodesics between each pair of vertices}
  \item{sigma}{A matrix containing the geodesic distances between each pair of vertices}
}
\references{ 
Brandes, U.  (2000).  ``Faster Evaluation of Shortest-Path Based Centrality Indices.'' \emph{Konstanzer Schriften in Mathematik und Informatik}, 120.

West, D.B.  (1996).  \emph{Introduction to Graph Theory.}  Upper Saddle River, N.J.: Prentice Hall. }
\author{ Carter T. Butts \email{ctb@andrew.cmu.edu} }
%\note{ }

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
\keyword{ math }%-- one or more ...
