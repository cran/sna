\name{gden}
\alias{gden}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Find the Density of a Graph }
\description{
  \code{gden} computes the density of graph \code{g} in stack \code{dat}, adjusting for the type of graph in question.
}
\usage{
gden(dat, g=NULL, diag=FALSE, mode="digraph")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ Data array to be analyzed.  By assumption, the first dimension of the array indexes the graph, with the next two indexing the actors.  If \code{dat} is a matrix, \code{g} will be ignored. }
  \item{g}{ Integer indicating the index of the graphs for which the density is to be calculated.  If \code{g==NULL} (the default), density is calculated for all graphs in \code{dat}. }
  \item{diag}{ Boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{mode}{ String indicating the type of graph being evaluated.  "digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  \code{mode} is set to "digraph" by default. }
}
\details{
  The density of a graph is here taken to be the sum of tie values divided by the number of possible ties (i.e., an unbiased estimator of the graph mean); hence, the result is interpretable for valued graphs as the mean tie value.  The number of possible ties is determined by the graph type (and by \code{diag}) in the usual fashion.
}
\value{
  The graph density
}
\references{ Wasserman, S., and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

%\seealso{  }

\examples{
#Draw three random graphs
dat<-rgraph(10,3)
#Find their densities
gden(dat)
}
\keyword{ univar }%-- one or more ...
\keyword{ math }
\keyword{graphs}
