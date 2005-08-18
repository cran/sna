\name{prestige}
\alias{prestige}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Calculate the Vertex Prestige Scores }
\description{
  \code{prestige} takes one or more graphs (\code{dat}) and returns the prestige scores of positions (selected by \code{nodes}) within the graphs indicated by \code{g}.  Depending on the specified mode, prestige based on any one of a number of different definitions will be returned. This function is compatible with \code{\link{centralization}}, and will return the theoretical maximum absolute deviation (from maximum) conditional on size (which is used by \code{\link{centralization}} to normalize the observed centralization score).
}
\usage{
prestige(dat, g=1, nodes=NULL, gmode="digraph", diag=FALSE, 
    cmode="indegree", tmaxdev=FALSE, rescale=FALSE, tol=1e-07)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{g}{ integer indicating the index of the graph for which centralities are to be calculated (or a vector thereof).  By default, \code{g==1}. }
  \item{nodes}{ vector indicating which nodes are to be included in the calculation.  By default, all nodes are included. }
  \item{gmode}{ string indicating the type of graph being evaluated.  "digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  \code{gmode} is set to "digraph" by default.}
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{cmode}{ one of "indegree", "indegree.rownorm", "indegree.rowcolnorm", "eigenvector", "eigenvector.rownorm", "eigenvector.colnorm", "eigenvector.rowcolnorm", "domain", or "domain.proximity". }
  \item{tmaxdev}{ boolean indicating whether or not the theoretical maximum absolute deviation from the maximum nodal centrality should be returned.  By default, \code{tmaxdev==FALSE}. }
  \item{rescale}{ if true, centrality scores are rescaled such that they sum to 1. }
  \item{tol}{ Currently ignored }
}
\details{
  "Prestige" is the name collectively given to a range of centrality scores which focus on the extent to which one is nominated by others.  The definitions supported here are as follows:
\enumerate{
\item indegree: indegree centrality

\item indegree.rownorm: indegree within the row-normalized graph

\item indegree.rowcolnorm: indegree within the row-column normalized graph

\item eigenvector: eigenvector centrality within the transposed graph (i.e., incoming ties recursively determine prestige)

\item eigenvector.rownorm: eigenvector centrality within the transposed row-normalized graph

\item eigenvector.colnorm: eigenvector centrality within the transposed column-normalized graph

\item eigenvector.rowcolnorm: eigenvector centrality within the transposed row/column-normalized graph

\item domain: indegree within the reachability graph (Lin's unweighted measure)

\item domain.proximity: Lin's proximity-weighted domain prestige
}

Note that the centralization of prestige is simply the extent to which one actor has substantially greater prestige than others; the underlying definition is the same.
}
\value{
  A vector, matrix, or list containing the prestige scores (depending on the number and size of the input graphs).
}
\references{ Lin, N.  (1976).  \emph{Foundations of Social Research}.  New York: McGraw Hill.

Wasserman, S., and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

 \section{Warning }{Making adjacency matrices doubly stochastic (row-column normalization) is not guaranteed to work.  In general, be wary of attempting to try normalizations on graphs with degenerate rows and columns.} 

\seealso{ \code{\link{centralization}}}

\examples{
g<-rgraph(10)                 #Draw a random graph with 10 members
prestige(g,cmode="domain")    #Compute domain prestige scores
}
\keyword{univar}
\keyword{ math }
\keyword{ graphs }

