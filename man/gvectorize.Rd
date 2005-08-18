\name{gvectorize}
\alias{gvectorize}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Vectorization of Adjacency Matrices}
\description{
\code{gvectorize} takes an input graph set and converts it into a corresponding number of vectors by row concatenation.
}
\usage{
gvectorize(mats, mode="digraph", diag=FALSE, censor.as.na=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{mats}{ one or more input graphs. }
  \item{mode}{ ``digraph'' if data is taken to be directed, else ``graph''. }
  \item{diag}{ boolean indicating whether diagonal entries (loops) are taken to contain meaningful data. }
  \item{censor.as.na}{ if \code{TRUE}, code unused parts of the adjacency matrix as \code{NA}s prior to vectorizing; otherwise, unused parts are simply removed. }
}
\details{
The output of \code{gvectorize} is a matrix in which each column corresponds to an input graph, and each row corresponds to an edge.  The columns of the output matrix are formed by simple row-concatenation of the original adjacency matrices, possibly after removing cells which are not meaningful (if \code{censor.as.na==FALSE}).  This is useful when preprocessing edge sets for use with \code{glm} or the like.
}
\value{
An nxk matrix, where n is the number of arcs and k is the number of graphs; if \code{censor.as.na==FALSE}, n will be reflect the relevant number of uncensored arcs.
}
%\references{ }
\author{ Carter T. Butts \email{buttsc@uci.edu}}
%\note{  }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

%\seealso{  }

\examples{
#Draw two random graphs
g<-rgraph(10,2)

#Examine the vectorized form of the adjacency structure
gvectorize(g)
}
\keyword{ manip }%-- one or more ...
\keyword{ math }
\keyword{ graphs}
