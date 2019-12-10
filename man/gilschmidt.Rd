\name{gilschmidt}
\alias{gilschmidt}
\alias{gilschmidt_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the Gil-Schmidt Power Index
}
\description{
\code{gilschmidt} computes the Gil-Schmidt Power Index for all nodes in \code{dat}, with or without normalization.
}
\usage{
gilschmidt(dat, g = 1, nodes = NULL, gmode = "digraph", diag = FALSE, 
    tmaxdev = FALSE, normalize = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
one or more input graphs (for best performance, sna edgelists or network objects are suggested).
  }
  \item{g}{ integer indicating the index of the graph for which centralities are to be calculated (or a vector thereof).  By default, \code{g}=1. }
  \item{nodes}{ list indicating which nodes are to be included in the calculation.  By default, all nodes are included. }
  \item{gmode}{ string indicating the type of graph being evaluated.  \code{"digraph"} indicates that edges should be interpreted as directed; \code{"graph"} indicates that edges are undirected.  \code{gmode} is set to \code{"digraph"} by default. }
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  (This has no effect on this index, but is included for compatibility with \code{\link{centralization}}. }
  \item{tmaxdev}{ boolean indicating whether or not the theoretical maximum absolute deviation from the maximum nodal centrality should be returned.  By default, \code{tmaxdev==FALSE}. }
  \item{normalize}{
logical; should the index scores be normalized?
  }
}
\details{
For graph \eqn{G=(V,E)}, let \eqn{R(v,G)} be the set of vertices reachable by \eqn{v} in \eqn{V\setminus v}{V \ v}.  Then the Gil-Schmidt power index is defined as
\deqn{C_{GS}(v) = \frac{\sum_{i \in R(v,G)} \frac{1}{d(v,i)}}{|R(v,G)|}.}{%
C_GS(v) = sum( 1/d(v,i), i in R(v,G) )/|R(v,G)|,}
where \eqn{d(v,i)} is the geodesic distance from \eqn{v} to \eqn{i} in \eqn{G}; the index is taken to be 0 for isolates.  The measure takes a value of 1 when \eqn{v} is adjacent to all reachable vertices, and approaches 0 as the distance from \eqn{v} to each vertex approaches infinity.  (For finite \eqn{N=|V|}, the minimum value is 0 if \eqn{v} is an isolate, and otherwise \eqn{1/(N-1)}.)

If \code{normalize=FALSE} is selected, then normalization by \eqn{|R(v,G)|} is not performed.  This measure has been proposed as a better-behaved alternative to closeness (to which it is closely related).

The \code{\link{closeness}} function in the sna library can also be used to compute this index.
}
\value{
A vector of centrality scores.
}
\references{
Gil, J. and Schmidt, S.  (1996).  \dQuote{The Origin of the Mexican Network of Power}.  Proceedings of the International Social Network Conference, Charleston, SC, 22-25.

Sinclair, P.A.  (2009).  \dQuote{Network Centralization with the Gil Schmidt Power Centrality Index}  \emph{Social Networks}, 29, 81-92.
}
\author{
Carter T. Butts, \email{buttsc@uci.edu}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{closeness}, \link{centralization}}
}
\examples{
data(coleman)  #Load Coleman friendship network
gs<-gilschmidt(coleman,g=1:2)  #Compute the Gil-Schmidt index

#Plot G-S values in the fall, versus spring
plot(gs,xlab="Fall",ylab="Spring",main="G-S Index")
abline(0,1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
\keyword{ graphs }% __ONLY ONE__ keyword per line
