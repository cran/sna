\name{centralization}
\alias{centralization}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Find the Centralization of a Given Network, for Some Measure of Centrality  }
\description{
  \code{Centralization} returns the centralization GLI (graph-level index) for a given graph in \code{dat}, given a (node) centrality measure \code{FUN}.  \code{Centralization} follows Freeman's (1979) generalized definition of network centralization, and can be used with any properly defined centrality measure.  This measure must be implemented separately; see the references below for examples.
}
\usage{
centralization(dat, FUN, g=NULL, mode="digraph", diag=FALSE, 
    normalize=TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{FUN}{ Function to return nodal centrality scores.}
  \item{g}{ Integer indicating the index of the graph for which centralization should be computed.  By default, all graphs are employed. }
  \item{mode}{ String indicating the type of graph being evaluated.  "digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  \code{mode} is set to "digraph" by default. }
  \item{diag}{ Boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{normalize}{ Boolean indicating whether or not the centralization score should be normalized to the theoretical maximum.  (Note that this function relies on \code{FUN} to return this value when called with \code{tmaxdev==TRUE}.)  By default, \code{tmaxdev==TRUE}. }
  \item{\dots}{ Additional arguments to \code{FUN}. }
}
\details{
  The centralization of a graph G for centrality measure \eqn{C(v)}{C(v)} is defined (as per Freeman (1979)) to be:

\deqn{C^*(G) = \sum_{i \in V(G)} \left|\max_{v \in V(G)}(C(v))-C(i)\right|}{%
C^*(G) = sum( |max(C(v))-C(i)|, i in V(G) )}

Or, equivalently, the absolute deviation from the maximum of C on G.  Generally, this value is normalized by the theoretical maximum centralization score, conditional on \eqn{|V(G)|}{|V(G)|}.  (Here, this functionality is activated by \code{normalize}.)  \code{Centralization} depends on the function specified by \code{FUN} to return the vector of nodal centralities when called with \code{dat} and \code{g}, and to return the theoretical maximum value when called with the above and \code{tmaxdev==TRUE}.  For an example of such a centrality routine, see \code{\link{degree}}.
}
\value{
The centralization of the specified graph.
}
\references{Freeman, L.C.  (1979).  ``Centrality in Social Networks I: Conceptual Clarification.'' \emph{Social Networks}, 1, 215-239.

Wasserman, S., and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press.
 }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ See \code{\link{cugtest}} for null hypothesis tests involving centralization scores. }

\seealso{ \code{\link{cugtest}} }

\examples{
#Generate some random graphs
dat<-rgraph(5,10)
#How centralized is the third one on indegree?
centralization(dat,g=3,degree,cmode="indegree")
#How about on total (Freeman) degree?
centralization(dat,g=3,degree)
}
\keyword{ univar }%-- one or more ...
\keyword{ math }
\keyword{ graphs }
