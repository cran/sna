\name{degree}
\alias{degree}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the Degree Centrality Scores of Network Positions }
\description{
  \code{Degree} takes a graph stack (\code{dat}) and returns the degree centralities of positions within one graph (indicated by \code{nodes} and \code{g}, respectively).  Depending on the specified mode, indegree, outdegree, or total (Freeman) degree will be returned; this function is compatible with \code{\link{centralization}}, and will return the theoretical maximum absolute deviation (from maximum) conditional on size (which is used by \code{\link{centralization}} to normalize the observed centralization score).
}
\usage{
degree(dat, g=1, nodes=c(1:dim(dat)[2]), gmode="digraph", diag=FALSE,
    tmaxdev=FALSE, cmode="freeman", rescale=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ Data array to be analyzed.  By assumption, the first dimension of the array indexes the graph, with the next two indexing the actors.  Provided that \code{FUN} is well-behaved, this can be an n x n matrix if only one graph is involved. }
  \item{g}{ Integer indicating the index of the graph for which centralities are to be calculated.  By default, \code{g==1}. }
  \item{nodes}{ List indicating which nodes are to be included in the calculation.  By default, all nodes are included. }
  \item{gmode}{ String indicating the type of graph being evaluated.  "digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  \code{gmode} is set to "digraph" by default. }
  \item{diag}{ Boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{tmaxdev}{ Boolean indicating whether or not the theoretical maximum absolute deviation from the maximum nodal centrality should be returned.  By default, \code{tmaxdev==FALSE}. }
  \item{cmode}{ String indicating the type of degree centrality being computed.  "indegree", "outdegree", and "freeman" refer to the indegree, outdegree, and total (Freeman) degree measures, respectively.  The default for \code{cmode} is "freeman". }
  \item{rescale}{ If true, centrality scores are rescaled such that they sum to 1. }
}
\details{
  Degree centrality is the social networker's term for various permutations of the graph theoretic notion of vertex degree: indegree of a vertex, v, corresponds to the cardinality of the vertex set \eqn{N^+(v)=\{i \in V(G) : (i,v) \in E(G)\}}{N^+(v) = \{i in V(G) : (i,v) in E(G)\}}; outdegree corresponds to the cardinality of the vertex set \eqn{N^-(v)=\{i \in V(G) : (v,i) \in E(G)\}}{N^-(v) = \{i in V(G) : (v,i) in E(G)\}}; and total (or "Freeman") degree corresponds to \eqn{\left|N^+(v)\right| + \left|N^-(v)\right|}{|N^+(v)|+|N^-(v)|}.  (Note that, for simple graphs, indegree=outdegree=total degree/2.)  Obviously, degree centrality can be interpreted in terms of the sizes of actors' neighborhoods within the larger structure.  See the references below for more details.
}
\value{
  A vector containing the degree centrality scores
}
\references{ Freeman, L.C.  (1979).  ``Centrality in Social Networks I: Conceptual Clarification.'' \emph{Social Networks}, 1, 215-239. }
\author{ Carter T. Butts \email{ctb@andrew.cmu.edu} }


\seealso{ \code{\link{centralization}} }

\examples{
#Create a random directed graph
dat<-rgraph(10)
#Find the indegrees, outdegrees, and total degrees
degree(dat,cmode="indegree")
degree(dat,cmode="outdegree")
degree(dat)
}
\keyword{univar}
\keyword{ network}
