\name{degree}
\alias{degree}
\alias{degree_R}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the Degree Centrality Scores of Network Positions }
\description{
  \code{Degree} takes one or more graphs (\code{dat}) and returns the degree centralities of positions (selected by \code{nodes}) within the graphs indicated by \code{g}.  Depending on the specified mode, indegree, outdegree, or total (Freeman) degree will be returned; this function is compatible with \code{\link{centralization}}, and will return the theoretical maximum absolute deviation (from maximum) conditional on size (which is used by \code{\link{centralization}} to normalize the observed centralization score).
}
\usage{
degree(dat, g=1, nodes=NULL, gmode="digraph", diag=FALSE,
    tmaxdev=FALSE, cmode="freeman", rescale=FALSE, ignore.eval=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{g}{ integer indicating the index of the graph for which centralities are to be calculated (or a vector thereof).  By default, \code{g}=1. }
  \item{nodes}{ vector indicating which nodes are to be included in the calculation.  By default, all nodes are included. }
  \item{gmode}{ string indicating the type of graph being evaluated.  \code{"digraph"} indicates that edges should be interpreted as directed; \code{"graph"} indicates that edges are undirected.  \code{gmode} is set to \code{"digraph"} by default. }
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{tmaxdev}{ boolean indicating whether or not the theoretical maximum absolute deviation from the maximum nodal centrality should be returned.  By default, \code{tmaxdev==FALSE}. }
  \item{cmode}{ string indicating the type of degree centrality being computed.  \code{"indegree"}, \code{"outdegree"}, and \code{"freeman"} refer to the indegree, outdegree, and total (Freeman) degree measures, respectively.  The default for \code{cmode} is \code{"freeman"}. }
  \item{rescale}{ if true, centrality scores are rescaled such that they sum to 1. }
  \item{ignore.eval}{ logical; should edge values be ignored when computing degree scores?}
}
\details{
  Degree centrality is the social networker's term for various permutations of the graph theoretic notion of vertex degree: for unvalued graphs, indegree of a vertex, \eqn{v}, corresponds to the cardinality of the vertex set \eqn{N^+(v)=\{i \in V(G) : (i,v) \in E(G)\}}{N^+(v) = \{i in V(G) : (i,v) in E(G)\}}; outdegree corresponds to the cardinality of the vertex set \eqn{N^-(v)=\{i \in V(G) : (v,i) \in E(G)\}}{N^-(v) = \{i in V(G) : (v,i) in E(G)\}}; and total (or \dQuote{Freeman}) degree corresponds to \eqn{\left|N^+(v)\right| + \left|N^-(v)\right|}{|N^+(v)|+|N^-(v)|}.  (Note that, for simple graphs, indegree=outdegree=total degree/2.)  Obviously, degree centrality can be interpreted in terms of the sizes of actors' neighborhoods within the larger structure.  See the references below for more details.

When \code{ignore.eval==FALSE}, \code{degree} weights edges by their values where supplied.  \code{ignore.eval==TRUE} ensures an unweighted degree score (independent of input).  Setting \code{gmode=="graph"} forces behavior equivalent to \code{cmode=="indegree"} (i.e., each edge is counted only once); to obtain a total degree score for an undirected graph in which both in- and out-neighborhoods are counted separately, simply use \code{gmode=="digraph"}.
}
\value{
  A vector, matrix, or list containing the degree scores (depending on the number and size of the input graphs).
}
\references{ Freeman, L.C.  (1979).  ``Centrality in Social Networks I: Conceptual Clarification.'' \emph{Social Networks}, 1, 215-239. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }


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
\keyword{ math }
