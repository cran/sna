\name{betweenness}
\alias{betweenness}
\alias{betweenness_R}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the Betweenness Centrality Scores of Network Positions }
\description{
   \code{betweenness} takes one or more graphs (\code{dat}) and returns the betweenness centralities of positions (selected by \code{nodes}) within the graphs indicated by \code{g}.  Depending on the specified mode, betweenness on directed or undirected geodesics will be returned; this function is compatible with \code{\link{centralization}}, and will return the theoretical maximum absolute deviation (from maximum) conditional on size (which is used by \code{\link{centralization}} to normalize the observed centralization score).
}
\usage{
betweenness(dat, g=1, nodes=NULL, gmode="digraph", diag=FALSE,
    tmaxdev=FALSE, cmode="directed", geodist.precomp=NULL, 
    rescale=FALSE, ignore.eval=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{g}{ integer indicating the index of the graph for which centralities are to be calculated (or a vector thereof).  By default, \code{g}=1. }
  \item{nodes}{ vector indicating which nodes are to be included in the calculation.  By default, all nodes are included. }
  \item{gmode}{ string indicating the type of graph being evaluated.  "digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  \code{gmode} is set to "digraph" by default. }
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{tmaxdev}{ boolean indicating whether or not the theoretical maximum absolute deviation from the maximum nodal centrality should be returned.  By default, \code{tmaxdev}==\code{FALSE}. }
  \item{cmode}{ string indicating the type of betweenness centrality being computed (directed or undirected geodesics, or a variant form -- see below). }
  \item{geodist.precomp}{ A \code{\link{geodist}} object precomputed for the graph to be analyzed (optional) }
  \item{rescale}{ if true, centrality scores are rescaled such that they sum to 1. }
  \item{ignore.eval}{ logical; ignore edge values when computing shortest paths?}
}
\details{
The shortest-path betweenness of a vertex, \eqn{v}, is given by

\deqn{
C_B(v) = \sum_{i,j : i \neq j, i \neq v, j \neq v} \frac{g_{ivj}}{g_{ij}}}{%
C_B(v) = sum( g_ivj / g_ij, i,j: i!=j,i!=v,j!=v )}

where \eqn{g_{ijk}}{g_ijk} is the number of geodesics from \eqn{i} to \eqn{k} through \eqn{j}.  Conceptually, high-betweenness vertices lie on a large number of non-redundant shortest paths between other vertices; they can thus be thought of as ``bridges'' or ``boundary spanners.''         

Several variant forms of shortest-path betweenness exist, and can be selected using the \code{cmode} argument.  Supported options are as follows:
\describe{
\item{\code{directed}}{ Standard betweenness (see above), calculated on directed pairs.  (This is the default option.)}
\item{\code{undirected}}{ Standard betweenness (as above), calculated on undirected pairs (undirected graphs only).}
\item{\code{endpoints}}{ Standard betweenness, with direct connections counted towards ego's score.  This expresses the intuition that individuals' control over their own direct contacts should be considered in their total score (e.g., when betweenness is interpreted as a measure of information control).}
\item{\code{proximalsrc}}{ Borgatti's \emph{proximal source betweenness,} given by
\deqn{
C_B(v) = \sum_{i,j : i \neq v, i\neq j, j \to v} \frac{g_{ivj}}{g_{ij}}.}{%
C_B(v) = sum( g_ivj / g_ij, i,j: i!=v,i!=j,j->v ).}
This variant allows betweenness to accumulate only for the last intermediating vertex in each incoming geodesic; this expresses the notion that, by serving as the \dQuote{proximal source} for the target, this particular intermediary will in some settings have greater influence or control than other intervening parties. 
}
\item{\code{proximaltar}}{ Borgatti's \emph{proximal target betweenness,} given by
\deqn{
C_B(v) = \sum_{i,j : i \neq v, i\to v, i\neq j} \frac{g_{ivj}}{g_{ij}}.}{%
C_B(v) = sum( g_ivj / g_ij, i,j: i!=j,i->v,j!=v ).}
This counterpart to proximal source betweenness (above) allows betweenness to accumulate only for the first intermediating vertex in each outgoing geodesic; this expresses the notion that, by serving as the \dQuote{proximal target} for the source, this particular intermediary will in some settings have greater influence or control than other intervening parties. }
\item{\code{proximalsum}}{The sum of Borgatti's proximal source and proximal target betweenness scores (above); this may be used when either role is regarded as relevant to the betweenness calculation.}
\item{\code{lengthscaled}}{Borgetti and Everett's \emph{length-scaled betweenness}, given by
\deqn{
C_B(v) = \sum_{i,j : i \neq j, i \neq v, j \neq v} \frac{1}{d_{ij}}\frac{g_{ivj}}{g_{ij}},}{%
C_B(v) = sum( (1/d_ij)*(g_ivj / g_ij), i,j: i!=j,i!=v,j!=v ),}
where \eqn{d_{ij}}{d_ij} is the geodesic distance from \eqn{i} to \eqn{j}.  This measure adjusts the standard betweenness score by downweighting long paths (e.g., as appropriate in circumstances for which such paths are less-often used).
}
\item{\code{linearscaled}}{Geisberger et al.'s \emph{linearly-scaled betweenness}:
\deqn{
C_B(v) = \sum_{i,j : i \neq j, i \neq v, j \neq v} \frac{1}{d_{ij}}\frac{g_{ivj}}{g_{ij}}.}{%
C_B(v) = sum( (d_iv/d_ij)*(g_ivj / g_ij), i,j: i!=j,i!=v,j!=v ).}
This variant modifies the standard betweenness score by giving more weight to intermediaries which are closer to their targets (much like proximal source betweenness, above).  This may be of use when those near the end of a path have greater direct control over the flow of influence or resources than those near its source.
}
}
See Brandes (2008) for details and additional references.  Geodesics for all of the above can be calculated using valued edges by setting \code{ignore.eval=TRUE}.  Edge values are interpreted as distances for this purpose; proximity data should be transformed accordingly before invoking this routine.
}
\value{
A vector, matrix, or list containing the betweenness scores (depending on the number and size of the input graphs).
}
\references{ 

Borgatti, S.P. and Everett, M.G.  (2006).  \dQuote{A Graph-Theoretic Perspective on Centrality.}  \emph{Social Networks}, 28, 466-484.

Brandes, U.  (2008).  \dQuote{On Variants of Shortest-Path Betweenness Centrality and their Generic Computation.}  \emph{Social Networks}, 30, 136--145.

Freeman, L.C.  (1979).  \dQuote{Centrality in Social Networks I: Conceptual Clarification.} \emph{Social Networks}, 1, 215-239.

Geisberger, R., Sanders, P., and Schultes, D.  (2008).  \dQuote{Better Approximation of Betweenness Centrality.}  In \emph{Proceedings of the 10th Workshop on Algorithm Engineering and Experimentation (ALENEX'08)}, 90-100.  SIAM.

}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ Judicious use of \code{geodist.precomp} can save a great deal of time when computing multiple path-based indices on the same network. }

\section{Warning }{Rescale may cause unexpected results if all actors have zero betweenness.}

\seealso{ \code{\link{centralization}}, \code{\link{stresscent}}, \code{\link{geodist}} }

\examples{
g<-rgraph(10)     #Draw a random graph with 10 members
betweenness(g)    #Compute betweenness scores
}

\keyword{univar}
\keyword{ graphs }
