\name{infocent}
\alias{infocent}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Find Information Centrality Scores of Network Positions }
\description{
  \code{infocent} takes one or more graphs (\code{dat}) and returns the information centralities of positions (selected by \code{nodes}) within the graphs indicated by \code{g}.  This function is compatible with \code{\link{centralization}}, and will return the theoretical maximum absolute deviation (from maximum) conditional on size (which is used by \code{\link{centralization}} to normalize the observed centralization score).
}
\usage{
infocent(dat, g=1, nodes=NULL, gmode="digraph", diag=FALSE, 
    cmode="weak", tmaxdev=FALSE, rescale=FALSE,tol=1e-20)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{g}{ integer indicating the index of the graph for which centralities are to be calculated (or a vector thereof).  By default, \code{g==1}. }
  \item{nodes}{ list indicating which nodes are to be included in the calculation.  By default, all nodes are included. }
  \item{gmode}{ string indicating the type of graph being evaluated.  \code{"digraph"} indicates that edges should be interpreted as directed; \code{"graph"} indicates that edges are undirected.  This is currently ignored. }
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{cmode}{ the rule to be used by \code{\link{symmetrize}} when symmetrizing dichotomous data; must be one of \code{"weak"} (for an \code{OR} rule), \code{"strong"} for an \code{AND} rule), \code{"upper"} (for a \code{max} rule), or \code{"lower"} (for a \code{min} rule).  Set to \code{"weak"} by default, this parameter obviously has no effect on symmetric data. }
  \item{tmaxdev}{ boolean indicating whether or not the theoretical maximum absolute deviation from the maximum nodal centrality should be returned.  By default, \code{tmaxdev==FALSE}. }
  \item{rescale}{ if true, centrality scores are rescaled such that they sum to 1. }
  \item{tol}{ tolerance for near-singularities during matrix inversion (see \code{\link{solve}}). }
}
\details{
  Actor information centrality is a hybrid measure which relates to both path-length indices (e.g., closeness, graph centrality) and to walk-based eigenmeasures (e.g., eigenvector centrality, Bonacich power).  In particular, the information centrality of a given actor can be understood to be the harmonic average of the ``bandwidth'' for all paths originating with said individual (where the bandwidth is taken to be inversely related to path length).  Formally, the index is constructed as follows.  First, we take \eqn{G} to be an undirected (but possibly valued) graph -- symmetrizing if necessary -- with (possibly valued) adjacency matrix \eqn{\mathbf{A}}{A}.  From this, we remove all isolates (whose information centralities are zero in any event) and proceed to create the weighted connection matrix

\deqn{
\mathbf{C} = \mathbf{B}^{-1}}{%
C = B^-1}

where \eqn{\mathbf{B}}{B} is a pseudo-adjacency matrix formed by replacing the diagonal of \eqn{1-\mathbf{A}}{1-A} with one plus each actor's degree.  Given the above, let \eqn{T} be the trace of \eqn{\mathbf{C}}{C} with sum \eqn{S_T}, and let \eqn{S_R} be an arbitrary row sum (all rows of \eqn{\mathbf{C}}{C} have the same sum).  The information centrality scores are then equal to

\deqn{
C_I = \frac{1}{T + \frac{S_T-2S_R}{|V(G)|}}}{%
C_I = ( T + (S_T-2S_R)/|V(G)| )^-1
}

(recalling that the scores for any omitted vertices are 0).

In general, actors with higher information centrality are predicted to have greater control over the flow of information within a network; highly information-central individuals tend to have a large number of short paths to many others within the social structure.  Because the raw centrality values can be difficult to interpret directly, rescaled values are sometimes preferred (see the \code{rescale} option).  Though the use of path weights suggest information centrality as a possible replacement for closeness, the problem of inverting the \eqn{\mathbf{B}}{B} matrix poses problems of its own; as with all such measures, caution is advised on disconnected or degenerate structures. 
}
\value{
A vector, matrix, or list containing the centrality scores (depending on the number and size of the input graphs).
}
\references{ Stephenson, K., and Zelen, M. (1989).  ``Rethinking Centrality: Methods and Applications.''  \emph{Social Networks}, 11, 1-37.

Wasserman, S., and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press.}
\author{ David Barron \email{david.barron@jesus.ox.ac.uk}

Carter T. Butts \email{buttsc@uci.edu} }
\note{ The theoretical maximum deviation used here is not obtained with the star network; rather, the maximum occurs for an empty graph with one complete dyad, which is the model used here. }

\seealso{ \code{\link{evcent}}, \code{\link{bonpow}}, \code{\link{closeness}}, \code{\link{graphcent}}, \code{\link{centralization}} }

\examples{
#Generate some test data
dat<-rgraph(10,mode="graph")
#Compute information centrality scores
infocent(dat)
}
\keyword{ math }% at least one, from doc/KEYWORDS
\keyword{ univar }% __ONLY ONE__ keyword per line
\keyword{ graphs }
