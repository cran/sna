\name{loadcent}
\Rdversion{1.1}
\alias{loadcent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the Load Centrality Scores of Network Positions
}
\description{
 \code{loadcent} takes one or more graphs (\code{dat}) and returns the load centralities of positions (selected by \code{nodes}) within the graphs indicated by \code{g}.  Depending on the specified mode, load on directed or undirected geodesics will be returned; this function is compatible with \code{\link{centralization}}, and will return the theoretical maximum absolute deviation (from maximum) conditional on size (which is used by \code{\link{centralization}} to normalize the observed centralization score).
}
\usage{
loadcent(dat, g = 1, nodes = NULL, gmode = "digraph", diag = FALSE, 
    tmaxdev = FALSE, cmode = "directed", geodist.precomp = NULL, 
    rescale = FALSE, ignore.eval = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
 one or more input graphs. 
}
  \item{g}{
integer indicating the index of the graph for which centralities are to be calculated (or a vector thereof).  By default, \code{g}=1.
}
  \item{nodes}{
vector indicating which nodes are to be included in the calculation.  By default, all nodes are included.
}
  \item{gmode}{
string indicating the type of graph being evaluated.  \code{digraph} indicates that edges should be interpreted as directed; \code{graph} indicates that edges are undirected.  \code{gmode} is set to \code{digraph} by default.
}
  \item{diag}{
logical; should self-ties be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default.
}
  \item{tmaxdev}{
logical; return the theoretical maximum absolute deviation from the maximum nodal centrality (instead of the observed centrality scores)?  By default, \code{tmaxdev}==\code{FALSE}.
}
  \item{cmode}{
string indicating the type of load centrality being computed (directed or undirected). 
}
  \item{geodist.precomp}{
a \code{\link{geodist}} object precomputed for the graph to be analyzed (optional).
}
  \item{rescale}{
logical; if true, centrality scores are rescaled such that they sum to 1.
}
  \item{ignore.eval}{
logical; ignore edge values when computing shortest paths?
}
}
\details{
Goh et al.'s \emph{load centrality} (as reformulated by Brandes (2008)) is a betweenness-like measure defined through a hypothetical flow process.  Specifically, it is assumed that each vertex sends a unit of some commodity to each other vertex to which it is connected (without edge or vertex capacity constraints), with routing based on a priority system: given an input of flow \eqn{x} arriving at vertex \eqn{v} with destination \eqn{v'}, \eqn{v} divides \eqn{x} equally among all neigbors of minumum geodesic distance to the target.  The total flow passing through a given \eqn{v} via this process is defined as \eqn{v}'s \emph{load}.  Load is a potential alternative to betweenness for the analysis of flow structures operating well below their capacity constraints.
}
\value{
A vector of centrality scores.
}
\references{
Brandes, U.  (2008).  \dQuote{On Variants of Shortest-Path Betweenness Centrality and their Generic Computation.}  \emph{Social Networks}, 30, 136-145.


Goh, K.-I.; Kahng, B.; and Kim, D.  (2001).  \dQuote{Universal Behavior of Load Distribution in Scale-free Networks.}  \emph{Physical Review Letters}, 87(27), 1-4.
}
\author{
Carter T. Butts \email{buttsc@uci.edu}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{betweenness}}
}
\examples{
g<-rgraph(10)     #Draw a random graph with 10 members
loadcent(g)       #Compute load scores
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ graphs }% __ONLY ONE__ keyword per line
