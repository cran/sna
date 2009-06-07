\name{flowbet}
\Rdversion{1.1}
\alias{flowbet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Flow Betweenness Scores of Network Positions
}
\description{
\code{flowbet} takes one or more graphs (\code{dat}) and returns the flow betweenness scores of positions (selected by \code{nodes}) within the graphs indicated by \code{g}.  Depending on the specified mode, flow betweenness on directed or undirected geodesics will be returned; this function is compatible with \code{\link{centralization}}, and will return the theoretical maximum absolute deviation (from maximum) conditional on size (which is used by \code{\link{centralization}} to normalize the observed centralization score).
}
\usage{
flowbet(dat, g = 1, nodes = NULL, gmode = "digraph", diag = FALSE,
    tmaxdev = FALSE, cmode = "rawflow", rescale = FALSE, 
    ignore.eval = FALSE)
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
string indicating the type of graph being evaluated.  \code{digraph} indicates that edges should be interpreted as directed (with flows summed over directed dyads); \code{graph} indicates that edges are undirected (with only undirected pairs considered).  \code{gmode} is set to \code{digraph} by default. 
}
  \item{diag}{
boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default.
}
  \item{tmaxdev}{
boolean indicating whether or not the theoretical maximum absolute deviation from the maximum nodal centrality should be returned.  By default, \code{tmaxdev}==\code{FALSE}.
}
  \item{cmode}{
one of \code{rawflow}, \code{normflow}, or \code{fracflow} (see below).
}
  \item{rescale}{
if true, centrality scores are rescaled such that they sum to 1.
}
  \item{ignore.eval}{
logical; ignore edge values when computing maximum flow (alternately, edge values will be assumed to carry capacity information)?
}
}
\details{
The (\dQuote{raw,} or unnormalized) flow betweenness of a vertex, \eqn{v \in V(G)}{v in V(G)}, is defined by Freeman et al. (1991) as
\deqn{
C_F(v) = \sum_{i,j : i \neq j, i \neq v, j \neq v} \left(f(i,j,G) - f(i,j,G\setminus v)\right),}{%
C_F(v) = sum( f(i,j,G) - f(i,j,G\v), i,j: i!=j,i!=v,j!=v ),}
where \eqn{f(i,j,G)} is the maximum flow from \eqn{i} to \eqn{j} within \eqn{G} (under the assumption of infinite vertex capacities, finite edge capacities, and non-simultaneity of pairwise flows).  Intuitively, unnormalized flow betweenness is simply the total maximum flow (aggregated across all pairs of third parties) mediated by \eqn{v}.

The above flow betweenness measure is computed by \code{flowbet} when \code{cmode=="rawflow"}.  In some cases, it may be desirable to normalize the raw flow betweenness by the total maximum flow among third parties (including \eqn{v}); this leads to the following normalized flow betweenness measure:
\deqn{
C'_F(v) = \frac{\sum_{i,j : i \neq j, i \neq v, j \neq v} \left(f(i,j,G) - f(i,j,G\setminus v)\right)}{\sum_{i,j : i \neq j, i \neq v, j \neq v} f(i,j,G)}.}{%
C'_F(v) = sum( f(i,j,G) - f(i,j,G\v), i,j: i!=j,i!=v,j!=v ) / sum( f(i,j,G), i,j: i!=j,i!=v,j!=v ).}
This variant can be selected by setting \code{cmode=="normflow"}.

Finally, it may be noted that the above normalization (from Freeman et al. (1991)) is rather different from that used in the definition of shortest-path betweenness, which normalizes within (rather than across) third-party dyads.  A third flow betweenness variant has been suggested by Koschutzki et al. (2005) based on a normalization of this type:
\deqn{
C''_F(v) = \sum_{i,j : i \neq j, i \neq v, j \neq v} \frac{ \left(f(i,j,G) - f(i,j,G\setminus v)\right)}{f(i,j,G)}}{%
C''_F(v) = sum( (f(i,j,G) - f(i,j,G\v)) / f(i,j,G), i,j: i!=j,i!=v,j!=v ),}
where 0/0 flow ratios are treated as 0 (as in shortest-path betweenness).  Setting \code{cmode=="fracflow"} selects this variant.
}
\value{
A vector of centrality scores.
}
\references{
Freeman, L.C.; Borgatti, S.P.; and White, D.R.  (1991).  \dQuote{Centrality in Valued Graphs: A Measure of Betweenness Based on Network Flow.}  \emph{Social Networks}, 13(2), 141-154.

Koschutzki, D.; Lehmann, K.A.; Peeters, L.; Richter, S.; Tenfelde-Podehl, D.; Zlotowski, O.  (2005).  \dQuote{Centrality Indices.}  In U. Brandes and T. Erlebach (eds.), \emph{Network Analysis: Methodological Foundations.}  Berlin: Springer.
}
\author{
Carter T. Butts \email{buttsc@uci.edu}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{betweenness}}, \code{\link{maxflow}}
}
\examples{
g<-rgraph(10)                                     #Draw a random graph
flowbet(g)                                        #Raw flow betweenness
flowbet(g,cmode="normflow")                       #Normalized flow betweenness

g<-g*matrix(rpois(100,4),10,10)                   #Add capacity constraints
flowbet(g)                                        #Note the difference!
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
\keyword{ univar }% __ONLY ONE__ keyword per line
