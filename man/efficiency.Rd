\name{efficiency}
\alias{efficiency}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute Graph Efficiency Scores }
\description{
  \code{efficiency} takes one or more graphs (\code{dat}) and returns the Krackhardt efficiency scores for the graphs selected by \code{g}.
}
\usage{
efficiency(dat, g=NULL, diag=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more graphs. }
  \item{g}{ index values for the graphs to be utilized; by default, all graphs are selected. }
  \item{diag}{ \code{TRUE} if the diagonal contains valid data; by default, \code{diag==FALSE}. }
}
\details{
Let \eqn{G=\cup_{i=1}^n G_i}{G= G_1 U ... U G_n} be a digraph with weak components \eqn{G_1,G_2,\dots,G_n}{G_1,G_2,...,G_n}.  For convenience, we denote the cardinalities of these components' vertex sets by \eqn{|V(G)|=N} and \eqn{|V(G_i)|=N_i}, \eqn{\forall i \in 1,\dots,n}{for i in 1,...,n}.  Then the Krackhardt efficiency of \eqn{G} is given by

\deqn{
1-\frac{|E(G)| - \sum_{i=1}^n \left(N_i-1\right)}{\sum_{i=1}^n \left(N_i \left(N_i-1\right)-\left(N_i-1\right)\right)}}{%
1 - ( |E(G)| - Sum(N_i-1,i=1,..,n) )/( Sum(N_i(N_i-1) - (N_i-1),i=1,..,n) )}

which can be interpreted as 1 minus the proportion of possible ``extra'' edges (above those needed to weakly connect the existing components) actually present in the graph.  A graph which an efficiency of 1 has precisely as many edges as are needed to connect its components; as additional edges are added, efficiency gradually falls towards 0.

Efficiency is one of four measures (\code{\link{connectedness}}, \code{\link{efficiency}}, \code{\link{hierarchy}}, and \code{\link{lubness}}) suggested by Krackhardt for summarizing hierarchical structures.  Each corresponds to one of four axioms which are necessary and sufficient for the structure in question to be an outtree; thus, the measures will be equal to 1 for a given graph iff that graph is an outtree.  Deviations from unity can be interpreted in terms of failure to satisfy one or more of the outtree conditions, information which may be useful in classifying its structural properties.
}
\value{
A vector of efficiency scores
}
\references{ Krackhardt, David.  (1994).  ``Graph Theoretical Dimensions of Informal Organizations.'' In K. M. Carley and M. J. Prietula (Eds.), \emph{Computational Organization Theory}, 89-111. Hillsdale, NJ: Lawrence Erlbaum and Associates. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ The four Krackhardt indices are, in general, nondegenerate for a relatively narrow band of size/density combinations (efficiency being the sole exception).  This is primarily due to their dependence on the reachability graph, which tends to become complete rapidly as size/density increase.  See Krackhardt (1994) for a useful simulation study. 

The violation normalization used before version 0.51 was \eqn{N\left(N-1\right) \sum_{i=1}^n \left(N_i-1\right)}{N(N-1) - Sum(N_i-1,i=1,..,n)}, based on a somewhat different interpretation of the definition in Krackhardt (1994).  The former version gave results which more closely matched those of the cited simulation study, but was less consistent with the textual definition.}


\seealso{ \code{\link{connectedness}}, \code{\link{efficiency}}, \code{\link{hierarchy}}, \code{\link{lubness}}, \code{\link{gden}} }

\examples{
#Get efficiency scores for graphs of varying densities
efficiency(rgraph(10,5,tprob=c(0.1,0.25,0.5,0.75,0.9)))
}
\keyword{ math }% at least one, from doc/KEYWORDS
\keyword{ univar }% __ONLY ONE__ keyword per line
\keyword{ graphs }
