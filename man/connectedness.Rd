\name{connectedness}
\alias{connectedness}
\alias{connectedness_R}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute Graph Connectedness Scores }
\description{
\code{connectedness} takes one or more graphs (\code{dat}) and returns the Krackhardt connectedness scores for the graphs selected by \code{g}. 
}
\usage{
connectedness(dat, g=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more graphs. }
  \item{g}{ index values for the graphs to be utilized; by default, all graphs are selected. }
}
\details{
Krackhardt's connectedness for a digraph \eqn{G} is equal to the fraction of all dyads, \eqn{\{i,j\}}, such that there exists an undirected path from \eqn{i} to \eqn{j} in \eqn{G}.  (This, in turn, is just the density of the weak \code{\link{reachability}} graph of \eqn{G}.)  Obviously, the connectedness score ranges from 0 (for the null graph) to 1 (for weakly connected graphs).

Connectedness is one of four measures (\code{\link{connectedness}}, \code{\link{efficiency}}, \code{\link{hierarchy}}, and \code{\link{lubness}}) suggested by Krackhardt for summarizing hierarchical structures.  Each corresponds to one of four axioms which are necessary and sufficient for the structure in question to be an outtree; thus, the measures will be equal to 1 for a given graph iff that graph is an outtree.  Deviations from unity can be interpreted in terms of failure to satisfy one or more of the outtree conditions, information which may be useful in classifying its structural properties.
}
\value{
A vector containing the connectedness scores
}
\references{ Krackhardt, David.  (1994).  ``Graph Theoretical Dimensions of Informal Organizations.'' In K. M. Carley and M. J. Prietula (Eds.), \emph{Computational Organization Theory}, 89-111. Hillsdale, NJ: Lawrence Erlbaum and Associates.}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{
The four Krackhardt indices are, in general, nondegenerate for a relatively narrow band of size/density combinations (efficiency being the sole exception).  This is primarily due to their dependence on the reachability graph, which tends to become complete rapidly as size/density increase.  See Krackhardt (1994) for a useful simulation study.
}
\seealso{ \code{\link{connectedness}}, \code{\link{efficiency}}, \code{\link{hierarchy}}, \code{\link{lubness}}, \code{\link{reachability}} }

\examples{
#Get connectedness scores for graphs of varying densities
connectedness(rgraph(10,5,tprob=c(0.1,0.25,0.5,0.75,0.9)))
}
\keyword{ math }% at least one, from doc/KEYWORDS
\keyword{ univar }% __ONLY ONE__ keyword per line
\keyword{ graphs }
