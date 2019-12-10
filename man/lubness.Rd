\name{lubness}
\alias{lubness}
\alias{lubness_con_R}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute Graph LUBness Scores }
\description{
  \code{lubness} takes a graph set (\code{dat}) and returns the Krackhardt LUBness scores for the graphs selected by \code{g}.
}
\usage{
lubness(dat, g=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{g}{ index values for the graphs to be utilized; by default, all graphs are selected. }
}
\details{
In the context of a directed graph \eqn{G}, two actors \eqn{i} and \eqn{j} may be said to have an \emph{upper bound} iff there exists some actor \eqn{k} such that directed \eqn{ki} and \eqn{kj} paths belong to \eqn{G}.  An upper bound \eqn{\ell}{l} is known as a \emph{least upper bound} for \eqn{i} and \eqn{j} iff it belongs to at least one \eqn{ki} and \eqn{kj} path (respectively) for all \eqn{i,j} upper bounds \eqn{k}; let \eqn{L(i,j)} be an indicator which returns 1 iff such an \eqn{\ell}{l} exists, otherwise returning 0.  Now, let \eqn{G_1,G_2,\dots,G_n}{G_1,G_2,...,G_n} represent the weak components of \eqn{G}.  For convenience, we denote the cardinalities of these graphs' vertex sets by \eqn{|V(G)|=N} and \eqn{|V(G_i)|=N_i}, \eqn{\forall i \in 1,\dots,n}{for i in 1,...,n}.  Given this, the Krackhardt LUBness of \eqn{G} is given by

\deqn{
1-\frac{\sum_{i=1}^n \sum_{v_j,v_k \in V(G_i)} \Bigl(1-L(v_j,v_k)\Bigr)}{\sum_{i=1}^n \frac{1}{2}(N_i-1)(N_i-2)}}{%
1-Sum(Sum(1-L(v_j,v_k),v_j,v_k in V(G_i)),i=1,...,n)/Sum((N_i-1)(N_i-2)/2,i=1,...,n)}


Where all vertex pairs possess a least upper bound, Krackhardt's LUBness is equal to 1; in general, it approaches 0 as this condition is broached.  (This convergence is problematic in certain cases due to the requirement that we sum violations across components; where a graph contains no components of size three or greater, Krackhardt's LUBness is not well-defined.  \code{lubness} returns a \code{NaN} in these cases.)  

LUBness is one of four measures (\code{\link{connectedness}}, \code{\link{efficiency}}, \code{\link{hierarchy}}, and \code{\link{lubness}}) suggested by Krackhardt for summarizing hierarchical structures.  Each corresponds to one of four axioms which are necessary and sufficient for the structure in question to be an outtree; thus, the measures will be equal to 1 for a given graph iff that graph is an outtree.  Deviations from unity can be interpreted in terms of failure to satisfy one or more of the outtree conditions, information which may be useful in classifying its structural properties.
}
\value{
A vector of LUBness scores
}
\references{ Krackhardt, David.  (1994).  ``Graph Theoretical Dimensions of Informal Organizations.'' In K. M. Carley and M. J. Prietula (Eds.), \emph{Computational Organization Theory}, 89-111. Hillsdale, NJ: Lawrence Erlbaum and Associates. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ The four Krackhardt indices are, in general, nondegenerate for a relatively narrow band of size/density combinations (efficiency being the sole exception).  This is primarily due to their dependence on the reachability graph, which tends to become complete rapidly as size/density increase.  See Krackhardt (1994) for a useful simulation study. }

\seealso{ \code{\link{connectedness}}, \code{\link{efficiency}}, \code{\link{hierarchy}}, \code{\link{lubness}}, \code{\link{reachability}} }

\examples{
#Get LUBness scores for graphs of varying densities
lubness(rgraph(10,5,tprob=c(0.1,0.25,0.5,0.75,0.9)))

}
\keyword{ math }% at least one, from doc/KEYWORDS
\keyword{ univar }% __ONLY ONE__ keyword per line
\keyword{ graphs }
