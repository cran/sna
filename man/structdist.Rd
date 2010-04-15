\name{structdist}
\alias{structdist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Find the Structural Distances Between Two or More Graphs }
\description{
  \code{structdist} returns the structural distance between the labeled graphs \code{g1} and \code{g2} in stack \code{dat} based on Hamming distance for dichotomous data, or else the absolute (manhattan) distance.  If \code{normalize} is true, this distance is divided by its dichotomous theoretical maximum (conditional on |V(G)|).
}
\usage{
structdist(dat, g1=NULL, g2=NULL, normalize=FALSE, diag=FALSE,
    mode="digraph", method="anneal",  reps=1000, prob.init=0.9,
    prob.decay=0.85, freeze.time=25, full.neighborhood=TRUE,
    mut=0.05, pop=20, trials=5, exchange.list=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{g1}{ a vector indicating which graphs to compare (by default, all elements of \code{dat}).  }
  \item{g2}{ a vector indicating against which the graphs of \code{g1} should be compared (by default, all graphs). }
  \item{normalize}{ divide by the number of available dyads? }
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{mode}{ string indicating the type of graph being evaluated.  \code{"digraph"} indicates that edges should be interpreted as directed; \code{"graph"} indicates that edges are undirected.  \code{mode} is set to \code{"digraph"} by default.}
  \item{method}{ method to be used to search the space of accessible permutations; must be one of \code{"none"}, \code{"exhaustive"}, \code{"anneal"}, \code{"hillclimb"}, or \code{"mc"}. }
  \item{reps}{ number of iterations for Monte Carlo method.}
  \item{prob.init}{ initial acceptance probability for the annealing routine. }
  \item{prob.decay}{ cooling multiplier for the annealing routine. }
  \item{freeze.time}{ freeze time for the annealing routine. }
  \item{full.neighborhood}{ should the annealer evaluate the full neighborhood of pair exchanges at each iteration? }
  \item{mut}{ GA Mutation rate (currently ignored). }
  \item{pop}{ GA population (currently ignored). }
  \item{trials}{ number of GA populations (currently ignored). }
  \item{exchange.list}{ information on which vertices are exchangeable (see below); this must be a single number, a vector of length n, or a nx2 matrix. }
}
\details{
The structural distance between two graphs G and H is defined as
\deqn{d_S\left(G,H \left| L_G,L_H\right.\right) = \min_{L_G,L_H} d\left(\ell\left(G\right),\ell\left(H\right)\right)}{%
d_S(G,H | L_G,L_H) = min_[L_G,L_H] d(l(G),l(H))}
where \eqn{L_G} is the set of accessible permutations/labelings of G, and \eqn{\ell(G)}{l(G)} is a permuation/relabeling of the vertices of G (\eqn{\ell(G) \in L_G}{l(G) in L_G}).  The set of accessible permutations on a given graph is determined by the \emph{theoretical exchangeability} of its vertices; in a  nutshell, two vertices are considered to be theoretically exchangeable for a given problem if all predictions under the conditioning theory are invariant to a relabeling of the vertices in question (see Butts and Carley (2001) for a more formal exposition).  Where no vertices are exchangeable, the structural distance becomes the its labeled counterpart (here, the Hamming distance).  Where \emph{all} vertices are exchangeable, the structural distance reflects the distance between unlabeled graphs; other cases correspond to distance under partial labeling.  

The accessible permutation set is determined by the \code{exchange.list} argument, which is dealt with in the following manner. First, \code{exchange.list} is expanded to fill an nx2 matrix.  If \code{exchange.list} is a single number, this is trivially accomplished by replication; if \code{exchange.list} is a vector of length n, the matrix is formed by \code{cbind}ing two copies together.  If \code{exchange.list} is already an nx2 matrix, it is left as-is.  Once the nx2 exchangeabiliy matrix has been formed, it is interpreted as follows: columns refer to graphs 1 and 2, respectively; rows refer to their corresponding vertices in the original adjacency matrices; and vertices are taken to be theoretically exchangeable iff their corresponding exchangeability matrix values are identical.  To obtain an unlabeled distance (the default), then, one could simply let \code{exchange.list} equal any single number.  To obtain the Hamming distance, one would use the vector \code{1:n}.

Because the set of accessible permutations is, in general, very large (\eqn{o(n!)}), searching the set for the minimum distance is a non-trivial affair.  Currently supported methods for estimating the structural distance are hill climbing, simulated annealing, blind monte carlo search, or exhaustive search (it is also possible to turn off searching entirely).  Exhaustive search is not recommended for graphs larger than size 8 or so, and even this may take days; still, this is a valid alternative for small graphs.  Blind monte carlo search and hill climbing tend to be suboptimal for this problem and are not, in general recommended, but they are available if desired.  The preferred (and default) option for permutation search is simulated annealing, which seems to work well on this problem (though some tinkering with the annealing parameters may be needed in order to get optimal performance).  See the help for \code{\link{lab.optimize}} for more information regarding these options.

Structural distance matrices may be used in the same manner as any other distance matrices (e.g., with multidimensional scaling, cluster analysis, etc.)  Classical null hypothesis tests should not be employed with structural distances, and QAP tests are almost never appropriate (save in the uniquely labeled case).  See \code{\link{cugtest}} for a more reasonable alternative.
}
\value{
  A structural distance matrix
}
\references{
Butts, C.T. and Carley, K.M.  (2005).  \dQuote{Some Simple Algorithms for Structural Comparison.}  \emph{Computational and Mathematical Organization Theory,} 11(4), 291-305.

Butts, C.T., and Carley, K.M.  (2001).  \dQuote{Multivariate Methods for Interstructural Analysis.}  CASOS Working Paper, Carnegie Mellon University. 
}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ Consult Butts and Carley (2001) for advice and examples on theoretical exchangeability. }

\section{Warning }{The search process can be \emph{very slow}, particularly for large graphs.  In particular, the \emph{exhaustive} method is order factorial, and will take approximately forever for unlabeled graphs of size greater than about 7-9.}

\seealso{ \code{\link{hdist}}, \code{\link{sdmat}} }

\examples{
#Generate two random graphs
g<-array(dim=c(3,5,5))
g[1,,]<-rgraph(5)
g[2,,]<-rgraph(5)

#Copy one of the graphs and permute it
g[3,,]<-rmperm(g[2,,])

#What are the structural distances between the labeled graphs?
structdist(g,exchange.list=1:5)

#What are the structural distances between the underlying unlabeled 
#graphs?
structdist(g,method="anneal", prob.init=0.9, prob.decay=0.85, 
    freeze.time=50, full.neighborhood=TRUE)
}
\keyword{ univar }
\keyword{ multivariate }%-- one or more ...
\keyword{ math }
\keyword{ graphs }
