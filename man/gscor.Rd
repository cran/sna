\name{gscor}
\alias{gscor}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Find the Structural Correlations Between Two or More Graphs }
\description{
  \code{gscor} finds the product-moment structural correlation between the adjacency matrices of graphs indicated by \code{g1} and \code{g2} in stack \code{dat} (or possibly \code{dat2}) given exchangeability list \code{exchange.list}.  Missing values are permitted. 
}
\usage{
gscor(dat, dat2=NULL, g1=NULL, g2=NULL,  diag=FALSE, 
    mode="digraph", method="anneal", reps=1000, prob.init=0.9,
    prob.decay=0.85, freeze.time=25, full.neighborhood=TRUE, 
    exchange.list=0)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ a stack of input graphs. }
  \item{dat2}{ optionally, a second graph stack. }
  \item{g1}{ the indices of \code{dat} reflecting the first set of graphs to be compared; by default, all members of \code{dat} are included. }
  \item{g2}{ the indices or \code{dat} (or \code{dat2}, if applicable) reflecting the second set of graphs to be compared; by default, all members of \code{dat} are included. }
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{mode}{ string indicating the type of graph being evaluated.  \code{"digraph"} indicates that edges should be interpreted as directed; \code{"graph"} indicates that edges are undirected.  \code{mode} is set to \code{"digraph"} by default. }
  \item{method}{ method to be used to search the space of accessible permutations; must be one of \code{"none"}, \code{"exhaustive"}, \code{"anneal"}, \code{"hillclimb"}, or \code{"mc"}.  }
  \item{reps}{ number of iterations for Monte Carlo method. }
  \item{prob.init}{ initial acceptance probability for the annealing routine. }
  \item{prob.decay}{ cooling multiplier for the annealing routine. }
  \item{freeze.time}{ freeze time for the annealing routine. }
  \item{full.neighborhood}{ should the annealer evaluate the full neighborhood of pair exchanges at each iteration? }
  \item{exchange.list}{ information on which vertices are exchangeable (see below); this must be a single number, a vector of length n, or a nx2 matrix. }
}
\details{

The structural correlation coefficient between two graphs G and H is defined as
\deqn{scor\left(G,H \left| L_G,L_H\right.\right) = \max_{L_G,L_H} cor(\ell(G),\ell(H))}{%
scor(G,H | L_G,L_H) = max_[L_G,L_H] cor(l(G),l(H))}
where \eqn{L_G} is the set of accessible permutations/labelings of G, \eqn{\ell(G)}{l(G)} is a permutation/relabeling of G, and \eqn{\ell(G) \in L_G}{l(G) in L_G}.  The set of accessible permutations on a given graph is determined by the \emph{theoretical exchangeability} of its vertices; in a  nutshell, two vertices are considered to be theoretically exchangeable for a given problem if all predictions under the conditioning theory are invariant to a relabeling of the vertices in question (see Butts and Carley (2001) for a more formal exposition).  Where no vertices are exchangeable, the structural correlation becomes the simple graph correlation.  Where \emph{all} vertices are exchangeable, the structural correlation reflects the correlation between unlabeled graphs; other cases correspond to correlation under partial labeling.  

The accessible permutation set is determined by the \code{exchange.list} argument, which is dealt with in the following manner. First, \code{exchange.list} is expanded to fill an nx2 matrix.  If \code{exchange.list} is a single number, this is trivially accomplished by replication; if \code{exchange.list} is a vector of length n, the matrix is formed by cbinding two copies together.  If \code{exchange.list} is already an nx2 matrix, it is left as-is.  Once the nx2 exchangeability matrix has been formed, it is interpreted as follows: columns refer to graphs 1 and 2, respectively; rows refer to their corresponding vertices in the original adjacency matrices; and vertices are taken to be theoretically exchangeable iff their corresponding exchangeability matrix values are identical.  To obtain an unlabeled graph correlation (the default), then, one could simply let \code{exchange.list} equal any single number.  To obtain the standard graph correlation, one would use the vector \code{1:n}.

Because the set of accessible permutations is, in general, very large (\eqn{o(n!)}), searching the set for the maximum correlation is a non-trivial affair.  Currently supported methods for estimating the structural correlation are hill climbing, simulated annealing, blind monte carlo search, or exhaustive search (it is also possible to turn off searching entirely).  Exhaustive search is not recommended for graphs larger than size 8 or so, and even this may take days; still, this is a valid alternative for small graphs.  Blind monte carlo search and hill climbing tend to be suboptimal for this problem and are not, in general recommended, but they are available if desired.  The preferred (and default) option for permutation search is simulated annealing, which seems to work well on this problem (though some tinkering with the annealing parameters may be needed in order to get optimal performance).  See the help for \code{\link{lab.optimize}} for more information regarding these options.

Structural correlation matrices are p.s.d., and are p.d. so long as no graph within the set is a linear combination of any other under any accessible permutation.  Their eigendecompositions are meaningful and they may be used in linear subspace analyses, so long as the researcher is careful to interpret the results in terms of the appropriate set of accessible labelings.  Classical null hypothesis tests should not be employed with structural correlations, and QAP tests are almost never appropriate (save in the uniquely labeled case).  See \code{\link{cugtest}} for a more reasonable alternative.
}
\value{
An estimate of the structural correlation matrix
}
\references{ Butts, C.T., and Carley, K.M.  (2001).  ``Multivariate Methods for Interstructural Analysis.''  CASOS Working Paper, Carnegie Mellon University.}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ Consult Butts and Carley (2001) for advice and examples on theoretical exchangeability. }

\section{Warning }{The search process can be \emph{very slow}, particularly for large graphs.  In particular, the \emph{exhaustive} method is order factorial, and will take approximately forever for unlabeled graphs of size greater than about 7-9.}

\seealso{ \code{\link{gscov}}, \code{\link{gcor}}, \code{\link{gcov}} }

\examples{
#Generate two random graphs
g.1<-rgraph(5)
g.2<-rgraph(5)

#Copy one of the graphs and permute it
perm<-sample(1:5)
g.3<-g.2[perm,perm]

#What are the structural correlations between the labeled graphs?
gscor(g.1,g.2,exchange.list=1:5)
gscor(g.1,g.3,exchange.list=1:5)
gscor(g.2,g.3,exchange.list=1:5)

#What are the structural correlations between the underlying 
#unlabeled graphs?
gscor(g.1,g.2)
gscor(g.1,g.3)
gscor(g.2,g.3)
}
\keyword{ univar }
\keyword{ multivariate }%-- one or more ...
\keyword{ graphs }
