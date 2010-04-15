\name{lab.optimize}
\alias{lab.optimize}
\alias{lab.optimize.anneal}
\alias{lab.optimize.exhaustive}
\alias{lab.optimize.gumbel}
\alias{lab.optimize.hillclimb}
\alias{lab.optimize.mc}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Optimize a Bivariate Graph Statistic Across a Set of Accessible Permutations}
\description{
  \code{lab.optimize} is the front-end to a series of heuristic optimization routines (see below), all of which seek to maximize/minimize some bivariate graph statistic (e.g., graph correlation) across a set of vertex relabelings.
}
\usage{
lab.optimize(d1, d2, FUN, exchange.list=0, seek="min", 
    opt.method=c("anneal", "exhaustive", "mc", "hillclimb", 
    "gumbel"), ...)
lab.optimize.anneal(d1, d2, FUN, exchange.list=0, seek="min", 
    prob.init=1, prob.decay=0.99, freeze.time=1000, 
    full.neighborhood=TRUE, ...)
lab.optimize.exhaustive(d1, d2, FUN, exchange.list=0, seek="min", ...)
lab.optimize.gumbel(d1, d2, FUN, exchange.list=0, seek="min", 
    draws=500, tol=1e-5, estimator="median", ...)
lab.optimize.hillclimb(d1, d2, FUN, exchange.list=0, seek="min", ...)
lab.optimize.mc(d1, d2, FUN, exchange.list=0, seek="min", 
    draws=1000, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{d1}{ a single graph. }
  \item{d2}{ another single graph. }
  \item{FUN}{ a function taking two graphs as its first two arguments, and returning a numeric value. }
  \item{exchange.list}{ information on which vertices are exchangeable (see below); this must be a single number, a vector of length n, or a nx2 matrix. }
  \item{seek}{ "min" if the optimizer should seek a minimum, or "max" if a maximum should be sought. }
  \item{opt.method}{ the particular optimization method to use. }
  \item{prob.init}{initial acceptance probability for a downhill move (\code{lab.optimize.anneal} only).}
  \item{prob.decay}{the decay (cooling) multiplier for the probability of accepting a downhill move (\code{lab.optimize.anneal} only).}
  \item{freeze.time}{number of iterations at which the annealer should be frozen (\code{lab.optimize.anneal} only).}
  \item{full.neighborhood}{should all moves in the binary-exchange neighborhood be evaluated at each iteration? (\code{lab.optimize.anneal} only).}
  \item{tol}{tolerance for estimation of gumbel distribution parameters (\code{lab.optimize.gumbel} only).}
  \item{estimator}{Gumbel distribution statistic to use as optimal value prediction; must be one of ``mean'', ``median'', or ``mode'' (\code{lab.optimize.gumbel} only).}
  \item{draws}{number of draws to take for gumbel and mc methods. }
  \item{\dots}{ additional arguments to \code{FUN}. }
}
\details{
\code{lab.optimize} is the front-end to a family of routines for optimizing a bivariate graph statistic over a set of permissible relabelings (or equivalently, permutations).  The accessible permutation set is determined by the \code{exchange.list} argument, which is dealt with in the following manner. First, \code{exchange.list} is expanded to fill an nx2 matrix.  If \code{exchange.list} is a single number, this is trivially accomplished by replication; if \code{exchange.list} is a vector of length n, the matrix is formed by cbinding two copies together.  If \code{exchange.list} is already an nx2 matrix, it is left as-is.  Once the nx2 exchangeabiliy matrix has been formed, it is interpreted as follows: columns refer to graphs 1 and 2, respectively; rows refer to their corresponding vertices in the original adjacency matrices; and vertices are taken to be theoretically exchangeable iff their corresponding exchangeability matrix values are identical.  To obtain an unlabeled graph statistic (the default), then, one could simply let \code{exchange.list} equal any single number.  To obtain the labeled statistic, one would use the vector \code{1:n}.

Assuming a non-degenerate set of accessible permutations/relabelings, optimization proceeds via the algorithm specified in \code{opt.method}. The optimization routines which are currently implemented use a variety of different techniques, each with certain advantages and disadvantages.  A brief summary of each is as follows:

\enumerate{
\item exhaustive search (``exhaustive''): Under exhaustive search, the entire space of accessible permutations is combed for the global optimum.  This guarantees a correct answer, but at a very high price: the set of all permutations grows with the factorial of the number of vertices, and even substantial exchangeability constraints are unlikely to keep the number of permutations from growing out of control.  While exhaustive search \emph{is} possible for small graphs, unlabeled structures of size approximately 10 or greater cannot be treated using this algorithm within a reasonable time frame.  

Approximate complexity: on the order of \eqn{\prod_{i \in L}|V_i|!}{prod( |V_i|!, i in L )}, where L is the set of exchangeability classes.

\item hill climbing (``hillclimb''): The hill climbing algorithm employed here searches, at each iteration, the set of all permissible binary exchanges of vertices.  If one or more exchanges are found which are superior to the current permutation, the best alternative is taken.  If no superior alternative is found, then the algorithm terminates.  As one would expect, this algorithm is guaranteed to terminate on a local optimum; unfortunately, however, it is quite prone to becoming ``stuck'' in suboptimal solutions.  In general, hill climbing is not recommended for permutation search, but the method may prove useful in certain circumstances.  

Approximate complexity: on the order of \eqn{|V(G)|^2} per iteration, total complexity dependent on the number of iterations.

\item simulated annealing (``anneal''): The (fairly simple) annealing procedure here employed proceeds as follows.  At each iteration, the set of all permissible binary exchanges (if \code{full.neighborhood==TRUE}) or a random selection from this set is evaluated.  If a superior option is identified, the best of these is chosen.  If no superior options are found, then the algorithm chooses randomly from the set of alternatives with probability equal to the current temperature, otherwise retaining its prior solution.  After each iteration, the current temperature is reduced by a factor equal to \code{prob.decay}; the initial temperature is set by \code{prob.init}.  When a number of iterations equal to \code{freeze.time} have been completed, the algorithm ``freezes.''  Once ``frozen,'' the annealer hillclimbs from its present location until no improvement is found, and terminates.  At termination, the best permutation identified so far is utilized; this need not be the most recent position (though it sometimes is).  

Simulated annealing is sometimes called ``noisy hill climbing'' because it uses the introduction of random variation to a hill climbing routine to avoid convergence to local optima; it works well on reasonably correlated search spaces with well-defined solution neighborhoods, and is far more robust than hill climbing algorithms.  As a general rule, simulated annealing is recommended here for most graphs up to size approximately 50.  At this point, computational complexity begins to become a serious barrier, and alternative methods may be more practical.  

Approximate complexity: on the order of \eqn{|V(G)|^2}*\code{freeze.time} if \code{full.neighborhood==TRUE}, otherwise complexity scales approximately linearly with \code{freeze.time}.  This can be misleading, however, since failing to search the full neighborhood generally requires that \code{freeze.time} be greatly increased.)

\item blind monte carlo search (``mc''): Blind monte carlo search, as the name implies, consists of randomly drawing a sample of permutations from the accessible permutation set and selecting the best.  Although this not such a bad option when A) a large fraction of points are optimal or nearly optimal and B) the search space is largely uncorrelated, these conditions do not seem to characterize most permutation search problems.  Blind monte carlo search is not generally recommended, but it is provided as an option should it be desired (e.g., when it is absolutely necessary to control the number of permutations examined).  

Approximate complexity: linear in \code{draws}.

\item extreme value estimation (``gumbel''): Extreme value estimation attempts to estimate a global optimum via stochastic modeling of the distribution of the graph statistic over the space of accessible permutations.  The algorithm currently proceeds as follows.  First, a random sample is taken from the accessible permutation set (as with monte carlo search, above).  Next, this sample is used to fit an extreme value (gumbel) model; the gumbel distribution is the limiting distribution of the extreme values from samples under a continuous, unbounded distribution, and we use it here as an approximation.  Having fit the model, an associated statistic (the mean, median, or mode as determined by \code{estimator}) is then used as an estimator of the global optimum.

Obviously, this approach has certain drawbacks.  First of all, our use of the gumbel model in particular assumes an unbounded, continuous underlying distribution, which may or may not be approximately true for any given problem.  Secondly, the inherent non-robustness of extremal problems makes the fact that our prediction rests on a string of approximations rather worrisome: our idea of the shape of the underlying distribution could be distorted by a bad sample, our parameter estimation could be somewhat off, etc., any of which could have serious consequences for our extremal prediction.  Finally, the prediction which is made by the extreme value model is \emph{nonconstructive}, in the sense that \emph{no permutation need have been found by the algorithm which induces the predicted value}.  On the bright side, this \emph{could} allow one to estimate the optimum without having to find it directly; on the dark side, this means that the reported optimum could be a numerical chimera.

At this time, extreme value estimation should be considered \emph{experimental}, and \emph{is not recommended for use on substantive problems.}  \code{lab.optimize.gumbel} is not guaranteed to work properly, or to produce intelligible results; this may eventually change in future revisions, or the routine may be scrapped altogether.

Approximate complexity: linear in \code{draws}.
}

This list of algorithms is itself somewhat unstable: some additional techniques (canonical labeling and genetic algorithms, for instance) may be added, and some existing methods (e.g., extreme value estimation) may be modified or removed.  Every attempt will be made to keep the command format as stable as possible for other routines (e.g., \code{\link{gscov}}, \code{\link{structdist}}) which depend on \code{lab.optimize} to do their heavy-lifting.  In general, it is not expected that the end-user will call \code{lab.optimize} directly; instead, most end-user interaction with these routines will be via the structural distance/covariance functions which used them.
}
\value{
The estimated global optimum of \code{FUN} over the set of relabelings permitted by \code{exchange.list}
}
\references{ 
Butts, C.T. and Carley, K.M.  (2005).  \dQuote{Some Simple Algorithms for Structural Comparison.}  \emph{Computational and Mathematical Organization Theory,} 11(4), 291-305.

Butts, C.T., and Carley, K.M.  (2001).  \dQuote{Multivariate Methods for Interstructural Analysis.}  CASOS Working Paper, Carnegie Mellon University. 
}
\author{ Carter T. Butts \email{buttsc@uci.edu}  }
%\note{  }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{gscov}}, \code{\link{gscor}}, \code{\link{structdist}}, \code{\link{sdmat}} }

\examples{
#Generate a random graph and copy it
g<-rgraph(10)
g2<-rmperm(g)  #Permute the copy randomly

#Seek the maximum correlation
lab.optimize(g,g2,gcor,seek="max",opt.method="anneal",freeze.time=50,
    prob.decay=0.9)

#These two don't do so well...
lab.optimize(g,g2,gcor,seek="max",opt.method="hillclimb")     
lab.optimize(g,g2,gcor,seek="max",opt.method="mc",draws=1000)

}
\keyword{ optimize }%-- one or more ...
\keyword{ math }
\keyword{ graphs }
