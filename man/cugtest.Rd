\name{cugtest}
\alias{cugtest}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Perform Conditional Uniform Graph (CUG) Hypothesis Tests for Graph-Level Indices}
\description{
 \code{cugtest} tests an arbitrary GLI (computed on \code{dat} by \code{FUN}) against a conditional uniform graph null hypothesis, via Monte Carlo simulation.  Some variation in the nature of the conditioning is available; currently, conditioning only on size, conditioning jointly on size and estimated tie probability (via expected density), and conditioning jointly on size and (bootstrapped) edge value distributions are implemented.  Note that fair amount of flexibility is possible regarding CUG tests on functions of GLIs (Anderson et al., 1999).  See below for more details. 
}
\usage{
cugtest(dat, FUN, reps=1000, gmode="digraph", cmode="density", 
    diag=FALSE, g1=1, g2=2, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ graph(s) to be analyzed. }
  \item{FUN}{ function to compute GLIs, or functions thereof.  \code{FUN} must accept \code{dat} and the specified \code{g} arguments, and should return a real number. }
  \item{reps}{ integer indicating the number of draws to use for quantile estimation.  Note that, as for all Monte Carlo procedures, convergence is slower for more extreme quantiles.  By default, \code{reps==1000}. }
  \item{gmode}{ string indicating the type of graph being evaluated.  "digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  \code{gmode} is set to "digraph" by default. }
  \item{cmode}{ string indicating the type of conditioning assumed by the null hypothesis.  If \code{cmode} is set to "density", then the density of the graph in question is used to determine the tie probabilities of the Bernoulli graph draws (which are also conditioned on |V(G)|).  If\code{cmode=="ties"}, then draws are bootstrapped from the distribution of edge values within the data matrices.  If \code{cmode="order"}, then draws are uniform over all graphs of the same order (size) as the graphs within the input stack.  By default, \code{cmode} is set to \code{"density"}.}
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{g1}{ integer indicating the index of the first graph input to the GLI.  By default, \code{g1==1}. }
  \item{g2}{ integer indicating the index of the second graph input to the GLI.  (\code{FUN} can ignore this, if one wishes to test the GLI value of a single graph, but it should recognize the argument.)  By default, \code{g2==2}. }
  \item{\dots}{ additional arguments to \code{FUN}. }
}
\details{
  The null hypothesis of the CUG test is that the observed GLI (or function thereof) was drawn from a distribution equivalent to that of said GLI evaluated (uniformly) on the space of all graphs conditional on one or more features.  The most common \dQuote{features} used for conditioning purposes are order (size) and density, both of which are known to have strong and nontrivial effects on other GLIs (Anderson et al., 1999) and which are, in many cases, exogenously determined.  (Note that maximum entropy distributions conditional on expected statistics are not in general correctly referred to as \dQuote{conditional uniform graphs}, but have been described as such for independent-dyad models; this is indeed the case for this function, although such terminology is not really proper.  See \code{\link{cug.test}} for CUG tests with exact conditioning.)  Since theoretical results regarding functions of arbitrary GLIs on the space of graphs are not available, the standard approach to CUG testing is to approximate the quantiles of the observed statistic associated with the null hypothesis using Monte Carlo methods.  This is the technique utilized by \code{cugtest}, which takes appropriately conditioned draws from the set of graphs and computes on them the GLI specified in \code{FUN}, thereby accumulating an approximation to the true quantiles.

The \code{cugtest} procedure returns a \code{cugtest} object containing the estimated distribution of the test GLI under the null hypothesis, the observed GLI value of the data, and the one-tailed p-values (estimated quantiles) associated with said observation.  As usual, the (upper tail) null hypothesis is rejected for significance level alpha if p>=observation is less than alpha (or p<=observation, for the lower tail).  Standard caveats regarding the use of null hypothesis testing procedures are relevant here: in particular, bear in mind that a significant result does not necessarily imply that the likelihood ratio of the null model and the alternative hypothesis favors the latter.

Informative and aesthetically pleasing portrayals of \code{cugtest} objects are available via the \code{\link{print.cugtest}} and \code{\link{summary.cugtest}} methods.  The \code{\link{plot.cugtest}} method displays the estimated distribution, with a reference line signifying the observed value.
}
\value{
An object of class \code{cugtest}, containing

\item{testval}{
The observed GLI value.
}
\item{dist}{
A vector containing the Monte Carlo draws. 
}
\item{pgreq}{
The proportion of draws which were greater than or equal to the observed GLI value.
}
\item{pleeq}{
The proportion of draws which were less than or equal to the observed GLI value.
}
}
\references{ 
Anderson, B.S.; Butts, C.T.; and Carley, K.M. (1999). \dQuote{The Interaction of Size and Density with Graph-Level Indices.} \emph{Social Networks}, 21(3), 239-267.
 }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{This function currently conditions only on expected statistics, and is somewhat cumbersome.  \code{\link{cug.test}} is now recommended for univariate CUG tests (and will eventually supplant this function).}



\seealso{ \code{\link{cug.test}}, \code{\link{qaptest}}, \code{\link{gliop}} }

\examples{
#Draw two random graphs, with different tie probabilities
dat<-rgraph(20,2,tprob=c(0.2,0.8))
#Is their correlation higher than would be expected, conditioning 
#only on size?
cug<-cugtest(dat,gcor,cmode="order")
summary(cug)
#Now, let's try conditioning on density as well.
cug<-cugtest(dat,gcor)
summary(cug)
}
\keyword{htest}
\keyword{ math }
\keyword{ graphs }

