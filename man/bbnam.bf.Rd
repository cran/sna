\name{bbnam.bf}
\alias{bbnam.bf}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Estimate Bayes Factors for the bbnam }
\description{
  This function uses monte carlo integration to estimate the BFs, and tests the fixed probability, pooled, and pooled by actor models. (See \code{\link{bbnam}} for details.)
}
\usage{
bbnam.bf(dat, nprior=matrix(rep(0.5, dim(dat)[1]^2), 
    nrow = dim(dat)[1], ncol = dim(dat)[1]), em.fp=0.5, ep.fp=0.5, 
    emprior.pooled=c(1, 1), epprior.pooled=c(1, 1), 
    emprior.actor=cbind(rep(1, dim(dat)[1]), rep(1, dim(dat)[1])), 
    epprior.actor=cbind(rep(1, dim(dat)[1]), rep(1, dim(dat)[1])), 
    diag=FALSE, mode="digraph", reps=1000)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ Data array to be analyzed.  This array must be of dimension m x n x n, where n is |V(G)|, the first dimension indexes the observer, the second indexes the sender of the relation, and the third dimension indexes the recipient of the relation.  (E.g., \code{dat[i,j,k]==1} implies that i observed j sending the relation in question to k.)  Note that only dichotomous data is supported at present, and missing values are permitted; the data collection pattern, however, is assumed to be ignorable, and hence the posterior inferences are implicitly conditional on the observation pattern. }
  \item{nprior}{ Network prior matrix.  This must be a matrix of dimension n x n, containing the arc/edge priors for the criterion network.  (E.g., \code{nprior[i,j]} gives the prior probability of i sending the relation to j in the criterion graph.)  If no network prior is provided, an uninformative prior on the space of networks will be assumed (i.e., p(i->j)=0.5).  Missing values are not allowed. }
  \item{em.fp}{ Probability of false negatives for the fixed probability model }
  \item{ep.fp}{ Probability of false positives for the fixed probability model }
  \item{emprior.pooled}{ (alpha,beta) pairs for the (beta) false negative prior under the pooled model }
  \item{epprior.pooled}{ (alpha,beta) pairs for the (beta) false positive prior under the pooled model }
  \item{emprior.actor}{ Matrix of per observer (alpha,beta) pairs for the (beta) false negative prior under the per observer/actor model }
  \item{epprior.actor}{ Matrix of per observer (alpha,beta) pairs for the (beta) false negative prior under the per observer/actor model }
  \item{diag}{ Boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the criterion graph can contain loops.  Diag is false by default. }
  \item{mode}{ String indicating the type of graph being evaluated.  "Digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  Mode is set to "digraph" by default. }
  \item{reps}{ Number of monte carlo draws to take }
}
\details{
  The bbnam model (detailed in the \code{\link{bbnam}} function help) is a fairly simple model for integrating informant reports regarding social network data.  \code{bbnam.bf} computes Bayes Factors (integrated likelihood ratios) for the three error submodels of the bbnam: fixed error probabilities, pooled error probabilities, and per observer/actor error probabilities.  
}
\value{
   An object of class \code{bayes.factor}.
}
\references{
Butts, C. T.  (2003).  ``Network Inference, Error, and Informant (In)Accuracy: A Bayesian Approach.''  \emph{Social Networks}, 25(2), 103-140.

Robert, C.  (1994).  \emph{The Bayesian Choice: A Decision-Theoretic Motivation.}  Springer. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ It is important to be aware that the model parameter priors are essential components of the models to be compared; inappropriate parameter priors will result in misleading Bayes Factors. }


\seealso{ \code{\link{bbnam}} }

\examples{
}

\keyword{models}
\keyword{ math }
