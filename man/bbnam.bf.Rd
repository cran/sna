\name{bbnam.bf}
\alias{bbnam.bf}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Estimate Bayes Factors for the bbnam }
\description{
  This function uses monte carlo integration to estimate the BFs, and tests the fixed probability, pooled, and pooled by actor models. (See \code{\link{bbnam}} for details.)
}
\usage{
bbnam.bf(dat, nprior=0.5, em.fp=0.5, ep.fp=0.5, emprior.pooled=c(1, 11),
    epprior.pooled=c(1, 11), emprior.actor=c(1, 11), epprior.actor=c(1, 11), 
    diag=FALSE, mode="digraph", reps=1000)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ Input networks to be analyzed.  This may be supplied in any reasonable form, but must be reducible to an array of dimension \eqn{m \times n \times n}{m x n x n}, where \eqn{n} is \eqn{|V(G)|}, the first dimension indexes the observer (or information source), the second indexes the sender of the relation, and the third dimension indexes the recipient of the relation.  (E.g., \code{dat[i,j,k]==1} implies that i observed j sending the relation in question to k.)  Note that only dichotomous data is supported at present, and missing values are permitted; the data collection pattern, however, is assumed to be ignorable, and hence the posterior draws are implicitly conditional on the observation pattern.}
  \item{nprior}{ Network prior matrix.  This must be a matrix of dimension \eqn{n} x \eqn{n}, containing the arc/edge priors for the criterion network.  (E.g., \code{nprior[i,j]} gives the prior probability of \code{i} sending the relation to \code{j} in the criterion graph.)  Non-matrix values will be coerced/expanded to matrix form as appropriate.  If no network prior is provided, an uninformative prior on the space of networks will be assumed (i.e., \eqn{\Pr(i\to j)=0.5}{Pr(i->j)=0.5}).  Missing values are not allowed. }
  \item{em.fp}{ Probability of false negatives for the fixed probability model }
  \item{ep.fp}{ Probability of false positives for the fixed probability model }
  \item{emprior.pooled}{ \eqn{(\alpha,\beta)}{(alpha,beta)} pairs for the (beta) false negative prior under the pooled model }
  \item{epprior.pooled}{ \eqn{(\alpha,\beta)}{(alpha,beta)} pairs for the (beta) false positive prior under the pooled model }
  \item{emprior.actor}{ Matrix of per observer \eqn{(\alpha,\beta)}{(alpha,beta)} pairs for the (beta) false negative prior under the per observer/actor model, or something that can be coerced to this form }
  \item{epprior.actor}{ Matrix of per observer (\eqn{(\alpha,\beta)}{(alpha,beta)} pairs for the (beta) false positive prior under the per observer/actor model, or something that can be coerced to this form }
  \item{diag}{ Boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the criterion graph can contain loops.  Diag is false by default. }
  \item{mode}{ String indicating the type of graph being evaluated.  \code{"digraph"} indicates that edges should be interpreted as directed; \code{"graph"} indicates that edges are undirected.  Mode is set to \code{"digraph"} by default. }
  \item{reps}{ Number of Monte Carlo draws to take }
}
\details{
  The bbnam model (detailed in the \code{\link{bbnam}} function help) is a fairly simple model for integrating informant reports regarding social network data.  \code{bbnam.bf} computes log Bayes Factors (integrated likelihood ratios) for the three error submodels of the bbnam: fixed error probabilities, pooled error probabilities, and per observer/actor error probabilities.
  
  By default, \code{bbnam.bf} uses weakly informative Beta(1,11) priors for false positive and false negative rates, which may not be appropriate for all cases.  (Likewise, the initial network prior is uniformative.)  Users are advised to consider adjusting the error rate priors when using this function in a practical context; for instance, it is often reasonable to expect higher false negative rates (on average) than false positive rates, and to expect the criterion graph density to be substantially less than 0.5.  See the reference below for a discussion of this issue.
}
\value{
   An object of class \code{bayes.factor}.
}
\references{
Butts, C. T.  (2003).  \dQuote{Network Inference, Error, and Informant (In)Accuracy: A Bayesian Approach.}  \emph{Social Networks}, 25(2), 103-140.

Robert, C.  (1994).  \emph{The Bayesian Choice: A Decision-Theoretic Motivation.}  Springer. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ It is important to be aware that the model parameter priors are essential components of the models to be compared; inappropriate parameter priors will result in misleading Bayes Factors. }


\seealso{ \code{\link{bbnam}} }

\examples{
#Create some random data from the "pooled" model
g<-rgraph(7)
g.p<-0.8*g+0.2*(1-g)
dat<-rgraph(7,7,tprob=g.p)

#Estimate the log Bayes Factors
b<-bbnam.bf(dat,emprior.pooled=c(3,5),epprior.pooled=c(3,5),
    emprior.actor=c(3,5),epprior.actor=c(3,5))
#Print the results
b
}

\keyword{models}
\keyword{ math }
