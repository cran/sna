\name{bn}
\alias{bn}
\alias{coef.bn}
\alias{bn.nlpl.dyad}
\alias{bn.nlpl.edge}
\alias{bn.nlpl.triad}
\alias{bn.nltl}
\alias{plot.bn}
\alias{print.bn}
\alias{print.summary.bn}
\alias{summary.bn}
\alias{bn_dyadstats_R}
\alias{bn_lpl_dyad_R}
\alias{bn_lpl_triad_R}
\alias{bn_ptriad_R}
\alias{bn_triadstats_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit a Biased Net Model }
\description{
  Fits a biased net model to an input graph, using moment-based or maximum pseudolikelihood techniques.
}
\usage{
bn(dat, method = c("mple.triad", "mple.dyad", "mple.edge", 
    "mtle"), param.seed = NULL, param.fixed = NULL, 
    optim.method = "BFGS", optim.control = list(), 
    epsilon = 1e-05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{ a single input graph. }
  \item{method}{ the fit method to use (see below). }
  \item{param.seed}{ seed values for the parameter estimates. }
  \item{param.fixed}{ parameter values to fix, if any. }
  \item{optim.method}{ method to be used by \code{\link{optim}}. }
  \item{optim.control}{ control parameter for \code{\link{optim}}. }
  \item{epsilon}{ tolerance for convergence to extreme parameter values (i.e., 0 or 1). }
}
\details{
The biased net model stems from early work by Rapoport, who attempted to model networks via a hypothetical "tracing" process.  This process may be described loosely as follows.  One begins with a small "seed" set of vertices, each member of which is assumed to nominate (generate ties to) other members of the population with some fixed probability.  These members, in turn, may nominate new members of the population, as well as members who have already been reached.  Such nominations may be "biased" in one fashion or another, leading to a non-uniform growth process.  Specifically, let \eqn{e_{ij}}{e_ij} be the random event that vertex \eqn{i} nominates vertex \eqn{j} when reached.  Then the conditional probability of \eqn{e_{ij}}{e_ij} is given by

\deqn{
\Pr(e_{ij}|T) = 1-\left(1-\Pr(B_e)\right) \prod_k \left(1-\Pr(B_k|T)\right)
}{%
Pr(e_ij | T) = 1 - (1-Pr(B_e)) prod_k (1 - Pr(B_k | T))
}

where \eqn{T} is the current state of the trace, \eqn{B_e} is the a Bernoulli event corresponding to the baseline probability of \eqn{e_{ij}}{e_ij}, and the \eqn{B_k} are "bias events."  Bias events are taken to be independent Bernoulli trials, given \eqn{T}, such that \eqn{e_{ij}}{e_ij} is observed with certainty if any bias event occurs.  The specification of a biased net model, then, involves defining the various bias events (which, in turn, influence the structure of the network).

Although other events have been proposed, the primary bias events employed in current biased net models are the "parent bias" (a tendency to return nominations); the "sibling bias" (a tendency to nominate alters who were nominated by the same third party); and the "double role bias" (a tendency to nominate alters who are both siblings and parents).  These bias events, together with the baseline edge events, are used to form the standard biased net model.  It is standard to assume homogeneity within bias class, leading to the four parameters \eqn{\pi}{pi} (probability of a parent bias event), \eqn{\sigma}{sigma} (probability of a sibling bias event), \eqn{\rho}{rho} (probability of a double role bias event), and \eqn{d} (probability of a baseline event).  

Unfortunately, there is no simple expression for the likelihood of a graph given these parameters (and hence, no basis for likelihood based inference).  However, Skvoretz et al. have derived a class of maximum pseudo-likelihood estimators for the the biased net model, based on local approximations to the likelihood at the edge, dyad, or triad level.  These estimators may be employed within \code{bn} by selecting the appropriate MPLE for the \emph{method} argument.  Alternately, it is also possible to derive expected triad census rates for the biased net model, allowing an estimator which maximizes the likelihood of the observed triad census (essentially, a method of moments procedure).  This last may be selected via the argument \code{mode="mtle"}.  In addition to estimating model parameters, \code{bn} generates predicted edge, dyad, and triad census statistics, as well as structure statistics (using the Fararo-Sunshine recurrence).  These can be used to evaluate goodness-of-fit.  

\code{print}, \code{summary}, and \code{plot} methods are available for \code{bn} objects.  See \code{\link{rgbn}} for simulation from biased net models.
}
\value{
  An object of class \code{bn}.
}
\references{ 
Fararo, T.J. and Sunshine, M.H. (1964).  ``A study of a biased friendship net.''  Syracuse, NY: Youth Development Center.

Rapoport, A.  (1957).  ``A contribution to the theory of random and biased nets.''  \emph{Bulletin of Mathematical Biophysics,} 15, 523-533.

Skvoretz, J.; Fararo, T.J.; and Agneessens, F.  (2004).  ``Advances in biased net theory: definitions, derivations, and estimations.''  \emph{Social Networks,} 26, 113-139.
}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ Asymptotic properties of the MPLE are not known for this model.  Caution is strongly advised. }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{rgbn}}, \code{\link{structure.statistics}} }
\examples{
#Generate a random graph
g<-rgraph(25)

#Fit a biased net model, using the triadic MPLE
gbn<-bn(g)

#Examine the results
summary(gbn)
plot(gbn)

#Now, fit a model containing only a density parameter
gbn<-bn(g,param.fixed=list(pi=0,sigma=0,rho=0))
summary(gbn)
plot(gbn)

}
\keyword{ graphs }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
