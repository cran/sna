\name{rgbn}
\alias{rgbn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Draw from a Skvoretz-Fararo Biased Net Process }
\description{
  Produces a series of draws from a Skvoretz-Fararo biased net process using a Gibbs sampler or exact sampling procedure.
}
\usage{
rgbn(n, nv, param = list(pi=0, sigma=0, rho=0, d=0.5),
    burn = nv*nv*5*100, thin = nv*nv*5, maxiter = Inf,
    method = c("mcmc","cftp"), return.as.edgelist = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ number of draws to take. }
  \item{nv}{ number of vertices in the graph to be simulated. }
  \item{param}{ a list containing the biased net parameters. }
  \item{burn}{ for the Gibbs sampler, the number of burn-in draws to take (and discard). }
  \item{thin}{ the thinning parameter for the Gibbs sampler. }
  \item{maxiter}{ for the CFTP method, the number of iterations to try before giving up.}
  \item{method}{ \code{"mcmc"} for the Gibbs sampler, or \code{"cftp"} for the exact sampling procedure.}
  \item{return.as.edgelist}{ logical; should the simulated draws be returned in edgelist format?} 
}
\details{
The biased net model stems from early work by Rapoport, who attempted to model networks via a hypothetical \dQuote{tracing} process.  This process may be described loosely as follows.  One begins with a small \dQuote{seed} set of vertices, each member of which is assumed to nominate (generate ties to) other members of the population with some fixed probability.  These members, in turn, may nominate new members of the population, as well as members who have already been reached.  Such nominations may be \dQuote{biased} in one fashion or another, leading to a non-uniform growth process. 

While the original biased net model depends upon the tracing process, a local interpretation has been put forward by Skvoretz and colleagues in recent years.  Using the standard four-parameter process, the conditional probability of an \eqn{(i,j)} edge given all other edges in a random graph \eqn{G} can be written as

\deqn{
\Pr(i \to j) = 1 - (1-\rho)^z (1-\sigma)^y (1-\pi)^x (1-d)
}{%
Pr(i->j) = 1 - (1-rho)^z (1-sigma)^y (1-pi)^x (1-d)
}

where \eqn{x=1} iff \eqn{j \to i} (and 0 otherwise), \eqn{y} is the number of vertices \eqn{k \neq i,j}{k!=i,j} such that \eqn{k \to i, k \to j}{k->i, k->j}, and \eqn{z=1} iff \eqn{x=1} and \eqn{y>0} (and 0 otherwise).  Thus, \eqn{x} is the number of \emph{parent bias} events, \eqn{y} is the number of \emph{sibling bias} events, and \eqn{z} is the number of \emph{double role bias} events.  \eqn{d} is the probability of the baseline edge event; note that an edge arises if the baseline event or any bias event occurs, and all events are assumed conditionally independent.  Written in this way, it is clear that the edges of \eqn{G} are conditionally independent if they share no endpoint.  Thus, the above model is a subfamily of the Markov graphs.

It should be noted that the above process is not entirely consistent with the tracing-based model, which is itself not uniformly well-specified in the literature.  For this reason, the local model is referred to here as a Skvoretz-Fararo graph process.  One significant advantage of this process is that it is well-defined, and easily simulated: the above equation can be used to form the basis of a Gibbs sampler, which is used by \eqn{rgbn} to take draws from the (local) biased net model.  Burn-in and thinning are controlled by the corresponding arguments; since degeneracy is common with models of this type, it is advisable to check for adequate mixing.  An alternative simulation strategy is the exact sampling procedure of Butts (2008), which employs a form of coupling from the past (CFTP).  The CFTP method generates exact, independent draws (up to numerical limits), but can be slow to attain coalescence.  Setting \code{maxiter<Inf} limits the search depth employed, at the possible cost of biasing the resulting sample.
}
\value{
  An adjacency array containing the simulated graphs.
}
\references{ 
Butts, C.T.  (2009).  \dQuote{A Perfect Sampling Method for Exponential Random Graph Models}.  Working paper, University of California, Irvine.

Rapoport, A.  (1957).  \dQuote{A Contribution to the Theory of Random and Biased Nets.}  \emph{Bulletin of Mathematical Biophysics,} 15, 523-533.

Skvoretz, J.; Fararo, T.J.; and Agneessens, F.  (2004).  \dQuote{Advances in Biased Net Theory: Definitions, Derivations, and Estimations.}  \emph{Social Networks,} 26, 113-139.
}
\author{ Carter T. Butts \email{buttsc@uci.edu}}
%\note{  }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{bn}} }
\examples{
#Generate draws with low density and no biases
g1<-rgbn(50,10,param=list(pi=0, sigma=0, rho=0, d=0.17))
apply(dyad.census(g1),2,mean) #Examine the dyad census

#Add a reciprocity bias
g2<-rgbn(50,10,param=list(pi=0.5, sigma=0, rho=0, d=0.17))
apply(dyad.census(g2),2,mean) #Compare with g1

#Alternately, add a sibling bias
g3<-rgbn(50,10,param=list(pi=0.0, sigma=0.3, rho=0, d=0.17))
mean(gtrans(g3))              #Compare transitivity scores
mean(gtrans(g1))

}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ graphs }% __ONLY ONE__ keyword per line
