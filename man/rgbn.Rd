\name{rgbn}
\alias{rgbn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Draw from a Skvoretz-Fararo Biased Net Process }
\description{
  Produces a series of draws from a Skvoretz-Fararo biased net process using a (pseudo) Gibbs sampler or exact sampling procedure.
}
\usage{
rgbn(n, nv, param = list(pi=0, sigma=0, rho=0, d=0.5, delta=0),
    burn = nv*nv*5*100, thin = nv*nv*5, maxiter = 1e7,
    method = c("mcmc","cftp"), dichotomize.sib.effects = FALSE,
    return.as.edgelist = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ number of draws to take. }
  \item{nv}{ number of vertices in the graph to be simulated. }
  \item{param}{ a list containing the biased net parameters (as described below); \eqn{d} may be given as a scalar or as an \code{nv x nv} matrix of edgewise baseline edge probabilities.}
  \item{burn}{ for the Gibbs sampler, the number of burn-in draws to take (and discard). }
  \item{thin}{ the thinning parameter for the Gibbs sampler. }
  \item{maxiter}{ for the CFTP method, the number of iterations to try before giving up.}
  \item{method}{ \code{"mcmc"} for the Gibbs sampler, or \code{"cftp"} for the exact sampling procedure.}
  \item{dichotomize.sib.effects}{ logical; should sibling and double role effects be dichotomized?} 
  \item{return.as.edgelist}{ logical; should the simulated draws be returned in edgelist format?} 
}
\details{
The biased net model stems from early work by Rapoport, who attempted to model networks via a hypothetical \dQuote{tracing} process.  This process may be described loosely as follows.  One begins with a small \dQuote{seed} set of vertices, each member of which is assumed to nominate (generate ties to) other members of the population with some fixed probability.  These members, in turn, may nominate new members of the population, as well as members who have already been reached.  Such nominations may be \dQuote{biased} in one fashion or another, leading to a non-uniform growth process. 

While the original biased net model depends upon the tracing process, a local interpretation has been put forward by Skvoretz and colleagues in recent years.  Using the standard four-parameter process, the conditional probability of an \eqn{(i,j)} edge given all other edges in a random graph \eqn{G} can be approximated as

\deqn{
\Pr(i \to j|G_{-ij}) \approx 1 - (1-\rho)^z (1-\sigma)^y (1-\pi)^x (1-d_{ij})
}{%
Pr(i->j | G_{-ij}) apx 1 - (1-rho)^z (1-sigma)^y (1-pi)^x (1-d_{ij})
}

where \eqn{x=1} iff \eqn{j \to i} (and 0 otherwise), \eqn{y} is the number of vertices \eqn{k \neq i,j}{k!=i,j} such that \eqn{k \to i, k \to j}{k->i, k->j}, and \eqn{z=1} iff \eqn{x=1} and \eqn{y>0} (and 0 otherwise).  Thus, \eqn{x} is the number of potential \emph{parent bias} events, \eqn{y} is the number of potential \emph{sibling bias} events, and \eqn{z} is the number of potential \emph{double role bias} events.  \eqn{d_{ij}} is the probability of the baseline edge event; note that an edge arises if the baseline event or any bias event occurs, and all events are assumed conditionally independent.  Written in this way, it is clear that the edges of \eqn{G} are conditionally independent if they share no endpoint.  Thus, a  model with the above structure should be a subfamily of the Markov graphs.

One potential problem with the above structure is that the hypothetical probabilities implied by the model are not guaranteed to be consistent - that is, the conditions under which there exists a joint pmf with the implied full conditionals are currently unknown (and may be restrictive).  The interpretation of the above as exact conditional probabilities is thus potentially problematic.  However, a well-defined process can be constructed by interpreting the above as transition probabilities for a Markov chain that evolves by updating a randomly selected edge variable at each time point; this is a Gibbs sampler for the implied joint pmf where it exists, and otherwise an irreducible and aperiodic Markov chain with a well-defined equilibrium distribution.

In the above process, all events act to promote the formation of edges; it is also possible to define events that inhibit them.  For instance, consider a \emph{satiation} event that, if it occurs, forbids the creation of an \eqn{i \to j}{i->j} edge; we assume that a potential satiation event occurs every time \eqn{i} emits an edge to some other vertex.  The associated approximate conditional (i.e., transition probability) is

\deqn{
\Pr(i \to j|G_{-ij}) \approx (1-\delta)^w\left(1 - (1-\rho)^z (1-\sigma)^y (1-\pi)^x (1-d_{ij})\right)
}{%
Pr(i->j | G_{-ij}) apx (1-delta)^w [1 - (1-rho)^z (1-sigma)^y (1-pi)^x (1-d_{ij})]
}

where \eqn{w} is the outdegree of \eqn{i} in \eqn{G_{-ij}} and \eqn{\delta}{delta} is the probability of the satiation event.  The net effect of satiation is to suppress edge formation (in roughly geometric fashion) on high degree nodes.  This may be useful in preventing degeneracy when using sigma and rho effects.  Degeneracy can also be reduced by employing the \code{dichotomize.sib.effects} argument, which counts only the first shared partner's contribution towards sibling and double role effects.

It should be noted that the above process is not entirely consistent with the tracing-based model, which is itself not uniformly well-specified in the literature.  For this reason, the local model is referred to here as a Skvoretz-Fararo graph process.  One significant advantage of this process is that it is well-defined, and easily simulated: the above equation can be used to form the basis of a (pseudo-) Gibbs sampler, which is used by \eqn{rgbn} to take draws from the (local) biased net model.  Burn-in and thinning are controlled by the corresponding arguments; since degeneracy is common with models of this type, it is advisable to check for adequate mixing.  An alternative simulation strategy is the exact sampling procedure of Butts (2009), which employs a form of coupling from the past (CFTP).  The CFTP method generates exact, independent draws from the equilibrium distribution of the biased net process (up to numerical limits), but can be slow to attain coalescence (and does not currently support satiation events).  Setting \code{maxiter} to smaller values limits the search depth employed, at the possible cost of biasing the resulting sample.
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
