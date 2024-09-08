\name{rgbn}
\alias{rgbn}
\alias{bn_cftp_R}
\alias{bn_mcmc_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Draw from a Skvoretz-Fararo Biased Net Process }
\description{
  Produces a series of draws from a Skvoretz-Fararo biased net process using a Markov chain Monte Carlo or exact sampling procedure.
}
\usage{
rgbn(n, nv, param = list(pi=0, sigma=0, rho=0, d=0.5, delta=0, 
    epsilon=0), burn = nv*nv*5*100, thin = nv*nv*5, maxiter = 1e7,
    method = c("mcmc","cftp"), dichotomize.sib.effects = FALSE,
    return.as.edgelist = FALSE, seed.graph = NULL, max.density = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ number of draws to take. }
  \item{nv}{ number of vertices in the graph to be simulated. }
  \item{param}{ a list containing the biased net parameters (as described below); \eqn{d} and \eqn{epsilon} may be given as scalars or as \code{nv x nv} matrices of edgewise event probabilities.}
  \item{burn}{ for the MCMC, the number of burn-in draws to take (and discard). }
  \item{thin}{ the thinning parameter for the MCMC algorithm. }
  \item{maxiter}{ for the CFTP method, the number of iterations to try before giving up.}
  \item{method}{ \code{"mcmc"} for the MCMC, or \code{"cftp"} for the exact sampling procedure.}
  \item{dichotomize.sib.effects}{ logical; should sibling and double role effects be dichotomized?} 
  \item{return.as.edgelist}{ logical; should the simulated draws be returned in edgelist format?} 
  \item{seed.graph}{ optionally, an initial state to use for MCMC.}
  \item{max.density}{ optional maximum density threshold for MCMC; if the chain encounters a graph of higher than max density, the chain is terminated (and the result flagged). }
}
\details{
The biased net model stems from early work by Rapoport, who attempted to model networks via a hypothetical \dQuote{tracing} process.  This process may be described loosely as follows.  One begins with a small \dQuote{seed} set of vertices, each member of which is assumed to nominate (generate ties to) other members of the population with some fixed probability.  These members, in turn, may nominate new members of the population, as well as members who have already been reached.  Such nominations may be \dQuote{biased} in one fashion or another, leading to a non-uniform growth process. 

While the original biased net model depends upon the tracing process, a local (conditional) interpretation was put forward by Skvoretz and colleagues (2004).  Using a four-parameter model, they propose approximating the conditional probability of an \eqn{(i,j)} edge given all other edges in a random graph \eqn{G} by

\deqn{
\Pr(i \to j|G_{-ij}) \approx 1 - (1-\rho)^z (1-\sigma)^y (1-\pi)^x (1-d_{ij})
}{%
Pr(i->j | G_{-ij}) apx 1 - (1-rho)^z (1-sigma)^y (1-pi)^x (1-d_{ij})
}

where \eqn{x=1} iff \eqn{j \to i} (and 0 otherwise), \eqn{y} is the number of vertices \eqn{k \neq i,j}{k!=i,j} such that \eqn{k \to i, k \to j}{k->i, k->j}, and \eqn{z=1} iff \eqn{x=1} and \eqn{y>0} (and 0 otherwise).  Thus, \eqn{x} is the number of potential \emph{parent bias} events, \eqn{y} is the number of potential \emph{sibling bias} events, and \eqn{z} is the number of potential \emph{double role bias} events.  \eqn{d_{ij}} is the probability of the baseline edge event; note that an edge arises if the baseline event or any bias event occurs, and all events are assumed conditionally independent.  Written in this way, it is clear that the edges of \eqn{G} are conditionally independent if they share no endpoint.  Thus, a  model with the above structure should be a subfamily of the Markov graphs.

One problem with the above structure is that the hypothetical probabilities implied by the model are not in general consistent - that is, there exist conditions under which there is no joint pmf for \eqn{G} with the implied full conditionals.  The interpretation of the above as exact conditional probabilities is thus potentially problematic.  However, a well-defined process can be constructed by interpreting the above as transition probabilities for a Markov chain that evolves by updating a randomly selected edge variable at each time point; this is a Gibbs sampler for the implied joint pmf where it exists, and otherwise an irreducible and aperiodic Markov chain with a well-defined equilibrium distribution (Butts, 2018).

In the above process, all events act to promote the formation of edges; it is also possible to define events that inhibit them (Butts, 2024).  Let an \emph{inhibition} event be one that, if it occurs, forbids the creation of an \eqn{i \to j}{i->j}.  As with \eqn{d}, we may specify a total probability \eqn{\epsilon_{ij}}{epsilon_{ij}} that such an event occurs exogenously for the \eqn{i,j} edge.  We may also specify endogenous inhibition events.  For instance, consider a \emph{satiation} event, which has the potential to occur every time \eqn{i} emits an edge to some other vertex; each existing edge has a chance of triggering \dQuote{satiation,} in which case the focal edge is inhibited. The associated approximate conditional (i.e., transition probability) with these effects is then

\deqn{
\Pr(i \to j|G_{-ij}) \approx (1-\epsilon_{ij}) (1-\delta)^w\left(1 - (1-\rho)^z (1-\sigma)^y (1-\pi)^x (1-d_{ij})\right)
}{%
Pr(i->j | G_{-ij}) apx (1-epsilon_{ij}) (1-delta)^w [1 - (1-rho)^z (1-sigma)^y (1-pi)^x (1-d_{ij})]
}

where \eqn{w} is the outdegree of \eqn{i} in \eqn{G_{-ij}} and \eqn{\delta}{delta} is the probability of the satiation event.  The net effect of satiation is to suppress edge formation (in roughly geometric fashion) on high degree nodes.  This may be useful in preventing degeneracy when using sigma and rho effects.  Degeneracy can also be reduced by employing the \code{dichotomize.sib.effects} argument, which counts only the first shared partner's contribution towards sibling and double role effects.

It should be noted that the above process is not entirely consistent with the tracing-based model, which is itself not uniformly well-specified in the literature.  For this reason, the local model is referred to here as a Skvoretz-Fararo or Markovian biased net graph process.  One significant advantage of this process is that it is well-defined, and easily simulated: the above equation can be used to form the transition rule for a Markov chain Monte Carlo algorithm, which is used by \eqn{rgbn} to take draws from the (local) biased net model.  (Note that while the underlying Markov chain is only a Gibbs sampler in the special cases for which the putative conditional distributions are jointly satisfiable, it \emph{always} can be interpreted as simulating draws from the equilibrium distribution of a SF/MBN graph process.)  Burn-in and thinning are controlled by the corresponding arguments; since degeneracy is common with models of this type, it is advisable to check for adequate mixing.  An alternative simulation strategy is the exact sampling procedure of Butts (2018), which employs a form of coupling from the past (CFTP).  The CFTP method generates exact, independent draws from the equilibrium distribution of the biased net process (up to numerical limits), but can be slow to attain coalescence (and does not currently support satiation events or other inhibition events).  Setting \code{maxiter} to smaller values limits the search depth employed, at the possible cost of biasing the resulting sample.  An initial condition may be specified for the MCMC using the \code{seed.graph}; if not specified, the empty graph is used.

For some applications (e.g., ABC rejection sampling), it can be useful to terminate simulation if the density is obviously too high for the draw to be useful.  (Compare to similar functionality in the \code{ergm} \dQuote{density guard} feature.)  This can be invoked for the MCMC algorithm by setting the \code{max.density} less than 1.  In this case, the chain is terminated as soon as the threshold density is reached.  The resulting object is marked with an attribute called \code{early.termination} with a value of \code{TRUE}, which should obviously be checked if this feature is used (since the terminated draws are not from the target distribution - especially if \code{n>1}!).  This feature cannot be used with CFTP, and is ignored when CFTP is selected.
}
\value{
  An adjacency array or list of sna edgelists containing the simulated graphs.
}
\references{ 
Butts, C.T.  (2018).  \dQuote{A Perfect Sampling Method for Exponential Family Random Graph Models.} \emph{Journal of Mathematical Sociology}, 42(1), 17-36.

Butts, C.T.  (2024).  \dQuote{A Return to Biased Nets: New Specifications and Approximate Bayesian Inference.}  \emph{Journal of Mathematical Sociology}.

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

#Create a two-group model with homophily
x<-rbinom(30,1,0.5)           #Generate group labels
d<-0.02+outer(x,x,"==")*0.2   #Set base tie probability
g4<-rgbn(1,30,param=list(pi=0.25, sigma=0.02, rho=0, d=d))
gplot(g4, vertex.col=1+x)     #Note the group structure

#Create a two-group model where cross-group ties are inhibited
x<-rbinom(30,1,0.5)           #Generate group labels
ep<-outer(x,x,"!=")*0.75      #Set inhibition probability
g5<-rgbn(1,30,param=list(pi=0.5, sigma=0.05, rho=0, d=0.1, 
    epsilon=ep))
gplot(g5, vertex.col=1+x)     #Note the group structure

}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ graphs }% __ONLY ONE__ keyword per line
