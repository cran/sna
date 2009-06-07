\name{bbnam}
\alias{bbnam}
\alias{bbnam.actor}
\alias{bbnam.pooled}
\alias{bbnam.fixed}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Butts' (Hierarchical) Bayesian Network Accuracy Model }
\description{
 Takes posterior draws from Butts' bayesian network accuracy/estimation model for multiple participant/observers (conditional on observed data and priors), using a Gibbs sampler.
}
\usage{
bbnam(dat, model="actor", ...)
bbnam.fixed(dat, nprior=0.5, em=0.25, ep=0.25, diag=FALSE,
    mode="digraph", draws=1500, outmode="draws", anames=NULL,
    onames=NULL)
bbnam.pooled(dat, nprior=0.5, emprior=c(1,11), epprior=c(1,11),
    diag=FALSE, mode="digraph", reps=5, draws=1500, burntime=500, 
    quiet=TRUE, anames=NULL, onames=NULL, compute.sqrtrhat=TRUE)
bbnam.actor(dat, nprior=0.5, emprior=c(1,11), epprior=c(1,11), 
    diag=FALSE, mode="digraph", reps=5, draws=1500, burntime=500,
    quiet=TRUE, anames=NULL, onames=NULL, compute.sqrtrhat=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ Input networks to be analyzed.  This may be supplied in any reasonable form, but must be reducible to an array of dimension \eqn{m \times n \times n}{m x n x n}, where \eqn{n} is \eqn{|V(G)|}, the first dimension indexes the observer (or information source), the second indexes the sender of the relation, and the third dimension indexes the recipient of the relation.  (E.g., \code{dat[i,j,k]==1} implies that i observed j sending the relation in question to k.)  Note that only dichotomous data is supported at present, and missing values are permitted; the data collection pattern, however, is assumed to be ignorable, and hence the posterior draws are implicitly conditional on the observation pattern. }
  \item{model}{ String containing the error model to use; options are \code{"actor"}, \code{"pooled"}, and \code{"fixed"}. }
  \item{...}{Arguments to be passed by \code{bbnam} to the particular model method.}
  \item{nprior}{ Network prior matrix.  This must be a matrix of dimension \eqn{n} x \eqn{n}, containing the arc/edge priors for the criterion network.  (E.g., \code{nprior[i,j]} gives the prior probability of \code{i} sending the relation to \code{j} in the criterion graph.)  Non-matrix values will be coerced/expanded to matrix form as appropriate.  If no network prior is provided, an uninformative prior on the space of networks will be assumed (i.e., \eqn{\Pr(i\to j)=0.5}{Pr(i->j)=0.5}).  Missing values are not allowed. }
  \item{em}{ Probability of a false negative; this may be in the form of a single number, one number per observation slice, one number per (directed) dyad, or one number per dyadic observation (fixed model only). }
  \item{ep}{ Probability of a false positive; this may be in the form of a single number, one number per observation slice, one number per (directed) dyad, or one number per dyadic observation (fixed model only). }
  \item{emprior}{ Parameters for the (Beta) false negative prior; these should be in the form of an \eqn{(\alpha,\beta)}{(alpha,beta)} pair for the pooled model, and of an \eqn{n \times 2}{n x 2} matrix of \eqn{(\alpha,\beta)}{(alpha,beta)} pairs for the actor model (or something which can be coerced to this form). If no \code{emprior} is given, a weakly informative prior (1,11) will be assumed; note that this may be inappropriate, as described below.  Missing values are not allowed. }
  \item{epprior}{  Parameters for the (Beta) false positive prior; these should be in the form of an \eqn{(\alpha,\beta)}{(alpha,beta)} pair for the pooled model, and of an \eqn{n \times 2}{n x 2} matrix of \eqn{(\alpha,\beta)}{(alpha,beta)} pairs for the actor model (or something which can be coerced to this form). If no \code{epprior} is given, a weakly informative prior (1,11) will be assumed; note that this may be inappropriate, as described below.  Missing values are not allowed. }
  \item{diag}{ Boolean indicating whether loops (matrix diagonals) should be counted as data. }
  \item{mode}{ A string indicating whether the data in question forms a \code{"graph"} or a \code{"digraph"} }
  \item{reps}{ Number of replicate chains for the Gibbs sampler (pooled and actor models only). }
  \item{draws}{ Integer indicating the total number of draws to take from the posterior distribution.  Draws are taken evenly from each replication (thus, the number of draws from a given chain is draws/reps). }
  \item{burntime}{ Integer indicating the burn-in time for the Markov Chain.  Each replication is iterated burntime times before taking draws (with these initial iterations being discarded); hence, one should realize that each increment to burntime increases execution time by a quantity proportional to reps. (pooled and actor models only) }
  \item{quiet}{ Boolean indicating whether MCMC diagnostics should be displayed (pooled and actor models only). }
  \item{outmode}{ \code{posterior} indicates that the exact posterior probability matrix for the criterion graph should be returned; otherwise draws from the joint posterior are returned instead (fixed model only). }
  \item{anames}{ A vector of names for the actors (vertices) in the graph. }
  \item{onames}{ A vector of names for the observers (possibly the actors themselves) whose reports are contained in the input data.}
  \item{compute.sqrtrhat}{ A boolean indicating whether or not Gelman et al.'s potential scale reduction measure (an MCMC convergence diagnostic) should be computed (pooled and actor models only). }
}
\details{
  The bbnam models a set of network data as reflecting a series of (noisy) observations by a set of participants/observers regarding an uncertain criterion structure.  Each observer is assumed to send false positives (i.e., reporting a tie when none exists in the criterion structure) with probability \eqn{e^+}{e^+}, and false negatives (i.e., reporting that no tie exists when one does in fact exist in the criterion structure) with probability \eqn{e^-}{e^-}.  The criterion network itself is taken to be a Bernoulli (di)graph.  Note that the present model includes three variants:

\enumerate{
\item Fixed error probabilities: Each edge is associated with a known pair of false negative/false positive error probabilities (provided by the researcher).  In this case, the posterior for the criterion graph takes the form of a matrix of Bernoulli parameters, with each edge being independent conditional on the parameter matrix.

\item Pooled error probabilities: One pair of (uncertain) false negative/false positive error probabilities is assumed to hold for all observations.  Here, we assume that the researcher's prior information regarding these parameters can be expressed as a pair of Beta distributions, with the additional assumption of independence in the prior distribution.  Note that error rates and edge probabilities are \emph{not} independent in the joint posterior, but the posterior marginals take the form of Beta mixtures and Bernoulli parameters, respectively. 

\item Per observer (``actor'') error probabilities: One pair of (uncertain) false negative/false positive error probabilities is assumed to hold for each observation slice.  Again, we assume that prior knowledge can be expressed in terms of independent Beta distributions (along with the Bernoulli prior for the criterion graph) and the resulting posterior marginals are Beta mixtures and a Bernoulli graph.  (Again, it should be noted that independence in the priors does \emph{not} imply independence in the joint posterior!)
}

By default, the \code{bbnam} routine returns (approximately) independent draws from the joint posterior distribution, each draw yielding one realization of the criterion network and one collection of accuracy parameters (i.e., probabilities of false positives/negatives).  This is accomplished via a Gibbs sampler in the case of the pooled/actor model, and by direct sampling for the fixed probability model. In the special case of the fixed probability model, it is also possible to obtain directly the posterior for the criterion graph (expressed as a matrix of Bernoulli parameters); this can be controlled by the \code{outmode} parameter.

As noted, the taking of posterior draws in the nontrivial case is accomplished via a Markov Chain Monte Carlo method, in particular the Gibbs sampler; the high dimensionality of the problem (\eqn{O(n^2+2n)}{O(n^2+2n)}) tends to preclude more direct approaches.  At present, chain burn-in is determined ex ante on a more or less arbitrary basis by specification of the \code{burntime} parameter.  Eventually, a more systematic approach will be utilized.  Note that insufficient burn-in will result in inaccurate posterior sampling, so it's not wise to skimp on burn time where otherwise possible.  Similarly, it is wise to employ more than one Markov Chain (set by \code{reps}), since it is possible for trajectories to become \dQuote{trapped} in metastable regions of the state space.  Number of draws per chain being equal, more replications are usually better than few; consult Gelman et al. for details.  A useful measure of chain convergence, Gelman and Rubin's potential scale reduction (\eqn{\sqrt{\hat{R}}}{\sqrt{\hat{R}}}), can be computed using the \code{compute.sqrtrhat} parameter.  The potential scale reduction measure is an ANOVA-like comparison of within-chain versus between-chain variance; it approaches 1 (from above) as the chain converges, and longer burn-in times are strongly recommended for chains with scale reductions in excess of 1.2 or thereabouts.  

Finally, a cautionary concerning prior distributions: it is important that the specified priors actually reflect the prior knowledge of the researcher; otherwise, the posterior will be inadequately informed.  In particular, note that an uninformative prior on the accuracy probabilities implies that it is a priori equally probable that any given actor's observations will be informative or \emph{negatively} informative (i.e., that \eqn{i} observing \eqn{j} sending a tie to \eqn{k} \emph{reduces} \eqn{\Pr(j\to k)}{Pr(j->k)}).  This is a highly unrealistic assumption, and it will tend to produce posteriors which are bimodal (one mode being related to the \dQuote{informative} solution, the other to the \dQuote{negatively informative} solution).  Currently, the default error parameter prior is Beta(1,11), which is both diffuse and which renders negatively informative observers extremely improbable (i.e., on the order of 1e-6).  Another plausible but still fairly diffuse prior would be Beta(3,5), which reduces the prior probability of an actor's being negatively informative to 0.16, and the prior probability of any given actor's being more than 50\% likely to make a particular error (on average) to around 0.22.  (This prior also puts substantial mass near the 0.5 point, which would seem consonant with the BKS studies.)  For network priors, a reasonable starting point can often be derived by considering the expected mean degree of the criterion graph: if \eqn{d} represents the user's prior expectation for the mean degree, then \eqn{d/(N-1)} is a natural starting point for the cell values of \code{nprior}.  Butts (2003) discusses a number of issues related to choice of priors for the \code{bbnam} model, and users should consult this reference if matters are unclear before defaulting to the uninformative solution.
}

\value{
An object of class bbnam, containing the posterior draws.  The components of the output are as follows:

\item{anames}{
A vector of actor names.
}
\item{draws}{
An integer containing the number of draws.
}
\item{em}{
A matrix containing the posterior draws for probability of producing false negatives, by actor.
}
\item{ep}{
A matrix containing the posterior draws for probability of producing false positives, by actor.
}
\item{nactors}{
An integer containing the number of actors.
}
\item{net}{
An array containing the posterior draws for the criterion network.
}
\item{reps}{
An integer indicating the number of replicate chains used by the Gibbs sampler.
}
}
\references{ 
Butts, C. T.  (2003).  \dQuote{Network Inference, Error, and Informant (In)Accuracy: A Bayesian Approach.}  \emph{Social Networks}, 25(2), 103-140.

Gelman, A.; Carlin, J.B.; Stern, H.S.; and Rubin, D.B.  (1995).  \emph{Bayesian Data Analysis.}  London: Chapman and Hall.

Gelman, A., and Rubin, D.B.  (1992).  \dQuote{Inference from Iterative Simulation Using Multiple Sequences.}  \emph{Statistical Science,} 7, 457-511.

Krackhardt, D.  (1987).  \dQuote{Cognitive Social Structures.} \emph{Social Networks,} 9, 109-134.
}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ As indicated, the posterior draws are conditional on the observed data, and hence on the data collection mechanism if the collection design is non-ignorable.  Complete data (e.g., a CSS) and random tie samples are examples of ignorable designs; see Gelman et al. for more information concerning ignorability.}
 
\seealso{\code{\link{npostpred}}, \code{\link{event2dichot}}, \code{\link{bbnam.bf}}}

\examples{
#Create some random data
g<-rgraph(5)
g.p<-0.8*g+0.2*(1-g)
dat<-rgraph(5,5,tprob=g.p)

#Define a network prior
pnet<-matrix(ncol=5,nrow=5)
pnet[,]<-0.5
#Define em and ep priors
pem<-matrix(nrow=5,ncol=2)
pem[,1]<-3
pem[,2]<-5
pep<-matrix(nrow=5,ncol=2)
pep[,1]<-3
pep[,2]<-5

#Draw from the posterior
b<-bbnam(dat,model="actor",nprior=pnet,emprior=pem,epprior=pep,
    burntime=100,draws=100)
#Print a summary of the posterior draws
summary(b)
}

\keyword{ models}
\keyword{ math }

