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
bbnam.fixed(dat, nprior=matrix(rep(0.5,dim(dat)[2]^2),
    nrow=dim(dat)[2],ncol=dim(dat)[2]), em=0.25, ep=0.25, diag=FALSE,
    mode="digraph", draws=1500, outmode="draws", anames=paste("a",
    1:dim(dat)[2],sep=""), onames=paste("o",1:dim(dat)[1], sep=""))
bbnam.pooled(dat, nprior=matrix(rep(0.5,dim(dat)[2]*dim(dat)[3]),
    nrow=dim(dat)[2],ncol=dim(dat)[3]), emprior=c(1,1), 
    epprior=c(1,1), diag=FALSE, mode="digraph", reps=5, draws=1500, 
    burntime=500, quiet=TRUE, anames=paste("a",1:dim(dat)[2],sep=""),
    onames=paste("o",1:dim(dat)[1],sep=""), compute.sqrtrhat=TRUE)
bbnam.actor(dat, nprior=matrix(rep(0.5,dim(dat)[2]*dim(dat)[3]),
    nrow=dim(dat)[2],ncol=dim(dat)[3]), 
    emprior=cbind(rep(1,dim(dat)[1]),rep(1,dim(dat)[1])), 
    epprior=cbind(rep(1,dim(dat)[1]),rep(1,dim(dat)[1])), diag=FALSE,
    mode="digraph", reps=5, draws=1500, burntime=500, quiet=TRUE, 
    anames=paste("a",1:dim(dat)[2],sep=""), 
    onames=paste("o",1:dim(dat)[1],sep=""), compute.sqrtrhat=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ Data array to be analyzed.  This array must be of dimension n x n x n, where n is |V(G)|, the first dimension indexes the observer, the second indexes the sender of the relation, and the third dimension indexes the recipient of the relation.  (E.g., \code{dat[i,j,k]==1} implies that i observed j sending the relation in question to k.)  Note that only dichotomous data is supported at present, and missing values are permitted; the data collection pattern, however, is assumed to be ignorable, and hence the posterior draws are implicitly conditional on the observation pattern. }
  \item{model}{ String containing the error model to use; options are ``actor,'' ``pooled,'' and ``fixed'' }
  \item{...}{Arguments to be passed by \code{bbnam} to the particular model method}
  \item{nprior}{ Network prior matrix.  This must be a matrix of dimension n x n, containing the arc/edge priors for the criterion network.  (E.g., \code{nprior[i,j]} gives the prior probability of i sending the relation to j in the criterion graph.)  If no network prior is provided, an uninformative prior on the space of networks will be assumed (i.e., p(i->j)=0.5).  Missing values are not allowed. }
  \item{em}{ Probability of a false negative; this may be in the form of a single number, one number per observation slice, one number per (directed) dyad, or one number per dyadic observation (fixed model only) }
  \item{ep}{ Probability of a false positive; this may be in the form of a single number, one number per observation slice, one number per (directed) dyad, or one number per dyadic observation (fixed model only) }
  \item{emprior}{ Parameters for the (beta) false negative prior; these should be in the form of an (alpha,beta) pair for the pooled model, and of an n x 2 matrix of (alpha,beta) pairs for the actor model. If no emprior is given, an uninformative prior (1,1) will be assumed; note that this is usually inappropriate, as described below.  Missing values are not allowed. }
  \item{epprior}{ Parameters for the (beta) false positive prior; these should be in the form of an (alpha,beta) pair for the pooled model, and of an n x 2 matrix of (alpha,beta) pairs for the actor model. If no epprior is given, an uninformative prior (1,1) will be assumed; note that this is usually inappropriate, as described below.  Missing values are not allowed. }
  \item{diag}{ Boolean indicating whether loops (matrix diagonals) should be counted as data }
  \item{mode}{ A string indicating whether the data in question forms a ``graph'' or a ``digraph'' }
  \item{reps}{ Number of replicate chains for the Gibbs sampler (pooled and actor models only) }
  \item{draws}{ Integer indicating the total number of draws to take from the posterior distribution.  Draws are taken evenly from each replication (thus, the number of draws from a given chain is draws/reps), and are randomly reordered to minimize dependence associated with position in the chain. }
  \item{burntime}{ Integer indicating the burn-in time for the Markov Chain.  Each replication is iterated burntime times before taking draws (with these initial iterations being discarded); hence, one should realize that each increment to burntime increases execution time by a quantity proportional to reps. (pooled and actor models only) }
  \item{quiet}{ Boolean indicating whether MCMC diagnostics should be displayed (pooled and actor models only) }
  \item{outmode}{ ``posterior'' indicates that the exact posterior probability matrix for the criterion graph should be returned, otherwise draws from the joint posterior are returned instead (fixed model only) }
  \item{anames}{ A vector of names for the actors (vertices) in the graph }
  \item{onames}{ A vector of names for the observers (possibly the actors themselves) whose reports are contained in the CSS}
  \item{compute.sqrtrhat}{ A boolean indicating whether or not Gelman et al.'s potential scale reduction measure (an MCMC convergence diagnostic) should be computed (pooled and actor models only) }
}
\details{
  The bbnam models a set of network data as reflecting a series of (noisy) observations by a set of participant/observers regarding an uncertain criterion structure.  Each observer is assumed to send false positives (i.e., reporting a tie when none exists in the criterion structure) with probability \eqn{e^+}{e^+}, and false negatives (i.e., reporting that no tie exists when one does in fact exist in the criterion structure) with probability \eqn{e^-}{e^-}.  The criterion network itself is taken to be a Bernoulli (di)graph.  Note that the present model includes three variants:

\enumerate{
\item Fixed error probabilities: Each edge is associated with a known pair of false negative/false positive error probabilities (provided by the researcher).  In this case, the posterior for the criterion graph takes the form of a matrix of Bernoulli parameters, with each edge being independent conditional on the parameter matrix.

\item Pooled error probabilities: One pair of (uncertain) false negative/false positive error probabilities is assumed to hold for all observations.  Here, we assume that the researcher's prior information regarding these parameters can be expressed as a pair of Beta distributions, with the additional assumption of independence in the prior distribution.  Note that error rates and edge probabilities are \emph{not} independent in the joint posterior, but the posterior marginals take the form of Beta mixtures and Bernoulli parameters, respectively. 

\item Per observer (``actor'') error probabilities: One pair of (uncertain) false negative/false positive error probabilities is assumed to hold for each observation slice.  Again, we assume that prior knowledge can be expressed in terms of independent Beta distributions (along with the Bernoulli prior for the criterion graph) and the resulting posterior marginals are Beta mixtures and a Bernoulli graph.  (Again, it should be noted that independence in the priors does \emph{not} imply independence in the joint posterior!)
}

By default, the bbnam routine returns (approximately) independent draws from the joint posterior distribution, each draw yielding one realization of the criterion network and one collection of accuracy parameters (i.e., probabilities of false positives/negatives).  This is accomplished via a Gibbs sampler in the case of the pooled/actor model, and by direct sampling for the fixed probability model. In the special case of the fixed probability model, it is also possible to obtain directly the posterior for the criterion graph (expressed as a matrix of Bernoulli parameters); this can be controlled by the \code{outmode} parameter.

As noted, the taking of posterior draws in the nontrivial case is accomplished via a Markov Chain Monte Carlo method, in particular the Gibbs sampler; the high dimensionality of the problem (\eqn{O(n^2+2n)}{O(n^2+2n)}) tends to preclude more direct approaches.  At present, chain burn-in is determined ex ante on a more or less arbitrary basis by specification of the burntime parameter.  Eventually, a more systematic approach will be utilized.  Note that insufficient burn-in will result in inaccurate posterior sampling, so it's not wise to skimp on burn time where otherwise possible.  Similarly, it is wise to employ more than one Markov Chain (set by reps), since it is possible for trajectories to become ``trapped'' in metastable regions of the state space.  Number of draws per chain being equal, more replications are usually better than few; consult Gelman et al. for details.  A useful measure of chain convergence, Gelman and Rubin's potential scale reduction (\eqn{\sqrt{\hat{R}}}{\sqrt{\hat{R}}}), can be computed using the \code{compute.sqrtrhat} parameter.  The potential scale reduction measure is an ANOVA-like comparison of within-chain versus between-chain variance; it approaches 1 (from above) as the chain converges, and longer burn-in times are strongly recommended for chains with scale reductions in excess of 1.1 or thereabouts.  

Finally, a cautionary concerning prior distributions: it is important that the specified priors actually reflect the prior knowledge of the researcher; otherwise, the posterior will be inadequately informed.  In particular, note that an uninformative prior on the accuracy probabilities implies that it is a priori equally probable that any given actor's observations will be informative or \emph{negatively} informative (i.e., that i observing j sending a tie to k \emph{reduces} p(j->k)).  This is a highly unrealistic assumption, and it will tend to produce posteriors which are bimodal (one mode being related to the ``informative'' solution, the other to the ``negatively informative'' solution).  A more plausible but still fairly diffuse prior would be Beta(3,5), which reduces the prior probability of an actor's being negatively informative to 0.16, and the prior probability of any given actor's being more than 50\% likely to make a particular error (on average) to around 0.22.  (This prior also puts substantial mass near the 0.5 point, which would seem consonant with the BKS studies.)  Butts(1999) discusses a number of issues related to choice of priors for the bbnam, and users should consult this reference if matters are unclear before defaulting to the uninformative solution.
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
\references{ Butts, C.T. (1999). ``Informant (In)Accuracy and Network Estimation: A Bayesian Approach.'' CASOS Working Paper, Carnegie Mellon University.

Gelman, A.; Carlin, J.B.; Stern, H.S.; and Rubin, D.B.  (1995).  \emph{Bayesian Data Analysis.}  London: Chapman and Hall.

Gelman, A., and Rubin, D.B.  (1992).  ``Inference from Iterative Simulation Using Multiple Sequences.''  \emph{Statistical Science,} 7, 457-511.

Krackhardt, D.  (1987).  ``Cognitive Social Structures.'' \emph{Social Networks,} 9, 109-134.
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

