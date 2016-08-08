\name{netlogit}
\alias{netlogit}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Logistic Regression for Network Data }
\description{
  \code{netlogit} performs a logistic regression of the network variable in \code{y} on the network variables in set \code{x}.  The resulting fits (and coefficients) are then tested against the indicated null hypothesis.
}
\usage{
netlogit(y, x, intercept=TRUE, mode="digraph", diag=FALSE,
    nullhyp=c("qap", "qapspp", "qapy", "qapx", "qapallx", 
    "cugtie", "cugden", "cuguman", "classical"), test.statistic = 
    c("z-value","beta"), tol=1e-7, reps=1000)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{ dependent network variable.  \code{NA}s are allowed, and the data should be dichotomous. }
  \item{x}{ the stack of independent network variables.  Note that \code{NA}s are permitted, as is dichotomous data. }
  \item{intercept}{logical; should an intercept term be fitted?}
  \item{mode}{ string indicating the type of graph being evaluated.  \code{"digraph"} indicates that edges should be interpreted as directed; \code{"graph"} indicates that edges are undirected.  \code{mode} is set to \code{"digraph"} by default. }
  \item{diag}{  boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{nullhyp}{ string indicating the particular null hypothesis against which to test the observed estimands. }
  \item{test.statistic}{string indicating the test statistic to be used for the Monte Carlo procedures.}
  \item{tol}{ tolerance parameter for \code{qr.solve}. }
  \item{reps}{ integer indicating the number of draws to use for quantile estimation.  (Relevant to the null hypothesis test only -- the analysis itself is unaffected by this parameter.)  Note that, as for all Monte Carlo procedures, convergence is slower for more extreme quantiles.  By default, \code{reps}=1000. }
}
\details{
  \code{netlogit} is primarily a front-end to the built-in \code{\link{glm.fit}} routine.  \code{netlogit} handles vectorization, sets up \code{\link{glm}} options, and deals with null hypothesis testing; the actual fitting is taken care of by \code{\link{glm.fit}}.  

Logistic network regression using is directly analogous to standard logistic regression elementwise on the appropriately vectorized adjacency matrices of the networks involved.  As such, it is often a more appropriate model for fitting dichotomous response networks than is linear network regression.  

Because of the frequent presence of row/column/block autocorrelation in network data, classical hull hypothesis tests (and associated standard errors) are generally suspect.  Further, it is sometimes of interest to compare fitted parameter values to those arising from various baseline models (e.g., uniform random graphs conditional on certain observed statistics).  The tests supported by \code{netlogit} are as follows:

\describe{
  \item{\code{classical}}{ tests based on classical asymptotics. }
  \item{\code{cug}}{conditional uniform graph test (see \code{\link{cugtest}}) controlling for order.}
  \item{\code{cugden}}{conditional uniform graph test, controlling for order and density.}
  \item{\code{cugtie}}{conditional uniform graph test, controlling for order and tie distribution.}
  \item{\code{qap}}{QAP permutation test (see \code{\link{qaptest}}); currently identical to \code{qapspp}. }
  \item{\code{qapallx}}{QAP permutation test, using independent x-permutations.}
  \item{\code{qapspp}}{QAP permutation test, using Dekker's \dQuote{semi-partialling plus} procedure. }
  \item{\code{qapx}}{QAP permutation test, using (single) x-permutations.}
  \item{\code{qapy}}{QAP permutation test, using y-permutations.}
}

Note that interpretation of quantiles for single coefficients can be complex in the presence of multicollinearity or third variable effects.  Although \code{qapspp} is known to be robust to these conditions in the OLS case, there are no equivalent results for logistic regression.  Caution is thus advised.  

The statistic to be employed in the above tests may be selected via \code{test.statistic}.  By default, the z-statistic (rather than estimated coefficient) is used, as this is more approximately pivotal; coefficient-based tests are not recommended for QAP null hypotheses, although they are provided here for legacy purposes.

Reasonable printing and summarizing of \code{netlogit} objects is provided by \code{\link{print.netlogit}} and \code{\link{summary.netlogit}}, respectively.  No plot methods exist at this time.


}
\value{
  An object of class \code{netlogit}
}
\references{Butts, C.T., and Carley, K.M.  (2001).  ``Multivariate Methods for Interstructural Analysis.''  CASOS working paper, Carnegie Mellon University.}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{}

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{glm}}, \code{\link{netlm}} }

\examples{
\dontrun{
#Create some input graphs
x<-rgraph(20,4)

#Create a response structure
y.l<-x[1,,]+4*x[2,,]+2*x[3,,]   #Note that the fourth graph is 
                                #unrelated
y.p<-apply(y.l,c(1,2),function(a){1/(1+exp(-a))})
y<-rgraph(20,tprob=y.p)

#Fit a netlogit model
nl<-netlogit(y,x,reps=100)

#Examine the results
summary(nl)
}
}
\keyword{ regression }%-- one or more ...
\keyword{ math }
\keyword{ graphs }
