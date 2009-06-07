\name{netlm}
\alias{netlm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Linear Regression for Network Data }
\description{
  \code{netlm} regresses the network variable in \code{y} on the network variables in stack \code{x} using ordinary least squares.  The resulting fits (and coefficients) are then tested against the indicated null hypothesis.
}
\usage{
netlm(y, x, intercept=TRUE, mode="digraph", diag=FALSE,
    nullhyp=c("qap", "qapspp", "qapy", "qapx", "qapallx", 
    "cugtie", "cugden", "cuguman", "classical"), 
    test.statistic = c("t-value", "beta"), tol=1e-7,
    reps=1000)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{ dependent network variable.  This should be a matrix, for obvious reasons; NAs are allowed, but dichotomous data is strongly discouraged due to the assumptions of the analysis. }
  \item{x}{ stack of independent network variables.  Note that NAs are permitted, as is dichotomous data. }
  \item{intercept}{ logical; should an intercept term be added? }
  \item{mode}{ string indicating the type of graph being evaluated.  \code{"digraph"} indicates that edges should be interpreted as directed; \code{"graph"} indicates that edges are undirected.  \code{mode} is set to \code{"digraph"} by default. }
  \item{diag}{ logical; should the diagonal be treated as valid data?  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{nullhyp}{ string indicating the particular null hypothesis against which to test the observed estimands. }
  \item{test.statistic}{ string indicating the test statistic to be used for the Monte Carlo procedures.}
  \item{tol}{ tolerance parameter for \code{\link{qr.solve}}. }
  \item{reps}{ integer indicating the number of draws to use for quantile estimation.  (Relevant to the null hypothesis test only - the analysis itself is unaffected by this parameter.)  Note that, as for all Monte Carlo procedures, convergence is slower for more extreme quantiles.  By default, \code{reps}=1000. }
}
\details{
\code{netlm} performs an OLS linear network regression of the graph \code{y} on the graphs in \code{x}.  Network regression using OLS is directly analogous to standard OLS regression elementwise on the appropriately vectorized adjacency matrices of the networks involved.  In particular, the network regression attempts to fit the model:

\deqn{\mathbf{A_y} = b_0 \mathbf{A_1} + b_1 \mathbf{A_{x_1}} + b_2 \mathbf{A_{x_2}} + \dots + \mathbf{Z}}{%
A_y = b_0 A_1 + b_1 A_x1 + b_2 A_x2 + \dots + Z}

where \eqn{\mathbf{A_y}}{A_y} is the dependent adjacency matrix, \eqn{\mathbf{A_{x_i}}}{A_xi} is the ith independent adjacency matrix, \eqn{\mathbf{A_1}}{A_1} is an n x n matrix of 1's, and \eqn{\mathbf{Z}}{Z} is an n x n matrix of independent normal random variables with mean 0 and variance \eqn{\sigma^2}{sigma^2}.  Clearly, this model is nonoptimal when \eqn{\mathbf{A_y}}{A_y} is dichotomous (or, for that matter, categorical in general); an alternative such as \code{netlogit} should be employed in such cases.  (Note that \code{netlm} will still attempt to fit such data...the user should consider him or herself to have been warned.)

Because of the frequent presence of row/column/block autocorrelation in network data, classical hull hypothesis tests (and associated standard errors) are generally suspect.  Further, it is sometimes of interest to compare fitted parameter values to those arising from various baseline models (e.g., uniform random graphs conditional on certain observed statistics).  The tests supported by \code{netlm} are as follows:

\describe{
  \item{\code{classical}}{ tests based on classical asymptotics. }
  \item{\code{cug}}{conditional uniform graph test (see \code{\link{cugtest}}) controlling for order.}
  \item{\code{cugden}}{conditional uniform graph test, controlling for order and density.}
  \item{\code{cugtie}}{conditional uniform graph test, controlling for order and tie distribution.}
  \item{\code{qap}}{QAP permutation test (see \code{\link{qaptest}}); currently identical to \code{qapspp}. }
  \item{\code{qapallx}}{QAP permutation test, using independent x-permutations.}
  \item{\code{qapspp}}{QAP permutation test, using Dekker's "semi-partialling plus" procedure. }
  \item{\code{qapx}}{QAP permutation test, using (single) x-permutations.}
  \item{\code{qapy}}{QAP permutation test, using y-permutations.}
}

The statistic to be employed in the above tests may be selected via \code{test.statistic}.  By default, the \eqn{t}-statistic (rather than estimated coefficient) is used, as this is more approximately pivotal; coefficient-based tests are not recommended for QAP null hypotheses, although they are provided here for legacy purposes.

Note that interpretation of quantiles for single coefficients can be complex in the presence of multicollinearity or third variable effects.  \code{qapspp}  is generally recommended for most multivariable analyses, as it is known to be fairly robust to these conditions.  Reasonable printing and summarizing of \code{netlm} objects is provided by \code{\link{print.netlm}} and \code{\link{summary.netlm}}, respectively.  No plot methods exist at this time, alas. 
}
\value{
  An object of class \code{netlm}
}
\references{
Dekker, D.; Krackhardt, D.; Snijders, T.A.B.  (2007).  \dQuote{Sensitivity of MRQAP Tests to Collinearity and Autocorrelation Conditions.}  \emph{Psychometrika}, 72(4), 563-581.

Dekker, D.; Krackhardt, D.; Snijders, T.A.B.  (2003).  \dQuote{Mulicollinearity Robust QAP for Multiple Regression.}  CASOS Working Paper, Carnegie Mellon University.

Krackhardt, D.  (1987).  \dQuote{QAP Partialling as a Test of Spuriousness.} \emph{Social Networks}, 9 171-186.

Krackhardt, D.  (1988).  \dQuote{Predicting With Networks: Nonparametric Multiple Regression Analyses of Dyadic Data.}  \emph{Social Networks}, 10, 359-382.
}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

 %~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{lm}}, \code{\link{netlogit}} }

\examples{
#Create some input graphs
x<-rgraph(20,4)

#Create a response structure
y<-x[1,,]+4*x[2,,]+2*x[3,,]   #Note that the fourth graph is unrelated

#Fit a netlm model
nl<-netlm(y,x,reps=100)

#Examine the results
summary(nl)
}
\keyword{ regression }%-- one or more ...
\keyword{ math }
\keyword{ graphs }
