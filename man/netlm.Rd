\name{netlm}
\alias{netlm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Linear Regression for Network Data }
\description{
  \code{netlm} regresses the network variable in \code{y} on the network variables in stack \code{x} using ordinary least squares.  The resulting fits (and coefficients) are then tested against the indicated null hypothesis.
}
\usage{
netlm(y, x, mode="digraph", diag=FALSE, nullhyp="cugtie", reps=1000)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{ Dependent network variable.  This should be a matrix, for obvious reasons; NAs are allowed, but dichotomous data is strongly discouraged due to the assumptions of the analysis. }
  \item{x}{ Data array containing the stack of independent network variables.  By assumption, the first dimension of the array indexes the graph, with the next two indexing the actors.  Note that NAs are permitted, as is dichotomous data. }
  \item{mode}{ String indicating the type of graph being evaluated.  "digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  \code{mode} is set to "digraph" by default. }
  \item{diag}{ Boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{Diag} is \code{FALSE} by default. }
  \item{nullhyp}{ String indicating the particular null hypothesis against which to test the observed estimands.  A value of "cug" implies a conditional uniform graph test (see \code{\link{cugtest}}) controlling for order \emph{only}; "cugden" controls for both order and tie probability; "cugtie" controls for order and tie distribution (via bootstrap); and "qap" implies that the QAP null hypothesis (see \code{\link{qaptest}}) should be used. }
  \item{reps}{ Integer indicating the number of draws to use for quantile estimation.  (Relevant to the null hypothesis test only - the analysis itself is unaffected by this parameter.)  Note that, as for all Monte Carlo procedures, convergence is slower for more extreme quantiles.  By default, \code{reps}=1000. }
}
\details{
 Unsurprisingly, \code{netlm} is really a front-end to the built-in \code{\link{lm}} routine.  \code{netlm} handles vectorization and null hypothesis testing; the actual fitting is taken care of by \code{\link{lm}}.  

Network regression using OLS is directly analogous to standard OLS regression elementwise on the appropriately vectorized adjacency matrices of the networks involved.  In particular, the network regression attempts to fit the model:

\deqn{\mathbf{A_y} = b_0 \mathbf{A_1} + b_1 \mathbf{A_{x_1}} + b_2 \mathbf{A_{x_2}} + \dots + \mathbf{Z}}{%
A_y = b_0 A_1 + b_1 A_x1 + b_2 A_x2 + \dots + Z}

where \eqn{\mathbf{A_y}}{A_y} is the dependent adjacency matrix, \eqn{\mathbf{A_{x_i}}}{A_xi} is the ith independent adjacency matrix, \eqn{\mathbf{A_1}}{A_1} is an n x n matrix of 1's, and \eqn{\mathbf{Z}}{Z} is an n x n matrix of independent normal random variables with mean 0 and variance \eqn{\sigma^2}{sigma^2}.  Clearly, this model is nonoptimal when \eqn{\mathbf{A_y}}{A_y} is dichotomous (or, for that matter, categorical in general); an alternative such as \code{netlogit} should be employed in such cases.  (Note that \code{netlm} will still attempt to fit such data...the user should consider him or herself to have been warned.)

Null hypothesis tests for the network regression are handled using either the conditional uniform graph hypothesis (the default) or QAP.  See the help pages for these tests for a fuller description of each.  Reasonable printing and summarizing of \code{netlm} objects is provided by \code{\link{print.netlm}} and \code{\link{summary.netlm}}, respectively.  No plot methods exist at this time, alas. 
}
\value{
  An object of class \code{netlm}
}
\references{Butts, C.T., and Carley, K.M.  (2001).  ``Multivariate Methods for Interstructural Analysis.''  CASOS working paper, Carnegie Mellon University. 

Krackhardt, D.  (1987).  ``QAP Partialling as a Test of Spuriousness.'' \emph{Social Networks}, 9 171-186.

Krackhardt, D.  (1988).  ``Predicting With Networks: Nonparametric Multiple Regression Analyses of Dyadic Data.''  \emph{Social Networks}, 10, 359-382.}
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
