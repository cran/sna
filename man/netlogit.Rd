\name{netlogit}
\alias{netlogit}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Logistic Regression for Network Data }
\description{
  \code{netlogit} performs a logistic regression of the network variable in \code{y} on the network variables in stack \code{x}.  The resulting fits (and coefficients) are then tested against the indicated null hypothesis.
}
\usage{
netlogit(y, x, mode="digraph", diag=FALSE, nullhyp="cugtie", 
    reps=1000)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{ Dependent network variable.  This should be a matrix, for obvious reasons; \code{NA}s are allowed, and the data should be dichotomous. }
  \item{x}{ Data array containing the stack of independent network variables.  By assumption, the first dimension of the array indexes the graph, with the next two indexing the actors.  Note that \code{NA}s are permitted, as is dichotomous data. }
  \item{mode}{ String indicating the type of graph being evaluated.  "Digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  \code{mode} is set to "digraph" by default. }
  \item{diag}{  Boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{nullhyp}{ String indicating the particular null hypothesis against which to test the observed estimands.  A value of "cug" implies a conditional uniform graph test (see \code{\link{cugtest}}) controlling for order \emph{only}; "cugden" controls for both order and tie probability; "cugtie" controls for order and tie distribution (via bootstrap); and "qap" implies that the QAP null hypothesis (see \code{\link{qaptest}}) should be used. }
  \item{reps}{ Integer indicating the number of draws to use for quantile estimation.  (Relevant to the null hypothesis test only -- the analysis itself is unaffected by this parameter.)  Note that, as for all Monte Carlo procedures, convergence is slower for more extreme quantiles.  By default, \code{reps}=1000. }
}
\details{
  \code{netlogit} is primarily a front-end to the built-in \code{\link{glm}} routine.  \code{netlogit} handles vectorization, sets up \code{\link{glm}} options, and deals with null hypothesis testing; the actual fitting is taken care of by \code{\link{glm}}.  

Logistic network regression using is directly analogous to standard logistic regression elementwise on the appropriately vectorized adjacency matrices of the networks involved.  As such, it is often a more appropriate model for fitting dichotomous response networks than is linear network regression.  

Null hypothesis tests for logistic network regression are handled using either the conditional uniform graph hypothesis (the default) or QAP.  See the help pages for these tests for a fuller description of each.  Reasonable printing and summarizing of \code{netlogit} objects is provided by \code{\link{print.netlogit}} and \code{\link{summary.netlogit}}, respectively.  No plot methods exist at this time.
}
\value{
  An object of class \code{netlogit}
}
\references{ \references{Butts, C.T., and Carley, K.M.  (2001).  ``Multivariate Methods for Interstructural Analysis.''  CASOS working paper, Carnegie Mellon University.}
\author{ Carter T. Butts \email{ctb@andrew.cmu.edu} }
%\note{}

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{glm}}, \code{\link{netlm}} }

\examples{
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
\keyword{ regression }%-- one or more ...
\keyword{ network }