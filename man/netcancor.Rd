\name{netcancor}
\alias{netcancor}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Canonical Correlation for Labeled Graphs }
\description{
 \code{netcancor} finds the canonical correlation(s) between the graph sets \code{x} and \code{y}, testing the result using either conditional uniform graph (CUG) or quadratic assignment procedure (QAP) null hypotheses.
}
\usage{
netcancor(y, x, mode="digraph", diag=FALSE, nullhyp="cugtie", 
    reps=1000)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{ one or more input graphs. }
  \item{x}{ one or more input graphs. }
  \item{mode}{ string indicating the type of graph being evaluated.  "digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  \code{mode} is set to "digraph" by default. }
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{nullhyp}{ string indicating the particular null hypothesis against which to test the observed estimands.  A value of "cug" implies a conditional uniform graph test (see \code{\link{cugtest}}) controlling for order \emph{only}; "cugden" controls for both order and tie probability; "cugtie" controls for order and tie distribution (via bootstrap); and "qap" implies that the QAP null hypothesis (see \code{\link{qaptest}}) should be used. }
  \item{reps}{integer indicating the number of draws to use for quantile estimation.  (Relevant to the null hypothesis test only - the analysis itself is unaffected by this parameter.)  Note that, as for all Monte Carlo procedures, convergence is slower for more extreme quantiles. }
}
\details{
The \code{netcancor} routine is actually a front-end to the \code{\link{cancor}} routine for computing canonical correlations between sets of vectors.  \code{netcancor} itself vectorizes the network variables (as per its graph type) and manages the appropriate null hypothesis tests; the actual canonical correlation is handled by \code{\link{cancor}}.  

Canonical correlation itself is a multivariate generalization of the product-moment correlation.  Specifically, the analysis seeks linear combinations of the variables in \code{y} which are well-explained by linear combinations of the variables in \code{x}.  The network version of this technique is performed elementwise on the adjacency matrices of the graphs in question; as usual, the result should be interpreted with an eye to the relationship between the type of data used and the assumptions of the underlying model.

Intelligent printing and summarizing of netcancor objects is provided by \code{\link{print.netcancor}} and \code{\link{summary.netcancor}}.  
}
\value{
An object of class \code{netcancor} with the following properties:

\item{xdist}{
Array containing the distribution of the X coefficients under the null hypothesis test.
}
\item{ydist}{
Array containing the distribution of the Y coefficients under the null hypothesis test.
}
\item{cdist}{
Array containing the distribution of the canonical correlation coefficients under the null hypothesis test.
}
\item{cor}{
Vector containing the observed canonical correlation coefficients.
}
\item{xcoef}{
Vector containing the observed X coefficients.
}
\item{ycoef}{
Vector containing the observed Y coefficients.
}
\item{cpgreq}{
Vector containing the estimated upper tail quantiles (p>=obs) for the observed canonical correlation coefficients under the null hypothesis.
}
\item{cpleeq}{
Vector containing the estimated lower tail quantiles (p<=obs) for the observed canonical correlation coefficients under the null hypothesis.
}
\item{xpgreq}{
Matrix containing the estimated upper tail quantiles (p>=obs) for the observed X coefficients under the null hypothesis.
}
\item{xpleeq}{
Matrix containing the estimated lower tail quantiles (p<=obs) for the observed X coefficients under the null hypothesis.
}
\item{ypgreq}{
Matrix containing the estimated upper tail quantiles (p>=obs) for the observed Y coefficients under the null hypothesis.
}
\item{ypleeq}{
Matrix containing the estimated lower tail quantiles (p<=obs) for the observed Y coefficients under the null hypothesis.
}
\item{cnames}{
Vector containing names for the canonical correlation coefficients.
}
\item{xnames}{
Vector containing names for the X vars.
}
\item{ynames}{
Vector containing names for the Y vars.
}
\item{xcenter}{
Values used to adjust the X variables.
}
\item{xcenter}{
Values used to adjust the Y variables.
}
\item{nullhyp}{
String indicating the null hypothesis employed.
}
}
\references{ Butts, C.T., and Carley, K.M.  (2001).  ``Multivariate Methods for Interstructural Analysis.''  CASOS working paper, Carnegie Mellon University. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ This will eventually be replaced with a superior cancor procedure with more interpretable output; the new version will handle arbitrary labeling as well. }

% ~Make other sections like WARNING with \section{WARNING }{....} ~
%\section{Requires }{\code{mva}}

\seealso{ \code{\link{gcor}}, \code{\link{cugtest}}, \code{\link{qaptest}}, \code{\link{cancor}} }

\examples{
#Generate a valued seed structure
cv<-matrix(rnorm(100),nrow=10,ncol=10)
#Produce two sets of valued graphs
x<-array(dim=c(3,10,10))
x[1,,]<-3*cv+matrix(rnorm(100,0,0.1),nrow=10,ncol=10)
x[2,,]<--1*cv+matrix(rnorm(100,0,0.1),nrow=10,ncol=10)
x[3,,]<-x[1,,]+2*x[2,,]+5*cv+matrix(rnorm(100,0,0.1),nrow=10,ncol=10)
y<-array(dim=c(2,10,10))
y[1,,]<--5*cv+matrix(rnorm(100,0,0.1),nrow=10,ncol=10)
y[2,,]<--2*cv+matrix(rnorm(100,0,0.1),nrow=10,ncol=10)
#Perform a canonical correlation analysis
nc<-netcancor(y,x,reps=100)
summary(nc)
}
\keyword{ multivariate }%-- one or more ...
\keyword{ math }
\keyword{ graphs }
