\name{eval.edgeperturbation}
\alias{eval.edgeperturbation}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the Effects of Single-Edge Perturbations on Structural Indices }
\description{
  Evaluates a given function on an input graph with and without a specified edge, returning the difference between the results in each case.
}
\usage{
eval.edgeperturbation(dat, i, j, FUN, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ A single adjacency matrix }
  \item{i}{ The row(s) of the edge(s) to be perturbed }
  \item{j}{ The column(s) of the edge(s) to be perturbed }
  \item{FUN}{ The function to be computed }
  \item{\dots}{ Additional arguments to \code{FUN} }
}
\details{
Although primarily a back-end utility for \code{\link{pstar}}, \code{eval.edgeperturbation} may be useful in any circumstance in which one wishes to assess the stability of a given structural index with respect to single edge perturbations. The function to be evaluated is calculated first on the input graph with all marked edges set to present, and then on the same graph with said edges absent.  (Obviously, this is sensible only for dichotomous data.)  The difference is then returned.  

In \code{\link{pstar}}, calls to \code{eval.edgeperturbation} are used to construct a perturbation effect matrix for the GLM. 
}
\value{
The difference in the values of \code{FUN} as computed on the perturbed graphs.
}
\references{ Anderson, C.; Wasserman, S.; and Crouch, B. (1999).  ``A p* Primer:  Logit Models for Social Networks.  \emph{Social Networks,} 21,37-66. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ \code{length(i)} and \code{length(j)} must be equal; where multiple edges are specified, the row and column listings are interpreted as pairs. }

\seealso{ \code{\link{pstar}} }

\examples{

#Create a random graph
g<-rgraph(5)

#How much does a one-edge change affect reciprocity?
eval.edgeperturbation(g,1,2,grecip)
}
\keyword{ math }%-- one or more ...
\keyword{ graphs } 
