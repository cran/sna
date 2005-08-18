\name{gcov}
\alias{gcov}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Find the Covariance(s) Between Two or More Labeled Graphs}
\description{
   \code{gcov} finds the covariances between the adjacency matrices of graphs indicated by \code{g1} and \code{g2} in stack \code{dat} (or possibly \code{dat2}).  Missing values are permitted.
}
\usage{
gcov(dat, dat2=NULL, g1=NULL, g2=NULL, diag=FALSE, mode="digraph")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{dat2}{ optionally, a second graph stack. }
  \item{g1}{ the indices of \code{dat} reflecting the first set of graphs to be compared; by default, all members of \code{dat} are included. }
  \item{g2}{ the indices or \code{dat} (or \code{dat2}, if applicable) reflecting the second set of graphs to be compared; by default, all members of \code{dat} are included. }
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{mode}{ string indicating the type of graph being evaluated.  "digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  \code{mode} is set to "digraph" by default. }
}
\details{
  The graph covariance between two labeled graphs is defined as
\deqn{cov(G,H) = \frac{1}{{|V| \choose 2}} \sum_{\{i,j\}} \left(A^G_{ij}-\mu_G\right)\left(A^H_{ij}-\mu_H\right)}{%
cov(G,H) = sum( (A^G_ij-mu_G)(A^H_ij-mu_H), \{i,j\} )/Choose(|V|,2)}
(with \eqn{A^G}{A^G} being the adjacency matrix of G).  The graph correlation/covariance is at the center of a number of graph comparison methods, including network variants of regression analysis, PCA, CCA, and the like.

Note that \code{gcov} computes only the covariance between \emph{uniquely labeled} graphs.  For the more general case, \code{\link{gscov}} is recommended.
}
\value{
  A graph covariance matrix
}
\references{ Butts, C.T., and Carley, K.M.  (2001).  ``Multivariate Methods for Interstructural Analysis.''  CASOS Working Paper, Carnegie Mellon University. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ The \code{gcov} routine is really just a front-end to the standard \code{\link{cov}} method; the primary value-added is the transparent vectorization of the input graphs (with intelligent handling of simple versus directed graphs, diagonals, etc.).  Classical null hypothesis testing procedures are not recommended for use with graph covariance; for nonparametric null hypothesis testing regarding graph covariance, see \code{\link{cugtest}} and \code{\link{qaptest}}.
}

 %~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{gscov}}, \code{\link{gcor}}, \code{\link{gscor}} }

\examples{
#Generate two random graphs each of low, medium, and high density
g<-rgraph(10,6,tprob=c(0.2,0.2,0.5,0.5,0.8,0.8))

#Examine the covariance matrix
gcov(g)
}
\keyword{ univar }
\keyword{ multivariate }%-- one or more ...
\keyword{ graphs }
