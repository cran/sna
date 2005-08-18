\name{gcor}
\alias{gcor}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ 
Find the (Product-Moment) Correlation Between Two or More Labeled Graphs }
\description{
 \code{gcor} finds the product-moment correlation between the adjacency matrices of graphs indicated by \code{g1} and \code{g2} in stack \code{dat} (or possibly \code{dat2}).  Missing values are permitted. 
}
\usage{
gcor(dat, dat2=NULL, g1=NULL, g2=NULL, diag=FALSE, mode="digraph")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{dat2}{ optionally, a second stack of graphs. }
  \item{g1}{ the indices of \code{dat} reflecting the first set of graphs to be compared; by default, all members of \code{dat} are included. }
  \item{g2}{ the indices or \code{dat} (or \code{dat2}, if applicable) reflecting the second set of graphs to be compared; by default, all members of \code{dat} are included. }
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{diag} is \code{FALSE} by default. }
  \item{mode}{ string indicating the type of graph being evaluated.  "Digraph" indicates that edges should be interpreted as directed; "graph" indicates that edges are undirected.  \code{mode} is set to "digraph" by default. }
}
\details{

The (product moment) graph correlation between labeled graphs G and H is given by 
\deqn{cor(G,H) = \frac{cov(G,H)}{\sqrt{cov(G,G) cov(H,H)}} }{%
cor(G,H) = cov(G,V)/sqrt(cov(G,G)cov(H,H))}
where the graph covariance is defined as
\deqn{cov(G,H) = \frac{1}{{|V| \choose 2}} \sum_{\{i,j\}} \left(A^G_{ij}-\mu_G\right)\left(A^H_{ij}-\mu_H\right)}{%
cov(G,H) = sum( (A^G_ij-mu_G)(A^H_ij-mu_H), \{i,j\} )/Choose(|V|,2)}
(with \eqn{A^G}{A^G} being the adjacency matrix of G).  The graph correlation/covariance is at the center of a number of graph comparison methods, including network variants of regression analysis, PCA, CCA, and the like.

Note that \code{gcor} computes only the correlation between \emph{uniquely labeled} graphs.  For the more general case, \code{\link{gscor}} is recommended.
}
\value{
A graph correlation matrix
}
\references{ Butts, C.T., and Carley, K.M.  (2001).  ``Multivariate Methods for Interstructural Analysis.''  CASOS Working Paper, Carnegie Mellon University.

Krackhardt, D.  (1987).  ``QAP Partialling as a Test of Spuriousness.''  \emph{Social Networks}, 9, 171-86
}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{  The \code{gcor} routine is really just a front-end to the standard \code{\link{cor}} method; the primary value-added is the transparent vectorization of the input graphs (with intelligent handling of simple versus directed graphs, diagonals, etc.).  As noted, the correlation coefficient returned is a standard Pearson's product-moment coefficient, and output should be interpreted accordingly.  Classical null hypothesis testing procedures are not recommended for use with graph correlations; for nonparametric null hypothesis testing regarding graph correlations, see \code{\link{cugtest}} and \code{\link{qaptest}}.  For multivariate correlations among graph sets, try \code{\link{netcancor}}. 
}

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{gscor}}, \code{\link{gcov}}, \code{\link{gscov}}  }

\examples{
#Generate two random graphs each of low, medium, and high density
g<-rgraph(10,6,tprob=c(0.2,0.2,0.5,0.5,0.8,0.8))

#Examine the correlation matrix
gcor(g)
}
\keyword{ univar }
\keyword{ multivariate }%-- one or more ...
\keyword{ graphs }
