\name{nacf}
\alias{nacf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Sample Network Covariance and Correlation Functions }
\description{
  \code{nacf} computes the sample network covariance/correlation function for a specified variable on a given input network.  Moran's \eqn{I} and Geary's \eqn{C} statistics at multiple orders may be computed as well.
}
\usage{
nacf(net, y, lag.max = NULL, type = c("correlation", "covariance",
    "moran", "geary"), neighborhood.type = c("in", "out", "total"),
    partial.neighborhood = TRUE, mode = "digraph", diag = FALSE,
    thresh = 0, demean = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{ one or more graphs. }
  \item{y}{ a numerical vector, of length equal to the order of \code{net}. }
  \item{lag.max}{ optionally, the maximum geodesic lag at which to compute dependence (defaults to order \code{net}-1). }
  \item{type}{ the type of dependence statistic to be computed. }
  \item{neighborhood.type}{ the type of neighborhood to be employed when assessing dependence (as per \code{\link{neighborhood}}). }
  \item{partial.neighborhood}{ logical; should partial (rather than cumulative) neighborhoods be employed at higher orders? }
  \item{mode}{ \code{"digraph"} for directed graphs, or \code{"graph"} if \code{net} is undirected. }
  \item{diag}{ logical; does the diagonal of \code{net} contain valid data?}
  \item{thresh}{ threshold at which to dichotomize \code{net}. }
  \item{demean}{ logical; demean \code{y} prior to analysis? }
}
\details{
  \code{nacf} computes dependence statistics for the vector \code{y} on network \code{net}, for neighborhoods of various orders.  Specifically, let \eqn{\mathbf{A}_i}{A_i} be the \eqn{i}th order adjacency matrix of \code{net}.  The sample network autocovariance of \eqn{\mathbf{y}}{y} on \eqn{\mathbf{A}_i}{A_i} is then given by
  \deqn{
    \sigma_i = \frac{\mathbf{y}^T \mathbf{A}_i \mathbf{y}}{E},
  }{%
    sigma_i = y^T A_i y / E,
  }
 where \eqn{E=\sum_{(j,k)}A_{ijk}}{E = sum(A_i)}.  Similarly, the sample network autocorrelation in the above case is \eqn{\rho_i=\sigma_i/\sigma_0}{sigma_i/sigma_0}, where \eqn{\sigma_0}{sigma_0} is the variance of \eqn{y}.  Moran's \eqn{I} and Geary's \eqn{C} statistics are defined in the usual fashion as
  \deqn{
  I_i = \frac{N \sum_{j=1}^N \sum_{k=1}^N (y_j-\bar{y}) (y_k-\bar{y}) A_{ijk}}{E \sum_{j=1}^N y_j^2},
  }{%
  I_i = N/E * ((y-ybar)^T A_i (y-ybar))/sum(y^2),
  }
and
  \deqn{
  C_i = \frac{(N-1) \sum_{j=1}^N \sum_{k=1}^N (y_j-y_k)^2 A_{ijk}}{2 E \sum_{j=1}^N (y-\bar{y})^2}
  }{%
  C_i = (N-1)/(2E) * (sum_jk (y_j-y_k)^2 A_ijk) / sum((y-ybar)^2)
  }
  respectively, where \eqn{N} is the order of \eqn{\mathbf{A}_i}{A_i} and \eqn{\bar{y}}{ybar} is the mean of \eqn{\mathbf{y}}{y}.
 
 The adjacency matrix associated with the \eqn{i}th order neighborhood is defined as the identity matrix for order 0, and otherwise depends on the type of neighborhood involved.  For input graph \eqn{G=(V,E)}, let the \emph{base relation}, \eqn{R}, be given by the underlying graph of \eqn{G} (i.e., \eqn{G \cup G^T}{G U G^T}) if total neighborhoods are sought, the transpose of \eqn{G} if incoming neighborhoods are sought, or \eqn{G} otherwise.  The partial neighborhood structure of order \eqn{i>0} on \eqn{R} is then defined to be the digraph on \eqn{V} whose edge set consists of the ordered pairs \eqn{(j,k)} having geodesic distance \eqn{i} in \eqn{R}.  The corresponding cumulative neighborhood is formed by the ordered pairs having geodesic distance less than or equal to \eqn{i} in \eqn{R}.  For purposes of \code{nacf}, these neighborhoods are calculated using \code{\link{neighborhood}}, with the specified parameters (including dichotomization at \code{thresh}).

The return value for \code{nacf} is the selected dependence statistic, calculated for each neighborhood structure from order 0 (the identity) through order \code{lag.max} (or \eqn{N-1}, if \code{lag.max==NULL}).  This vector can be used much like the conventional autocorrelation function, to identify dependencies at various lags.  This may, in turn, suggest a starting point for modeling via routines such as \code{\link{lnam}}.
  
}
\value{
  A vector containing the dependence statistics (ascending from order 0).
}
\references{
  Geary, R.C. (1954). \dQuote{The Contiguity Ratio and Statistical Mapping.} \emph{The Incorporated Statistician,} 5: 115-145.

  Moran, P.A.P.  (1950).  \dQuote{Notes on Continuous Stochastic Phenomena.}  \emph{Biometrika,} 37: 17-23.
}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{geodist}}, \code{\link{gapply}}, \code{\link{neighborhood}}, \code{\link{lnam}}, \code{\link{acf}}}
\examples{
#Create a random graph, and an autocorrelated variable
g<-rgraph(50,tp=4/49)
y<-qr.solve(diag(50)-0.8*g,rnorm(50,0,0.05))

#Examine the network autocorrelation function
nacf(g,y)                             #Partial neighborhoods
nacf(g,y,partial.neighborhood=FALSE)  #Cumulative neighborhoods

#Repeat, using Moran's I on the underlying graph
nacf(g,y,type="moran") 
nacf(g,y,partial.neighborhood=FALSE,type="moran")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
