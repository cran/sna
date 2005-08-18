\name{bonpow}
\alias{bonpow}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Find Bonacich Power Centrality Scores of Network Positions }
\description{
  \code{bonpow} takes one or more graphs (\code{dat}) and returns the Boncich power centralities of positions (selected by \code{nodes}) within the graphs indicated by \code{g}.  The decay rate for power contributions is specified by \code{exponent} (1 by default).  This function is compatible with \code{\link{centralization}}, and will return the theoretical maximum absolute deviation (from maximum) conditional on size (which is used by \code{\link{centralization}} to normalize the observed centralization score).
}
\usage{
bonpow(dat, g=1, nodes=NULL, gmode="digraph", diag=FALSE,
    tmaxdev=FALSE, exponent=1, rescale=FALSE, tol=1e-07)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{g}{ integer indicating the index of the graph for which centralities are to be calculated (or a vector thereof).  By default, \code{g}=1. }
  \item{nodes}{ vector indicating which nodes are to be included in the calculation.  By default, all nodes are included. }
  \item{gmode}{ string indicating the type of graph being evaluated.  \code{"digraph"} indicates that edges should be interpreted as directed; \code{"graph"} indicates that edges are undirected.  This is currently ignored. }
  \item{diag}{ boolean indicating whether or not the diagonal should be treated as valid data.  Set this true if and only if the data can contain loops.  \code{Diag} is \code{FALSE} by default. }
  \item{tmaxdev}{ boolean indicating whether or not the theoretical maximum absolute deviation from the maximum nodal centrality should be returned.  By default, \code{tmaxdev}=\code{FALSE}. }
  \item{exponent}{ exponent (decay rate) for the Bonacich power centrality score; can be negative }
  \item{rescale}{ if true, centrality scores are rescaled such that they sum to 1. }
  \item{tol}{ tolerance for near-singularities during matrix inversion (see \code{\link{solve}}) }
}
\details{
Bonacich's power centrality measure is defined by \eqn{C_{BP}\left(\alpha,\beta\right)=\alpha\left(\mathbf{I}-\beta\mathbf{A}\right)^{-1}\mathbf{A}\mathbf{1}}{C_BP(alpha,beta)=alpha (I-A)^-1 A 1}, where \eqn{\beta}{beta} is an attenuation parameter (set here by \code{exponent}) and \eqn{\mathbf{A}}{A} is the graph adjacency matrix.  (The coefficient \eqn{\alpha}{alpha} acts as a scaling parameter, and is set here (following Bonacich (1987)) such that the sum of squared scores is equal to the number of vertices.  This allows 1 to be used as a reference value for the ``middle'' of the centrality range.)  When \eqn{\beta \rightarrow 1/\lambda_{\mathbf{A}1}}{beta->1/lambda_A1} (the reciprocal of the largest eigenvalue of \eqn{\mathbf{A}}{A}), this is to within a constant multiple of the familiar eigenvector centrality score; for other values of \eqn{\beta}, the behavior of the measure is quite different.  In particular, \eqn{\beta} gives positive and negative weight to even and odd walks, respectively, as can be seen from the series expansion \eqn{C_{BP}\left(\alpha,\beta\right)=\alpha \sum_{k=0}^\infty \beta^k \mathbf{A}^{k+1} \mathbf{1}}{C_BP(alpha,beta) = alpha sum( beta^k A^(k+1) 1, k in 0..infinity )} which converges so long as \eqn{|\beta| < 1/\lambda_{\mathbf{A}1}}{|beta|<1/lambda_A1}.  The magnitude of \eqn{\beta}{beta} controls the influence of distant actors on ego's centrality score, with larger magnitudes indicating slower rates of decay.  (High rates, hence, imply a greater sensitivity to edge effects.)

Interpretively, the Bonacich power measure corresponds to the notion that the power of a vertex is recursively defined by the sum of the power of its alters.  The nature of the recursion involved is then controlled by the power exponent: positive values imply that vertices become more powerful as their alters become more powerful (as occurs in cooperative relations), while negative values imply that vertices become more powerful only as their alters become \emph{weaker} (as occurs in competitive or antagonistic relations).  The magnitude of the exponent indicates the tendency of the effect to decay across long walks; higher magnitudes imply slower decay.  One interesting feature of this measure is its relative instability to changes in exponent magnitude (particularly in the negative case).  If your theory motivates use of this measure, you should be very careful to choose a decay parameter on a non-ad hoc basis.
}
\value{
A vector, matrix, or list containing the centrality scores (depending on the number and size of the input graphs).
}
\references{ Bonacich, P.  (1972).  ``Factoring and Weighting Approaches to Status Scores and Clique Identification.''  \emph{Journal of Mathematical Sociology}, 2, 113-120.

Bonacich, P.  (1987).  ``Power and Centrality: A Family of Measures.'' \emph{American Journal of Sociology}, 92, 1170-1182.
 }
\author{ Carter T. Butts \email{buttsc@uci.edu}}
\note{ The theoretical maximum deviation used here is not obtained with the star network, in general.  For positive exponents, at least, the symmetric maximum occurs for an empty graph with one complete dyad (the asymmetric maximum is generated by the outstar).  UCINET V seems not to adjust for this fact, which can cause some oddities in their centralization scores (thus, don't expect to get the same numbers with both packages).}

\section{Warning }{Singular adjacency matrices cause no end of headaches for this algorithm; thus, the routine may fail in certain cases.  This will be fixed when I get a better algorithm.  \code{bonpow} will not symmetrize your data before extracting eigenvectors; don't send this routine asymmetric matrices unless you really mean to do so.}

\seealso{ \code{\link{centralization}}, \code{\link{evcent}} }

\examples{
#Generate some test data
dat<-rgraph(10,mode="graph")
#Compute Bonpow scores
bonpow(dat,exponent=1,tol=1e-20)
bonpow(dat,exponent=-1,tol=1e-20)
}
\keyword{ univar }%-- one or more ...
\keyword{ math }
\keyword{ graphs}
