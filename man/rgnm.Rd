\name{rgnm}
\alias{rgnm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Draw Density-Conditioned Random Graphs }
\description{
  \code{rgnm} generates random draws from a density-conditioned uniform random graph distribution. 
}
\usage{
rgnm(n, nv, m, mode = "digraph", diag = FALSE, 
    return.as.edgelist = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ the number of graphs to generate. }
  \item{nv}{ the size of the vertex set (\eqn{|V(G)|}) for the random graphs. }
  \item{m}{ the number of edges on which to condition. }
  \item{mode}{ \code{"digraph"} for directed graphs, or \code{"graph"} for undirected graphs. }
  \item{diag}{ logical; should loops be allowed? }
  \item{return.as.edgelist}{ logical; should the resulting graphs be returned in edgelist form?}
}
\details{
  \code{rgnm} returns draws from the density-conditioned uniform random graph first popularized by the famous work of Erdos and Renyi (the \eqn{G(N,M)} process).  In particular, the pmf of a \eqn{G(N,M)} process is given by

\deqn{
p(G=g|N,M) = \left( {E_m \atop M} \right)^{-1}
}{%
p(G=g|N,M) = 1/Choose(E_m,M)
}

where \eqn{E_m} is the maximum number of edges in the graph.  (\eqn{E_m} is equal to \code{nv*(nv-diag)/(1+(mode=="graph"))}.)

The \eqn{G(N,M)} process is one of several process which are used as baseline models of social structure.  Other well-known baseline models include the Bernoulli graph (the \eqn{G(N,p)} model of Erdos and Renyi) and the U|MAN model of dyadic independence.  These are implemented within \code{sna} as \code{\link{rgraph}} and \code{\link{rgnm}}, respectively.
}
\value{
  A matrix or array containing the drawn adjacency matrices
}
\references{
Erdos, P. and Renyi, A.  (1960).  \dQuote{On the Evolution of Random Graphs.}  \emph{Public Mathematical Institute of Hungary Academy of Sciences,} 5:17-61. 
}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ The famous mathematicians referenced in this man page now have misspelled names, due to R's difficulty with accent marks.   }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{rgraph}}, \code{\link{rguman}} }
\examples{
#Draw 5 random graphs of order 10 
all(gden(rgnm(5,10,9,mode="graph"))==0.2) #Density 0.2
all(gden(rgnm(5,10,9))==0.1)              #Density 0.1

#Plot a random graph
gplot(rgnm(1,10,20))
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ graphs }% __ONLY ONE__ keyword per line
