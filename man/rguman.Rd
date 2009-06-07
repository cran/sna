\name{rguman}
\alias{rguman}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Draw Dyad Census-Conditioned Random Graphs }
\description{
  \code{rguman} generates random draws from a dyad census-conditioned uniform random graph distribution.
}
\usage{
rguman(n, nv, mut = 0.25, asym = 0.5, null = 0.25, 
    method = c("probability", "exact"), return.as.edgelist = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ the number of graphs to generate. }
  \item{nv}{ the size of the vertex set (\eqn{|V(G)|}) for the random graphs. }
  \item{mut}{ if \code{method=="probability"}, the probability of obtaining a mutual dyad; otherwise, the number of mutual dyads. }
  \item{asym}{ if \code{method=="probability"}, the probability of obtaining an asymmetric dyad; otherwise, the number of asymmetric dyads. }
  \item{null}{ if \code{method=="probability"}, the probability of obtaining a null dyad; otherwise, the number of null dyads. }
  \item{method}{ the generation method to use.  \code{"probability"} results in a multinomial dyad distribution (conditional on the underlying rates), while \code{"exact"} results in a uniform draw conditional on the exact dyad distribution. }
  \item{return.as.edgelist}{ logical; should the resulting graphs be returned in edgelist form?}
}
\details{
A simple generalization of the Erdos-Renyi family, the U|MAN distributions are uniform on the set of graphs, conditional on order (size) and the dyad census.  As with the E-R case, there are two U|MAN variants.  The first (corresponding to \code{method=="probability"}) takes dyad states as independent multinomials with parameters \eqn{m} (for mutuals), \eqn{a} (for asymmetrics), and \eqn{n} (for nulls).  The resulting pmf is then
\deqn{
p(G=g|m,a,n) = \frac{(M+A+N)!}{M!A!N!} m^M a^A n^N,
}{%
p(G=g|m,a,n) = (M+A+N)!/(M!A!N!) m^M a^A n^N,
}
where \eqn{M}, \eqn{A}, and \eqn{N} are realized counts of mutual, asymmetric, and null dyads, respectively.  (See \code{\link{dyad.census}} for an explication of dyad types.) 

The second U|MAN variant is selected by \code{method=="exact"}, and places equal mass on all graphs having the specified (exact) dyad census.  The corresponding pmf is
\deqn{
p(G=g|M,A,N) = \frac{M!A!N!}{(M+A+N)!}.
}{%
p(G=g|M,A,N) = M!A!N!/(M+A+N)!.
}

U|MAN graphs provide a natural baseline model for networks which are constrained by size, density, and reciprocity.  In this way, they provide a bridge between edgewise models (e.g., the E-R family) and models with higher order dependence (e.g., the Markov graphs).
}
\value{
   A matrix or array containing the drawn adjacency matrices
}
\references{ 
Holland, P.W. and Leinhardt, S.  (1976).  \dQuote{Local Structure in Social Networks.}  In D. Heise (Ed.), \emph{Sociological Methodology}, pp 1-45.  San Francisco: Jossey-Bass.

Wasserman, S. and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press.
}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ The famous mathematicians referenced in this man page now have misspelled names, due to R's difficulty with accent marks.  }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{rgraph}}, \code{\link{rgnm}}, \code{\link{dyad.census}} }
\examples{
#Show some examples of extreme U|MAN graphs
gplot(rguman(1,10,mut=45,asym=0,null=0,method="exact")) #Clique
gplot(rguman(1,10,mut=0,asym=45,null=0,method="exact")) #Tournament
gplot(rguman(1,10,mut=0,asym=0,null=45,method="exact")) #Empty

#Draw a sample of multinomial U|MAN graphs
g<-rguman(5,10,mut=0.15,asym=0.05,null=0.8)

#Examine the dyad census
dyad.census(g)
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ graphs}% __ONLY ONE__ keyword per line
