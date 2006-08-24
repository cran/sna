\name{ego.extract}
\alias{ego.extract}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract Egocentric Networks from Complete Network Data }
\description{
  \code{ego.extract} takes one or more input graphs (\code{dat}) and returns a list containing the egocentric networks centered on vertices named in \code{ego}, using adjacency rule \emph{neighborhood} to define inclusion.
}
\usage{
ego.extract(dat, ego = NULL, neighborhood = c("combined", "in",
    "out"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more graphs. }
  \item{ego}{ a vector of vertex IDs, or \code{NULL} if all are to be selected. }
  \item{neighborhood}{ the neighborhood to use. }
}
\details{
The egocentric network (or \dQuote{ego net}) of vertex \eqn{v} in graph \eqn{G} is defined as \eqn{G[v \cup N(v)]}{G[v U N(v)]} (i.e., the subgraph of \eqn{G} induced by \eqn{v} and its neighborhood).  The neighborhood employed by \code{ego.extract} is selected by the eponymous argument: \code{"in"} selects in-neighbors, \code{"out"} selects out-neighbors, and \code{"combined"} selects all neighbors.  In the event that one of the vertices selected by \code{ego} has no qualifying neighbors, \code{ego.extract} will return a degenerate (1 by 1) adjacency matrix containing that individual's diagonal entry.

Vertices within the returned matrices are maintained in their original order, save for ego (who is always listed first).  The ego nets themselves are returned in the order specified in the \code{ego} parameter (or their vertex order, if no value was specified).

\code{ego.extract} is useful for finding local properties associated with particular vertices.  To compute functions of neighbors' covariates, see \code{\link{gapply}}.
}
\value{
  A list containing the adjacency matrices for the ego nets of each vertex in \code{ego}.
}
\references{ Wasserman, S. and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press.}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} %~
%}
\seealso{ \code{\link{gapply}} }
\examples{
#Generate a sample network
g<-rgraph(10,tp=1.5/9)

#Extract some ego nets
g.in<-ego.extract(g,neighborhood="in")
g.out<-ego.extract(g,neighborhood="out")
g.comb<-ego.extract(g,neighborhood="in")

#View some networks
g.comb

#Compare ego net size with degree
all(sapply(g.in,NROW)==degree(g,cmode="indegree")+1)    #TRUE
all(sapply(g.out,NROW)==degree(g,cmode="outdegree")+1)  #TRUE
all(sapply(g.comb,NROW)==degree(g)/2+1)        #Usually FALSE!

#Calculate egocentric network density
ego.size<-sapply(g.comb,NROW)
if(any(ego.size>2))
  sapply(g.comb[ego.size>2],function(x){gden(x[-1,-1])})
}
\keyword{ math }% at least one, from doc/KEYWORDS
\keyword{ graphs }% __ONLY ONE__ keyword per line
