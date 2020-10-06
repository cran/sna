\name{simmelian}
\alias{simmelian}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Find the Simmelian Tie Structure of a Graph }
\description{
  \code{simmelian} takes one or more (possibly directed) graphs as input, producing the associated Simmelian tie structures.
}
\usage{
simmelian(dat, dichotomize=TRUE, return.as.edgelist=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more graphs (directed or otherwise). }
  \item{dichotomize}{ logical; should the presence or absence of Simmelian edges be returned?  If \code{FALSE}, returned edges are valued by the number of 3-clique co-memberships for each dyad. }
  \item{return.as.edgelist}{ logical; return the result as an sna edgelist? }
}
\details{
For a digraph \eqn{G=(V,E)} with vertices \eqn{i} and \eqn{j}, then \eqn{i} and \eqn{j} are said to have a \emph{Simmelian tie} iff \eqn{i} and \eqn{j} belong to a 3-clique of \eqn{G}.  (Note that, in the undirected case, we simply treat \eqn{G} as a fully mutual digraph.)  Because they have both a mutual dyad and mutual ties to/from at least one third party, vertex pairs with Simmelian ties in interpersonal networks are often expected to have strong relationships; Simmelian ties may also be more stable than other relationships, due to reinforcement from the mutual shared partner.  In other settings, the derived network of Simmelian ties (which is simply the co-membership network of non-trivial cliques) may be useful for identifying cohesively connected elements in a larger graph, or for finding \dQuote{backbone} structures in networks with large numbers of unreciprocated and/or bridging ties.

Currently, Simmelian tie calculation is performed using \code{\link{kcycle.census}}.  While the bulk of the calculations and data handling are performed using edgelists, \code{\link{kcycle.census}} currently returns co-memberships in adjacency form.  The implication for the end user is that performance for \code{simmelian} will begin to degrade for networks on the order of ten thousand vertices or so (due to the cost of allocating the adjacency structure), irrespective of the content of the network or other settings.  This bottleneck will likely be removed in future versions.
}
\value{
An adjacency matrix containing the Simmelian ties, or the equivalent edgelist representation
}
\references{ Krackhardt, David.  (1999).  \dQuote{The Ties That Torture: Simmelian Tie Analysis in Organizations.}  \emph{Research in the Sociology of Organizations}, 16:183-210.}
\author{ Carter T. Butts \email{buttsc@uci.edu} }

\seealso{ \code{\link{kcycle.census}}, \code{\link{clique.census}} }

\examples{
#Contrast the Simmelian ties in the Coleman friendship network with the "raw" ties
data(coleman)
fall<-coleman[1,,]                   #Fall ties
spring<-coleman[2,,]                 #Spring ties
sim.fall<-simmelian(coleman[1,,])    #Fall Simmelian ties
sim.spring<-simmelian(coleman[2,,])  #Spring Simmelian ties

par(mfrow=c(2,2))
gplot(fall,main="Nominations in Fall")
gplot(spring,main="Nominations in Spring")
gplot(sim.fall,main="Simmelian Ties in Fall")
gplot(sim.spring,main="Simmelian Ties in Spring")

#Which ties shall survive?
table(fall=gvectorize(fall),spring=gvectorize(spring))   #Fall vs. spring
table(sim.fall=gvectorize(sim.fall),spring=gvectorize(spring))
sum(fall&spring)/sum(fall)                #About 58% of ties survive, overall...
sum(sim.fall&spring)/sum(sim.fall)        #...but 74% of Simmelian ties survive!
sum(sim.fall&sim.spring)/sum(sim.fall)    #(About 44% stay Simmelian.)
sum(sim.fall&sim.spring)/sum(sim.spring)  #39% of spring Simmelian ties were so in fall
sum(fall&sim.spring)/sum(sim.spring)      #and 67% had at least some tie in fall
}
\keyword{ algebra }% __ONLY ONE__ keyword per line
\keyword{ graphs }
