\name{rgws}
\alias{rgws}
\alias{rewire.ws}
\alias{rewire.ud}
\alias{udrewire_R}
\alias{wsrewire_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Draw From the Watts-Strogatz Rewiring Model }
\description{
  \code{rgws} generates draws from the Watts-Strogatz rewired lattice model.  Given a set of input graphs, \code{rewire.ws} performs a (dyadic) rewiring of those graphs.
}
\usage{
rgws(n, nv, d, z, p, return.as.edgelist = FALSE)
rewire.ud(g, p, return.as.edgelist = FALSE)
rewire.ws(g, p, return.as.edgelist = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ the number of draws to take. }
  \item{nv}{ the number of vertices per lattice dimension. }
  \item{d}{ the dimensionality of the underlying lattice. }
  \item{z}{ the nearest-neighbor threshold for local ties. }
  \item{p}{ the dyadic rewiring probability. }
  \item{g}{ a graph or graph stack. }
  \item{return.as.edgelist}{ logical; should the resulting graphs be returned in edgelist form?}
}
\details{
  A Watts-Strogatz graph process generates a random graph via the following procedure.  First, a \code{d}-dimensional uniform lattice is generated, here with \code{nv} vertices per dimension (i.e., \code{nv^d} vertices total).  Next, all \code{z} neighbors are connected, based on geodesics of the underlying lattice.  Finally, each non-null dyad in the resulting augmented lattice is "rewired" with probability \code{p}, where the rewiring operation exchanges the initial dyad state with the state of a uniformly selected null dyad sharing exactly one endpoint with the original dyad.  (In the standard case, this is equivalent to choosing an endpoint of the dyad at random, and then transferring the dyadic edges to/from that endpoint to another randomly chosen vertex.  Hence the "rewiring" metaphor.)  For \code{p==0}, the W-S process generates (deterministic) uniform lattices, approximating a uniform G(N,M) process as \code{p} approaches 1.  Thus, \code{p} can be used to tune overall entropy of the process.  A well-known property of the W-S process is that (for large \code{nv^d} and small \code{p}) it generates draws with short expected mean geodesic distances (approaching those found in uniform graphs) while maintaining high levels of local "clustering" (i.e., transitivity).  It has thus been proposed as one potential mechanism for obtaining "small world" structures.
  
\code{rgws} produces independent draws from the above process, returning them as an adjacency matrix (if \code{n==1}) or array (otherwise).  \code{rewire.ws}, on the other hand, applies the rewiring phase of the W-S process to one or more input graphs.  This can be used to explore local perturbations of the original graphs, conditioning on the dyad census.  \code{rewire.ud} is similar to \code{rewire.ws}, save in that all dyads are eligible for rewiring (not just non-null dyads), and exchanges with non-null dyads are permitted.  This process may be easier to work with than standard W-S rewiring in some cases.
}
\value{
A graph or graph stack containing draws from the appropriate W-S process.
}
\references{ 
Watts, D. and Strogatz, S. (1998).  ``Collective Dynamics of Small-world Networks.''  \emph{Nature,} 393:440-442.
 }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ \code{rgws} generates non-toroidal lattices; some published work in this area utilizes underlying toroids, so users should check for this prior to comparing simulations against published results. }

\section{Warning }{Remember that the total number of vertices in the graph is \code{nv^d}.  This can get out of hand \emph{very} quickly.}

\seealso{ \code{\link{rgnm}}, \code{\link{rgraph}} }
\examples{

#Generate Watts-Strogatz graphs, w/increasing levels of rewiring
gplot(rgws(1,100,1,2,0))     #No rewiring
gplot(rgws(1,100,1,2,0.01))  #1% rewiring
gplot(rgws(1,100,1,2,0.05))  #5% rewiring
gplot(rgws(1,100,1,2,0.1))   #10% rewiring
gplot(rgws(1,100,1,2,1))     #100% rewiring 

#Start with a simple graph, then rewire it
g<-matrix(0,50,50)
g[1,]<-1; g[,1]<-1    #Create a star
gplot(g)
gplot(rewire.ws(g,0.05))  #5% rewiring

}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ graphs }% __ONLY ONE__ keyword per line
