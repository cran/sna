\name{component.dist}
\alias{component.dist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Calculate the Component Size Distribution of a Graph }
\description{
  \code{component.dist} returns a data frame containing a vector of length n such that the ith element contains the number of components of G having size i, and a vector of length n giving component membership.  Component strength is determined by the \code{connected} parameter; see below for details.
}
\usage{
component.dist(dat, connected=c("strong","weak","unilateral",
     "recursive"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ A single nxn adjacency matrix }
  \item{connected}{ A string selecting strong, weak, unilateral or recursively connected components; by default, \code{"strong"} components are used.}
}
\details{
Components are maximal sets of mutually connected vertices; depending on the definition of ``connected'' one employs, one can arrive at several types of components.  Those supported here are as follows (in increasing order of restrictiveness):

\enumerate{
\item \code{weak}:  \eqn{v_1} is connected to \eqn{v_2} iff there exists a semi-path from \eqn{v_1} to \eqn{v_2} (i.e., a path in the weakly symmetrized graph)

\item \code{unilateral}:  \eqn{v_1} is connected to \eqn{v_2} iff there exists a directed path from \eqn{v_1} to \eqn{v_2} \emph{or} a directed path from \eqn{v_2} to \eqn{v_1}

\item \code{strong}:  \eqn{v_1} is connected to \eqn{v_2} iff there exists a directed path from \eqn{v_1} to \eqn{v_2} \emph{and} a directed path from \eqn{v_2} to \eqn{v_1}

\item \code{recursive}:  \eqn{v_1} is connected to \eqn{v_2} iff there exists a vertex sequence \eqn{v_a,\ldots,v_z}{v_a,...,v_z} such that \eqn{v_1,v_a,\ldots,v_z,v_2}{v_1,v_a,...,v_z,v_2} and \eqn{v_2,v_z,\ldots,v_a,v_1}{v_2,v_z,...,v_a,v_1} are directed paths
}

Note that the above definitions are distinct for directed graphs only; if \code{dat} is symmetric, then the \code{connected} parameter has no effect.  
}
\value{
A data frame containing:

  \item{membership }{A vector of component memberships, by vertex}
  \item{csize}{A vector of component sizes, by component}
  \item{cdist}{A vector of length |V(G)| with the (unnormalized) empirical distribution function of component sizes}
}
\references{ West, D.B.  (1996).  \emph{Introduction to Graph Theory.}  Upper Saddle River, N.J.: Prentice Hall.}
\author{ Carter T. Butts \email{buttsc@uci.edu} }

\section{Note}{Unilaterally connected component partitions may not be well-defined, since it is possible for a given vertex to be unilaterally connected to two vertices which are not unilaterally connected with one another.  Consider, for instance, the graph \eqn{a \rightarrow b \leftarrow c \rightarrow d}{a->b<-c<-d}.  In this case, the maximal unilateral components are \eqn{ab} and \eqn{bcd}, with vertex \eqn{b} properly belonging to both components.  For such graphs, a unique partition of vertices by component does not exist, and we ``solve'' the problem by allocating each ``problem vertex'' to one of its components on an essentially arbitrary basis.  (\code{component.dist} generates a warning when this occurs.)  It is recommended that the \code{unilateral} option be avoided where possible.}

\seealso{ \code{\link{components}}, \code{\link{symmetrize}}, \code{\link{reachability}} \code{\link{geodist}} }

\examples{
g<-rgraph(20,tprob=0.075)   #Generate a sparse random graph

#Find weak components
cd<-component.dist(g,connected="weak")
cd$membership              #Who's in what component?
cd$csize                   #What are the component sizes?
                           #Plot the size distribution
plot(1:length(cd$cdist),cd$cdist/sum(cd$cdist),ylim=c(0,1),type="h")  

#Find strong components
cd<-component.dist(g,connected="strong")
cd$membership              #Who's in what component?
cd$csize                   #What are the component sizes?
                           #Plot the size distribution
plot(1:length(cd$cdist),cd$cdist/sum(cd$cdist),ylim=c(0,1),type="h")

}
\keyword{ math }%-- one or more ...
