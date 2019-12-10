\name{component.dist}
\alias{component.dist}
\alias{component.largest}
\alias{undirComponents_R}
\alias{component_dist_R}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Calculate the Component Size Distribution of a Graph }
\description{
  \code{component.dist} returns a list containing a vector of length \code{n} such that the \code{i}th element contains the number of components of graph \eqn{G} having size \code{i}, and a vector of length \code{n} giving component membership (where \code{n} is the graph order).  Component strength is determined by the \code{connected} parameter; see below for details.

  \code{component.largest} identifies the component(s) of maximum order within graph \code{G}.  It returns either a \code{logical} vector indicating membership in a maximum component or the adjacency matrix of the subgraph of \eqn{G} induced by the maximum component(s), as determined by \code{result}.  Component strength is determined as per \code{component.dist}.
}
\usage{
component.dist(dat, connected=c("strong","weak","unilateral",
     "recursive"))

component.largest(dat, connected=c("strong","weak","unilateral",
     "recursive"), result = c("membership", "graph"), return.as.edgelist =
     FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{connected}{ a string selecting strong, weak, unilateral or recursively connected components; by default, \code{"strong"} components are used.}
  \item{result}{ a string indicating whether a vector of membership indicators or the induced subgraph of the component should be returned. }
  \item{return.as.edgelist}{ logical; if \code{result=="graph"}, should the resulting structure be returned in edgelist form?}
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
For \code{component.dist}, a list containing:

  \item{membership }{A vector of component memberships, by vertex}
  \item{csize}{A vector of component sizes, by component}
  \item{cdist}{A vector of length |V(G)| with the (unnormalized) empirical distribution function of component sizes}

If multiple input graphs are given, the return value is a list of lists.

For \code{component.largest}, either a \code{logical} vector of component membership indicators or the adjacency matrix/edgelist of the subgraph induced by the largest component(s) is returned.  If multiple graphs were given as input, a list of results is returned.
}
\references{ West, D.B.  (1996).  \emph{Introduction to Graph Theory.}  Upper Saddle River, N.J.: Prentice Hall.}
\author{ Carter T. Butts \email{buttsc@uci.edu} }

\section{Note}{Unilaterally connected component partitions may not be well-defined, since it is possible for a given vertex to be unilaterally connected to two vertices that are not unilaterally connected with one another.  Consider, for instance, the graph \eqn{a \rightarrow b \leftarrow c \rightarrow d}{a->b<-c<-d}.  In this case, the maximal unilateral components are \eqn{ab} and \eqn{bcd}, with vertex \eqn{b} properly belonging to both components.  For such graphs, a unique partition of vertices by component does not exist, and we ``solve'' the problem by allocating each ``problem vertex'' to one of its components on an essentially arbitrary basis.  (\code{component.dist} generates a warning when this occurs.)  It is recommended that the \code{unilateral} option be avoided where possible.

Do not make the mistake of assuming that the subgraphs returned by \code{component.largest} are necessarily connected.  This is \emph{usually} the case, but depends upon the uniqueness of the largest component.
}

\seealso{ \code{\link{components}}, \code{\link{symmetrize}}, \code{\link{reachability}} \code{\link{geodist}} }

\examples{
g<-rgraph(20,tprob=0.06)   #Generate a sparse random graph

#Find weak components
cd<-component.dist(g,connected="weak")
cd$membership              #Who's in what component?
cd$csize                   #What are the component sizes?
                           #Plot the size distribution
plot(1:length(cd$cdist),cd$cdist/sum(cd$cdist),ylim=c(0,1),type="h")  
lgc<-component.largest(g,connected="weak")  #Get largest component
gplot(g,vertex.col=2+lgc)  #Plot g, with component membership
                           #Plot largest component itself 
gplot(component.largest(g,connected="weak",result="graph"))

#Find strong components
cd<-component.dist(g,connected="strong")
cd$membership              #Who's in what component?
cd$csize                   #What are the component sizes?
                           #Plot the size distribution
plot(1:length(cd$cdist),cd$cdist/sum(cd$cdist),ylim=c(0,1),type="h")
lgc<-component.largest(g,connected="strong")  #Get largest component
gplot(g,vertex.col=2+lgc)  #Plot g, with component membership
                           #Plot largest component itself 
gplot(component.largest(g,connected="strong",result="graph"))
}
\keyword{ math }%-- one or more ...
\keyword{ graphs } 
