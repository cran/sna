\name{component.dist}
\alias{component.dist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Calculate the Component Size Distribution of a Graph }
\description{
  \code{component.dist} returns a data frame containing a vector of length n such that the ith element contains the number of components of G having size i, and a vector of length n giving component membership.  Component strength is determined by the rule used to symmetrize the adjacency matrix; this is controlled by the eponymous parameter given to the \code{\link{symmetrize}} command.
}
\usage{
component.dist(dat, connected="strong")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ A single nxn adjacency matrix }
  \item{connected}{ A string corresponding to the \code{rule} parameter of the \code{\link{symmetrize}} command; by default, ``strong'' components are used.}
}
\details{
  Use ``strong'' to find strongly connected components, ``weak'' for weakly connected components, etc.  If \code{dat} is already symmetric, then the \code{connected} parameter has no effect.
}
\value{
A data frame containing:

  \item{membership }{A vector of component memberships, by vertex}
  \item{csize}{A vector of component sizes, by component}
  \item{cdist}{A vector of length |V(G)| with the (unnormalized) empirical distribution function of component sizes}
}
\references{ West, D.B.  (1996).  \emph{Introduction to Graph Theory.}  Upper Saddle River, N.J.: Prentice Hall.}
\author{ Carter T. Butts \email{ctb@andrew.cmu.edu} }


\seealso{ \code{\link{components}}, \code{\link{symmetrize}}, \code{\link{geodist}} }

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
\keyword{ network }%-- one or more ...
