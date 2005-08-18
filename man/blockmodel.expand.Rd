\name{blockmodel.expand}
\alias{blockmodel.expand}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Generate a Graph (or Stack) from a Given Blockmodel Using Particular Expansion Rules }
\description{
  \code{blockmodel.expand} takes a blockmodel and an expansion vector, and expands the former by making copies of the vertices.
}
\usage{
blockmodel.expand(b, ev, mode="digraph", diag=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{b}{ blockmodel object. }
  \item{ev}{ a vector indicating the number of copies to make of each class (respectively). }
  \item{mode}{ a string indicating whether the result should be a ``graph'' or ``digraph''. }
  \item{diag}{ a boolean indicating whether or not loops should be permitted. }
}
\details{
  The primary use of blockmodel expansion is in generating test data from a blockmodeling hypothesis.  Expansion is performed depending on the content type of the blockmodel; at present, only density is supported.  For the density content type, expansion is performed by interpreting the interclass density as an edge probability, and by drawing random graphs from the Bernoulli parameter matrix formed by expanding the density model.  Thus, repeated calls to \code{blockmodel.expand} can be used to generate a sample for monte carlo null hypothesis tests under a Bernoulli graph model.
}
\value{
An adjacency matrix, or stack thereof.
}
\references{ 
Doreian, P.; Batagelj, V.; and Ferligoj, A.  (2005).  \emph{Generalized Blockmodeling.}  Cambridge: Cambridge University Press.

White, H.C.; Boorman, S.A.; and Breiger, R.L.  (1976).  ``Social Structure from Multiple Networks I: Blockmodels of Roles and Positions.''  \emph{American Journal of Sociology}, 81, 730-779. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ Eventually, other content types will be supported. }


\seealso{ \code{\link{blockmodel}} }

\examples{
#Create a random graph with _some_ edge structure
g.p<-sapply(runif(20,0,1),rep,20)  #Create a matrix of edge 
                                   #probabilities
g<-rgraph(20,tprob=g.p)            #Draw from a Bernoulli graph 
                                   #distribution

#Cluster based on structural equivalence
eq<-equiv.clust(g)

#Form a blockmodel with distance relaxation of 15
b<-blockmodel(g,eq,h=15)

#Draw from an expanded density blockmodel
g.e<-blockmodel.expand(b,rep(2,length(b$rlabels)))  #Two of each class
g.e
}
\keyword{ manip }%-- one or more ...
\keyword{ math }
\keyword{ graphs}
