\name{plot.blockmodel}
\alias{plot.blockmodel}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Plotting for blockmodel Objects }
\description{
  Displays a plot of the blocked data matrix, given a blockmodel object.
}
\usage{
\method{plot}{blockmodel}(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ An object of class \code{blockmodel} }
  \item{\dots}{ Further arguments passed to or from other methods }
}
\details{
Plots of the blocked data matrix (i.e., the data matrix with rows and columns permuted to match block membership) can be useful in assessing the strength of the block solution (particularly for clique detection and/or regular equivalence).
}
\value{
  None
}
\references{ White, H.C.; Boorman, S.A.; and Breiger, R.L.  (1976).  ``Social Structure from Multiple Networks I: Blockmodels of Roles and Positions.''  \emph{American Journal of Sociology}, 81, 730-779.

Wasserman, S., and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press.}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

 %~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{blockmodel}}, \code{\link{plot.sociomatrix}} }

\examples{
#Create a random graph with _some_ edge structure
g.p<-sapply(runif(20,0,1),rep,20)  #Create a matrix of edge 
                                   #probabilities
g<-rgraph(20,tprob=g.p)            #Draw from a Bernoulli graph 
                                   #distribution

#Cluster based on structural equivalence
eq<-equiv.clust(g)

#Form a blockmodel with distance relaxation of 10
b<-blockmodel(g,eq,h=10)
plot(b)                            #Plot it
}
\keyword{ hplot }%-- one or more ...
