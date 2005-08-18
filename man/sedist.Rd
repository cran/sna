\name{sedist}
\alias{sedist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Find a Matrix of Distances Between Positions Based on Structural Equivalence }
\description{
\code{sedist} uses the graphs indicated by \code{g} in \code{dat} to assess the extent to which each vertex is structurally equivalent; \code{joint.analysis} determines whether this analysis is simultaneous, and \code{method} determines the measure of approximate equivalence which is used.
}
\usage{
sedist(dat, g=c(1:dim(dat)[1]), method="hamming", 
    joint.analysis=FALSE, mode="digraph", diag=FALSE, code.diss=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ a graph or set thereof. }
  \item{g}{ a vector indicating which elements of \code{dat} should be examined. }
  \item{method}{ one of \code{"correlation"}, \code{"euclidean"}, \code{"hamming"}, or \code{"gamma"}.}
  \item{joint.analysis}{ should equivalence be assessed across all networks jointly (\code{TRUE}), or individually within each (\code{FALSE})? }
  \item{mode}{ \code{"digraph"} for directed data, otherwise \code{"graph"}. }
  \item{diag}{ boolean indicating whether diagonal entries (loops) should be treated as meaningful data. }
  \item{code.diss}{ reverse-code the raw comparison values.}
}
\details{
\code{sedist} provides a basic tool for assessing the (approximate) structural equivalence of actors.  (Two vertices i and j are said to be structurally equivalent if i->k iff j->k for all k.)  SE similarity/difference scores are computed by comparing vertex rows and columns using the measure indicated by \code{method}:
\enumerate{
\item correlation: the product-moment correlation

\item euclidean: the euclidean distance

\item hamming: the Hamming distance

\item gamma: the gamma correlation
}
Once these similarities/differences are calculated, the results can be used with a clustering routine (such as \code{\link{equiv.clust}}) or an MDS (such as \code{\link{cmdscale}}).
}
\value{
A matrix of similarity/difference scores
}
\references{Breiger, R.L.; Boorman, S.A.; and Arabie, P.  (1975).  ``An Algorithm for Clustering Relational Data with Applications to Social Network Analysis and Comparison with Multidimensional Scaling.''  \emph{Journal of Mathematical Psychology}, 12, 328-383.

Burt, R.S.  (1976).  ``Positions in Networks.''  \emph{Social Forces}, 55, 93-122.

Wasserman, S., and Faust, K.  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press.   }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ Be careful to verify that you have computed what you meant to compute, with respect to similarities/differences. Also, note that (despite its popularity) the product-moment correlation can give rather strange results in some cases.
}

\seealso{ \code{\link{equiv.clust}}, \code{\link{blockmodel}} }

\examples{
#Create a random graph with _some_ edge structure
g.p<-sapply(runif(20,0,1),rep,20)  #Create a matrix of edge 
                                   #probabilities
g<-rgraph(20,tprob=g.p)            #Draw from a Bernoulli graph 
                                   #distribution

#Get SE distances
g.se<-sedist(g)

#Plot a metric MDS of vertex positions in two dimensions
plot(cmdscale(as.dist(g.se)))
}
\keyword{ classif }
\keyword{ cluster }%-- one or more ...
\keyword{ math }
\keyword{ graphs} 
