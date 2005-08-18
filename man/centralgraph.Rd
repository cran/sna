\name{centralgraph}
\alias{centralgraph}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Find the Central Graph of a Labeled Graph Stack }
\description{
  Returns the central graph of a set of labeled graphs, i.e. that graph in which i->j iff i->j in >=50\% of the graphs within the set.  If \code{normalize==TRUE}, then the value of the i,jth edge is given as the proportion of graphs in which i->j.
}
\usage{
centralgraph(dat, normalize=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{normalize}{ boolean indicating whether the results should be normalized.  The result of this is the "mean matrix".  By default, \code{normalize==FALSE}.}
}
\details{
  The central graph of a set of graphs S is that graph C which minimizes the sum of Hamming distances between C and G in S.  As such, it turns out (for the dichotomous case, at least), to be analogous to both the mean and median for sets of graphs.  The central graph is useful in a variety of contexts; see the references below for more details.
}
\value{
A matrix containing the central graph (or mean matrix)
}
\references{ 
Banks, D.L., and Carley, K.M.  (1994).  ``Metric Inference for Social Networks.''  \emph{Journal of Classification}, 11(1), 121-49. 
 }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ 0.5 is used as the cutoff value regardless of whether or not the data is dichotomous (as is tacitly assumed).  The routine is unaffected by data type when \code{normalize==TRUE}. }

\seealso{ \code{\link{hdist}} }

\examples{
#Generate some random graphs
dat<-rgraph(10,5)
#Find the central graph
cg<-centralgraph(dat)
#Plot the central graph
gplot(cg)
#Now, look at the mean matrix
cg<-centralgraph(dat,normalize=TRUE)
print(cg)
}
\keyword{ math }%-- one or more ...
\keyword{ graphs}
