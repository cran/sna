\name{rgraph}
\alias{rgraph}
\alias{rgbern_R}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Generate Bernoulli Random Graphs }
\description{
  \code{rgraph} generates random draws from a Bernoulli graph distribution, with various parameters for controlling the nature of the data so generated.
}
\usage{
rgraph(n, m=1, tprob=0.5, mode="digraph", diag=FALSE, replace=FALSE,
    tielist=NULL, return.as.edgelist=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n}{ The size of the vertex set (|V(G)|) for the random graphs }
  \item{m}{ The number of graphs to generate }
  \item{tprob}{ Information regarding tie (edge) probabilities; see below }
  \item{mode}{ ``digraph'' for directed data, ``graph'' for undirected data }
  \item{diag}{ Should the diagonal entries (loops) be set to zero? }
  \item{replace}{ Sample with or without replacement from a tie list (ignored if \code{tielist==NULL} }
  \item{tielist}{ A vector of edge values, from which the new graphs should be bootstrapped }
  \item{return.as.edgelist}{ logical; should the resulting graphs be returned in edgelist form?}
}
\details{
   \code{rgraph} is a reasonably versatile routine for generating random network data.  The graphs so generated are either Bernoulli graphs (graphs in which each edge is a Bernoulli trial, independent conditional on the Bernoulli parameters), or are bootstrapped from a user-provided edge distribution (very handy for CUG tests).  In the latter case, edge data should be provided using the \code{tielist} argument; the exact form taken by the data is irrelevant, so long as it can be coerced to a vector.  In the former case, Bernoulli graph probabilities are set by the \code{tprob} argument as follows:
\enumerate{
\item If \code{tprob} contains a single number, this number is used as the probability of all edges.

\item If \code{tprob} contains a vector, each entry is assumed to correspond to a separate graph (in order).  Thus, each entry is used as the probability of all edges within its corresponding graph.

\item If \code{tprob} contains a matrix, then each entry is assumed to correspond to a separate edge.  Thus, each entry is used as the probability of its associated edge in each graph which is generated.

\item Finally, if \code{tprob} contains a three-dimensional array, then each entry is assumed to correspond to a particular edge in a particular graph, and is used as the associated probability parameter.
}

Finally, note that \code{rgraph} will symmetrize all generated networks if \code{mode} is set to ``graph'' by copying down the upper triangle.  The lower half of \code{tprob}, where applicable, must still be specified, however.
}
\value{
  A graph stack
}
\references{ 
Erdos, P. and Renyi, A.  (1960).  \dQuote{On the Evolution of Random Graphs.}  \emph{Public Mathematical Institute of Hungary Academy of Sciences,} 5:17-61.

Wasserman, S., and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications}.  Cambridge: Cambridge University Press. 
}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ The famous mathematicians referenced in this man page now have misspelled names, due to R's difficulty with accent marks.  }

 %~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{rmperm}}, \code{\link{rgnm}}, \code{\link{rguman}} }

\examples{

#Generate three graphs with different densities
g<-rgraph(10,3,tprob=c(0.1,0.9,0.5))

#Generate from a matrix of Bernoulli parameters
g.p<-matrix(runif(25,0,1),nrow=5)
g<-rgraph(5,2,tprob=g.p)
}
\keyword{ distribution }%-- one or more ...
\keyword{ math }
