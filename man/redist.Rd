\name{redist}
\Rdversion{1.1}
\alias{redist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find a Matrix of Distances Between Positions Based on Regular Equivalence
}
\description{
\code{redist} uses the graphs indicated by \code{g} in \code{dat} to assess the extent to which each vertex is regularly equivalent; \code{method} determines the measure of approximate equivalence which is used (currently, only CATREGE).
}
\usage{
redist(dat, g = NULL, method = c("catrege"), mode = "digraph", 
    diag = FALSE, seed.partition = NULL, code.diss = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
a graph or set thereof.
}
  \item{g}{
a vector indicating which elements of \code{dat} should be examined (by default, all are used).
}
  \item{method}{
method to use when assessing regular equivalence (currently, only \code{"catrege"}).
}
  \item{mode}{
\code{"digraph"} for directed data, otherwise \code{"graph"}. 
}
  \item{diag}{
 logical; should diagonal entries (loops) should be treated as meaningful data?
 }
  \item{seed.partition}{
optionally, an initial equivalence partition to \dQuote{seed} the CATREGE algorithm.
}
  \item{code.diss}{
logical; return as dissimilarities (rather than similarities)?
}
  \item{\dots}{
additional parameters (currently ignored).
}
}
\details{
\code{redist} provides a basic tool for assessing the (approximate) regular equivalence of actors.  Two vertices \eqn{i} and \eqn{j} are said to be regularly equivalent with respect to role assignment \code{r} if \eqn{\{r(u): u\in N^+(i)\}=\{r(u): u\in N^+(j)\}}{{r(u): u in N^+(i)} = {r(u): u in N^+(j)}} and \eqn{\{r(u): u\in N^-(i)\}=\{r(u): u\in N^-(j)\}}{{r(u): u in N^-(i)} = {r(u): u in N^-(j)}}, where \eqn{N^+} and \eqn{N^-} denote out- and in-neighborhoods (respectively).  RE similarity/difference scores are computed by \code{method}, currently Borgatti and Everett's CATREGE algorithm (which is based on the multiplex maximal regular equivalence on \eqn{G} and its transpose).  The \dQuote{distance} between positions in this case is the inverse of the number of iterative refinements of the initial equivalence (i.e., role) structure required to allocate the positions to regularly equivalent roles (with 0 indicating positions which ultimately belong in the same role).  By default, the initial equivalence structure is one in which all vertices are treated as occupying the same role; the \code{seed.partition} option can be used to impose alternative constraints.  From this initial structure, vertices within the same role having non-identical mixes of neighbor types are re-allocated to different roles (where \dQuote{neighbor type} is initially due to the pattern of (possibly valued) in- and out-ties, cross-classified by current alter type).  This procedure is then iterated until no further division of roles is necessary to satisfy the regularity condition.

Once the similarities/differences are calculated, the results can be used with a clustering routine (such as \code{\link{equiv.clust}}) or an MDS (such as \code{\link{cmdscale}}) to identify the underlying role structure.
}
\value{
A matrix of similarity/difference scores.
}
\references{
Borgatti, S.P. and Everett, M.G.  (1993).  \dQuote{Two Algorithms for Computing Regular Equivalence.}  \emph{Social Networks}, 15, 361-376.
}
\author{
Carter T. Butts \email{buttsc@uci.edu}
}
\note{
The maximal regular equivalence is often very uninteresting (i.e., degenerate) for unvalued, undirected graphs.  An exogenous constraint (e.g., via the \code{seed.partition}) may be required to uncover a more useful refinement of the unconstrained maximal equivalence.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sedist}}, \code{\link{equiv.clust}}
}
\examples{
#Create a random graph with _some_ edge structure
g.p<-sapply(runif(20,0,1),rep,20)  #Create a matrix of edge 
                                   #probabilities
g<-rgraph(20,tprob=g.p)            #Draw from a Bernoulli graph 
                                   #distribution

#Get RE distances
g.re<-redist(g)

#Plot a metric MDS of vertex positions in two dimensions
plot(cmdscale(as.dist(g.re)))

#What if there were already something known to be different about
#the first five vertices?
sp<-rep(1:2,times=c(5,15))            #Create "seed" partition
g.spre<-redist(g,seed.partition=sp)   #Get new RE distances
g.spre
plot.sociomatrix(g.spre)              #Note the blocking! 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
\keyword{ graphs }% __ONLY ONE__ keyword per line
\keyword{ classif }
\keyword{ cluster }
