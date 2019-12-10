\name{maxflow}
\Rdversion{1.1}
\alias{maxflow}
\alias{maxflow_EK_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Maximum Flows Between Vertices
}
\description{
\code{maxflow} calculates a matrix of maximum pairwise flows within a (possibly valued) input network.
}
\usage{
maxflow(dat, src = NULL, sink = NULL, ignore.eval = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
one or more input graphs.
}
  \item{src}{
optionally, a vector of source vertices; by default, all vertices are selected.
}
  \item{sink}{
optionally, a vector of sink (or target) vertices; by default, all vertices are selected.
}
  \item{ignore.eval}{
logical; ignore edge values (i.e., assume unit capacities) when computing flow?
}
}
\details{
\code{maxflow} computes the maximum flow from each source vertex to each sink vertex, assuming infinite vertex capacities and limited edge capacities.  If \code{ignore.eval==FALSE}, supplied edge values are assumed to contain capacity information; otherwise, all non-zero edges are assumed to have unit capacity.

Note that all flows computed here are pairwise -- i.e., when computing the flow from \eqn{v} to \eqn{v'}, we ignore any other flows which could also be taking place within the network.  As a result, it should not be assumed that these flows can be realized \emph{simultaneously}.  (For the latter purpose, the values returned by \code{maxflow} can be treated as upper bounds.)
}
\value{
A matrix of pairwise maximum flows (if multiple sources/sinks selected), or a single maximum flow value (otherwise).
}
\references{
Edmonds, J. and Karp, R.M. (1972). \dQuote{Theoretical Improvements in Algorithmic Efficiency for Network Flow Problems.}  \emph{Journal of the ACM,} 19(2), 248-264.
}
\author{
Carter T. Butts \email{buttsc@uci.edu}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{flowbet}}, \code{\link{geodist}}
}
\examples{
g<-rgraph(10,tp=2/9)                     #Generate a sparse random graph
maxflow(g)                               #Compute all-pairs max flow
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
\keyword{ graphs }% __ONLY ONE__ keyword per line
