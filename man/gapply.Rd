\name{gapply}
\alias{gapply}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Apply Functions Over Vertex Neighborhoods }
\description{
  Returns a vector or array or list of values obtained by applying a function to vertex neighborhoods of a given order.
}
\usage{
gapply(X, MARGIN, STATS, FUN, ..., mode = "digraph", diag = FALSE, 
    distance = 1, thresh = 0, simplify = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ one or more input graphs. }
  \item{MARGIN}{ a vector giving the ``margin'' of \code{X} to be used in calculating neighborhoods.  1 indicates rows (out-neighbors), 2 indicates columns (in-neighbors), and c(1,2) indicates rows and columns (total neighborhood). }
  \item{STATS}{ the vector or matrix of vertex statistics to be used. }
  \item{FUN}{ the function to be applied.  In the case of operators, the function name must be quoted. }
  \item{\dots}{ additional arguments to \code{FUN}. }
  \item{mode}{ \code{"graph"} if \code{X} is a simple graph, else \code{"digraph"}. }
  \item{diag}{ boolean; are the diagonals of \code{X} meaningful? }
  \item{distance}{ the maximum geodesic distance at which neighborhoods are to be taken.  1 signifies first-order neighborhoods, 2 signifies second-order neighborhoods, etc. }
  \item{thresh}{ the threshold to be used in dichotomizing \code{X}. }
  \item{simplify}{ boolean; should we attempt to coerce output to a vector if possible? }
}
\details{
  For each vertex in \code{X}, \code{gapply} first identifies all members of the relevant neighborhood (as determined by \code{MARGIN} and \code{distance}) and pulls the rows of \code{STATS} associated with each.  \code{FUN} is then applied to this collection of values.  This provides a very quick and easy way to answer questions like:
\itemize{
\item How many persons are in each ego's 3rd-order neighborhood?
\item What fraction of each ego's alters are female?
\item What is the mean income for each ego's trading partners?
\item etc.
}
With clever use of \code{FUN} and \code{STATS}, a wide range of functionality can be obtained.
}
\value{
The result of the iterated application of \code{FUN} to each vertex neighborhood's \code{STATS}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Carter T. Butts \email{buttsc@uci.edu}}
%\note{ ~~further notes~~ }
%
% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{apply}}, \code{\link{sapply}} }
\examples{
#Generate a random graph
g<-rgraph(6)

#Calculate the degree of g using gapply
all(gapply(g,1,rep(1,6),sum)==degree(g,cmode="outdegree"))
all(gapply(g,2,rep(1,6),sum)==degree(g,cmode="indegree"))
all(gapply(g,c(1,2),rep(1,6),sum)==degree(symmetrize(g),cmode="freeman")/2)

#Find first and second order neighborhood means on some variable
gapply(g,c(1,2),1:6,mean)
gapply(g,c(1,2),1:6,mean,distance=2)

}
\keyword{ iteration }% at least one, from doc/KEYWORDS
\keyword{ manip }% __ONLY ONE__ keyword per line
\keyword{ graphs }
