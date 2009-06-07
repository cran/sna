\name{add.isolates}
\alias{add.isolates}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Add Isolates to a Graph }
\description{
  Adds \code{n} isolates to the graph (or graphs) in \code{dat}.
}
\usage{
add.isolates(dat, n, return.as.edgelist = FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{n}{ the number of isolates to add. }
  \item{return.as.edgelist}{ logical; should the input graph be returned as an edgelist (rather than an adjacency matrix)?}
}
\details{
  If \code{dat} contains more than one graph, the \code{n} isolates are added to each member of \code{dat}.
}
\value{
   The updated graph(s).
}
\references{ Butts, C.T., and Carley, K.M.  (2001).  ``Multivariate Methods for Inter-Structural Analysis.''  CASOS Working Paper, Carnegie Mellon University. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ Isolate addition is particularly useful when computing structural distances between graphs of different orders; see the above reference for details. }

\seealso{ \code{\link{isolates}} }

\examples{

g<-rgraph(10,5)		#Produce some random graphs

dim(g)			#Get the dimensions of g

g<-add.isolates(g,2)	#Add 2 isolates to each graph in g

dim(g)			#Now examine g
g

}
\keyword{ manip }
\keyword{ math }
\keyword{ graphs }

