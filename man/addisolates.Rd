\name{addisolates}
\alias{addisolates}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Add Isolates to a Graph }
\description{
  Adds \code{n} isolates to the adjacency matrix (or matrices) in \code{dat}.
}
\usage{
addisolates(dat, n)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ One or more adjacency matrices }
  \item{n}{ The number of isolates to add }
}
\details{
  If \code{dat} contains more than one adjacency matrix, the \code{n} isolates are added to each member of \code{dat}.
}
\value{
   The updated graph stack.
}
\references{ Butts, C.T., and Carley, K.M.  (2001).  ``Multivariate Methods for Inter-Structural Analysis.''  CASOS Working Paper, Carnegie Mellon University. }
\author{ Carter T. Butts \email{ctb@andrew.cmu.edu} }
\note{ Isolate addition is particularly useful when computing structural distances between graphs of different orders; see the above reference for details. }

\seealso{ \code{\link{isolates}} }

\examples{

g<-rgraph(10,5)		#Produce some random graphs

dim(g)			#Get the dimensions of g

g<-addisolates(g,2)	#Add 2 isolates to each graph in g

dim(g)			#Now examine g
g

}
\keyword{ manip }
\keyword{ network }
