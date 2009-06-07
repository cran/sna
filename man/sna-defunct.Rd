\name{sna-defunct}
\alias{addisolates}
\alias{addisolates-defunct}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Defunct sna Objects }
\description{
  These objects have been removed from \code{sna}, and should no longer be used.
}
\usage{
addisolates(dat, n)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ One or more adjacency matrices }
  \item{n}{ The number of isolates to add }
}
%\details{
%  If \code{dat} contains more than one adjacency matrix, the \code{n} isolates %are added to each member of \code{dat}.
%}
%\value{
%   The updated graph stack.
%}
%\references{ Butts, C.T., and Carley, K.M.  (2001).  \dQuote{Multivariate Methods for Inter-Structural Analysis.}  CASOS Working Paper, Carnegie Mellon University. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ Isolate addition is particularly useful when computing structural distances between graphs of different orders; see the above reference for details. }

%\seealso{ \code{\link{isolates}} }
%
%\examples{
%
%g<-rgraph(10,5)		#Produce some random graphs
%
%dim(g)			#Get the dimensions of g
%
%\dontrun{g<-addisolates(g,2)}	#Add 2 isolates to each graph in g
%
%dim(g)			#Now examine g
%g
%
%}
\keyword{ internal }
%\keyword{ math }
