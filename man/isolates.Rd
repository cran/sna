\name{isolates}
\alias{isolates}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ List the Isolates in a Graph or Graph Stack }
\description{
  Returns a list of the isolates in the graph or graph set given by \code{dat}.
}
\usage{
isolates(dat, diag=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{diag}{ boolean indicating whether adjacency matrix diagonals (i.e., loops) contain meaningful data. }
}
%\details{
%}
\value{
  A vector containing the isolates, or a list of vectors if more than one graph was specified
}
\references{ Wasserman, S., and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press.

West, D.B.  (1996).  \emph{Introduction to Graph Theory}.  Upper Saddle River, NJ: Prentice Hall.}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{is.isolate}}, \code{\link{add.isolates}} }

\examples{
#Generate a test graph
g<-rgraph(20)
g[,4]<-0          #Create an isolate
g[4,]<-0

#List the isolates
isolates(g)
}
\keyword{ math }%-- one or more ...
\keyword{ graphs }
 
