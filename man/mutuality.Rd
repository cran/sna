\name{mutuality}
\alias{mutuality}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Find the Mutuality of a Graph }
\description{
  Returns the mutuality scores of the graphs indicated by \code{g} in \code{dat}.
}
\usage{
mutuality(dat, g=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{g}{ a vector indicating which elements of \code{dat} should be analyzed; by default, all graphs are included. }
}
\details{
  The mutuality of a digraph G is defined as the number of complete dyads (i.e., i<->j) within G.  (Compare this to dyadic reciprocity, the fraction of dyads within G which are symmetric.)  Mutuality is commonly employed as a measure of reciprocal tendency within the p* literature; although mutuality can be very hard to interpret in practice, it is much better behaved than many alternative measures.
}
\value{
One or more mutuality scores
}
\references{ Moreno, J.L., and Jennings, H.H.  (1938).  ``Statistics of Social Configurations.''  \emph{Sociometry}, 1, 342-374. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{grecip}} }

\examples{
#Create some random graphs
g<-rgraph(15,3)

#Get mutuality and reciprocity scores
mutuality(g)
grecip(g)         #Compare with mutuality
}
\keyword{ univar }%-- one or more ...
\keyword{ math }
\keyword{ graphs }
