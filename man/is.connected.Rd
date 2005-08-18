\name{is.connected}
\alias{is.connected}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Is a Given Graph Connected? }
\description{
  Returns \code{TRUE} iff the specified graphs are connected.
}
\usage{
is.connected(g, connected = "strong", comp.dist.precomp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{ one or more input graphs. }
  \item{connected}{ definition of connectedness to use; must be one of \code{"strong"}, \code{"weak"}, \code{"unilateral"}, or \code{"recursive"}. }
  \item{comp.dist.precomp}{ a \code{\link{component.dist}} object precomputed for the graph to be analyzed (optional). }
}
\details{
  \code{is.connected} determines whether the elements of \code{g} are connected under the definition specified in \code{connected}.  (See \code{\link{component.dist}} for details.)  Since \code{is.connected} is really just a wrapper for \code{\link{component.dist}}, an object created with the latter can be supplied (via \code{comp.dist.precomp}) to speed computation.
}
\value{
\code{TRUE} iff \code{g} is connected, otherwise \code{FALSE}
}
\references{ 
West, D.B.  (1996).  \emph{Introduction to Graph Theory.}  Upper Saddle River, N.J.: Prentice Hall.
}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{component.dist}}, \code{\link{components}} }
\examples{
#Generate two graphs:
g1<-rgraph(10,tp=0.1)
g2<-rgraph(10)

#Check for connectedness
is.connected(g1)  #Probably not
is.connected(g2)  #Probably so
}
\keyword{ graphs }% at least one, from doc/KEYWORDS
\keyword{ logic }% __ONLY ONE__ keyword per line
