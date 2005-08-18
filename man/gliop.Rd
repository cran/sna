\name{gliop}
\alias{gliop}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Return a Binary Operation on GLI Values Computed on Two Graphs }
\description{
  \code{gliop} is a wrapper which allows for an arbitrary binary operation on GLIs to be treated as a single call.  This is particularly useful for test routines such as \code{\link{cugtest}} and \code{\link{qaptest}}.
}
\usage{
gliop(dat, GFUN, OP="-", g1=1, g2=2, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ a collection of graphs. }
  \item{GFUN}{ a function taking single graphs as input. }
  \item{OP}{ the operator to use on the output of \code{GFUN}. }
  \item{g1}{ the index of the first input graph. }
  \item{g2}{ the index of the second input graph. }
  \item{\dots}{ Additional arguments to \code{GFUN} }
}
\details{
  \code{gliop} operates by evaluating \code{GFUN} on the graphs indexed by \code{g1} and \code{g2} and returning the result of \code{OP} as applied to the \code{GFUN} output.
}
\value{
\code{OP(GFUN(dat[g1, , ],...),GFUN(dat[g2, , ],...))}
}
\references{ Anderson, B.S.; Butts, C.T.; and Carley, K.M. (1999). ``The Interaction of Size and Density with Graph-Level Indices.'' \emph{Social Networks}, 21(3), 239-267. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ If the output of \code{GFUN} is not sufficiently well-behaved, undefined behavior may occur.  Common sense is advised. }

 %~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{cugtest}}, \code{\link{qaptest}} }

\examples{
#Draw two random graphs
g<-rgraph(10,2,tprob=c(0.2,0.5))

#What is their difference in density?
gliop(g,gden,"-",1,2)
}
\keyword{ utilities }%-- one or more ...
\keyword{ math }
\keyword{ graphs }
