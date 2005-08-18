\name{nties}
\alias{nties}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Find the Number of Possible Ties in a Given Graph or Graph Stack }
\description{
\code{nties} returns the number of possible edges in each element of \code{dat}, given \code{mode} and \code{diag}.
}
\usage{
nties(dat, mode="digraph", diag=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ a graph or set thereof. }
  \item{mode}{ one of ``digraph'', ``graph'', and ``hgraph''. }
  \item{diag}{ a boolean indicating whether or not diagonal entries (loops) should be treated as valid data; ignored for hypergraphic (``hgraph'') data. }
}
\details{
\code{nties} is used primarily to automate maximum edge counts for use with normalization routines.
}
\value{
The number of possible edges, or a vector of the same
}
%\references{  }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ For two-mode (hypergraphic) data, the value returned isn't technically the number of edges per se, but rather the number of edge memberships. }

 %~Make other sections like WARNING with \section{WARNING }{....} ~

%\seealso{ }

\examples{
#How many possible edges in a loopless digraph of order 15?
nties(rgraph(15),diag=FALSE)
}
\keyword{ univar }%-- one or more ...
\keyword{ math }
\keyword{ graphs }
