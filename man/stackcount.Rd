\name{stackcount}
\alias{stackcount}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ How Many Graphs are in a Graph Stack? }
\description{
Returns the number of graphs in the stack provided by \code{d}.
}
\usage{
stackcount(d)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{d}{ a graph or graph stack. }
}
%\details{
%}
\value{
The number of graphs in \code{d}
}
%\references{  }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

 %~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{nties}} }

\examples{
stackcount(rgraph(4,8))==8
}
\keyword{ utilities }%-- one or more ...
\keyword{ math }
