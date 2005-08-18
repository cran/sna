\name{sr2css}
\alias{sr2css}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Convert a Row-wise Self-Report Matrix to a CSS Matrix with Missing Observations }
\description{
  Given a matrix in which the ith row corresponds to i's reported relations, \code{sr2css} creates a graph stack in which each element represents a CSS slice with missing observations.
}
\usage{
sr2css(net)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{net}{ an adjacency matrix. }
}
\details{
 A cognitive social structure (CSS) is an nxnxn array in which the ith matrix corresponds to the ith actor's perception of the entire network.  Here, we take a conventional self-report data structure and put it in CSS format for routines (such as \code{\link{bbnam}}) which require this.
}
\value{
An array (graph stack) containing the CSS
}
\references{ Krackhardt, D.  (1987).  \emph{Cognitive Social Structures}, 9, 109-134.  }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ A row-wise self-report matrix doesn't contain a great deal of data, and the data in question is certainly not an ignorable sample of the individual's CSS for most purposes.  The provision of this routine should not be perceived as license to substitute SR for CSS data at will. }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

%\seealso{ }

\examples{
#Start with some random reports
g<-rgraph(10)

#Transform to CSS format
c<-sr2css(g)
}
\keyword{ manip }%-- one or more ...
\keyword{ math }
\keyword{ graphs }
