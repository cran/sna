\name{grecip}
\alias{grecip}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the Reciprocity of an Input Graph or Graph Stack }
\description{
  \code{grecip} calculates the dyadic reciprocity of the elements of \code{dat} selected by \code{g}.
}
\usage{
grecip(dat, g = NULL, measure = c("dyadic", "edgewise"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ An adjacency matrix, or stack thereof }
  \item{g}{ A vector indicating which matrices to evaluate (optional) }
  \item{measure}{ One of \code{"dyadic"} (default) or \code{"edgewise"} }
}
\details{
  The dyadic reciprocity of a graph is the proportion of dyads which are symmetric; this is computed and returned by \code{grecip} for the graphs indicated.  Note that the dyadic reciprocity is distinct from the \emph{edgewise} or \emph{tie reciprocity}, which is the proportion of \emph{edges} which are reciprocated.  This latter form may be obtained by setting \code{measure="edgewise"}.
}
\value{
  The graph reciprocity value(s)
}
\references{ Wasserman, S., and Faust, K.  (1994).  ``Social Network Analysis: Methods and Applications.''  Cambridge: Cambridge University Press. }
\author{ Carter T. Butts \email{ctb@andrew.cmu.edu} }
%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{mutuality}}, \code{\link{symmetrize}} }

\examples{
#Calculate the dyadic reciprocity scores for some random graphs
grecip(rgraph(10,5))
}
\keyword{ math }%-- one or more ...
\keyword{ univar }%-- one or more ...
