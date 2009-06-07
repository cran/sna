\name{grecip}
\alias{grecip}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the Reciprocity of an Input Graph or Graph Stack }
\description{
  \code{grecip} calculates the dyadic reciprocity of the elements of \code{dat} selected by \code{g}.
}
\usage{
grecip(dat, g = NULL, measure = c("dyadic", "dyadic.nonnull",
    "edgewise", "edgewise.lrr"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{g}{ a vector indicating which graphs to evaluate (optional). }
  \item{measure}{ one of \code{"dyadic"} (default), \code{"dyadic.nonnull"}, \code{"edgewise"}, or \code{"edgewise.lrr"}. }
}
\details{
  The dyadic reciprocity of a graph is the proportion of dyads which are symmetric; this is computed and returned by \code{grecip} for the graphs indicated.  (\code{dyadic.nonnull} returns the ratio of mutuals to non-null dyads.)  Note that the dyadic reciprocity is distinct from the \emph{edgewise} or \emph{tie reciprocity}, which is the proportion of \emph{edges} which are reciprocated.  This latter form may be obtained by setting \code{measure="edgewise"}.  Setting \code{measure="edgewise.lrr"} returns the log of the ratio of the edgewise reciprocity to the density; this is measure (called \eqn{r_4} by Butts (2008)) can be interpreted as the relative log-odds of an edge given a reciprocation, versus the baseline probability of an edge.
}
\value{
  The graph reciprocity value(s)
}
\references{ Wasserman, S., and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press. 

Butts, C.T.  (2008).  \dQuote{Social Networks: A Methodological Introduction.} \emph{Asian Journal of Social Psychology,} 11(1), 13-41. 
}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{mutuality}}, \code{\link{symmetrize}} }

\examples{
#Calculate the dyadic reciprocity scores for some random graphs
grecip(rgraph(10,5))
}
\keyword{ math }%-- one or more ...
\keyword{ univar }%-- one or more ...
\keyword{ graphs }
