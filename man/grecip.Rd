\name{grecip}
\alias{grecip}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the Reciprocity of an Input Graph or Graph Stack }
\description{
  \code{grecip} calculates the dyadic reciprocity of the elements of \code{dat} selected by \code{g}.
}
\usage{
grecip(dat, g = NULL, measure = c("dyadic", "dyadic.nonnull",
    "edgewise", "edgewise.lrr", "correlation"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{g}{ a vector indicating which graphs to evaluate (optional). }
  \item{measure}{ one of \code{"dyadic"} (default), \code{"dyadic.nonnull"}, \code{"edgewise"}, \code{"edgewise.lrr"}, or \code{"correlation"}. }
}
\details{
The dyadic reciprocity of a graph is the proportion of dyads which are symmetric; this is computed and returned by \code{grecip} for the graphs indicated.  (\code{dyadic.nonnull} returns the ratio of mutuals to non-null dyads.)  Note that the dyadic reciprocity is distinct from the \emph{edgewise} or \emph{tie reciprocity}, which is the proportion of \emph{edges} which are reciprocated.  This latter form may be obtained by setting \code{measure="edgewise"}.  Setting \code{measure="edgewise.lrr"} returns the log of the ratio of the edgewise reciprocity to the density; this is measure (called \eqn{r_4} by Butts (2008)) can be interpreted as the relative log-odds of an edge given a reciprocation, versus the baseline probability of an edge.  Finally, \code{measure="correlation"} returns the correlation between within-dyad edge values, where this is defined by
  \deqn{\frac{2\sum_{\{i,j\}} (Y_{ij}-\mu_G)(Y_{ji}-\mu_G)}{(2N_d-1) \sigma^2_G}}{%
    2*(sum_{i,j} (Y_ij-mu_G)*(Y_ji-mu_G))/((2*N_d-1)*sigma^2_G)
  }
with \eqn{Y} being the graph adjacency matrix, \eqn{\mu_G}{mu_G} being the mean non-loop edge value, \eqn{\sigma^2_G}{sigma^2_G} being the variance of non-loop edge values, and \eqn{N_d} being the number of dyads.  (Note that this quantity is unaffected by dyad orientation.)  The correlation measure may be interpreted as the net tendency for edges of similar relative value (with respect to the mean edge value) to occur within the same dyads.  For dichotomous data, adjacencies are interpreted as having values of 0 (no edge present) or 1 (edge present), but edge values are used where supplied.  In cases where all edge values are identical (e.g., the complete or empty graph), the correlation reciprocity is taken to be 1 by definition.

Note that \code{grecip} calculates values based on non-missing data; dyads containing missing data are removed from consideration when calculating reciprocity scores (except for the correlation measure, which uses non-missing edges within missing dyads when calculating the graph mean and variance).
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
