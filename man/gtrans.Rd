\name{gtrans}
\alias{gtrans}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the Transitivity of an Input Graph or Graph Stack }
\description{
  \code{gtrans} returns the fraction of triads which are transitive, i.e., for which \eqn{A^2=A}{A^2=A}.  Triads involving missing values are omitted from the analysis.
}
\usage{
gtrans(dat, g=NULL, diag=FALSE, mode="digraph")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ A graph or graph stack }
  \item{g}{ A vector indicating the graphs which are to be analyzed; by default, all graphs are analyzed }
  \item{diag}{ A boolean indicating whether or not diagonal entries (loops) are to be taken as valid data }
  \item{mode}{ ``digraph'' if directed triads are sought, or else ``graph'' }
}
\details{

}
\value{
A vector of transitivity scores
}
\references{ Holland, P.W., and Leinhardt, S.  (1972).  ``Some Evidence on the Transitivity of Positive Interpersonal Sentiment.'' \emph{American Journal of Sociology,} 72, 1205-1209.

Wasserman, S., and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications}.  Cambridge: Cambridge University Press. }
\author{ Carter T. Butts \email{ctb@andrew.cmu.edu} }
%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{cugtest}} }

\examples{
#Draw some random graphs
g<-rgraph(5,10)

#Find transitivity scores
gtrans(g)
}
\keyword{ univar }%-- one or more ...
\keyword{ network }