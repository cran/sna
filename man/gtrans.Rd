\name{gtrans}
\alias{gtrans}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the Transitivity of an Input Graph or Graph Stack }
\description{
  \code{gtrans} returns the transitivity of the elements of \code{dat} selected by \code{g}, using the definition of \code{measure}.  Triads involving missing values are omitted from the analysis.
}
\usage{
gtrans(dat, g=NULL, diag=FALSE, mode="digraph", measure = c("weak", 
    "strong", "weakcensus", "strongcensus"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ A graph or graph stack }
  \item{g}{ A vector indicating the graphs which are to be analyzed; by default, all graphs are analyzed }
  \item{diag}{ A boolean indicating whether or not diagonal entries (loops) are to be taken as valid data }
  \item{mode}{ \code{"digraph"} if directed triads are sought, or else \code{"graph"} }
  \item{measure}{ One of \code{"weak"} (default), \code{"strong"}, \code{"weakcensus"}, or \code{"strongcensus"} }
}
\details{
Transitivity is a triadic, algebraic structural constraint.  In its weak form, the transitive constraint corresponds to \eqn{a \rightarrow b \rightarrow c \Rightarrow a \rightarrow c}{a->b->c => a->c}.  In the corresponding strong form, the constraint is \eqn{a \rightarrow b \rightarrow c \Leftrightarrow a \rightarrow c}{a->b->c <=> a->c}.  (Note that the weak form is that most commonly employed.)  Where \code{measure=="weak"}, the fraction of potentially intransitive triads obeying the weak condition is returned.  With the \code{measure=="weakcensus"} setting, by contrast, the total \emph{number} of transitive triads is computed.  The \code{strong} versions of the measures are similar to the above, save in that the set of all triads is considered (since all are ``at risk'' for intransitivity).
}
\value{
A vector of transitivity scores
}
\references{ Holland, P.W., and Leinhardt, S.  (1972).  ``Some Evidence on the Transitivity of Positive Interpersonal Sentiment.'' \emph{American Journal of Sociology,} 72, 1205-1209.

Wasserman, S., and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications}.  Cambridge: Cambridge University Press. }
\author{ Carter T. Butts \email{ctb@andrew.cmu.edu} }
\note{ In version 0.3, the strong form was the default (and only) option for \code{gtrans}. }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{triad.classify}}, \code{\link{cugtest}} }

\examples{
#Draw some random graphs
g<-rgraph(5,10)

#Find transitivity scores
gtrans(g)
}
\keyword{ network}
\keyword{ algebra }%-- one or more ...
