\name{gtrans}
\alias{gtrans}
\alias{transitivity_R}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the Transitivity of an Input Graph or Graph Stack }
\description{
  \code{gtrans} returns the transitivity of the elements of \code{dat} selected by \code{g}, using the definition of \code{measure}.  Triads involving missing values are omitted from the analysis.
}
\usage{
gtrans(dat, g=NULL, diag=FALSE, mode="digraph", measure = c("weak", 
    "strong", "weakcensus", "strongcensus", "rank", "correlation"), 
    use.adjacency = TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ a collection of input graphs. }
  \item{g}{ a vector indicating the graphs which are to be analyzed; by default, all graphs are analyzed. }
  \item{diag}{ a boolean indicating whether or not diagonal entries (loops) are to be taken as valid data. }
  \item{mode}{ \code{"digraph"} if directed triads are sought, or else \code{"graph"}. }
  \item{measure}{ one of \code{"weak"} (default), \code{"strong"}, \code{"weakcensus"}, \code{"strongcensus"}, \code{"rank"}, or \code{"correlation"}. }
  \item{use.adjacency}{ logical; should adjacency matrices (versus sparse graph methods) be used in the transitivity computation?}
}
\details{
Transitivity is a triadic, algebraic structural constraint.  In its weak form, the transitive constraint corresponds to \eqn{a \rightarrow b \rightarrow c \Rightarrow a \rightarrow c}{a->b->c => a->c}.  In the corresponding strong form, the constraint is \eqn{a \rightarrow b \rightarrow c \Leftrightarrow a \rightarrow c}{a->b->c <=> a->c}.  (Note that the weak form is that most commonly employed.)  Where \code{measure=="weak"}, the fraction of potentially intransitive triads obeying the weak condition is returned.  With the \code{measure=="weakcensus"} setting, by contrast, the total \emph{number} of transitive triads is computed.  The \code{strong} versions of the measures are similar to the above, save in that the set of all triads is considered (since all are ``at risk'' for intransitivity).  

Note that where missing values prevent the assessment of whether a triple is transitive, that triple is omitted.

Generalizations of transitivity to valued graphs are numerous.  The above strong and weak forms ignore edge values, treating any non-zero edge as present.  Two additional notions of transitivity are also supported valued data.  The \code{"rank"} condition treads an \eqn{(i, j, k)} triple as transitive if the value of the \eqn{(i,k)} directed dyad is greater than or equal to the minimum of the values of the \eqn{(i,j)} and \eqn{(j,k)} dyads.  The \code{"correlation"} option implements the correlation transitivity of David Dekker, which is defined as the matrix correlation of the valued adjacency matrix \eqn{A} with its second power (i.e., \eqn{A^2}), omitting diagonal entries where inapplicable.

Note that the base forms of transitivity can be calculated using either matrix multiplication or sparse graph methods.  For very large, sparse graphs, the sparse graph method (which can be forced by \code{use.adjacency=FALSE}) may be preferred.  The latter provides much better scaling, but is significantly slower for networks of typical size due to the overhead involved (and R's highly optimized matrix operations).  Where \code{use.adjacency} is set to \code{TRUE}, \code{gtrans} will attempt some simple heuristics to determine if the edgelist method should be used instead (and will do so if indicated).  These heuristics depend on recognition of the input data type, and hence may behave slightly differently depending on the form in which \code{dat} is given.  Note that the rank measure can at present be calculated only via sparse graph methods, and the correlation measure only by adjacency matrices.  For these measures, the \code{use.adjacency} argument is ignored.
}
\value{
A vector of transitivity scores
}
\references{ Holland, P.W., and Leinhardt, S.  (1972).  ``Some Evidence on the Transitivity of Positive Interpersonal Sentiment.'' \emph{American Journal of Sociology,} 72, 1205-1209.

Wasserman, S., and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications}.  Cambridge: Cambridge University Press. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ In version 0.3, the strong form was the default (and only) option for \code{gtrans}. }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{triad.classify}}, \code{\link{cugtest}} }

\examples{
#Draw some random graphs
g<-rgraph(5,10)

#Find transitivity scores
gtrans(g)
}
\keyword{ algebra }%-- one or more ...
\keyword{ graphs }
