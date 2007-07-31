\name{dyad.census}
\alias{dyad.census}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute a Holland and Leinhardt MAN Dyad Census }
\description{
  \code{dyad.census} computes a Holland and Leinhardt dyad census on the graphs of \code{dat} selected by \code{g}.
}
\usage{
dyad.census(dat, g=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more graphs. }
  \item{g}{ the elements of \code{dat} to be included; by default, all graphs are processed. }
}
\details{
Each dyad in a directed graph may be in one of four states: the null state (\eqn{a \not\leftrightarrow b}{a <-!-> b}), the complete or mutual state (\eqn{a \leftrightarrow b}{a <--> b}), and either of two asymmetric states (\eqn{a \leftarrow b}{a --> b} or \eqn{a \rightarrow b}{a <-- b}).  Holland and Leinhardt's dyad census classifies each dyad into the mutual, asymmetric, or null categories, counting the number of each within the digraph.  These counts can be used as the basis for null hypothesis tests (since their distributions are known under assumptions such as constant edge probability), or for the generation of random graphs (e.g., via the U|MAN distribution, which conditions on the numbers of mutual, asymmetric, and null dyads in each graph).  
}
\value{
A matrix whose three columns contain the counts of mutual, asymmetric, and null dyads (respectively) for each graph
}
\references{ Holland, P.W. and Leinhardt, S.  (1970).  ``A Method for Detecting Structure in Sociometric Data.''  \emph{American Journal of Sociology}, 76, 492-513. 

Wasserman, S., and Faust, K.  (1994).  ``Social Network Analysis: Methods and Applications.''  Cambridge: Cambridge University Press.}
\author{ Carter T. Butts \email{buttsc@uci.edu} }

\seealso{ \code{\link{mutuality}}, \code{\link{grecip}}, \code{\link{rguman}} \code{\link{triad.census}}, \code{\link{kcycle.census}}, \code{\link{kpath.census}} }

\examples{
#Generate a dyad census of random data with varying densities
dyad.census(rgraph(15,5,tprob=c(0.1,0.25,0.5,0.75,0.9)))
}
\keyword{ math }% __ONLY ONE__ keyword per line
\keyword{ graphs }
