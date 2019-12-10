\name{triad.census}
\alias{triad.census}
\alias{triad_census_R}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the Davis and Leinhardt Triad Census }
\description{
  \code{triad.census} returns the Davis and Leinhardt triad census of the elements of \code{dat} indicated by \code{g}.
}
\usage{
triad.census(dat, g=NULL, mode = c("digraph", "graph"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ a graph or graph stack. }
  \item{g}{ the elements of \code{dat} to process. }
  \item{mode}{string indicating the directedness of edges; \code{"digraph"} implies a directed structure, whereas \code{"graph"} implies an undirected structure.}
}
\details{
The Davis and Leinhardt triad census consists of a classification of all directed triads into one of 16 different categories; the resulting distribution can be compared against various null models to test for the presence of configural biases (e.g., transitivity bias).  \code{triad.census} is a front end for the \code{\link{triad.classify}} routine, performing the classification for all triads within the selected graphs.  The results are placed in the order indicated by the column names; this is the same order as presented in the \code{\link{triad.classify}} documentation, to which the reader is referred for additional details.

In the undirected case, the triad census reduces to four states (based on the number of edges in each triad.  Where \code{mode=="graph"}, this is returned instead.

Compare \code{\link{triad.census}} to \code{\link{dyad.census}}, the dyadic equivalent.
}
\value{
A matrix whose 16 columns contain the counts of triads by class for each graph, in the directed case.  In the undirected case, only 4 columns are used.
}
\references{ Davis, J.A. and Leinhardt, S.  (1972).  ``The Structure of Positive Interpersonal Relations in Small Groups.''  In J. Berger (Ed.), \emph{Sociological Theories in Progress, Volume 2}, 218-251.  Boston: Houghton Mifflin.

Wasserman, S., and Faust, K.  (1994).  ``Social Network Analysis: Methods and Applications.''  Cambridge: Cambridge University Press. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ }

\seealso{ \code{\link{triad.classify}}, \code{\link{dyad.census}}, \code{\link{kcycle.census}}, \code{\link{kpath.census}}, \code{\link{gtrans}} }

\section{Warning }{Valued data may cause strange behavior with this routine.  Dichotomize the data first.}

\examples{
#Generate a triad census of random data with varying densities
triad.census(rgraph(15,5,tprob=c(0.1,0.25,0.5,0.75,0.9)))
}
\keyword{ math }% __ONLY ONE__ keyword per line
\keyword{ graphs }
