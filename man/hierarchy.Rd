\name{hierarchy}
\alias{hierarchy}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute Graph Hierarchy Scores }
\description{
  \code{hierarchy} takes a graph set (\code{dat}) and returns reciprocity or Krackhardt hierarchy scores for the graphs selected by \code{g}.
}
\usage{
hierarchy(dat, g=NULL, measure=c("reciprocity", "krackhardt"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{a stack of input graphs. }
  \item{g}{ index values for the graphs to be utilized; by default, all graphs are selected. }
  \item{measure}{ one of \code{"reciprocity"} or \code{"krackhardt"}. }
}
\details{
Hierarchy measures quantify the extent of asymmetry in a structure; the greater the extent of asymmetry, the more hierarchical the structure is said to be.  (This should not be confused with how \emph{centralized} the structure is, i.e., the extent to which centralities of vertex positions are highly concentrated.)  \code{hierarchy} provides two measures (selected by the \code{measure} argument) as follows:

\enumerate{
\item \code{reciprocity}: This setting returns one minus the dyadic reciprocity for each input graph (see \code{\link{grecip}})

\item \code{krackhardt}: This setting returns the Krackhardt hierarchy score for each input graph.  The Krackhardt hierarchy is defined as the fraction of non-null dyads in the \code{\link{reachability}} graph which are asymmetric.  Thus, when no directed paths are reciprocated (e.g., in an in/outtree), Krackhardt hierarchy is equal to 1; when all such paths are reciprocated, by contrast (e.g., in a cycle or clique), the measure falls to 0. 

Hierarchy is one of four measures (\code{\link{connectedness}}, \code{\link{efficiency}}, \code{\link{hierarchy}}, and \code{\link{lubness}}) suggested by Krackhardt for summarizing hierarchical structures.  Each corresponds to one of four axioms which are necessary and sufficient for the structure in question to be an outtree; thus, the measures will be equal to 1 for a given graph iff that graph is an outtree.  Deviations from unity can be interpreted in terms of failure to satisfy one or more of the outtree conditions, information which may be useful in classifying its structural properties.
}

Note that hierarchy is inherently density-constrained: as densities climb above 0.5, the proportion of mutual dyads must (by the pigeonhole principle) increase rapidly, thereby reducing possibilities for asymmetry.  Thus, the interpretation of hierarchy scores should take density into account, particularly if density is artifactual (e.g., due to a particular dichotomization procedure).  
}
\value{
A vector of hierarchy scores
}
\references{ Krackhardt, David.  (1994).  ``Graph Theoretical Dimensions of Informal Organizations.'' In K. M. Carley and M. J. Prietula (Eds.), \emph{Computational Organization Theory}, 89-111. Hillsdale, NJ: Lawrence Erlbaum and Associates. 

Wasserman, S., and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press.}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ The four Krackhardt indices are, in general, nondegenerate for a relatively narrow band of size/density combinations (efficiency being the sole exception).  This is primarily due to their dependence on the reachability graph, which tends to become complete rapidly as size/density increase.  See Krackhardt (1994) for a useful simulation study. }

\seealso{ \code{\link{connectedness}}, \code{\link{efficiency}}, \code{\link{hierarchy}}, \code{\link{lubness}}, \code{\link{grecip}}, \code{\link{mutuality}}, \code{\link{dyad.census}} }

\examples{
#Get hierarchy scores for graphs of varying densities
hierarchy(rgraph(10,5,tprob=c(0.1,0.25,0.5,0.75,0.9)),
    measure="reciprocity")
hierarchy(rgraph(10,5,tprob=c(0.1,0.25,0.5,0.75,0.9)),
    measure="krackhardt")
}
\keyword{ math }% at least one, from doc/KEYWORDS
\keyword{ univar }% __ONLY ONE__ keyword per line
\keyword{ graphs }
