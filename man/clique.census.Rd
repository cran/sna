\name{clique.census}
\Rdversion{1.1}
\alias{clique.census}
\alias{cliques_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Cycle Census Information
}
\description{
\code{clique.census} computes clique census statistics on one or more input graphs.  In addition to aggregate counts of maximal cliques, results may be disaggregated by vertex and co-membership information may be computed.
}
\usage{
clique.census(dat, mode = "digraph", tabulate.by.vertex = TRUE,
    clique.comembership = c("none", "sum", "bysize"), enumerate = TRUE,
    na.omit = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
one or more input graphs.
}
  \item{mode}{
\code{"digraph"} for directed graphs, or \code{"graph"} for undirected graphs.
}
  \item{tabulate.by.vertex}{
logical; should maximal clique counts be tabulated by vertex? 
}
  \item{clique.comembership}{
the type of clique co-membership information to be tabulated, if any.  \code{"sum"} returns a vertex by vertex matrix of clique co-membership counts; these are disaggregated by clique size if \code{"bysize"} is used.  If \code{"none"} is given, no co-membership information is computed.
}
  \item{enumerate}{
logical; should an enumeration of all maximal cliques be returned?
}
  \item{na.omit}{
logical; should missing edges be omitted?
}
}
\details{
A (maximal) \emph{clique} is a maximal set of mutually adjacenct vertices.  Cliques are important for their role as cohesive subgroups, but show up in many other contexts as well.

A \emph{subgraph census statistic} is a function which, for any given graph and subgraph, gives the number of copies of the latter contained in the former.  A collection of subgraph census statistics is referred to as a \emph{subgraph census}; widely used examples include the dyad and triad censuses, implemented in \code{sna} by the \code{\link{dyad.census}} and \code{\link{triad.census}} functions (respectively).  Likewise, \code{kpath.census} and \code{kcycle.census} compute a range of census statistics related to \eqn{k}-paths and \eqn{k}-cycles.   \code{clique.census} provides similar functionality for the census of maximal cliques, including:

\itemize{
\item Aggregate counts of maximal cliques by size. 
\item Counts of cliques to which each vertex belongs (when \code{tabulate.byvertex==TRUE}). 
\item Counts of clique co-memberships, potentially disaggregated by size (when the appropriate co-membership argument is set to \code{bylength}). 
}

These calculations are intrinsically expensive (clique enumeration is NP hard in the general case), and users should be aware that computing the census can be impractical on large graphs (unless they are very sparse).  On the other hand, the algorithm employed here (a variant of Makino and Uno (2004)) is generally fast enough to suport enumeration for even dense graphs of several hundred vertices on a typical desktop computer.

Calling this function with \code{mode=="digraph"}, forces and initial symmetrization step, which can be avoided with \code{mode=="graph"}.  While incorrectly employing the default is harmless (except for the relatively small cost of verifying symmetry), setting \code{mode=="graph"} incorrectly may result in problematic behavior.  When in doubt, stick with the default.
}
\value{
  A list with the following elements:
  \item{clique.count }{If \code{tabulate.byvertex==FALSE}, a vector of aggregate counts by clique size.  Otherwise, a matrix whose first column is a vector of aggregate clique counts, and whose succeeding columns contain vectors of clique counts for each vertex.}
  \item{clique.comemb }{If \code{clique.comembership!="none"}, a matrix or array containing co-membership in cliques by vertex pairs.  If \code{clique.comembership=="sum"}, only a matrix of co-memberships is returned; if \code{bysize} is used, however, co-memberships are returned in a \code{maxsize} by \eqn{n} by \eqn{n} array whose \eqn{i,j,k}th cell is the number of cliques of size \eqn{i} containing \code{j} and \code{k} (with \code{maxsize} being the size of the largest maximal clique).}
  \item{cliques }{If \code{enumerate=TRUE}, a list of length equal to the maximum clique size, each element of which is in turn a list of all cliques of corresponding size (given as vectors of vertices).}
}
\references{
Wasserman, S. and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press.

Makino, K. and Uno, T.  (2004.)  \dQuote{New Algorithms for Enumerating All Maximal Cliques.}  In T. Hagerup and J. Katajainen (eds.), \emph{SWAT 2004}, LNCS 3111, 260-272.  Berlin: Springer-Verlag.
}
\author{
Carter T. Butts \email{buttsc@uci.edu}
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning }{ The computational cost of calculating cliques grows very sharply in size and network density.  It is possible that the expected completion time for your calculation may exceed your life expectancy (and those of subsequent generations). }

\seealso{
\code{\link{dyad.census}}, \code{\link{triad.census}}, \code{\link{kcycle.census}}, \code{\link{kpath.census}}
}
\examples{
#Generate a fairly dense random graph
g<-rgraph(25)

#Obtain cliques by vertex, with co-membership by size
cc<-clique.census(g,clique.comembership="bysize")
cc$clique.count                             #Examine clique counts
cc$clique.comemb[1,,]                       #Isolate co-membership is trivial
cc$clique.comemb[2,,]                       #Co-membership for 2-cliques
cc$clique.comemb[3,,]                       #Co-membership for 3-cliques
cc$cliques                                  #Enumerate the cliques
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
\keyword{ math }% __ONLY ONE__ keyword per line
