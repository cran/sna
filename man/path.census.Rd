\name{kpath.census}
\alias{kpath.census}
\alias{kcycle.census}
\alias{cycleCensus_R}
\alias{pathCensus_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute Path or Cycle Census Information }
\description{
  \code{kpath.census} and \code{kcycle.census} compute \eqn{k}-path or \eqn{k}-cycle census statistics (respectively) on one or more input graphs.  In addition to aggregate counts of paths or cycles, results may be disaggregated by vertex and co-membership information may be computed.
}
\usage{
kcycle.census(dat, maxlen = 3, mode = "digraph", 
    tabulate.by.vertex = TRUE, cycle.comembership = c("none", "sum",
    "bylength"))

kpath.census(dat, maxlen = 3, mode = "digraph", 
    tabulate.by.vertex = TRUE, path.comembership = c("none", "sum",
    "bylength"), dyadic.tabulation = c("none", "sum", "bylength"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cycle.comembership}{ the type of cycle co-membership information to be tabulated, if any.  \code{"sum"} returns a vertex by vertex matrix of cycle co-membership counts; these are disaggregated by cycle length if \code{"bylength"} is used.  If \code{"none"} is given, no co-membership information is computed.}
  \item{dat}{ one or more input graphs. }
  \item{maxlen}{ the maximum path/cycle length to evaluate. }
  \item{mode}{ \code{"digraph"} for directed graphs, or \code{"graph"} for undirected graphs.}
  \item{tabulate.by.vertex}{ logical; should path or cycle incidence counts be tabulated by vertex? }
  \item{path.comembership}{ as per \code{cycle.comembership}, for paths rather than cycles.  }
  \item{dyadic.tabulation}{ the type of dyadic path count information to be tabulated, if any.  \code{"sum"} returns a vertex by vertex matrix of source/destination path counts, while \code{"bylength"} disaggregates these counts by path length.  Selecting \code{"none"} disables this computation. }
}
\details{
There are several equivalent characterizations of paths and cycles, of which the following is one example.  For an arbitrary graph \eqn{G}, a \emph{path} is a sequence of distinct vertices \eqn{v_1, v_2, \ldots, v_n}{v_1, v_2, .... v_n} and included edges such that \eqn{v_i} is adjacent to \eqn{v_{i+1}} for all \eqn{i \in 1, 2, \ldots, n-1}{i in 1, 2, ... k-1} via the pair's included edge. (Contrast this with a \emph{walk}, in which edges and/or vertices may be repeated.)  A \emph{cycle} is the union of a path and an edge making \eqn{v_n} adjacent to \eqn{v_i}.  \eqn{k}-paths and \eqn{k}-cycles are respective paths and cycles having \eqn{k} edges (in the former case) or \eqn{k} vertices (in the latter).  The above definitions may be applied in both directed and undirected contexts, by substituting the appropriate notion of adjacency.  (Note that authors do not always employ the same terminology for these concepts, especially in older texts -- it is wise to verify the definitions being used in any particular context.)

A \emph{subgraph census statistic} is a function which, for any given graph and subgraph, gives the number of copies of the latter contained in the former.  A collection of subgraph census statistics is referred to as a \emph{subgraph census}; widely used examples include the dyad and triad censuses, implemented in \code{sna} by the \code{\link{dyad.census}} and \code{\link{triad.census}} functions (respectively).  \code{kpath.census} and \code{kcycle.census} compute a range of census statistics related to \eqn{k}-paths and \eqn{k}-cycles, including:

\itemize{
\item Aggregate counts of paths/cycles by length (i.e., \eqn{k}). 
\item Counts of paths/cycles to which each vertex belongs (when \code{tabulate.byvertex==TRUE}). 
\item Counts of path/cycle co-memberships, potentially disaggregated by length (when the appropriate co-membership argument is set to \code{bylength}). 
\item For \code{path.census}, counts of the total number of paths from each vertex to each other vertex, possibly disaggregated by length (if \code{dyadic.tabulation=="bylength"}).
}

The length of the maximum-length path/cycle to compute is given by \code{maxlen}.  These calculations are intrinsically expensive (path/cycle computation is NP complete in the general case), and users should hence be wary when increasing \code{maxlen}.  On the other hand, it may be possible to enumerate even long paths or cycles on a very sparse graph; scaling is approximately \eqn{c^k}, where \eqn{k} is given by \code{maxlen} and \eqn{c} is the size of the largest dense cluster.  

The paths or cycles computed by this function are directed if \code{mode=="digraph"}, or undirected if \code{mode=="graph"}.  Failing to set \code{mode} correctly may result in problematic behavior.
}
\value{
  For \code{kpath.census}, a list with the following elements:
  \item{path.count }{If \code{tabulate.byvertex==FALSE}, a vector of aggregate counts by path length.  Otherwise, a matrix whose first column is a vector of aggregate path counts, and whose succeeding columns contain vectors of path counts for each vertex.}
  \item{path.comemb }{If \code{path.comembership!="none"}, a matrix or array containing co-membership in paths by vertex pairs.  If \code{path.comembership=="sum"}, only a matrix of co-memberships is returned; if \code{bylength} is used, however, co-memberships are returned in a \code{maxlen} by \eqn{n} by \eqn{n} array whose \eqn{i,j,k}th cell is the number of paths of length \eqn{i} containing \code{j} and \code{k}.}
  \item{paths.bydyad }{If \code{dyadic.tabulation!="none"}, a matrix or array containing the number of paths originating at a particular vertex and terminating.  If \code{bylength} is used, dyadic path counts are supplied via a \code{maxlen} by \eqn{n} by \eqn{n} array whose \eqn{i,j,k}th cell is the number of paths of length \eqn{i} starting at \code{j} and ending with \code{k}.  If \code{sum} is used instead, only a matrix whose \eqn{i,j} cell contains the total number of paths from \eqn{i} to \eqn{j} is returned.}

  For \code{kcycle.census}, a similar list:
  \item{cycle.count }{If \code{tabulate.byvertex==FALSE}, a vector of aggregate counts by cycle length.  Otherwise, a matrix whose first column is a vector of aggregate cycle counts, and whose succeeding columns contain vectors of cycle counts for each vertex.}
  \item{cycle.comemb }{If \code{cycle.comembership!="none"}, a matrix or array containing co-membership in cycles by vertex pairs.  If \code{cycle.comembership=="sum"}, only a matrix of co-memberships is returned; if \code{bylength} is used, however, co-memberships are returned in a \code{maxlen} by \eqn{n} by \eqn{n} array whose \eqn{i,j,k}th cell is the number of cycles of length \eqn{i} containing \code{j} and \code{k}.}
  
}
\references{ 
Butts, C.T.  (2006).  \dQuote{Cycle Census Statistics for Exponential Random Graph Models.}  IMBS Technical Report MBS 06-05, University of California, Irvine.

West, D.B.  (1996).  \emph{Introduction to Graph Theory.}  Upper Saddle River, N.J.: Prentice Hall.
}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ 
%
\section{Warning }{ The computational cost of calculating paths and cycles grows very sharply in both \code{maxlen} and network density.  Be wary of setting \code{maxlen} greater than 5-6, unless you know what you are doing.  Otherwise, the expected completion time for your calculation may exceed your life expectancy (and those of subsequent generations). }
%}
\seealso{ \code{\link{dyad.census}}, \code{\link{triad.census}}, \code{\link{clique.census}}, \code{\link{geodist}} }
\examples{
g<-rgraph(20,tp=1.5/19)

#Obtain paths by vertex, with dyadic path counts
pc<-kpath.census(g,maxlen=5,dyadic.tabulation="sum")
pc$path.count                                 #Examine path counts
pc$paths.bydyad                               #Examine dyadic paths

#Obtain aggregate cycle counts, with co-membership by length
cc<-kcycle.census(g,maxlen=5,tabulate.by.vertex=FALSE,
    cycle.comembership="bylength")
cc$cycle.count                             #Examine cycle counts
cc$cycle.comemb[1,,]                       #Co-membership for 2-cycles
cc$cycle.comemb[2,,]                       #Co-membership for 3-cycles
cc$cycle.comemb[3,,]                       #Co-membership for 4-cycles
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
\keyword{ math }% __ONLY ONE__ keyword per line
