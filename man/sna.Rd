\name{sna}
\alias{sna}
%\docType{data}
\title{Tools for Social Network Analysis}
\description{
 \code{sna} is a package containing a range of tools for social network analysis.  Supported functionality includes node and graph-level indices, structural distance and covariance methods, structural equivalence detection, p* modeling, random graph generation, and 2D/3D network visualization (among other things).  
}
%\usage{data(sna)}
%\format{
%  The format is:
% logi NA
%}
\details{
Network data for \code{sna} routines can (except as noted otherwise) appear in any of the following forms:

\itemize{
\item adjacency matrices (dimension N x N);
\item arrays of adjacency matrices, aka \dQuote{graph stacks} (dimension m x N x N);
\item sna edge lists (see below);
\item sparse matrix objects (from the SparseM package);
\item \code{network} objects (from the \link[network]{network} package); or
\item lists of adjacency matrices/arrays, sparse matrices, and/or \code{network} objects.
}

Within the package documentation, the term \dQuote{graph} is used generically to refer to any or all of the above (with multiple graphs being referred to as a \dQuote{graph stack}).  Note that usage of sparse matrix objects requires that the SparseM package be installed.  (No additional packages are required for use of adjacency matrices/arrays or lists thereof, though the network package, on which sna depends as of 2.4, is used for network objects.)  In general, \code{sna} routines attempt to make intelligent decisions regarding the processing of multiple graphs, but common sense is always advised; certain functions, in particular, have more specific data requirements.  Calling \code{sna} functions with inappropriate input data can produce \dQuote{interesting} results.

One special data type supported by the sna package (as of version 2.0) is the \emph{sna edgelist}.  This is a simple data format that is well-suited to representing large, sparse graphs.  (As of version 2.0, many - now most - package routines also process data in this form natively, so using it can produce significant savings of time and/or memory.  Prior to 2.0, all package functions coerced input data to adjacency matrix form.)  An sna edgelist is a three-column matrix, containing (respectively) senders, receivers, and values for each edge in the graph.  (Unvalued edges should have a value of 1.)  Note that this form is invariant to the number of edges in the graph: if there are no edges, then the edgelist is a degenerate matrix of dimension 0 by 3.  Edgelists for undirected graphs should be coded as fully mutual digraphs (as would be the case with an adjacency matrix), with two edges per dyad (one (i,j) edge, and one (j,i) edge).  Graph size for an sna edgelist matrix is indicated by a mandatory numeric attribute, named \code{"n"}.  Vertex names may be optionally specified by a vector-valued attribute named \code{"vnames"}.  In the case of two-mode data (i.e., data with an enforced bipartition), it is possible to indicate this status via the optional \code{"bipartite"} attribute.  Vertices in a two-mode edgelist should be grouped in mode order, with \code{"n"} equal to the total number of vertices (across both modes) and \code{"bipartite"} equal to the number of vertices in the first mode.

Direct creation of sna edgelists can be performed by creating a three-column matrix and using the \code{\link{attr}} function to create the required \code{"n"} attribute.  Alternately, the function \code{\link{as.edgelist.sna}} can be used to coerce data in any of the above forms to an sna edgelist.  By turns, the function \code{\link{as.sociomatrix.sna}} can be used to convert any of these data types to adjacency matrix form.

To get started with \code{sna}, try obtaining viewing the list of available functions.  This can be accomplished via the command \code{library(help=sna)}.
}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{
If you use this package and/or software manual in your work, a citation would be appreciated.  The \code{link{citation}} function has helpful information in this regard.  See also the following paper, which explores the package in some detail:

Butts, Carter T.  (2008).  \dQuote{Social Network Analysis with sna.}  \emph{Journal of Statistical Software}, 24(6).

If utilizing a contributed routine, please also consider recognizing the author(s) of that specific function.  Contributing authors, if any, are listed on the relevant manual pages.  Your support helps to encourage the growth of the \code{sna} package, and is greatly valued!
}
%\source{
%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
%\examples{
%data(sna)
%## maybe str(sna) ; plot(sna) ...
%}
\keyword{graphs}
\keyword{misc}

