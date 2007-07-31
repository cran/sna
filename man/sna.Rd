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
\item arrays of adjacency matrices, aka "graph stacks" (dimension m x N x N);
\item sparse matrix objects (from the \code{\link[SparseM]{SparseM}} package);
\item \code{network} objects (from the \code{\link[network]{network}} package); or
\item lists of adjacency matrices/arrays, sparse matrices, and/or \code{network} objects.
}

Within the package documentation, the term "graph" is used generically to refer to any or all of the above (with multiple graphs being referred to as a "graph stack").  Note that usage of sparse matrix or \code{network} objects requires that the \link[SparseM]{SparseM} or \link[network]{network} package be installed, respectively.  (No additional packages are required for use of adjacency matrices/arrays or lists thereof.)  In general, \code{sna} routines attempt to make intelligent decisions regarding the processing of multiple graphs, but common sense is always advised; certain functions, in particular, have more specific data requirements.  Calling \code{sna} functions with inappropriate input data can produce "interesting" results.

To get started with \code{sna}, try obtaining viewing the list of available functions.  This can be accomplished via the command \code{library(help=sna)}.
}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{
If you use this package and/or software manual in your work, a citation would be appreciated.  References to the current versions are:

Butts, Carter T.  (2007).  ``sna: Tools for Social Network Analysis.''  R package version 1.5.

Butts, Carter T.  (2007).  ``Software Manual for the R sna Package.''  R package version 1.5.

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

