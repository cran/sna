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
\item \code{network} objects; or
\item lists of adjacency matrices/arrays and/or \code{network} objects.
}

Within the package documentation, the term "graph" is used generically to refer to any or all of the above (with multiple graphs being referred to as a "graph stack").  Note that usage of \code{network} objects requires that the \link[network]{network} package be installed.  (No additional packages are required for use of adjacency matrices/arrays or lists thereof.)  In general, \code{sna} routines attempt to make intelligent decisions regarding the processing of multiple graphs, but common sense is always advised; certain functions, in particular, have more specific data requirements.  Calling \code{sna} functions with inappropriate input data can produce "interesting" results.

To get started with \code{sna}, try obtaining viewing the list of available functions.  This can be accomplished via the command \code{library(help=sna)}.
}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
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

