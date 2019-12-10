\name{component.size.byvertex}
\alias{component.size.byvertex}
\alias{compsizes_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get Component Sizes, by Vertex
}
\description{
This function computes the component structure of the input network, and returns a vector whose \eqn{i}th entry is the size of the component to which \eqn{i} belongs.  This is useful e.g. for studies of diffusion or similar applications.
}
\usage{
component.size.byvertex(dat, connected = c("strong", "weak", 
    "unilateral", "recursive"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
one or more input graphs (for best performance, sna edgelists or network objects are suggested).
}
  \item{connected}{
a string selecting the connectedness definition to use; by default, \code{"strong"} components are used.
}
}
\details{
Component sizes are here computed using \code{\link{component.dist}}; see this function for additional information.

In an undirected graph, the size of \eqn{v}'s component represents the maximum number of nodes that can be reached by a diffusion process along the edges of the graph originating with node \eqn{v}; the expectation of component sizes by vertex (rather than the mean component size) is thus one measure of the maximum average diffusion potential of a graph.  Because this quantity is monotone with respect to edge addition, it can be bounded using Bernoulli graphs (see Butts (2011)).  In the directed case, multiple types of components are possible; see \code{\link{component.dist}} for details.
}
\value{
A vector of length equal to the number of vertices in \code{dat}, whose \eqn{i}th element is the number of vertices in the component to which the \eqn{i}th vertex belongs.
}
\references{
West, D.B.  (1996).  \emph{Introduction to Graph Theory.} Upper Saddle River, N.J.: Prentice Hall.

Butts, C.T.  (2011).  \dQuote{Bernoulli Bounds for General Random Graphs.}  \emph{Sociological Methodology}, 41, 299-345.
}
\author{
Carter T. Butts \email{buttsc@uci.edu}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{component.dist}}
}
\examples{
#Generate a random undirected graph
g<-rgraph(100,tprob=1.5/99,mode="graph",return.as.edgelist=TRUE)

#Get the component sizes for each vertex
cs<-component.size.byvertex(g)
cs
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
\keyword{ graphs }% __ONLY ONE__ keyword per line
