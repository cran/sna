\name{gt}
\Rdversion{1.1}
\alias{gt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transpose an Input Graph
}
\description{
\code{gt} returns the graph transpose of its input.  For an adjacency matrix, this is the same as using \code{\link{t}}; however, this function is also applicable to sna edgelists (which cannot be transposed in the usual fashion).  Code written using \code{gt} instead of \code{t} is thus guaranteed to be safe for either form of input.
}
\usage{
gt(x, return.as.edgelist = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
one or more graphs.
}
  \item{return.as.edgelist}{
logical; should the result be returned in sna edgelist form?
}
}
\details{
The transpose of a (di)graph, \eqn{G=(V,E)}, is the graph \eqn{G=(V,E')} where \eqn{E'=\{(j,i): (i,j) \in E\}}{E'={(j,i): (i,j) in E}}.  This is simply the graph formed by reversing the sense of the edges.
}
\value{
The transposed graph(s).
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Carter T. Butts \email{buttsc@uci.edu}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{symmetrize}}, \code{\link{t}}
}
\examples{
#Create a graph....
g<-rgraph(5)
g

#Transpose it
gt(g)
gt(g)==t(g)                  #For adjacency matrices, same as t(g)

#Now, see both versions in edgelist form
as.edgelist.sna(g)
gt(g,return.as.edgelist=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
\keyword{ manip }% __ONLY ONE__ keyword per line
