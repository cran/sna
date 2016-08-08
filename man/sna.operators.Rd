\name{sna.operators}
\alias{sna.operators}
\alias{\%c\%.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Graphical Operators }
\description{
 These operators allow for algebraic manupulation of graph adjacency matrices.
}
\usage{
\method{\%c\%}{matrix}(e1,e2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{e1}{ an (unvalued) adjacency matrix. }
  \item{e2}{ another (unvalued) adjacency matrix. }
}
\details{
  Currently, only one operator is supported.  \code{x \%c\% y} returns the adjacency matrix of the composition of graphs with adjacency matrices \code{x} and \code{y} (respectively).  (Note that this may contain loops.)
}
\value{
   The resulting adjacency matrix.
}
\references{ Wasserman, S. and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: University of Cambridge Press. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
#Create an in-star
g<-matrix(0,6,6)
g[2:6,1]<-1
gplot(g)

#Compose g with its transpose
gcgt<-g\%c\%t(g)
gplot(gcgt,diag=TRUE)
gcgt
}
\keyword{ math }% at least one, from doc/KEYWORDS
\keyword{ graphs }% __ONLY ONE__ keyword per line
