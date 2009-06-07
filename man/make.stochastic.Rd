\name{make.stochastic}
\alias{make.stochastic}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Make a Graph Stack Row, Column, or Row-column Stochastic }
\description{
  Returns a graph stack in which each adjacency matrix in \code{dat} has been normalized to row stochastic, column stochastic, or row-column stochastic form, as specified by \code{mode}.
}
\usage{
make.stochastic(dat, mode="rowcol", tol=0.005, 
    maxiter=prod(dim(dat)) * 100, anneal.decay=0.01, errpow=1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ a collection of input graphs. }
  \item{mode}{ one of ``row,'' ``col,'' or ``rowcol''. }
  \item{tol}{ tolerance parameter for the row-column normalization algorithm. }
  \item{maxiter}{ maximum iterations for the rwo-column normalization algorithm. }
  \item{anneal.decay}{ probability decay factor for the row-column annealer. }
  \item{errpow}{ power to which absolute row-column normalization errors should be raised for the annealer (i.e., the penalty function). }
}
\details{
  Row and column stochastic matrices are those whose rows and columns sum to 1 (respectively).  These are quite straightforwardly produced here by dividing each row (or column) by its sum.  Row-column stochastic matrices, by contrast, are those in which each row \emph{and} each column sums to 1.  Here, we try to produce row-column stochastic matrices whose values are as close in proportion to the original data as possible by means of an annealing algorithm.  This is probably not optimal in the long term, but the results seem to be consistent where row-column stochasticization of the original data is possible (which it is not in all cases).
}
\value{
The stochasticized adjacency matrices
}
%\references{  }
\author{ Carter T. Butts \email{buttsc@uci.edu}}
%\note{  }

\section{Warning }{Rows or columns which sum to 0 in the original data will generate undefined results.  This can happen if, for instance, your input graphs contain in- or out-isolates.} 

%\seealso{  }

\examples{
#Generate a test matrix
g<-rgraph(15)

#Make it row stochastic
make.stochastic(g,mode="row")

#Make it column stochastic
make.stochastic(g,mode="col")

#(Try to) make it row-column stochastic
make.stochastic(g,mode="rowcol")
}
\keyword{ manip }
\keyword{ array }
\keyword{ algebra }%-- one or more ...
