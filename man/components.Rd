\name{components}
\alias{components}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Find the Number of (Maximal) Components Within a Given Graph }
\description{
  Returns the number of components within \code{dat}, using the connectedness rule given in \code{connected}.
}
\usage{
components(dat, connected="strong", comp.dist.precomp=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more input graphs. }
  \item{connected}{ the the component definition to be used by \code{\link{component.dist}} during component extraction. }
  \item{comp.dist.precomp}{ a component size distribution object from \code{\link{component.dist}} (optional). }
}
\details{
  The \code{connected} parameter corresponds to the \code{rule} parameter of \code{\link{component.dist}}.  By default, \code{components} returns the number of strong components, but other component types can be returned if so desired.  (See \code{\link{component.dist}} for details.)  For symmetric matrices, this is obviously a moot point.
}
\value{
  A vector containing the number of components for each graph in \code{dat}
}
\references{ West, D.B.  (1996).  \emph{Introduction to Graph Theory}.  Upper Saddle River, NJ: Prentice Hall. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }

\seealso{ \code{\link{component.dist}}, \code{\link{symmetrize}} }

\examples{
g<-rgraph(20,tprob=0.05)   #Generate a sparse random graph

#Find weak components
components(g,connected="weak")

#Find strong components
components(g,connected="strong")

}
\keyword{ math }%-- one or more ...
\keyword{ graphs }
