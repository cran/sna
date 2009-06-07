\name{sna-internal}
\alias{sna-internal}
\alias{bbnam.jntlik}
\alias{bbnam.jntlik.slice} 
\alias{bbnam.probtie}
\alias{as.edgelist.sna}
\alias{as.sociomatrix.sna}
\alias{logSum}
\alias{logMean}
\alias{logSub}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Internal sna Functions }
\description{
  Internal \code{sna} functions.
}
\usage{
as.sociomatrix.sna(x, attrname=NULL, simplify=TRUE, force.bipartite=FALSE)
as.edgelist.sna(x, attrname = NULL, as.digraph = TRUE, 
    suppress.diag = FALSE, force.bipartite = FALSE) 
bbnam.jntlik(dat, log=FALSE, ...)
bbnam.jntlik.slice(s, dat, a, em, ep, log=FALSE)
bbnam.probtie(dat, i, j, npriorij, em, ep)
logMean(x)
logSub(x, y)
logSum(x)
}
%- maybe also `usage' for other objects documented here.
\details{
These are not to be called by the end user.
}
\keyword{internal}
