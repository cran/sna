\name{symmetrize}
\alias{symmetrize}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Symmetrize an Adjacency Matrix }
\description{
 Symmetrizes the elements of \code{mats} according to the rule in \code{rule}.
}
\usage{
symmetrize(mats, rule="weak", return.as.edgelist=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{mats}{ a graph or graph stack.}
  \item{rule}{ one of ``upper'', ``lower'', ``strong'' or ``weak''. }
  \item{return.as.edgelist}{ logical; should the symmetrized graphs be returned in edgelist form?}
}
\details{
 The rules used by \code{symmetrize} are as follows:
\enumerate{
\item upper: Copy the upper triangle over the lower triangle

\item lower: Copy the lower triangle over the upper triangle

\item strong: i<->j iff i->j and i<-j  (AND rule)

\item weak: i<->j iff i->j or i<-j  (OR rule)
}
}
\value{
The symmetrized graph stack
}
\references{ Wasserman, S., and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications}.  Cambridge: Cambridge University Press. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

 %~Make other sections like WARNING with \section{WARNING }{....} ~

%\seealso{ }

\examples{
#Generate a graph
g<-rgraph(5)

#Weak symmetrization
symmetrize(g)

#Strong symmetrization
symmetrize(g,rule="strong")
}
\keyword{ manip }
\keyword{ array }%-- one or more ...
\keyword{ math }
\keyword{ graphs }
