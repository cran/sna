\name{neighborhood}
\alias{neighborhood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute Neighborhood Structures of Specified Order }
\description{
  For a given graph, returns the specified neighborhood structure at the selected order(s).
}
\usage{
neighborhood(dat, order, neighborhood.type = c("in", "out", "total"),
    mode = "digraph", diag = FALSE, thresh = 0, return.all = FALSE,
    partial = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{ one or more graphs. }
  \item{order}{ order of the neighborhood to extract. }
  \item{neighborhood.type}{ neighborhood type to employ. }
  \item{mode}{ \code{"digraph"} if \code{dat} is directed, otherwise \code{"graph"}. }
  \item{diag}{ logical; do the diagonal entries of \code{dat} contain valid data? }
  \item{thresh}{ dichotomization threshold to use for \code{dat}; edges whose values are greater than \code{thresh} are treated as present. }
  \item{return.all}{ logical; return neighborhoods for all orders up to \code{order}? }
  \item{partial}{ logical; return partial (rather than cumulative) neighborhoods? }
}
\details{
The adjacency matrix associated with the \eqn{i}th order neighborhood is defined as the identity matrix for order 0, and otherwise depends on the type of neighborhood involved.  For input graph \eqn{G=(V,E)}, let the \emph{base relation}, \eqn{R}, be given by the underlying graph of \eqn{G} (i.e., \eqn{G \cup G^T}{G U G^T}) if total neighborhoods are sought, the transpose of \eqn{G} if incoming neighborhoods are sought, or \eqn{G} otherwise.  The partial neighborhood structure of order \eqn{i>0} on \eqn{R} is then defined to be the digraph on \eqn{V} whose edge set consists of the ordered pairs \eqn{(j,k)} having geodesic distance \eqn{i} in \eqn{R}.  The corresponding cumulative neighborhood is formed by the ordered pairs having geodesic distance less than or equal to \eqn{i} in \eqn{R}.

Neighborhood structures are commonly used to parameterize various types of network autocorrelation models.  They may also be used in the calculation of certain types of local structural indices; \code{\link{gapply}} provides an alternative function which can be used for this purpose.
}
\value{
  An array or adjacency matrix containing the neighborhood structures (if \code{dat} is a single graph); if \code{dat} contains multiple graphs, then a list of such structures is returned.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{gapply}}, \code{\link{nacf}}}
\examples{
#Draw a random graph
g<-rgraph(10,tp=2/9)

#Show the total partial out-neighborhoods
neigh<-neighborhood(g,9,neighborhood.type="out",return.all=TRUE)
par(mfrow=c(3,3))
for(i in 1:9)
  gplot(neigh[i,,],main=paste("Partial Neighborhood of Order",i))

#Show the total cumulative out-neighborhoods
neigh<-neighborhood(g,9,neighborhood.type="out",return.all=TRUE,
    partial=FALSE)
par(mfrow=c(3,3))
for(i in 1:9)
  gplot(neigh[i,,],main=paste("Cumulative Neighborhood of Order",i))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
\keyword{ manip }% __ONLY ONE__ keyword per line
