\name{event2dichot}
\alias{event2dichot}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Convert an Observed Event Matrix to a Dichotomous matrix }
\description{
  Given a valued adjacency matrix (possibly derived from observed interaction ``events''), \code{event2dichot} returns a dichotomous adjacency matrix.
}
\usage{
event2dichot(m, method="quantile", thresh=0.5, leq=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{m}{ An adjacency matrix }
  \item{method}{ One of ``quantile,'' ``rquantile,'' ``cquantile,'' ``mean,'' ``rmean,'' ``cmean,'' ``absolute,'' ``rank,'' ``rrank,'' or ``crank''  }
  \item{thresh}{ Dichotomization thresholds for ranks or quantiles }
  \item{leq}{ Boolean indicating whether values less than or equal to the threshold should be taken as existing edges; the alternative is to use values strictly greater than the threshold }
}
\details{
The methods used for choosing dichotomization thresholds are as follows:

\enumerate{
\item quantile: specified quantile over the distribution of all edge values

\item rquantile: specified quantile by row

\item cquantile: specified quantile by column

\item mean: grand mean

\item rmean: row mean

\item cmean: column mean

\item absolute: the value of \code{thresh} itself

\item rank: specified rank over the distribution of all edge values

\item rrank: specified rank by row

\item crank: specified rank by column
}

Note that when a quantile, rank, or value is said to be ``specified,'' this refers to the value of \code{thresh}.  
}
\value{
The dichotomized data matrix
}
\references{ ~put references to the literature/web site here ~ }
\author{ Carter T. Butts \email{ctb@andrew.cmu.edu} }
%\note{ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

%\seealso{ }

\examples{
#Draw a matrix of normal values
n<-matrix(rnorm(25),nrow=5,ncol=5)

#Dichotomize by the mean value
event2dichot(n,"mean")

#Dichotomize by the 0.95 quantile
event2dichot(n,"quantile",0.95)

}
\keyword{ manip }%-- one or more ...
\keyword{ math }
