\name{read.dot}
\alias{read.dot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read Graphviz DOT Files }
\description{
  Reads network information in Graphviz's DOT format, returning an adjacency matrix.
}
\usage{
read.dot(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{ The name of the file whence to import the data, or else a connection object (suitable for processing by \code{\link{readLines}}. }
}
\details{
  The Graphviz project's DOT language is a simple but flexible tool for describing graphs.  See the included reference for details.
}
\value{
  The imported graph, in adjacency matrix form.
}
\references{ Graphviz Project.  "The DOT Language."  http://www.graphviz.org/doc/info/lang.html }
\author{ Matthijs den Besten \email{matthijs.denbesten@gmail.com} }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{read.nos}}, \code{\link{write.nos}}, \code{\link{write.dl}} }
%\examples{
%}
\keyword{ graphs }% at least one, from doc/KEYWORDS
\keyword{ file }% __ONLY ONE__ keyword per line
