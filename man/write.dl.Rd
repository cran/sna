\name{write.dl}
\alias{write.dl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Write Output Graphs in DL Format }
\description{
  Writes a graph stack to an output file in DL format.
}
\usage{
write.dl(x, file, vertex.lab = NULL, matrix.lab = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a graph or graph stack, of common order. }
  \item{file}{ a string containing the filename to which the data should be written. }
  \item{vertex.lab}{ an optional vector of vertex labels. }
  \item{matrix.lab}{ an optional vector of matrix labels. }
}
\details{
DL format is used by a number of software packages (including UCINET and Pajek) to store network data.  \code{write.dl} saves one or more (possibly valued) graphs in DL edgelist format, along with vertex and graph labels (if desired).  These files can, in turn, be used to import data into other software packages.
}
\value{
  None.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{write.nos}}  }
\examples{
\dontrun{
#Generate a random graph stack
g<-rgraph(5,10)

#This would save the graphs in DL format
write.dl(g,file="testfile.dl")
}
}
\keyword{ file }% at least one, from doc/KEYWORDS
\keyword{ graphs }% __ONLY ONE__ keyword per line
