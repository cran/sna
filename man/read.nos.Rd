\name{read.nos}
\alias{read.nos}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Read (N)eo-(O)rg(S)tat Input Files}
\description{
  Reads an input file in NOS format, returning the result as a graph set.
}
\usage{
read.nos(file, return.as.edgelist = FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{file}{ the file to be imported }
  \item{return.as.edgelist}{ logical; should the resulting graphs be returned in sna edgelist format?}
}
\details{
  NOS format consists of three header lines, followed by a whitespace delimited stack of raw adjacency matrices; the format is not particularly elegant, but turns up in certain legacy applications (mostly at CMU).  \code{read.nos} provides a quick and dirty way of reading in these files, without the headache of messing with \code{\link{read.table}} settings.

The content of the NOS format is as follows:

<m>

<n> <o>

<kr1> <kr2> ... <krn> <kc1> <kc2> ... <kcn>

<a111> <a112> ... <a11o>

<a121> <a122> ... <a12o>

...

<a1n1> <a1n2> ... <a1no>

<a211> <a212> ... <a21o>

...

<a2n1> <a2n2> ... <a2no>

...

<amn1> <amn2> ... <amno>


where <abcd> is understood to be the value of the c->d edge in the bth graph of the file.  (As one might expect, m, n, and o are the numbers of graphs (matrices), rows, and columns for the data, respectively.)  The "k" line contains a list of row and column "colors", categorical variables associated with each row and column, respectively.  Although originally intended to communicate exchangability information, these can be used for other purposes (though there are easier ways to deal with attribute data these days).
}
\value{
  The imported graph set (in adjacency array or edgelist form).
}
%\references{  }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ \code{read.nos} currently ignores the coloring information. }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{write.nos}}, \code{\link{scan}}, \code{\link{read.table}} }

%\examples{
%
%}
\keyword{ file }%-- one or more ...
\keyword{graphs}

