\name{write.nos}
\alias{write.nos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Write Output Graphs in (N)eo-(O)rg(S)tat Format }
\description{
  Writes a graph stack to an output file in NOS format.
}
\usage{
write.nos(x, file, row.col = NULL, col.col = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a graph or graph stack (all graphs must be of common order). }
  \item{file}{ string containing the output file name. }
  \item{row.col}{ vector of row labels (or "row colors"). }
  \item{col.col}{ vector of column labels ("column colors"). }
}
\details{
NOS format consists of three header lines, followed by a whitespace delimited stack of raw adjacency matrices; the format is not particularly elegant, but turns up in certain legacy applications (mostly at CMU).  \code{write.nos} provides a quick and dirty way of writing files NOS, which can later be retrieved using \code{\link{read.nos}}.

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

Note that NOS format only supports graph stacks of common order; graphs of different sizes cannot be stored within the same file.
}
\value{
  None.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{read.nos}}, \code{\link{write.dl}}, \code{\link{write.table}} }
\examples{
\dontrun{
#Generate a random graph stack
g<-rgraph(5,10)

#This would save the graphs in NOS format
write.nos(g,file="testfile.nos")

#We can also read them back, like so:
g2<-read.nos("testfile.nos")
}
}
\keyword{ graphs }% at least one, from doc/KEYWORDS
\keyword{ file }% __ONLY ONE__ keyword per line
