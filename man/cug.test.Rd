\name{cug.test}
\Rdversion{1.1}
\alias{cug.test}
\alias{plot.cug.test}
\alias{print.cug.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Univariate Conditional Uniform Graph Tests
}
\description{
\code{cug.test} takes an input network and conducts a conditional uniform graph (CUG) test of the statistic in \code{FUN}, using the conditioning statistics in \code{cmode}.  The resulting test object has custom print and plot methods.
}
\usage{
cug.test(dat, FUN, mode = c("digraph", "graph"), cmode = c("size", 
    "edges", "dyad.census"), diag = FALSE, reps = 1000, 
    ignore.eval = TRUE, FUN.args = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
one or more input graphs.
}
  \item{FUN}{
the function generating the test statistic; note that this must take a graph as its first argument, and return a single numerical value.
}
  \item{mode}{
\code{graph} if \code{dat} is an undirected graph, else \code{digraph}.
}
  \item{cmode}{
string indicating the type of conditioning to be applied.
}
  \item{diag}{
logical; should self-ties be treated as valid data?
}
  \item{reps}{
number of Monte Carlo replications to use.
}
  \item{ignore.eval}{
logical; should edge values be ignored?  (Note: \code{TRUE} is usually more efficient.)
}
  \item{FUN.args}{
a list containing any additional arguments to \code{FUN}.
}
}
\details{
\code{cug.test} is an improved version of \code{cugtest}, for use only with univariate CUG hypotheses.  Depending on \code{cmode}, conditioning on the realized size, edge count (or exact edge value distribution), or dyad census (or dyad value distribution) can be selected.  Edges are treated as unvalued unless \code{ignore.eval=FALSE}; since the latter setting is less efficient for sparse graphs, it should be used only when necessary.

A brief summary of the theory and goals of conditional uniform graph testing can be found in the reference below.  See also \code{\link{cugtest}} for a somewhat informal description.
}
\value{
An object of class \code{cug.test}.
}
\references{
Butts, Carter T.  (2008).  \dQuote{Social Networks: A Methodological Introduction.}  \emph{Asian Journal of Social Psychology,} 11(1), 13--41.
}
\author{
Carter T. Butts \email{buttsc@uci.edu}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{cugtest}}
}
\examples{
#Draw a highly reciprocal network
g<-rguman(1,15,mut=0.25,asym=0.05,null=0.7)

#Test transitivity against size, density, and the dyad census
cug.test(g,gtrans,cmode="size")
cug.test(g,gtrans,cmode="edges")
cug.test(g,gtrans,cmode="dyad.census")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{ math }
\keyword{ graphs }
