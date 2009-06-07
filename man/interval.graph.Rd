\name{interval.graph}
\alias{interval.graph}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Convert Spell Data to Interval Graphs }
\description{
  Constructs one or more interval graphs (and exchangeability vectors) from a set of spells.
}
\usage{
interval.graph(slist, type="simple", diag=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{slist}{ A spell list.  This must consist of an nxmx3 array, with n being the number of actors, m being the maximum number of spells (one per row) and with the three columns of the last dimension containing a (categorical) spell type code, the time of spell onset (any units), and the time of spell termination (same units), respectively. }
  \item{type}{ One of ``simple'', ``overlap'', ``fracxy'', ``fracyx'', or ``jntfrac''. }
  \item{diag}{ Include the dyadic entries? }
}
\details{
  Given some ordering dimension T (usually time), a ``spell'' is defined as the interval between a specified onset and a specified termination (with onset preceding the termination).  An interval graph, then, on spell set V, is \eqn{G=\{V,E\}}, where \eqn{\{i,j\} \in E}{\{i,j\} in E} iff there exists some point \eqn{t \in T}{e \in T} such that \eqn{t \in i}{t \in i} and \eqn{t \in j}{t \in j}.  In more prosaic terms, an interval graph on a given spell set has each spell as a vertex, with vertices adjacent iff they overlap.  Such structures are useful for quantifying life history data (where spells might represent marriages, periods of child custody/co-residence, periods of employment, etc.), organizational history data (where spells might reflect periods of strategic alliances, participation in a particular product market, etc.), task scheduling (with spells representing the dedication of a particular resource to a given task), etc.  By giving complex historical data a graphic representation, it is possible to easily perform a range of analyses which would otherwise be difficult and/or impossible (see Butts and Pixley (2004) for examples).

In addition to the simple interval graph (described above), \code{interval.graph} can also generate valued interval graphs using a number of different edge definitions.  This is controlled by the \code{type} argument, with edge values as follows:
\enumerate{
\item simple: dichotomous coding based on simple overlap (i.e., (x,y)=1 iff x overlaps y)

\item overlap: edge value equals the total magnitude of the overlap between spells

\item fracxy: the (x,y) edge value equals the fraction of the duration of y which is covered by x

\item fracyx: the (x,y) edge value equals the fraction of the duration of x which is covered by y

\item jntfrac: edge value equals the total magnitude of the overlap between spells divided by the mean of the spells' lengths
}
Note that ``simple,'' ``overlap,'' and ``jntfrac'' are symmetric relations, while ``fracxy'' and ``fracyx'' are directed.  As always, the specific edge type used should reflect the application to which the interval graph is being put.
}
\value{
  A data frame containing:

  \item{graph }{A graph stack containing the interval graphs}
  \item{exchange.list }{Matrix containing the vector of spell types associated with each interval graph}
}
\references{ 
Butts, C.T. and Pixley, J.E.  (2004).  ``A Structural Approach to the Representation of Life History Data.''  \emph{Journal of Mathematical Sociology}, 28(2), 81-124.

West, D.B.  (1996).  \emph{Introduction to Graph Theory}.  Upper Saddle River, NJ: Prentice Hall. }
\author{ Carter T. Butts \email{buttsc@uci.edu}}
%\note{ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

%\seealso{ }

%\examples{
%}
\keyword{ manip }%-- one or more ...
\keyword{ math }
\keyword{ graphs }
