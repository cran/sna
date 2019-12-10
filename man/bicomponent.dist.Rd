\name{bicomponent.dist}
\Rdversion{1.1}
\alias{bicomponent.dist}
\alias{bicomponents_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the Bicomponents of a Graph
}
\description{
\code{bicomponent.dist} returns the bicomponents of an input graph, along with size distribution and membership information.
}
\usage{
bicomponent.dist(dat, symmetrize = c("strong", "weak"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
a graph or graph stack.
}
  \item{symmetrize}{
symmetrization rule to apply when pre-processing the input (see \code{\link{symmetrize}}).
}
}
\details{
The bicomponents of undirected graph \code{G} are its maximal 2-connected vertex sets.  \code{bicomponent.dist} calculates the bicomponents of \eqn{G}, after first coercing to undirected form using the symmetrization rule in \code{symmetrize}.  In addition to bicomponent memberships, various summary statistics regarding the bicomponent distribution are returned; see below.
}
\value{
A list containing
\item{members }{A list, with one entry per bicomponent, containing component members.}
\item{memberships }{A vector of component memberships, by vertex.  (Note: memberships may not be unique.)  Vertices not belonging to any bicomponent have membership values of \code{NA}.}
\item{csize }{A vector of component sizes, by bicomponent.}
\item{cdist }{A vector of length \eqn{|V(G)|} with the (unnormalized) empirical distribution function of bicomponent sizes.}
}
\references{
Brandes, U. and Erlebach, T.  (2005).  \emph{Network Analysis: Methodological Foundations.}  Berlin: Springer.
}
\author{
Carter T. Butts \email{buttsc@uci.edu}
}
\note{
Remember that bicomponents can intersect; when this occurs, the relevant vertices' entries in the membership vector are assigned to one of the overlapping bicomponents on an arbitrary basis.   The \code{members} element of the return list is the safe way to recover membership information.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{component.dist}}, \code{\link{cutpoints}}, \code{\link{symmetrize}}
}
\examples{
#Draw a moderately sparse graph
g<-rgraph(25,tp=2/24,mode="graph")

#Compute the bicomponents
bicomponent.dist(g)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
\keyword{ graphs }% __ONLY ONE__ keyword per line
