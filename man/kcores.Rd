\name{kcores}
\Rdversion{1.1}
\alias{kcores}
\alias{kcores_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the k-Core Structure of a Graph
}
\description{
\code{kcores} calculates the k-core structure of the input network, using the centrality measure indicated in \code{cmode}.
}
\usage{
kcores(dat, mode = "digraph", diag = FALSE, cmode = "freeman",
    ignore.eval = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
one or more (possibly valued) graphs.
}
  \item{mode}{
\code{"digraph"} for directed data, otherwise \code{"graph"}.
}
  \item{diag}{
logical; should self-ties be included in the degree calculations?
}
  \item{cmode}{
the \code{\link{degree}} centrality mode to use when constructing the cores.
}
  \item{ignore.eval}{
logical; should edge values be ignored when computing degree?
}
}
\details{
Let \eqn{G=(V,E)} be a graph, and let \eqn{f(v,S,G)} for \eqn{v \in V, S\subseteq V}{v in V, S subseteq V} be a real-valued \emph{vertex property function} (in the language of Batagelj and Zaversnik).  Then some set \eqn{H \subseteq V}{V subseteq H} is a \emph{generalized k-core} for \eqn{f} if \eqn{H} is a maximal set such that \eqn{f(v,H,G)\ge k}{f(v,H,G)>=k} for all \eqn{v \in H}{v in H}.  Typically, \eqn{f} is chosen to be a degree measure with respect to \eqn{S} (e.g., the number of ties to vertices in \eqn{S}).  In this case, the resulting k-cores have the intuitive property of being maximal sets such that every set member is tied (in the appropriate manner) to at least k others within the set.

Degree-based k-cores are a simple tool for identifying well-connected structures within large graphs.  Let the \emph{core number} of vertex \eqn{v} be the value of the highest-value core containing \eqn{v}.  Then, intuitively, vertices with high core numbers belong to relatively well-connected sets (in the sense of sets with high minimum internal degree).  It is important to note that, while a given k-core need not be connected, it is composed of subsets which are themselves well-connected; thus, the k-cores can be thought of as unions of relatively cohesive subgroups.  As k-cores are nested, it is also natural to think of each k-core as representing a \dQuote{slice} through a hypothetical \dQuote{cohesion surface} on \eqn{G}.  (Indeed, k-cores are often visualized in exactly this manner.)

The \code{kcores} function produces degree-based k-cores, for various degree measures (with or without edge values).  The return value is the vector of core numbers for \eqn{V}, based on the selected degree measure.  Missing (i.e., \code{NA}) edge are removed for purposes of the degree calculation.
}
\value{
A vector containing the maximum core membership for each vertex.
}
\references{
 Batagelj, V. and Zaversnik, M.  (2002).  \dQuote{An \eqn{O(m)} Algorithm for Cores Decomposition of Networks.}  arXiv:cs/0310049v1
 
 Batagelj, V. and Zaversnik, M.  (2002).  \dQuote{Generalized Cores.} 	arXiv:cs/0202039v1

 Wasserman, S. and Faust,K.  (1994).  \emph{Social Network Analysis: Methods and Applications}.  Cambridge: Cambridge University Press.
}
\author{
Carter T. Butts \email{buttsc@uci.edu}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{degree}}
}
\examples{
#Generate a graph with core-periphery structure
cv<-runif(30)
g<-rgraph(30,tp=cv\%o\%cv)

#Compute the k-cores based on total degree
kc<-kcores(g)
kc

#Plot the result
gplot(g,vertex.col=kc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
\keyword{ graphs }% __ONLY ONE__ keyword per line
\keyword{ cluster }
