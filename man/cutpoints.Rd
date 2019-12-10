\name{cutpoints}
\Rdversion{1.1}
\alias{cutpoints}
\alias{cutpointsDir_R}
\alias{cutpointsUndir_R}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify the Cutpoints of a Graph or Digraph
}
\description{
\code{cutpoints} identifies the cutpoints of an input graph.  Depending on \code{mode}, either a directed or undirected notion of \dQuote{cutpoint} can be used.
}
\usage{
cutpoints(dat, mode = "digraph", connected = c("strong","weak","recursive"),
    return.indicator = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
one or more input graphs.
}
  \item{mode}{
\code{"digraph"} for directed graphs, or \code{"graph"} for undirected graphs.
}
  \item{connected}{
string indicating the type of connectedness rule to apply (only relevant where \code{mode=="digraph"}).
}
  \item{return.indicator}{
logical; should the results be returned as a logical (\code{TRUE/FALSE}) vector of indicators, rather than as a vector of vertex IDs?
}
}
\details{
A \emph{cutpoint} (also known as an \emph{articulation point} or \emph{cut-vertex}) of an undirected graph, \eqn{G} is a vertex whose removal increases the number of components of \eqn{G}.  Several generalizations to the directed case exist.  Here, we define a \emph{strong cutpoint} of directed graph \eqn{G} to be a vertex whose removal increases the number of strongly connected components of \eqn{G} (see \code{\link{component.dist}}).  Likewise, \emph{weak} and \emph{recursive} cutpoints of \emph{G} are those vertices whose removal increases the number of weak or recursive cutpoints (respectively).  By default, strong cutpoints are used; alternatives may be selected via the \code{connected} argument.

Cutpoints are of particular interest when seeking to identify critical positions in flow networks, since their removal by definition alters the connectivity properties of the graph.  In this context, cutpoint status can be thought of as a primitive form of centrality (with some similarities to \code{\link{betweenness}}).

Cutpoint computation is significantly faster for the undirected case (and for the weak/recursive cases) than for the strong directed case.  While calling \code{cutpoints} with \code{mode="digraph"} on an undirected graph will give the same answer as \code{mode="graph"}, it is thus to one's advantage to use the latter form.  Do not, however, employ \code{mode="graph"} with directed data, unless you enjoy unpredictable behavior.
}
\value{
A vector of cutpoints (if \code{return.indicator==FALSE}), or else a logical vector indicating cutpoint status for each vertex.
}
\references{
Berge, Claude.  (1966).  \emph{The Theory of Graphs.}  New York: John Wiley and Sons.
}
\author{
Carter T. Butts \email{buttsc@uci.edu}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{component.dist}}, \code{\link{bicomponent.dist}}, \code{\link{betweenness}}
}
\examples{
#Generate some sparse random graph
gd<-rgraph(25,tp=1.5/24)                               #Directed
gu<-rgraph(25,tp=1.5/24,mode="graph")                  #Undirected

#Calculate the cutpoints (as an indicator vector)
cpu<-cutpoints(gu,mode="graph",return.indicator=TRUE)
cpd<-cutpoints(gd,return.indicator=TRUE)

#Plot the result
gplot(gu,gmode="graph",vertex.col=2+cpu)
gplot(gd,vertex.col=2+cpd)

#Repeat with alternate connectivity modes
cpdw<-cutpoints(gd,connected="weak",return.indicator=TRUE)
cpdr<-cutpoints(gd,connected="recursive",return.indicator=TRUE)

#Visualize the difference
gplot(gd,vertex.col=2+cpdw)
gplot(gd,vertex.col=2+cpdr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
\keyword{ graphs }% __ONLY ONE__ keyword per line
