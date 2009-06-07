\name{rgnmix}
\Rdversion{1.1}
\alias{rgnmix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw Mixing-Conditioned Random Graphs
}
\description{
 \code{rgnmix} generates random draws from a mixing-conditioned uniform random graph distribution.
}
\usage{
rgnmix(n, tv, mix, mode = "digraph", diag = FALSE, 
    method = c("probability", "exact"), return.as.edgelist = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
the number of graphs to generate.
}
  \item{tv}{
a vector of types or classes (one entry per vertex), corresponding to the rows and columns of \code{mix}.  (Note that the total number of vertices generated will be \code{length(tv)}.)
}
  \item{mix}{
a class-by-class mixing matrix, containing either mixing rates (for \code{method=="probability"}) or edge counts (for \code{method=="exact"}).
}
  \item{mode}{
\code{"digraph"} for directed graphs, or \code{"graph"} for undirected graphs.
}
  \item{diag}{
logical; should loops be allowed?
}
  \item{method}{
the generation method to use.  \code{"probability"} results in a Bernoulli edge distribution (conditional on the underlying rates), while \code{"exact"} results in a uniform draw conditional on the exact per-block edge distribution.
}
  \item{return.as.edgelist}{
logical; should the resulting graphs be returned in sna edgelist form?
}
}
\details{
The generated graphs (in either adjacency or edgelist form).
}
\value{
\code{rgnmix} draws from a simple generalization of the Erdos-Renyi N,M family (and the related N,p family), generating graphs with fixed expected or realized mixing rates.  Mixing is determined by the \code{mix} argument, which must contain a class-by-class matrix of mixing rates (either edge probabilities or number of realized edges).  The class for each vertex is specified in \code{tv}, whose entries must correspond to the rows and columns of \code{mix}.  The resulting functionality is much like \code{\link{blockmodel.expand}}, although more general (and in some cases more efficient).
}
\references{
Wasserman, S. and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press.
}
\author{
Carter T. Butts \email{buttsc@uci.edu}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rguman}}, \code{\link{rgnm}}, \code{\link{blockmodel.expand}}
}
\examples{
#Draw a random mixing matrix
mix<-matrix(runif(9),3,3)

#Generate a graph with 4 members per class
g<-rgnmix(1,rep(1:3,each=4),mix)
plot.sociomatrix(g)                                 #Visualize the result

#Repeat the exercise, using the exact method
mix2<-round(mix*8)                                  #Draw an exact matrix
g<-rgnmix(1,rep(1:3,each=4),mix2,method="exact")
plot.sociomatrix(g)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
\keyword{ graphs }% __ONLY ONE__ keyword per line
