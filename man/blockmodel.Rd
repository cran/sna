\name{blockmodel}
\alias{blockmodel}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Generate Blockmodels Based on Partitions of Network Positions }
\description{
  Given an equivalence clustering object and a graph stack, \code{blockmodel} will form a blockmodel of the graph stack based on the classes in question, using the specified block content type.
}
\usage{
blockmodel(dat, ec, k=NULL, h=NULL, block.content="density", 
    plabels=ec$plabels, glabels=ec$glabels, rlabels=NULL, 
    mode="digraph", diag=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ A graph stack }
  \item{ec}{ An object of class \code{equiv.clust} }
  \item{k}{ The number of classes to form (using \code{\link{cutree}}) }
  \item{h}{ The height at which to split classes (using \code{\link{cutree}}) }
  \item{block.content}{ String indicating block content type (see below)}
  \item{plabels}{ A vector of labels to be applied to the individual nodes}
  \item{glabels}{ A vector of labels to be applied to the graphs being modeled }
  \item{rlabels}{ A vector of labels to be applied to the (reduced) roles }
  \item{mode}{ A string indicating whether we are dealing with graphs or digraphs }
  \item{diag}{ A boolean indicating whether loops are permitted }
}
\details{
\code{blockmodel} forms its eponymous models by using \code{\link{cutree}} to cut an equivalence clustering in the fashion specified by \code{k} and \code{h}.  After forming clusters (roles), the input graphs are reordered and blockmodel reduction is applied.  Currently supported reductions are:

\enumerate{
\item \code{density}: block density, computed as the mean value of the block

\item \code{meanrowsum}: mean row sums for the block

\item \code{meancolsum}: mean column sums for the block

\item \code{sum}: total block sum

\item \code{median}: median block value

\item \code{min}: minimum block value

\item \code{max}: maximum block value

\item \code{types}: semi-intelligent coding of blocks by ``type.''  Currently recognized types are (in order of precedence) ``\code{NA}'' (i.e., blocks with no valid data), ``null'' (i.e., all values equal to zero), ``complete'' (i.e., all values equal to 1), ``1 covered'' (i.e., all rows/cols contain a 1), ``1 row-covered'' (i.e., all rows contain a 1), ``1 col-covered'' (i.e., all cols contain a 1), and ``other'' (i.e., none of the above). 
}

Density or median-based reductions are probably the most interpretable for most conventional analyses, though type-based reduction can be useful in examining certain equivalence class hypotheses (e.g., 1 covered and null blocks can be used to infer regular equivalence classes).  Once a given reduction is performed, the model can be analyzed and/or expansion can be used to generate new graphs based on the inferred role structure.
}
\value{
An object of class \code{blockmodel}.
}
\references{White, H.C.; Boorman, S.A.; and Breiger, R.L.  (1976).  ``Social Structure from Multiple Networks I: Blockmodels of Roles and Positions.''  \emph{American Journal of Sociology}, 81, 730-779.}

\author{ Carter T. Butts \email{buttsc@uci.edu} }

%\note{ The \code{mva} library is required for \code{blockmodel} to function. }


\seealso{ \code{\link{equiv.clust}}, \code{\link{blockmodel.expand}} }

\examples{
#Create a random graph with _some_ edge structure
g.p<-sapply(runif(20,0,1),rep,20)  #Create a matrix of edge 
                                   #probabilities
g<-rgraph(20,tprob=g.p)            #Draw from a Bernoulli graph 
                                   #distribution

#Cluster based on structural equivalence
eq<-equiv.clust(g)

#Form a blockmodel with distance relaxation of 10
b<-blockmodel(g,eq,h=10)
plot(b)                            #Plot it
}
\keyword{ classif }%-- one or more ...
\keyword{ math }
