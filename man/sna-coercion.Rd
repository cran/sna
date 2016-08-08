\name{sna-coercion}
\alias{as.edgelist.sna}
\alias{as.sociomatrix.sna}
\alias{is.edgelist.sna}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ sna Coercion Functions }
\description{
  Functions to coerce network data into one form or another; these are generally internal, but may in some cases be helpful to the end user.
}
\usage{
as.sociomatrix.sna(x, attrname=NULL, simplify=TRUE, force.bipartite=FALSE)
\method{as.edgelist}{sna}(x, attrname = NULL, as.digraph = TRUE, 
    suppress.diag = FALSE, force.bipartite = FALSE, ...)
is.edgelist.sna(x)
}
\arguments{
  \item{x}{network data in any of several acceptable forms (see below).}
  \item{attrname}{if \code{x} is a \code{\link[network]{network}} object, the (optional) edge attribute to be used to obtain edge values.}
  \item{simplify}{logical; should output be simplified by collapsing adjacency matrices of identical dimension into adjacency arrays?}
  \item{force.bipartite}{logical; should the data be interpreted as bipartite (with rows and columns representing different data modes)?}
  \item{as.digraph}{logical; should \code{\link[network]{network}} objects be coded as digraphs, regardless of object properties?  (Recommended)}
  \item{suppress.diag}{logical; should loops be suppressed?}
  \item{...}{ additional arguments to \code{sna.edgelist} (currently ignored). }
}
%- maybe also `usage' for other objects documented here.
\details{
The \code{\link{sna}} coercion functions are normally called internally within user-level \code{\link{sna}} functions to convert network data from various supported forms into a format usable by the function in question.  With few (if any) exceptions, formats acceptable by these functions should be usable with any user-level function in the \code{\link{sna}} library.

\code{as.sociomatrix.sna} takes one or more input graphs, and returns them in adjacency matrix (and/or array) form.  If \code{simplify==TRUE}, consolidation of matrices having the same dimensions into adjacency arrays is attempted; otherwise, elements are returned as lists of matrices/arrays.

\code{as.edgelist.sna} takes one or more input graphs, and returns them in \code{sna} edgelist form -- i.e., a three-column matrix whose rows represent edges, and whose columns contain (respectively) the sender, receiver, and value of each edge.  (Undirected graphs are generally assumed to be coded as fully mutual digraphs; edges may be listed in any order.)  \code{sna} edgelists must also carry an attribute named \code{n} indicating the number of vertices in the graph, and may optionally contain the attributes \code{vnames} (carrying a vector of vertex names, in order) and/or \code{bipartite} (optionally, containing the number of row vertices in a two-mode network).  If the bipartite attribute is present and non-false, vertices whose numbers are less than or equal to the attribute value are taken to belong to the first mode (i.e., row vertices), and those of value greater than the attribute are taken to belong to the second mode (i.e., column vertices).  Note that the \code{bipartite} attribute is not strictly necessary to represent two-mode data, and may not be utilized by all \code{\link{sna}} functions.

\code{is.edgelist.sna} returns \code{TRUE} if its argument is a \code{sna} edgelist, or \code{FALSE} otherwise; if called with a list, this check is performed (recursively) on the list elements.

Data for \code{sna} coercion routines may currently consist of any combination of standard or sparse (via \code{SparseM}) adjacency matrices or arrays, \code{\link[network]{network}} objects, or \code{sna} edgelists.  If multiple items are given, they must be contained within a \code{\link{list}}.  Where adjacency arrays are specified, they must be in three-dimensional form, with dimensions given in graph/sender/receiver order.  Matrices or arrays having different numbers of rows and columns are taken to be two-mode adjacency structures, and are treated accordingly; setting \code{force.bipartite} will cause square matrices to be treated in similar fashion.  In the case of \code{\link[network]{network}} or \code{sna} edgelist matrices, bipartition information is normally read from the object's internal properties.
}
\value{
An adjacency or edgelist structure, or a list thereof.
}
%\references{
%}
\author{ Carter T. Butts \email{buttsc@uci.edu}
}
\note{
For large, sparse graphs, edgelists can be dramatically more efficient than adjacency matrices.  Where such savings can be realized, \code{\link{sna}} package functions usually employ \code{sna} edgelists as their \dQuote{native} format (coercing input data with \code{as.edgelist.sna} as needed).  For this reason, users of large graphs can often obtain considerable savings by storing data in edgelist form, and passing edgelists (rather than adjacency matrices) to \code{\link{sna}} functions.

The maximum size of adjacency matrices and edgelists depends upon \code{R}'s vector allocation limits.  On a 64-bit platform, these limits are currently around 4.6e4 vertices (adjacency case) or 7.1e8 edges (edgelist case).  The number of vertices in the edgelist case is effectively unlimited (and can technically be infinite), although not all functions will handle such objects gracefully.  (Use of vertex names will limit the number of edgelist vertices to around 2e9.)
}
%Just FYI, precision problems for vertices don't arise until around 1e210 vertices, which is far larger than the estimated volume of the universe in cubic Planck lengths (about 1e186, according to Wikipedia).

\seealso{
\code{\link{sna}}, \code{\link[network]{network}}
}
\examples{
#Produce some random data, and transform it
g<-rgraph(5)
g
all(g==as.sociomatrix.sna(g))     #TRUE
as.edgelist.sna(g)                #View in edgelist form
as.edgelist.sna(list(g,g))        #Double the fun
g2<-as.sociomatrix.sna(list(g,g)) #Will simplify to an array
dim(g2)
g3<-as.sociomatrix.sna(list(g,g),simplify=FALSE)  #Do not simplify
g3                                                #Now a list

#We can also build edgelists from scratch...
n<-6
edges<-rbind(
c(1,2,1),
c(2,1,2),
c(1,3,1),
c(1,5,2),
c(4,5,1),
c(5,4,1)
)
attr(edges,"n")<-n
attr(edges,"vnames")<-letters[1:n]
gplot(edges,displaylabels=TRUE)               #Plot the graph
as.sociomatrix.sna(edges)                     #Show in matrix form

#Two-mode data works similarly
n<-6
edges<-rbind(
c(1,4,1),
c(1,5,2),
c(4,1,1),
c(5,1,2),
c(2,5,1),
c(5,2,1),
c(3,5,1),
c(3,6,2),
c(6,3,2)
)
attr(edges,"n")<-n
attr(edges,"vnames")<-c(letters[1:3],LETTERS[4:6])
attr(edges,"bipartite")<-3
edges
gplot(edges,displaylabels=TRUE,gmode="twomode")  #Plot
as.sociomatrix.sna(edges)                        #Convert to matrix

}
\keyword{manip}
\keyword{array}
\keyword{math}
\keyword{graphs}

