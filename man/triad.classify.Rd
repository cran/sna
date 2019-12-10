\name{triad.classify}
\alias{triad.classify}
\alias{triad_classify_R}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the Davis and Leinhardt Classification of a Given Triad }
\description{
  \code{triad.classify} returns the Davis and Leinhardt classification of the triad indicated by \code{tri} in the \code{g}th graph of stack \code{dat}.
}
\usage{
triad.classify(dat, g=1, tri=c(1, 2, 3), mode=c("digraph", 
    "graph"))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ a graph or graph stack. }
  \item{g}{ the index of the graph to be analyzed. }
  \item{tri}{ a triple containing the indices of the triad to be classified.}
  \item{mode}{string indicating the directedness of edges; \code{"digraph"} implies a directed structure, whereas \code{"graph"} implies an undirected structure.}
}
\details{
Every unoriented directed triad may occupy one of 16 distinct states.  These states were used by Davis and Leinhardt as a basis for classifying triads within a larger structure; the distribution of triads within a graph (see \code{\link{triad.census}}), for instance, is linked to a range of substantive hypotheses (e.g., concerning structural balance).  The Davis and Leinhardt classification scheme describes each triad by a string of four elements: the number of mutual (complete) dyads within the triad; the number of asymmetric dyads within the triad; the number of null (empty) dyads within the triad; and a configuration code for the triads which are not uniquely distinguished by the first three distinctions.  The complete list of classes is as follows.

\describe{
\item{\code{003}}{ \eqn{a \not\leftrightarrow b \not\leftrightarrow c, a \not\leftrightarrow c}{a<-!->b<-!->c, a<-!->c}}
\item{\code{012}}{ \eqn{a \rightarrow b \not\leftrightarrow c, a \not\leftrightarrow c}{a-->b<-!->c, a<-!->c}}
\item{\code{102}}{ \eqn{a \leftrightarrow b \not\leftrightarrow c, a \not\leftrightarrow c}{a<-->b<-!->c, a<-!->c}}
\item{\code{021D}}{ \eqn{a \leftarrow b \rightarrow c, a \not\leftrightarrow c}{a<--b-->c, a<-!->c}}
\item{\code{021U}}{ \eqn{a \rightarrow b \leftarrow c, a \not\leftrightarrow c}{a-->b<--c, a<-!->c}}
\item{\code{021C}}{ \eqn{a \rightarrow b \rightarrow c, a \not\leftrightarrow c}{a-->b-->c, a<-!->c}}
\item{\code{111D}}{ \eqn{a \not\leftrightarrow b \rightarrow c, a \leftrightarrow c}{a<-!->b-->c, a<-->c}}
\item{\code{111U}}{ \eqn{a \not\leftrightarrow b \leftarrow c, a \leftrightarrow c}{a<-!->b-->c, a<-->c}}
\item{\code{030T}}{ \eqn{a \rightarrow b \leftarrow c, a \rightarrow c}{a-->b<--c, a-->c}}
\item{\code{030C}}{ \eqn{a \leftarrow b \leftarrow c, a \rightarrow c}{a<--b<--c, a-->c}}
\item{\code{201}}{ \eqn{a \leftrightarrow b \not\leftrightarrow c, a \leftrightarrow c}{a<-->b<-!->c, a<-->c}}
\item{\code{120D}}{ \eqn{a \leftarrow b \rightarrow c, a \leftrightarrow c}{a<--b-->c, a<-->c}}
\item{\code{120U}}{ \eqn{a \rightarrow b \leftarrow c, a \leftrightarrow c}{a-->b<--c, a<-->c}}
\item{\code{120C}}{ \eqn{a \rightarrow b \rightarrow c, a \leftrightarrow c}{a-->b-->c, a<-->c}}
\item{\code{210}}{ \eqn{a \rightarrow b \leftrightarrow c, a \leftrightarrow c}{a-->b<-->c, a<-->c}}
\item{\code{300}}{ \eqn{a \leftrightarrow b \leftrightarrow c, a \leftrightarrow c}{a<-->b<-->c, a<-->c}}
}

These codes are returned by \code{triad.classify} as strings.  In the undirected case, only four triad states are possible (corresponding to the number of edges in the triad).  These are evaluated for \code{mode=="graph"}, with the return value being the number of edges.
}
\value{
A string containing the triad classification, or \code{NA} if one or more edges were missing
}
\references{ Davis, J.A. and Leinhardt, S.  (1972).  ``The Structure of Positive Interpersonal Relations in Small Groups.''  In J. Berger (Ed.), \emph{Sociological Theories in Progress, Volume 2}, 218-251.  Boston: Houghton Mifflin.

Wasserman, S., and Faust, K.  (1994).  \emph{Social Network Analysis: Methods and Applications.}  Cambridge: Cambridge University Press.}
\author{ Carter T. Butts  \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ }


\seealso{ \code{\link{triad.census}}, \code{\link{gtrans}} }

\section{Warning }{Valued data and/or loops may cause strange behavior with this routine.  Dichotomize/remove loops first.}

\examples{
#Generate a random graph
g<-rgraph(10)

#Classify the triads (1,2,3) and (2,3,4)
triad.classify(g,tri=c(1,2,3))
triad.classify(g,tri=c(1,2,3))

#Plot the triads in question
gplot(g[1:3,1:3])
gplot(g[2:4,2:4])
}
\keyword{ math }% __ONLY ONE__ keyword per line
\keyword{ graphs }
