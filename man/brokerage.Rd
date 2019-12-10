\name{brokerage}
\alias{brokerage}
\alias{summary.brokerage}
\alias{print.summary.brokerage}
\alias{brokerage_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Perform a Gould-Fernandez Brokerage Analysis }
\description{
  Performs the brokerage analysis of Gould and Fernandez on one or more input graphs, given a class membership vector.
}
\usage{
brokerage(g, cl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{ one or more input graphs. }
  \item{cl}{ a vector of class memberships. }
}
\details{
Gould and Fernandez (following Marsden and others) describe \emph{brokerage} as the role played by a social actor who mediates contact between two alters.  More formally, vertex \eqn{v} is a broker for distinct vertices \eqn{a} and \eqn{b} iff \eqn{a \to v \to b}{a -> v -> b} and \eqn{a \not\to b}{a -!> b}.  Where actors belong to a priori distinct groups, group membership may be used to segment brokerage roles into particular types.  Let \eqn{A \to B \to C}{A -> B -> C} denote the two-path associated with a brokerage structure, such that some vertex from group \eqn{B} brokers the connection from some vertex from group \eqn{A} to a vertex in group \eqn{C}.  The types of brokerage roles defined by Gould and Fernandez (and their accompanying two-path structures) are then defined in terms of group membership as follows:

\itemize{
  \item \eqn{w_I}: Coordinator role; the broker mediates contact between two individuals from his or her own group.  Two-path structure: \eqn{A \to A \to A}{A -> A -> A}  
  \item \eqn{w_O}: Itinerant broker role; the broker mediates contact between two individuals from a single group to which he or she does not belong.  Two-path structure: \eqn{A \to B \to A}{A -> B -> A}  
  \item \eqn{b_{OI}}: Gatekeeper role; the broker mediates an incoming contact from an out-group member to an in-group member.  Two-path structure: \eqn{A \to B \to B}{A -> B -> B}  
  \item \eqn{b_{IO}}: Representative role; the broker mediates an outgoing contact from an in-group member to an out-group member.  Two-path structure: \eqn{A \to A \to B}{A -> A -> B}  
  \item \eqn{b_O}: Liaison role; the broker mediates contact between two individuals from different groups, neither of which is the group to which he or she belongs.  Two-path structure: \eqn{A \to B \to C}{A -> B -> C}  
  \item \eqn{t}: Total (cumulative) brokerage role occupancy.  (Any of the above two-paths.)   
}

The \emph{brokerage score} for a given vertex with respect to a given role is the number of ordered pairs having the appropriate group membership(s) brokered by said vertex.  \code{brokerage} computes the brokerage scores for each vertex, given an input graph and vector of class memberships.  Aggregate scores are also computed at the graph level, which correspond to the total frequency of each role type within the network structure.  Expectations and variances of the brokerage scores conditional on size and density are computed, along with approximate \eqn{z}-tests for incidence of brokerage.  (Note that the accuracy of the normality assumption is not known in the general case; see Gould and Fernandez (1989) for details.  Simulation-based tests may be desirable as an alternative.)

}
\value{
  An object of class \code{brokerage}, containing the following elements:
  \item{raw.nli }{The matrix of observed brokerage scores, by vertex}
  \item{exp.nli }{The matrix of expected brokerage scores, by vertex}
  \item{sd.nli }{The matrix of predicted brokerage score standard deviations, by vertex}
  \item{z.nli }{The matrix of standardized brokerage scores, by vertex}
  \item{raw.gli }{The vector of observed aggregate brokerage scores}
  \item{exp.gli }{The vector of expected aggregate brokerage scores}
  \item{sd.gli }{The vector of predicted aggregate brokerage score standard deviations}
  \item{z.gli }{The vector of standardized aggregate brokerage scores}
  \item{exp.grp }{The matrix of expected brokerage scores, by group}
  \item{sd.grp }{The matrix of predicted brokerage score standard deviations, by group}
  \item{cl }{The vector of class memberships}
  \item{clid }{The original class names}
  \item{n }{The input class sizes}
  \item{N }{The order of the input network}
}
\references{ Gould, R.V. and Fernandez, R.M.  1989.  \dQuote{Structures of Mediation: A Formal Approach to Brokerage in Transaction Networks.}  \emph{Sociological Methodology,} 19: 89-126. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{triad.census}}, \code{\link{gtrans}} }
\examples{
#Draw a random network with 3 groups
g<-rgraph(15)
cl<-rep(1:3,5)

#Compute a brokerage object
b<-brokerage(g,cl)
summary(b)
}
\keyword{ graphs }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
