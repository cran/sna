\name{pstar}
\alias{pstar}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Fit a p*/ERG Model Using a Logistic Approximation }
\description{
  Fits a p*/ERG model to the graph in \code{dat} containing the effects listed in \code{effects}.  The result is returned as a \code{\link{glm}} object.
}
\usage{
pstar(dat, effects=c("choice", "mutuality", "density", "reciprocity",
    "stransitivity", "wtransitivity", "stranstri",  "wtranstri", 
    "outdegree", "indegree", "betweenness", "closeness", 
    "degcentralization", "betcentralization", "clocentralization",
    "connectedness", "hierarchy", "lubness", "efficiency"), 
    attr=NULL, memb=NULL, diag=FALSE, mode="digraph")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ a single graph }
  \item{effects}{ a vector of strings indicating which effects should be fit. }
  \item{attr}{ a matrix whose columns contain individual attributes (one row per vertex) whose differences should be used as supplemental predictors. }
  \item{memb}{ a matrix whose columns contain group memberships whose categorical similarities (same group/not same group) should be used as supplemental predictors.}
  \item{diag}{ a boolean indicating whether or not diagonal entries (loops) should be counted as meaningful data. }
  \item{mode}{ \code{"digraph"} if \code{dat} is directed, else \code{"graph"} }
}
\details{
The Exponential Family-Random Graph Model (ERGM) family, referred to as \dQuote{p*} in older literature, is an exponential family specification for network data.  In this specification, it is assumed that 
\deqn{p(G=g) \propto \exp(\beta_0 \gamma_0(g) + \beta_1 \gamma_1(g) + \dots)}{%
p(G=g) propto exp(beta_0 gamma_0(g) + beta_1 gamma_1(g) + \dots)}
for all g, where the betas represent real coefficients and the gammas represent functions of g.  Unfortunately, the unknown normalizing factor in the above expression makes evaluation difficult in the general case.  One solution to this problem is to operate instead on the edgewise log odds; in this case, the ERGM/p* MLE can be approximated by a logistic regression of each edge on the \emph{differences} in the gamma scores induced by the presence and absence of said edge in the graph (conditional on all other edges).  It is this approximation (known as autologistic regression, or maximum pseudo-likelihood estimation) that is employed here.  

Note that ERGM modeling is considerably more advanced than it was when this function was created, and estimation by MPLE is now used only in special cases.  Guidelines for model specification and assessment have also evolved.  The \code{ergm} package within the \code{statnet} library reflects the current state of the art, and use of the \code{ergm()} function in said library is highly recommended.  This function is retained primarily as a legacy tool, for users who are nostalgic for 2000-vintage ERGM (\dQuote{p*}) modeling experience.  Caveat emptor.

Using the \code{effects} argument, a range of different potential parameters can be estimated.  The network measure associated with each is, in turn, the edge-perturbed difference in:
\enumerate{
\item \code{choice}: the number of edges in the graph (acts as a constant)

\item \code{mutuality}: the number of reciprocated dyads in the graph

\item \code{density}: the density of the graph

\item \code{reciprocity}: the edgewise reciprocity of the graph

\item \code{stransitivity}: the strong transitivity of the graph

\item \code{wtransitivity}: the weak transitivity of the graph

\item \code{stranstri}: the number of strongly transitive triads in the graph

\item \code{wtranstri}: the number of weakly transitive triads in the graph

\item \code{outdegree}: the outdegree of each actor (|V| parameters)

\item \code{indegree}: the indegree of each actor (|V| parameters)

\item \code{betweenness}: the betweenness of each actor (|V| parameters)

\item \code{closeness}: the closeness of each actor (|V| parameters)

\item \code{degcentralization}: the Freeman degree centralization of the graph

\item \code{betcentralization}: the betweenness centralization of the graph

\item \code{clocentralization}: the closeness centralization of the graph

\item \code{connectedness}: the Krackhardt connectedness of the graph

\item \code{hierarchy}: the Krackhardt hierarchy of the graph

\item \code{efficiency}: the Krackhardt efficiency of the graph

\item \code{lubness}: the Krackhardt LUBness of the graph
}
(Note that some of these do differ somewhat from the common specifications employed in the older p* literature, e.g. quantities such as density and reciprocity are computed as per the \code{\link{gden}} and \code{\link{grecip}} functions rather than via the unnormalized "choice" and "mutual" quantities that were generally used.)  \emph{Please do not attempt to use all effects simultaneously!!!}  In addition to the above, the user may specify a matrix of individual attributes whose absolute dyadic differences are to be used as predictors, as well as a matrix of individual memberships whose dyadic categorical similarities (same/different) are used in the same manner.

Although the ERGM framework is quite versatile in its ability to accommodate a range of structural predictors, it should be noted that the \emph{substantial} collinearity of many of the terms provided here can lead to very unstable model fits.  Measurement and specification errors compound this problem, as does the use of the MPLE; thus, it is somewhat risky to use \code{pstar} in an exploratory capacity (i.e., when there is little prior knowledge to constrain choice of parameters).  While raw instability due to multicollinearity should decline with graph size, improper specification will still result in biased coefficient estimates so long as an omitted predictor correlates with an included predictor.  Moreover, many models created using these effects are at risk of degeneracy, which is difficult to assess without simulation-based model assessment.  Caution is advised - or, better, use of the \code{ergm} package.
}
\value{
A \code{\link{glm}} object
}
\references{ Anderson, C.; Wasserman, S.; and Crouch, B. (1999).  ``A p* Primer:  Logit Models for Social Networks.  \emph{Social Networks,} 21,37-66.

Holland, P.W., and Leinhardt, S. (1981).  ``An Exponential Family of Probability Distributions for Directed Graphs.'' \emph{Journal of the American statistical Association}, 81, 51-67.

Wasserman, S., and Pattison, P. (1996).  ``Logit Models and Logistic Regressions for Social Networks:  I.  An introduction to Markov Graphs and p*.''  \emph{Psychometrika,} 60, 401-426.
 }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ This is a legacy function - use of the \code{ergm} package is now strongly advised.}

\section{WARNING }{Estimation of p* models by maximum pseudo-likelihood is now known to be a dangerous practice.  Use at your own risk.} 

\seealso{ \code{\link{eval.edgeperturbation}} }

\examples{
\dontrun{
#Create a graph with expansiveness and popularity effects
in.str<-rnorm(20,0,3)
out.str<-rnorm(20,0,3)
tie.str<-outer(out.str,in.str,"+")
tie.p<-apply(tie.str,c(1,2),function(a){1/(1+exp(-a))})
g<-rgraph(20,tprob=tie.p)

#Fit a model with expansiveness only
p1<-pstar(g,effects="outdegree")
#Fit a model with expansiveness and popularity
p2<-pstar(g,effects=c("outdegree","indegree"))
#Fit a model with expansiveness, popularity, and mutuality
p3<-pstar(g,effects=c("outdegree","indegree","mutuality"))

#Compare the model AICs -- use ONLY as heuristics!!!
extractAIC(p1)
extractAIC(p2)
extractAIC(p3)
}
}
\keyword{ models }
\keyword{ multivariate }
\keyword{ regression }%-- one or more ...
\keyword{ graphs }
