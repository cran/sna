\name{qaptest}
\alias{qaptest}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Perform Quadratic Assignment Procedure (QAP) Hypothesis Tests for Graph-Level Statistics }
\description{
   \code{qaptest} tests an arbitrary graph-level statistic (computed on \code{dat} by \code{FUN}) against a QAP null hypothesis, via Monte Carlo simulation of likelihood quantiles.  Note that fair amount of flexibility is possible regarding QAP tests on functions of such statistics (see an equivalent discussion with respect to CUG null hypothesis tests in Anderson et al. (1999)).  See below for more details. 
}
\usage{
qaptest(dat, FUN, reps=1000, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dat}{ graphs to be analyzed.  Though one could in principle use a single graph, this is rarely if ever sensible in a QAP-test context.}
  \item{FUN}{ function to generate the test statistic.  \code{FUN} must accept \code{dat} and the specified \code{g} arguments, and should return a real number. }
  \item{reps}{ integer indicating the number of draws to use for quantile estimation.  Note that, as for all Monte Carlo procedures, convergence is slower for more extreme quantiles.  By default, \code{reps}=1000. }
  \item{\dots}{ additional arguments to \code{FUN}. }
}
\details{
  The null hypothesis of the QAP test is that the observed graph-level statistic on graphs \eqn{G_1,G_2,\dots} was drawn from the distribution of said statistic evaluated (uniformly) on the set of all relabelings of \eqn{G_1,G_2,\dots}.  Pragmatically, this test is performed by repeatedly (randomly) relabeling the input graphs, recalculating the test statistic, and then evaluating the fraction of draws greater than or equal to (and less than or equal to) the observed value.  This accumulated fraction approximates the integral of the distribution of the test statistic over the set of unlabeled input graphs.

The \code{qaptest} procedure returns a \code{qaptest} object containing the estimated likelihood (distribution of the test statistic under the null hypothesis), the observed value of the test statistic on the input data, and the one-tailed p-values (estimated quantiles) associated with said observation.  As usual, the (upper tail) null hypothesis is rejected for significance level alpha if p>=observation is less than alpha (or p<=observation, for the lower tail); if the hypothesis is undirected, then one rejects if either p<=observation or p>=observation is less then alpha/2.  Standard caveats regarding the use of null hypothesis testing procedures are relevant here: in particular, bear in mind that a significant result does not necessarily imply that the likelihood ratio of the null model and the alternative hypothesis favors the latter.

In interpreting a QAP test, it is important to bear in mind the nature of the QAP null hypothesis.  The QAP test should \emph{not} be interpreted as evaluating underlying structural differences; indeed, QAP is more accurately understood as testing differences induced by a particular vertex labeling \emph{controlling for} underlying structure.  Where there is substantial automorphism in the underling structures, QAP will tend to given non-significant results.  (In fact, it is \emph{impossible} to obtain a one-tailed significance level in excess of \eqn{\max_{g \in \{G,H\}} \frac{|Aut(g)|}{|Perm(g)|}}{max_[g in \{G,H\}] |Aut(g)|/|Perm(g)|} when using a QAP test on a bivariate graph statistic \eqn{f(G,H)}, where Aut(g) and Perm(g) are the automorphism and permutation groups on g, respectively.  This follows from the fact that all members of Aut(g) will induce the same values of \eqn{f()}.)  By turns, significance under QAP does not necessarily imply that the observed structural relationship is unusual relative to what one would expect from typical structures with (for instance) the sizes and densities of the graphs in question.  In contexts in which one's research question implies a particular labeling of vertices (e.g., "within this group of individuals, do friends also tend to give advice to one another"), QAP can be a very useful way of ruling out spurious structural influences (e.g., some respondents tend to indiscriminately nominate many people (without regard to whom), resulting in a structural similarity which has nothing to do with the identities of those involved).  Where one's question does not imply a labeled relationship (e.g., is the \emph{shape} of this group's friendship network similar to that of its advice network), the QAP null hypothesis is inappropriate.
}
\value{
  An object of class \code{qaptest}, containing

\item{testval}{
The observed value of the test statistic.
}
\item{dist}{
A vector containing the Monte Carlo draws. 
}
\item{pgreq}{
The proportion of draws which were greater than or equal to the observed value.
}
\item{pleeq}{
The proportion of draws which were less than or equal to the observed value.
}
}
\references{ 
Anderson, B.S.; Butts, C.T.; and Carley, K.M. (1999). ``The Interaction of Size and Density with Graph-Level Indices.'' \emph{Social Networks}, 21(3), 239-267.

Hubert, L.J., and Arabie, P.  (1989).  ``Combinatorial Data Analysis: Confirmatory Comparisons Between Sets of Matrices.''  \emph{Applied Stochastic Models and Data Analysis}, 5, 273-325.

Krackhardt, D.  (1987).  ``QAP Partialling as a Test of Spuriousness.'' \emph{Social Networks}, 9 171-186.

Krackhardt, D.  (1988).  ``Predicting With Networks: Nonparametric Multiple Regression Analyses of Dyadic Data.''  \emph{Social Networks}, 10, 359-382. 
}
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{cugtest}} }

\examples{

#Generate three graphs
g<-array(dim=c(3,10,10))
g[1,,]<-rgraph(10)
g[2,,]<-rgraph(10,tprob=g[1,,]*0.8)
g[3,,]<-1; g[3,1,2]<-0              #This is nearly a clique

#Perform qap tests of graph correlation
q.12<-qaptest(g,gcor,g1=1,g2=2)
q.13<-qaptest(g,gcor,g1=1,g2=3)

#Examine the results
summary(q.12)
plot(q.12)
summary(q.13)
plot(q.13)
}
\keyword{ htest }%-- one or more ...
\keyword{ math }
\keyword{ graphs }
