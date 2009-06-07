\name{potscalered.mcmc}
\alias{potscalered.mcmc}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute Gelman and Rubin's Potential Scale Reduction Measure for a Markov Chain Monte Carlo Simulation}
\description{
  Computes Gelman and Rubin's (simplified) measure of scale reduction for draws of a single scalar estimand from parallel MCMC chains. 
}
\usage{
potscalered.mcmc(psi)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{psi}{ An nxm matrix, with columns corresponding to chains and rows corresponding to iterations. }
}
\details{
  The Gelman and Rubin potential scale reduction (\eqn{\sqrt{\hat{R}}}{sqrt(hat(R))}) provides an ANOVA-like comparison of the between-chain to within-chain variance on a given scalar estimand; the disparity between these gives an indication of the extent to which the scale of the simulated distribution can be reduced via further sampling.  As the parallel chains converge \eqn{\sqrt{\hat{R}}}{sqrt(hat(R))} approaches 1 (from above), and it is generally recommended that values of 1.2 or less be obtained before a series of draws can be considered well-mixed.  (Even so, one should ideally examine other indicators of chain mixing, and verify that the properties of the draws are as they should be.  There is currently no fool-proof way to verify burn-in of an MCMC, but using multiple indicators should help one avoid falling prey to the idiosyncrasies of any one index.)

Note that the particular estimators used in the \eqn{\sqrt{\hat{R}}}{sqrt(hat(R))} formulation are based on normal-theory results, and as such have been criticized vis a vis their behavior on other distributions.  Where simulating distributions whose properties differ greatly from the normal, an alternative form of the measure using robust measures of scale (e.g., the IQR) may be preferable.
}
\value{
The potential scale reduction measure
}
\references{ Gelman, A.; Carlin, J.B.; Stern, H.S.; and Rubin, D.B.  (1995).  \emph{Bayesian Data Analysis.}  London: Chapman and Hall.

Gelman, A., and Rubin, D.B.  (1992).  ``Inference from Iterative Simulation Using Multiple Sequences.''  \emph{Statistical Science,} 7, 457-511. }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

 %~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{ \code{\link{bbnam}}}

%\examples{
%
%}
\keyword{ models }%-- one or more ...
\keyword{univar}
