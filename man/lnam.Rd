\name{lnam}
\alias{lnam}
\alias{coef.lnam}
\alias{se.lnam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit a Linear Network Autocorrelation Model }
\description{
  \code{lnam} is used to fit linear network autocorrelation models.  These include standard OLS as a special case, although \code{\link{lm}} is to be preferred for such analyses.
}
\usage{
lnam(y, x = NULL, W1 = NULL, W2 = NULL, theta.seed = NULL, 
    null.model = c("meanstd", "mean", "std", "none"), method = "BFGS", 
    control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ a vector of responses. }
  \item{x}{ a vector or matrix of covariates; if the latter, each column should contain a single covariate. }
  \item{W1}{ a (possibly valued) adjacency matrix on the elements of \code{y}. }
  \item{W2}{ another (possibly valued) adjacency matrix on the elements of \code{y}. }
  \item{theta.seed}{ an optional seed value for \code{\link{optim}}. }
  \item{null.model}{ the null model to be fit; must be one of \code{"meanstd"}, \code{"mean"}, \code{"std"}, or \code{"none"}. }
  \item{method}{ method to be used with \code{\link{optim}}. }
  \item{control}{ optional control parameters for \code{\link{optim}}. }
}
\details{
  \code{lnam} fits the linear network autocorrelation model given by
  
  \deqn{y = \rho_1 W_1 y + X \beta + e, \quad e=\rho_2 W_2 e + \nu}{%
  y = rho1 * W1 \%*\% y + X \%*\% beta + e, e = rho2 * W2 \%*\% e + nu}

  where \eqn{y} is a vector of responses, \eqn{X} is a covariate matrix, \eqn{W_1}{W1} and \eqn{W_2}{W2} are (possibly valued) adjacency matrices, and \eqn{\nu \sim N(0,\sigma^2)}{nu ~ Norm(0,sigma^2)}.  Intuitively, \eqn{\rho_1}{rho1} is an ``AR''-like parameter (parameterizing the autoregression of each \eqn{y} value on its neighbors in \eqn{W_1}{W1}) while \eqn{\rho_2}{rho2} is an ``MA''-like parameter (parameterizing the autocorrelation of each \emph{disturbance} in \eqn{y} on its neighbors in \eqn{W_2}{W2}).  In general, the two models are distinct, and either or both effects may be selected by including the appropriate matrix arguments.

Model parameters are estimated by maximum likelihood, and asymptotic standard errors are provided as well; all of the above (and more) can be obtained by means of the appropriate \code{print} and \code{summary} methods.  A plotting method is also provided, which supplies fit basic diagnostics for the estimated model.  For purposes of comparison, fits may be evaluated against one of four null models:

\enumerate{
\item \code{meanstd}: mean and standard deviation estimated (default).
\item \code{mean}: mean estimated; standard deviation assumed equal to 1.
\item \code{std}: standard deviation estimated; mean assumed equal to 0.
\item \code{none}: no parameters estimated; data assumed to be drawn from a standard normal density.
}

The default setting should be appropriate for the vast majority of cases, although the others may have use when fitting ``pure'' autoregressive models (e.g., without covariates).  Although a major use of the \code{lnam} is in controlling for network autocorrelation within a regression context, the model is subtle and has a variety of uses.  (See the references below for suggestions.)
}
\value{
  An object of class \code{"lnam"} containing the following elements:

  \item{y}{the response vector used.}
  \item{x}{if supplied, the coefficient matrix.}
  \item{W1}{if supplied, the W1 matrix.}
  \item{W2}{if supplied, the W2 matrix.}
  \item{model}{a code indicating the model terms fit.}
  \item{infomat}{the estimated Fisher information matrix for the fitted model.}
  \item{acvm}{the estimated asymptotic covariance matrix for the model parameters.}
  \item{null.model}{a string indicating the null model fit.}
  \item{lnlik.null}{the log-likelihood of y under the null model.}
  \item{df.null.resid}{the residual degrees of freedom under the null model.}
  \item{df.null}{the model degrees of freedom under the null model.}
  \item{null.param}{parameter estimates for the null model.}
  \item{lnlik.model}{the log-likelihood of y under the fitted model.}
  \item{df.model}{the model degrees of freedom.}
  \item{df.residual}{the residual degrees of freedom.}
  \item{df.total}{the total degrees of freedom.}
  \item{rho1}{if applicable, the MLE for rho1.}
  \item{rho1.se}{if applicable, the asymptotic standard error for rho1.}
  \item{rho2}{if applicable, the MLE for rho2.}
  \item{rho2.se}{if applicable, the asymptotic standard error for rho2.}
  \item{sigma}{the MLE for sigma.}
  \item{sigma.se}{the standard error for sigma}
  \item{beta}{if applicable, the MLE for beta.}
  \item{beta.se}{if applicable, the asymptotic standard errors for beta.}
  \item{fitted.values}{the fitted mean values.}
  \item{residuals}{the residuals (response minus fitted); note that these correspond to \eqn{\hat{e}}{e-hat} in the model equation, not \eqn{\hat{\nu}}{nu-hat}.}
  \item{disturbances}{the estimated disturbances, i.e., \eqn{\hat{\nu}}{nu-hat}.}
  \item{call}{the matched call.}
}
\references{ 
Leenders, T.Th.A.J.  (2002)  ``Modeling Social Influence Through Network Autocorrelation: Constructing the Weight Matrix''  \emph{Social Networks}, 24(1), 21-47. 

Anselin, L.  (1988)  \emph{Spatial Econometrics: Methods and Models.}  Norwell, MA: Kluwer.
 }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ Actual optimization is performed by calls to \code{\link{optim}}.  Information on algorithms and control parameters can be found via the appropriate man pages. }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{lm}}, \code{\link{optim}} }
\examples{
#Construct a simple, random example:
w1<-rgraph(100)               #Draw the AR matrix
w2<-rgraph(100)               #Draw the MA matrix
x<-matrix(rnorm(100*5),100,5) #Draw some covariates
r1<-0.2                       #Set the model parameters
r2<-0.1
sigma<-0.1
beta<-rnorm(5)
#Assemble y from its components:
nu<-rnorm(100,0,sigma)          #Draw the disturbances
e<-qr.solve(diag(100)-r2*w2,nu) #Draw the effective errors
y<-qr.solve(diag(100)-r1*w1,x\%*\%beta+e)  #Compute y

#Now, fit the autocorrelation model:
fit<-lnam(y,x,w1,w2)
summary(fit)
plot(fit)
}
\keyword{ regression }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
