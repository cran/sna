\name{coleman}
\Rdversion{1.1}
\alias{coleman}
\docType{data}
\title{
Coleman's High School Friendship Data
}
\description{
James Coleman (1964) reports research on self-reported friendship ties among 73 boys in a small high school in Illinois over the 1957-1958 academic year.  Networks of reported ties for all 73 informants are provided for two time points (fall and spring). 
}
\usage{data(coleman)}
\format{
  An adjacency array containing two directed, unvalued 73-node networks:

   \tabular{rlll}{
    [1,,] \tab Fall  \tab binary matrix \tab Friendship for Fall, 1957\cr
    [2,,] \tab Spring  \tab binary matrix \tab Friendship for Spring, 1958\cr
   } 
}
\details{
Both networks reflect answers to the question, \dQuote{What fellows here in school do you go around with most often?} with the presence of an \eqn{(i,j,k)} edge indicating that \eqn{j} nominated \eqn{k} in time period \eqn{i}.  The data are unvalued and directed; although the self-reported ties are highly reciprocal, unreciprocated nominations are possible.

It should be noted that, although this data is usually described as \dQuote{friendship,} the sociometric item employed might be more accurately characterized as eliciting \dQuote{frequent elective interaction.}  This should be borne in mind when interpreting this data.
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
Coleman, J. S.  (1964).  \emph{Introduction to Mathermatical Sociology.} New York: Free Press.
}
\examples{
data(coleman)

#Plot showing edges by time point
gplot(coleman[1,,]|coleman[2,,],edge.col=2*coleman[1,,]+3*coleman[2,,])
}
\keyword{datasets}
