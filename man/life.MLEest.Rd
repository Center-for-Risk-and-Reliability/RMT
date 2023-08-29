\name{life.MLEest}
\alias{life.MLEest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Likelihood Life-Model Estimator
}
\description{
Computes maximum likelihood estimate of Life model parameters.
}
\usage{
life.MLEest(LSQest, dist, TTF, Tc, confid, sided)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LSQest}{
The least-squares or another initial estimate for the maximum likelihood process to begin with.  Make sure the order matches the model being entered and the order ouput by the \code{probplot} code.  That is, if the Weibull paramters are the form c(\eqn{\alpha, \beta}), then the input must also be some form of c(\eqn{\alpha, \beta}).
}
  \item{dist}{
  The distribution that represents time or life.
  \enumerate{
  \item \strong{Weibull} (Enter as "Weibull".  See more on \code{\link{probplot.wbl}}.)
  \item \strong{Three Parameter Weibull} (Enter as "3PWeibull".  See more on \code{\link{probplot.wbl3P}}.)
  \item \strong{Exponential} (Enter as "Exponential".  See more on \code{\link{probplot.exp}}.)
  \item \strong{Two Parameter Exponential} (Enter as "2PExponential".  See more on \code{\link{probplot.exp2P}}.)
  \item \strong{Normal} (Enter as "Normal".  See more on \code{\link{probplot.nor}}.)
  \item \strong{Lognormal} (Enter as "Lognormal".  See more on \code{\link{probplot.logn}}.)
  }
}
  \item{TTF}{
A vector consisting of the failure or primary event data in a given set.  Output[[1]] of \code{sort.xircstressdata}.
}
  \item{Tc}{
A vector consisting of the right censored data of a given set of data (NOTE: If no censored data exisits, this will return \code{NULL}).  Output[[2]] of \code{sort.xircdata}.
}
  \item{confid}{
The confidence in the maximum likelihood parameter estimates entered as a probability between \eqn{0} and \eqn{1}.  Default is \eqn{0.95} (or \eqn{95\%})
}
  \item{sided}{
Indicates whether the limit calculations are based on a two-sided confidence or one-sided confidence.  The three options for this heading are \code{"twosided"} for two-sided confidence (default), \code{"onesidedlow"} for one-sided (low) confidence, and \code{"onesidedhigh"} for one-sided (high) confidence.
}
}
\details{
The maximum likelihood estimation is performed on the negated log-likelihood of the entered life model using the built-in \code{\link{nlm}} function, a non-linear minimization script that is based on a Newton-type algorithm.  Certain output will not respond well to a poor initial estimate \code{LSQest}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{inv.fish}{
  The variance-covariance matrix as computed by non-linear minimization (in order of \code{LSQest}).  Output[[1]]
}
\item{conflim}{
  Vector-by-vector output of the confidence bounds for each parameter (in order of \code{LSQest}).  Output[[2]]
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Dr. Reuel Smith (2021-2022)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1: Single data set.  NOTE: Data must include a stress column because the group stress definition is one of the output.
mat0 <- matrix.failcen(xi = c(70,71,75,78,89), rc = c(80,80,84), nx = 8)
temp <- rep(373,length(mat0[,1]))
mat1 <- matrix(c(mat0[,1], mat0[,2], temp), nrow = length(mat0[,1]), ncol = 3, byrow = FALSE)
est1 <- c(probplot.logn(mat1, "Blom", "Time (hours)")[[2]])
life.MLEest(est1,"Weibull",c(70,71,75,78,89),c(80,80,84),0.90)

## Example 2: Same example but with the one-sided upper confidence computed
life.MLEest(est1,"Weibull",c(70,71,75,78,89),c(80,80,84),0.80,"onesidedhigh")

## Example 3: Dual Stress values (Temperature in Kelvin and Relative Humidity (RH))
time <- c(102, 115, 151, 196, 210, 184, 193, 244, 265, 281, 298, 311, 337, 347, 405)
cens <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
temp <- c(393, 393, 393, 393, 393, 353, 353, 353, 353, 353, 353, 353, 353, 353, 353)
RH <- c(0.6, 0.6, 0.6, 0.6, 0.6, 0.85, 0.85, 0.85, 0.85, 0.85, 0.6, 0.6, 0.6, 0.6, 0.6)
mat3 <- matrix(c(time, cens, temp, RH), nrow = length(time), ncol = 4, byrow = FALSE)
probplot.logn(mat3, "KaplanMeier", "Time (hours)")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
