\name{probplotparam.exp2P}
\alias{probplotparam.exp2P}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two-Parameter Exponential Probability Plot Parameters
}
\description{
Computes least-squares parameter estimates for two-parameter exponential probability distribution for preparation for two-parameter exponential probability plotting.
}
\usage{
probplotparam.exp2P(xi, F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xi}{
  A vector consisting of the failure or primary event data in a given set.  Data must consist of real and finite numbers.
  }
  \item{F}{
  A vector of non-parametric estimates for failure probability \eqn{F} at data \code{xi}.
  }
}
\details{
  Two-parameter exponential distribution is positively bounded (\eqn{x \geq 0}) consisting of a two-parameter exponential rate parameter \eqn{\sigma} and a two-parameter exponential location parameter \eqn{\theta}.

  PDF:
  \deqn{f(x) = \frac{1}{\sigma} \exp[- \frac{x - \theta}{\sigma}]}

  CDF:

  \deqn{F(x) = 1 - \exp[- \frac{x - \theta}{\sigma}]}
}
\value{
Returns several parameters for a single probability plot from a set of data \code{xi}.

\item{x_p at (0.1\% and 99.9\%)}{
  Returns the \eqn{x_{p\;at\;0.1\%}} and \eqn{x_{p\;at\;99.9\%}} based on a two-parameter exponential distribution.  \code{probplot.exp2P} uses this output to plot the best-fit line on the two-parameter exponential probability plots.
}
\item{fcB}{
Returns the transformed y-axis parameter for probabilities \eqn{0.1}\% and \eqn{99.9}\%.  For the two-parameter exponential plot the y-axis is \eqn{- ln(R)}.  \code{probplot.exp2P} uses this output to plot the best-fit line on the two-parameter exponential probability plots.
}
\item{exp2Presults}{
  The least squares estimates for the two-parameter exponential rate parameter \eqn{\sigma} and the two-parameter exponential location parameter \eqn{\theta}.
}
\item{R2}{
  Coeefficient of determination \eqn{R^2}.
}
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
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
## Example 1
Xi1 <- c(70,71,75,78,89)
Rc1 <- c(80,80,84)
x_R <- plotposit.blom(i = rankcalc(xi = Xi1, rc = Rc1), xi = Xi1, rc = Rc1)
probplotparam.exp2P(x_R[,1], x_R[,3])

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
