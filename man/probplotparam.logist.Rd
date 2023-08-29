\name{probplotparam.logist}
\alias{probplotparam.logist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Logistic Probability Plot Parameters
}
\description{
Computes least-squares parameter estimates for logistic probability distribution for preparation for logistic  probability plotting.
}
\usage{
probplotparam.logist(xi,F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xi}{
  A vector consisting of the failure or primary event data in a given set.  Data must consist of real and finite numbers.
  }
  \item{F}{
  A vector of non-parametric estimates for unreliability \eqn{F} at data \code{xi}.
  }
}
\details{
The Logistic distribution is unbounded (\eqn{-\infty < x < \infty}) and consists of a logistic location parameter \eqn{\mu} and a logistic scale parameter \eqn{\sigma}.

  PDF:
  \deqn{f(x) = \frac{\exp(\frac{x - \mu}{\sigma})}{\sigma [1+\exp(\frac{x - \mu}{\sigma})]^2}}

  CDF:

  \deqn{F(x) = \frac{\exp(\frac{x - \mu}{\sigma})}{1+\exp(\frac{x - \mu}{\sigma})}}
}
\value{
Returns several parameters for a single probability plot from a set of data \code{xi}.

\item{x_p at (0.1\% and 99.9\%)}{
  Returns the \eqn{x_{p\;at\;0.1\%}} and \eqn{x_{p\;at\;99.9\%}} based on a logistic distribution.  \code{probplot.logist} uses this output to plot the best-fit line on the logistic probability plots.
}
\item{fcB}{
Returns the transformed y-axis parameter for probabilities \eqn{0.1}\% and \eqn{99.9}\%.  For the logistic plot the y-axis is \eqn{\ln (\frac{1}{F} - 1)}.  \code{probplot.logist} uses this output to plot the best-fit line on the logistic probability plots.
}
\item{logistresults}{
  The least squares estimates for the logistic location parameter \eqn{\mu} and a logistic scale parameter \eqn{\sigma}.
}
\item{R2}{
  Coefficient of determination \eqn{R^2}.
}
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2023)
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
probplotparam.logist(x_R[,1], x_R[,2])
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
