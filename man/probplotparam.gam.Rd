\name{probplotparam.gam}
\alias{probplotparam.gam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gamma Probability Plot Parameters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
probplotparam.gam(xi, F)
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
  The gamma distribution is positively bounded (\eqn{x \geq 0}) consisting of a Gamma shape parameter \eqn{\alpha} and a Gamma scale parameter \eqn{\beta} (sometimes defined as a rate or frequency where rate \eqn{\beta = \frac{1}{\theta}} and \eqn{\theta} is the average time between events).

  PDF:
  \deqn{f(x) = \frac{1}{\Gamma(\alpha)} \beta^\alpha x^{\alpha-1} \exp(-\beta x)}

  CDF:

  \deqn{F(x) = \frac{1}{\Gamma(\alpha)} \int_{0}^{\beta x} t^{\alpha-1} \exp(-t) dt}
}
\value{
Returns several parameters for a single probability plot from a set of data \code{xi}.

\item{x_p at (0.1\% and 99.9\%)}{
  Returns the \eqn{x_{p\;at\;0.1\%}} and \eqn{x_{p\;at\;99.9\%}} based on a Gamma distribution.  \code{probplot.gam} uses this output to plot the best-fit line on the Gamma probability plots.
}
\item{fcB}{
Returns the transformed y-axis parameter for probabilities \eqn{0.1}\% and \eqn{99.9}\%.  For the Gamma plot the y-axis is \eqn{\Phi_{gamma}^{-1}(\alpha,F)} where \eqn{\Phi_{gamma}^{-1}(\alpha,F)} is the inverse of the incomplete gamma function for a given shape parameter \eqn{\alpha}.

\deqn{\Phi_{gamma} = \frac{1}{\Gamma(\alpha)} \int_{0}^{x} t^{\alpha-1} \exp(-t) dt}

\code{probplot.gam} uses this output to plot the best-fit line on the Weibull probability plots.
}
\item{gamresults}{
  The least squares estimates for the Gamma shape parameter \eqn{\alpha} and the Gamma scale parameter \eqn{\beta}.
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
probplotparam.gam(x_R[,1], x_R[,3])
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
