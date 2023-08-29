\name{probplotparam.nor}
\alias{probplotparam.nor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normal Probability Plot Parameters
}
\description{
Computes least-squares parameter estimates for normal probability distribution for preparation for normal probability plotting.
}
\usage{
probplotparam.nor(xi, F)
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
Normal distribution is unbounded (\eqn{-\infty < x < \infty}) consisting of a normal mean parameter \eqn{\mu} and a normal standard deviation parameter \eqn{\sigma}.

  PDF:
  \deqn{f(x) = \frac{1}{\sigma \sqrt{2 \pi}} \exp[-\frac{1}{2} (\frac{x - \mu}{\sigma})^2]}

  CDF:

  \deqn{F(x) = \Phi[\frac{x - \mu}{\sigma}]
  }
}
\value{
Returns several parameters for a single probability plot from a set of data \code{xi}.

\item{x_p at (0.1\% and 99.9\%)}{
  Returns the \eqn{x_{p\;at\;0.1\%}} and \eqn{x_{p\;at\;99.9\%}} based on a normal distribution.  \code{probplot.nor} uses this output to plot the best-fit line on the normal probability plots.
}
\item{fcB}{
Returns the transformed y-axis parameter for probabilities \eqn{0.1}\% and \eqn{99.9}\%.  For the normal plot the y-axis is \eqn{\Phi^{-1} (p)}.  \code{probplot.nor} uses this output to plot the best-fit line on the normal probability plots.
}
\item{normresults}{
  The least squares estimates for the normal mean parameter \eqn{\mu} and the normal standard deviation parameter \eqn{\sigma}.
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
probplotparam.nor(x_R[,1], x_R[,3])

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
