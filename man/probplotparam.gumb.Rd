\name{probplotparam.gumb}
\alias{probplotparam.gumb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gumbel Probability Plot Parameters
}
\description{
Computes least-squares parameter estimates for Gumbel probability distribution for preparation for Gumbel probability plotting.
}
\usage{
probplotparam.gumb(xi, R)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xi}{
  A vector consisting of the failure or primary event data in a given set.  Data must consist of real and finite numbers.
  }
  \item{R}{
  A vector of non-parametric estimates for reliability \eqn{R} at data \code{xi}.
  }
}

\details{
The Gumbel distribution is unbounded (\eqn{-\infty < x < \infty}) consisting of a Gumbel location parameter \eqn{\mu} and a Gumbel scale parameter \eqn{\sigma}.

  PDF:
  \deqn{f(x) = \frac{1}{\sigma}\exp[\frac{x - \mu}{\sigma} - \exp(\frac{x - \mu}{\sigma})]}

  CDF:

  \deqn{F(x) = 1 - \exp[-\exp(\frac{x - \mu}{\sigma})]}
}

\value{
Returns several parameters for a single probability plot from a set of data \code{xi}.

\item{x_p at (0.1\% and 99.9\%)}{
  Returns the \eqn{x_{p\;at\;0.1\%}} and \eqn{x_{p\;at\;99.9\%}} based on a Gumbel distribution.  \code{probplot.gumb} uses this output to plot the best-fit line on the Gumbel probability plots.
}
\item{fcB}{
Returns the transformed y-axis parameter for probabilities \eqn{0.1}\% and \eqn{99.9}\%.  For the Gumbel plot the y-axis is \eqn{\ln[-\ln(R)]}.  \code{probplot.gumb} uses this output to plot the best-fit line on the Gumbel probability plots.
}
\item{gumbresults}{
  The least squares estimates for the Gumbel location parameter \eqn{\mu} and the Gumbel scale parameter \eqn{\sigma}.
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
probplotparam.gumb(x_R[,1], x_R[,3])
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
