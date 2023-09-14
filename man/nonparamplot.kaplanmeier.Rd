\name{nonparamplot.kaplanmeier}
\alias{nonparamplot.kaplanmeier}
\title{
Kaplan-Meier Non-Parametric Plot
}
\description{
Nonparametric plot based on Kaplan-Meier Plottng Position of either failure probability or unreliability, reliability, hazard function, cumulative hazard, or probability density.
}
\usage{
nonparamplot.kaplanmeier(xi, rc, relfcn, confbnd, xlabel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xi}{
  A vector consisting of the failure or primary event data in a given set.  Data must consist of real and finite numbers.
  }
  \item{rc}{
  A vector consisting of the right censored data of a given set of data.  Data must consist of real and finite numbers.
  }
  \item{relfcn}{
Name of the non-parametric curve to plot.  For (1) \emph{failure probability} or \emph{unreliability} \eqn{F(x)} enter \strong{"unreliability"}.  For (2) \emph{reliability} \eqn{R(x)} enter \strong{"reliability"}.  For (3) \emph{hazard function} \eqn{h(x)} enter \strong{"hazard"}.  For (4) \emph{cumulative hazard function} \eqn{H(x)} enter \strong{"cumulativehazard"}.  Finally for (5) \emph{probability density function} \eqn{f(x)} enter \strong{"probabilitydensity"}.
}
  \item{confbnd}{
Confidence interval set to the nonparametric output (enter \eqn{95 \%} confidence interval as "95")
}
  \item{xlabel}{
Label for the x-axis for the plot.  If left blank, the default label is "X".
}
}
\details{
Definition of \emph{Kaplan-Meier} plotting position equation for failure probability (or unreliability) \eqn{F(x)}:
\deqn{F_i(x_i) = 1 - \Pi_{x_j \le x} (1 - \frac{d}{j})} (\emph{NOTE}: \eqn{d} is the number of failed units and \eqn{j} is the reverse rank)
For survival probability (or reliability) \eqn{R(x)}:
\deqn{R_i(x_i) = \prod_{x_j \leq x} \left(1 - \frac{d}{j}\right)}
For cumulative hazard \eqn{H(x)}:
\deqn{H_i(x_i) = -\ln\left[\prod_{x_j \leq x} \left(1 - \frac{d}{j}\right)\right]}
\deqn{\text{or }H_i(x_i) = -\sum_{x_j \leq x} \ln\left(1 - \frac{d}{j}\right)}
For hazard rate \eqn{h(x)}:
\deqn{h_i(x_i) = -\ln\left(1 - \frac{d}{j}\right)}
For probability density \eqn{f(x)}:
\deqn{f_i(x_i) = -\ln\left(1 - \frac{d}{j}\right)\times\prod_{x_j \leq x} \left(1 - \frac{d}{j}\right)}
}
\value{
Outputs the matrix of non-parametric mean point estimates of probability density, unreliability, reliability, hazard function, and cumulative hazard function versus the random variable \eqn{x} (he upper and lower confidence limits are also included for the unreliability, reliability, and cumulative hazard function). Output also includes the nonparametric plot stated by \code{relfcn} input.
}
\references{
Modarres, M., Kaminskiy, M., Krivstov, V.  (2022)	Reliability Engineering and Risk Analysis: A Practical Guide (4th Edition).  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2020-2023)
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
nonparamplot.kaplanmeier(xi = Xi1, rc = Rc1, relfcn = "reliability", confbnd = 95, xlabel = "Time (hours)")

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
