\name{nonparamplot.bernbos}
\alias{nonparamplot.bernbos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bernard Bos-Levenbach Non-Parametric Plot
}
\description{
Nonparametric plot based on Benard and Bos-Levenbach Plottng Position of either failure probability or unreliability, reliability, hazard function, or cumulative hazard
}
\usage{
nonparamplot.bernbos(xi, rc, relfcn, confbnd, xlabel)
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
Name of the non-parametric curve to plot.  For (1) \emph{failure probability} or \emph{unreliability} \eqn{F(x)} enter \strong{"unreliability"}.  For (2) \emph{reliability} \eqn{R(x)} enter \strong{"reliability"}.  For (3) \emph{hazard function} \eqn{h(x)} enter \strong{"hazard"}.  Finally for (4) \emph{cumulative hazard function} \eqn{H(x)} enter \strong{"cumulativehazard"}.
}
  \item{confbnd}{
Confidence interval set to the nonparametric output (enter \eqn{95 \%} confidence interval as "95")
}
  \item{xlabel}{
Label for the x-axis for the plot.  If left blank, the default label is "X".
}
}
\details{
Definition of \emph{Benard and Bos-Levenbach} plotting position equation for failure probability (or unreliability) \eqn{F(x)}:
\deqn{F_i(x_i) = \frac{i - 0.3}{n + 0.2}}
}
\value{
Outputs the matrix of non-parametric mean point estimates and upper and lower confidence limits of unreliability, reliability, hazard function, and cumulative hazard function versus the random variable \eqn{x}.  Output also includes the nonparametric plot stated by \code{relfcn} input.
}
\references{
Benard, A. and Bos-Levenbach, E.C. (1953). The Plotting of Observations on Probability. Statististica Neerlandica. Vol. 7 : pp. 163--173.
}
\author{
Dr. Reuel Smith (2020-2022)
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
nonparamplot.bernbos(xi = Xi1, rc = Rc1, relfcn = "reliability", confbnd = 95, xlabel = "Time (hours)")

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
