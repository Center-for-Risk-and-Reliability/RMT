\name{fit.gumb}
\alias{fit.gumb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gumbel Fit
}
\description{
Computes the maximum likelihood estimates for the point estimates for the Gumbel parameters.
}
\usage{
fit.gumb(xi, rc = NULL, pp = "Blom", confid = 0.95, sided = "twosided")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xi}{
    A vector consisting of the failure or primary event data in a given set.  Data must consist of real and finite numbers.
  }
  \item{rc}{
    A vector consisting of the right censored data of a given set of data.  Data must consist of real and finite numbers.
  }
  \item{pp}{
    Name of the plotting position to select for computation.  See \code{\link{plotposit.select}} for plotting position definitions.  The default position is \code{"Blom"}.
  }
  \item{confid}{
    The confidence in the maximum likelihood parameter estimates entered as a probability between \eqn{0} and \eqn{1}.  Default is \eqn{0.95} (or \eqn{95\%})
  }
  \item{sided}{
    Indicates whether the limit calculations are based on a two-sided confidence or one-sided confidence.  The three options for this heading are \code{"twosided"} for two-sided confidence (default), \code{"onesidedlow"} for one-sided (low) confidence, and \code{"onesidedhigh"} for one-sided (high) confidence.
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
Returns a table made up of the MLE point estimates and confidence limits for the Gumbel location parameter \eqn{\mu} and the Gumbel scale parameter \eqn{\sigma}.  Additionally generates a list made up of the following:
\item{gumbparams}{
  Maximum likelihood mean estimate of the Gumbel location parameter \eqn{\mu} and the Gumbel scale parameter \eqn{\sigma}.
}
\item{CImu}{
  The confidence interval for the Gumbel location parameter \eqn{\mu}.
}
\item{CIsigma}{
  The confidence interval for the Gumbel scale parameter \eqn{\sigma}.
}
\item{varcovmat}{
  The variance-covariance matrix as computed by non-linear minimization.
}
\item{loglikelihood}{
  The log-likelihood of the point estimate.
}
\item{AIC}{
  The Akaikie Information Criteria (AIC) score for the point estimate.
}
\item{BIC}{
  The Bayesian Information Criteria (BIC) score for the point estimate.
}
}
\references{
Modarres, M., Kaminskiy, M., Krivstov, V.  (2022)	Reliability Engineering and Risk Analysis: A Practical Guide (4th Edition).  Center for Risk and Reliability
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
Xi1 <- c(70,71,75,78,89)
Rc1 <- c(80,80,84)
fit.gumb(Xi1,Rc1)

## Example 2
Xi1 <- c(70,71,75,78,89)
Rc1 <- c(80,80,84)
fit.gumb(Xi1)
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
