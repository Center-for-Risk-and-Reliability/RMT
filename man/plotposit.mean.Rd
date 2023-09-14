\name{plotposit.mean}
\alias{plotposit.mean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean Non-Parametric Output Tabulation
}
\description{
Tabulation of nonparametric computations based on Mean (also known as "Expected" or "Weibull") Plotting Position (failure probability or unreliability, reliability, hazard, cumulative hazard, and probability distribution).
}
\usage{
plotposit.mean(i, xi, rc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{
Rank of the failure data set \code{xi}.
}
  \item{xi}{
  A vector consisting of the failure or primary event data in a given set.  Data must consist of real and finite numbers.
  }
  \item{rc}{
  A vector consisting of the right censored data of a given set of data.  Data must consist of real and finite numbers.
  }
}
\details{
Definition of \emph{Mean} (or \emph{Expected}) plotting position equation for failure probability (or unreliability) \eqn{F(x)}:
\deqn{F_i(x_i) = \frac{i}{n + 1}}
For survival probability (or reliability) \eqn{R(x)}:
\deqn{R_i(x_i) = \frac{n - i + 1}{n + 1}}
For cumulative hazard \eqn{H(x)}:
\deqn{H_i(x_i) = -\ln\left(\frac{n - i + 1}{n + 1}\right)}
\deqn{\text{or }H_i(x_i) = \ln\left(n + 1\right) - \ln\left(n - i + 1\right)}
For hazard rate \eqn{h(x)}:
\deqn{h_i(x_i) = \frac{1}{(n - i + 1)(x_{i+1}-x_i)}}
For probability density \eqn{f(x)}:
\deqn{f_i(x_i) = \frac{1}{(n + 1)(x_{i+1}-x_i)}}
}
\value{
Returns a matrix of nonparameter estimates for each failure data xi: Column 2 is the \emph{failure probability}, Column 3 is the \emph{reliability}, Column 4 is the \emph{hazard rate}, Column 5 is the \emph{cumulative hazard}, and Column 6 is the \emph{probability density}.
}
\references{
Weibull, W. (1939) A statistical Theory of Strength of Materials. Ingeniors Vetenskaps Academy Handlingar.  (Stockholm), 151, pp 1–45.

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
plotposit.mean(i = rankcalc(xi = Xi1, rc = Rc1), xi = Xi1, rc = Rc1)
plotposit.mean(rankcalc(Xi1, Rc1), Xi1, Rc1)

## Example 2
Xi1 <- c(70,71,75,78,89)
plotposit.mean(i = rankcalc(xi = Xi1), xi = Xi1)
plotposit.mean(rankcalc(Xi1, NULL), Xi1, NULL)
plotposit.mean(rankcalc(Xi1), Xi1)
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
