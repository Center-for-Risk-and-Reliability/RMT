\name{plotposit.grigorten}
\alias{plotposit.grigorten}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Grigorten Non-Parametric Output Tabulation
}
\description{
Tabulation of nonparametric computations based on Grigorten Plotting Position (failure probability or unreliability, reliability, hazard, cumulative hazard).
}
\usage{
plotposit.grigorten(i, xi, rc)
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
Definition of \emph{Grigorten} plotting position equation for failure probability (or unreliability) \eqn{F(x)}:
\deqn{F_i(x_i) = \frac{i - 0.44}{n + 0.12}}
}
\value{
Returns a matrix of nonparameter estimates for each failure data xi: Column 2 is the \emph{failure probability}, Column 3 is the \emph{reliability}, Column 4 is the \emph{hazard rate}, and Column 5 is the \emph{cumulative hazard}.
}
\references{
Gringorten, I. I. (1963) A Plotting Rule for Extreme Probability. Journal of Geophysical Research. vol. 68, pp. 813–-814
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
plotposit.grigorten(i = rankcalc(xi = Xi1, rc = Rc1), xi = Xi1, rc = Rc1)
plotposit.grigorten(rankcalc(Xi1, Rc1), Xi1, Rc1)

## Example 2
Xi1 <- c(70,71,75,78,89)
plotposit.grigorten(i = rankcalc(xi = Xi1), xi = Xi1)
plotposit.grigorten(rankcalc(Xi1, NULL), Xi1, NULL)
plotposit.grigorten(rankcalc(Xi1), Xi1)
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
