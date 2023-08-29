\name{plot.nonparam}
\alias{plot.nonparam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Non-Parametric Plot Fitter
}
\description{
Fits the y-output for nonparametric plotting to a square-like curve.  Used for nonparametric failure probability and Reliablity curves.
}
\usage{
plot.nonparam(xi, yi)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xi}{
vector of x-values (failure data).
}
  \item{yi}{
Vector of y-values (failure probabilty or reliability data).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{xout}{
  Post processed x-values (rectangularized) for nonparametric plotting.  Output[[1]]
}
\item{yout}{
  Post processed y-values (rectangularized) for nonparametric plotting.  Output[[2]]
}
}
\references{
%% ~put references to the literature/web site here ~
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
## Example 1: Failure probability
tab1 <- plotposit.select(xi = c(70,71,75,78,89), rc = c(80,80,84), pp="Blom")
plot.nonparam(tab1[,1],tab1[,2])

## Example 2: Reliability
tab1 <- plotposit.select(xi = c(70,71,75,78,89), rc = c(80,80,84), pp="Blom")
plot.nonparam(tab1[,1],tab1[,3])

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
