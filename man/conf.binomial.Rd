\name{conf.binomial}
\alias{conf.binomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Non-Parametric Binomial Confidence
}
\description{
Computes the non-parametric confidence bounds of a non-parametric mean curve for unreliability, relibaility, hazard, and/or cumulative hazard based on Binomial approximation for variance.  This limit is primarily used for all non-parametric estimates \emph{except} by Kaplan-Meier.
}
\usage{
conf.binomial(alp, n, R)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alp}{
The \eqn{\alpha} value for the confidence interval \eqn{100(1 - \alpha)}.
}
  \item{n}{
Sample size of data.
}
  \item{R}{
The non-parametric estimate for reliability.
}
}
\details{
The binomial definition of the upper and lower limits of reliability and the cumulative hazard function are as follows:
\enumerate{
  \item \strong{Reliability Bounds}:
  \deqn{R(x) +/- t_(1 - \alpha/2, n - 1) \sqrt [R(1 - R)]}

  \item \strong{Cumulative Hazard Bounds}:
  \deqn{- ln(R(x)) +/- t_(1 - \alpha/2, n - 1) \sqrt [(1 - R)/R]}
}
}
\value{
A list of y-curves output consisting of the following:

\item{Rlow}{
  Lower percentile of non-parametric reliability \code{[[1]]}.
}
\item{Rhi}{
  Upper percentile of non-parametric reliability \code{[[2]]}.
}
\item{Flow}{
  Lower percentile of non-parametric failure probability or unreliability \code{[[3]]}.
}
\item{Fhi}{
  Upper percentile of non-parametric failure probability or unreliability \code{[[4]]}.
}
\item{Hlow}{
  Lower percentile of non-parametric cumulative hazard \code{[[5]]}.
}
\item{Hhi}{
  Upper percentile of non-parametric cumulative hazard \code{[[6]]}.
}
\item{hlow}{
  Lower percentile of non-parametric hazard function \code{[[7]]}.
}
\item{hhi}{
  Upper percentile of non-parametric hazard function \code{[[8]]}.
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
## Example:
Xi <- c(65,75,75.2,87.5,88.3,94.2,101.7,109.2,130)
RC <- c(31.7,39.2,57.2,65.8,70,105.8,110)
tab1 <- plotposit.select(xi = Xi, rc = RC, pp="Blom")
N <- length(Xi)+length(RC)
conf.binomial(alp = 0.05, n = N, R = tab1[,3])

## The function is currently defined as
function (alp, n, R)
{
    Rlow <- R - qt(1 - alp/2, df = n - 1) * sqrt(R * (1 - R))
    Rhi <- R + qt(1 - alp/2, df = n - 1) * sqrt(R * (1 - R))
    Hlow <- -log(R) - qt(1 - alp/2, df = n - 1) * sqrt((1 - R)/R)
    Hhi <- -log(R) + qt(1 - alp/2, df = n - 1) * sqrt((1 - R)/R)
    Rlow[which(Rlow < 0)] <- 0
    Rlow[which(is.nan(Rlow))] <- 0
    Rhi[which(Rhi > 1)] <- 1
    Rhi[which(is.nan(Rhi))] <- 0
    Flow <- 1 - Rhi
    Fhi <- 1 - Rlow
    Hlow[which(Hlow < 0)] <- 0
    hlow <- c(Hlow[1], diff(Hlow))
    hlow[which(hlow < 0)] <- 0
    hhi <- c(Hhi[1], diff(Hhi))
    return(list(Rlow, Rhi, Flow, Fhi, Hlow, Hhi, hlow, hhi))
  }
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
