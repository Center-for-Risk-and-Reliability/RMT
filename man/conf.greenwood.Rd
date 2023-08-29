\name{conf.greenwood}
\alias{conf.greenwood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Non-Parametric Greenwood Confidence
}
\description{
Computes the non-parametric confidence bounds of a non-parametric mean curve for unreliability, relibaility, hazard, and/or cumulative hazard based on Greenwood approximation for variance.
}
\usage{
conf.greenwood(alp, n, R, nandd)
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
The non-parametric estimate for reliability.  (NOTE: Pull from output column 3 of \code{plotposit.select(xi, rc, pp)}))
}
  \item{nandd}{
Matrix representing the \emph{remaining units} \eqn{n_i} at \eqn{x_i} and the \emph{number of units that failed} \eqn{d_i} at \eqn{x_i} (i.e. The output from \code{count.failcen(mat, n)} where \code{mat} is the output from \code{matrix.failcen(xi, rc, nx)}).
}
}
\details{
The Greenwood definition of the upper and lower limits of reliability \eqn{R(x)} and the cumulative hazard function \eqn{H(x)} are as follows:
\enumerate{
  \item \strong{Reliability Bounds}:
  \deqn{R(x) +/- t_(\alpha/2, n - 1) \sqrt [R(x)^2 \Sigma _(x_j <= x) d_j/[n_j (n_j - d_j)]]}

  \item \strong{Cumulative Hazard Bounds}:
  \deqn{- ln(R(x)) +/- t_(\alpha/2, n - 1) \sqrt [\Sigma _(x_j <= x) d_j/[n_j (n_j - d_j)]]}
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
tab1 <- plotposit.select(xi = Xi, rc = RC, pp="KaplanMeier")
N <- length(c(Xi,RC))
mat1 <- matrix.failcen(Xi, RC, N)
n_i.and.d_i <- count.failcen(mat1, N)

conf.greenwood(alp = 0.05, n = N, R = tab1[,3], nandd <- n_i.and.d_i)

## The function is currently defined as
function (alp, n, R, nandd)
{
    di <- nandd[, 2][which(nandd[, 2] > 0)]
    ni <- nandd[, 4][which(nandd[, 2] > 0)]
    Rlow <- R - qt(1 - alp/2, df = n - 1) * sqrt((R^2) * cumsum(di/(ni *
        (ni - di))))
    Rhi <- R + qt(1 - alp/2, df = n - 1) * sqrt((R^2) * cumsum(di/(ni *
        (ni - di))))
    Hlow <- -log(R) - qt(1 - alp/2, df = n - 1) * sqrt(cumsum(di/(ni *
        (ni - di))))
    Hhi <- -log(R) + qt(1 - alp/2, df = n - 1) * sqrt(cumsum(di/(ni *
        (ni - di))))
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
