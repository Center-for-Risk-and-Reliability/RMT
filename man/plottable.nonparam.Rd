\name{plottable.nonparam}
\alias{plottable.nonparam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Non-Parametric Output Tabulation
}
\description{
Tabulation and plotting of non-parametric computations (failure probability or unreliability, reliability, hazard function, and cumulative hazard function).
}
\usage{
plottable.nonparam(xi, rc, FRhH, relfcn, alpha, xlabel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xi}{
  A vector consisting of the failure or primary event data in a given set.  Data must consist of real and finite numbers.
  }
  \item{rc}{
  A vector consisting of the right censored data of a given set of data.  Data must consist of real and finite numbers.
  }
  \item{FRhH}{
The output from tabulation of nonparametric computations (CDF, reliability, hazard, cumulative hazard).  (i.e. The output from \code{plotposit.select(xi, rc, pp)})
}
  \item{relfcn}{
Name of the non-parametric curve to plot.  For (1) \emph{failure probability} or \emph{unreliability} \eqn{F(x)} enter \strong{"unreliability"}.  For (2) \emph{reliability} \eqn{R(x)} enter \strong{"reliability"}.  For (3) \emph{hazard function} \eqn{h(x)} enter \strong{"hazard"}.  Finally for (4) \emph{cumulative hazard function} \eqn{H(x)} enter \strong{"cumulativehazard"}.
}
  \item{alpha}{
The \eqn{\alpha} value for the confidence interval \eqn{100(1 - \alpha)}.
}
  \item{xlabel}{
Label for the x-axis for the plot.  If left blank, the default label is "X".
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Outputs the matrix of non-parametric mean point estimates and upper and lower confidence limits of unreliability, reliability, hazard function, and cumulative hazard function versus the random variable \eqn{x}.
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
## Example 1a: Plot reliability mean and confidence intervals
Xi <- c(65,75,75.2,87.5,88.3,94.2,101.7,109.2,130)
RC <- c(31.7,39.2,57.2,65.8,70,105.8,110)
nonparamtable <- plotposit.select(xi = Xi, rc = RC, pp = "Blom")
plottable.nonparam(xi = Xi, rc = RC, FRhH = nonparamtable, relfcn = "reliability", alpha = 0.05, xlabel = "time (hours")
##         1b: Plot cumulative hazard funtion mean and confidence intervals
plottable.nonparam(xi = Xi, rc = RC, FRhH = nonparamtable, relfcn = "cumulativehazard", alpha = 0.05, xlabel = "time (hours")
##         1c: Plot failure probability mean and confidence
plottable.nonparam(xi = Xi, rc = RC, FRhH = nonparamtable, relfcn = "unreliability", alpha = 0.05, xlabel = "time (hours")
##         1d: Plot hazard function mean and confidence
plottable.nonparam(xi = Xi, rc = RC, FRhH = nonparamtable, relfcn = "hazard", alpha = 0.05, xlabel = "time (hours")

## Example 1a: Plot reliability mean and confidence intervals
Xi <- c(65,75,75.2,87.5,88.3,94.2,101.7,109.2,130)
RC <- c(31.7,39.2,57.2,65.8,70,105.8,110)
nonparamtable2 <- plotposit.select(xi = Xi, rc = RC, pp = "KaplanMeier")
plottable.nonparam(xi = Xi, rc = RC, FRhH = nonparamtable2, relfcn = "reliability", alpha = 0.05, xlabel = "time (hours")
##         1b: Plot cumulative hazard funtion mean and confidence intervals
plottable.nonparam(xi = Xi, rc = RC, FRhH = nonparamtable2, relfcn = "cumulativehazard", alpha = 0.05, xlabel = "time (hours")
##         1c: Plot failure probability mean and confidence
plottable.nonparam(xi = Xi, rc = RC, FRhH = nonparamtable2, relfcn = "unreliability", alpha = 0.05, xlabel = "time (hours")
##         1d: Plot hazard function mean and confidence
plottable.nonparam(xi = Xi, rc = RC, FRhH = nonparamtable2, relfcn = "hazard", alpha = 0.05, xlabel = "time (hours")

## The function is currently defined as
function (xi, rc, FRhH, relfcn, alpha, xlabel)
{
    if (missing(rc)) {
        nxi <- length(xi)
        rawdat <- matrix.failcen(xi, nx = nxi)
    }
    else {
        nxi <- length(c(xi, rc))
        rawdat <- matrix.failcen(xi, rc, nx = nxi)
    }
    if (missing(xlabel)){
        xlabel<-"X"
    }
    bounds <- conf.greenwood(alpha, nxi, FRhH[, 3], count.failcen(rawdat,
        nxi))
    if (relfcn == "unreliability") {
        medplot <- plot.nonparam(FRhH[, 1], FRhH[, 2])
        lowplot <- plot.nonparam(FRhH[, 1], bounds[[3]])
        hiplot <- plot.nonparam(FRhH[, 1], bounds[[4]])
        plot(medplot[[1]], medplot[[2]], type = "l", xlab = xlabel,
            ylab = "Unreliability", ylim = c(0, 1), col = "blue")
        lines(medplot[[1]], lowplot[[2]], type = "l", lty = 2,
            col = "blue")
        lines(medplot[[1]], hiplot[[2]], type = "l", lty = 2,
            col = "blue")
    }
    if (relfcn == "reliability") {
        medplot <- plot.nonparam(FRhH[, 1], FRhH[, 3])
        lowplot <- plot.nonparam(FRhH[, 1], bounds[[1]])
        hiplot <- plot.nonparam(FRhH[, 1], bounds[[2]])
        plot(medplot[[1]], medplot[[2]], type = "l", xlab = xlabel,
            ylab = "Reliability", ylim = c(0, 1), col = "blue")
        lines(medplot[[1]], lowplot[[2]], type = "l", lty = 2,
            col = "blue")
        lines(medplot[[1]], hiplot[[2]], type = "l", lty = 2,
            col = "blue")
    }
    if (relfcn == "hazard") {
        medplot <- plot.nonparamhaz(FRhH[, 1], FRhH[, 4])
        lowplot <- plot.nonparamhaz(FRhH[, 1], bounds[[7]])
        hiplot <- plot.nonparamhaz(FRhH[, 1], bounds[[8]])
        if (is.finite(max(hiplot[[2]]))) {
            ylims <- c(0, max(hiplot[[2]]))
        }
        else {
            ylims <- c(0, hiplot[[2]][length(hiplot[[2]]) - 2])
        }
        plot(medplot[[1]], medplot[[2]], type = "l", xlab = xlabel,
            ylab = "Hazard", ylim = ylims, col = "blue")
        lines(medplot[[1]], lowplot[[2]], type = "l", lty = 2,
            col = "blue")
        lines(medplot[[1]], hiplot[[2]], type = "l", lty = 2,
            col = "blue")
    }
    if (relfcn == "cumulativehazard") {
        medplot <- plot.nonparamhaz(FRhH[, 1], FRhH[, 5])
        lowplot <- plot.nonparamhaz(FRhH[, 1], bounds[[5]])
        hiplot <- plot.nonparamhaz(FRhH[, 1], bounds[[6]])
        if (is.finite(max(hiplot[[2]]))) {
            ylims <- c(0, max(hiplot[[2]]))
        }
        else {
            ylims <- c(0, hiplot[[2]][length(hiplot[[2]]) - 2])
        }
        plot(medplot[[1]], medplot[[2]], type = "l", xlab = xlabel,
            ylab = "Cumulative Hazard", ylim = ylims, col = "blue")
        lines(medplot[[1]], lowplot[[2]], type = "l", lty = 2,
            col = "blue")
        lines(medplot[[1]], hiplot[[2]], type = "l", lty = 2,
            col = "blue")
    }
    nonparamset <- matrix(c(FRhH[, 1], FRhH[, 2], bounds[[3]],
        bounds[[4]], FRhH[, 3], bounds[[1]], bounds[[2]], FRhH[,
            4], bounds[[7]], bounds[[8]], FRhH[, 5], bounds[[5]],
        bounds[[6]]), nrow = length(xi), ncol = 13)
    colnames(nonparamset) <- c(xlabel, "Unreliability", "Low bound",
        "High bound", "Reliability", "Low bound", "High bound",
        "Hazard", "Low bound", "High bound", "Cumulative Hazard",
        "Low bound", "High bound")
    return(nonparamset)
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
