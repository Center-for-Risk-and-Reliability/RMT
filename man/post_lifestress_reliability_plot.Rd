\name{post_lifestress_reliability_plot}
\alias{post_lifestress_reliability_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior Life-Stress Reliability Plot
}
\description{
Post-processes the reliability plot for a given stress-life model and output parameters (Bayesian parameter estimate).
}
\usage{
post_lifestress_reliability_plot(post_params, ls, dist, confid, Xmin, Xmax, Xlab, S0, t0, R0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{post_params}{
  The mean and standard deviation of the posterior parameter markov chains as well as the log-likelihood (Output[[3]] from \code{lifestress.BAYESest}).
}
  \item{ls}{
The life-stress model definition for analysis.  Choose from the list of life-models outlined in \code{\link{lifestress.select}}.
}
  \item{dist}{
%%     ~~Describe \code{dist} here~~
}
  \item{confid}{
The confidence in the Bayesian parameter estimates entered as a probability between \eqn{0} and \eqn{1}.  Default is \eqn{0.95} (or \eqn{95\%})
}
  \item{Xmin}{
%%     ~~Describe \code{Xmin} here~~
}
  \item{Xmax}{
%%     ~~Describe \code{Xmax} here~~
}
  \item{Xlab}{
Written label used for the charateristic stress axis or axes.  When the Arrhenius or Eyring life-stress model is selected, this label defaults to \eqn{K^-1} where \eqn{K} is the temperature notation for Kelvin.  Dual-stresses must be entered as vectors (Example: For "Temperature" and "Humidity" enter \code{c("Temperature", "Humidity")}).
}
  \item{S0}{
The use-level stress to be used for the relationship plot (if dual stress enter as a vector \code{c(Stress1,Stress2)}).  This is an optional input but if left out, then the output values will vary (see below).
}
  \item{t0}{
%%     ~~Describe \code{t0} here~~
}
  \item{R0}{
%%     ~~Describe \code{R0} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (post_params, ls, dist, confid, Xmin, Xmax, Xlab, S0,
    t0, R0)
{
    if (missing(confid)) {
        conf.level <- 0.95
    }
    else {
        conf.level <- confid
    }
    if (missing(Xlab)) {
        Xlab <- "X"
    }
    else {
        Xlab <- Xlab
    }
    N <- length(post_params[[1]])
    M <- 100
    X_v <- linspace(Xmin, Xmax, n = M)
    X_m <- matrix(rep(X_v, N), nrow = N, ncol = M, byrow = TRUE)
    if (ls == "Linear") {
        a_m <- matrix(rep(post_params$a, M), nrow = N, ncol = M,
            byrow = FALSE)
        b_m <- matrix(rep(post_params$b, M), nrow = N, ncol = M,
            byrow = FALSE)
        life_m <- b_m + S0 * a_m
        life_v <- post_params$b + S0 * post_params$a
    }
    if (ls == "Exponential") {
        a_m <- matrix(rep(post_params$a, M), nrow = N, ncol = M,
            byrow = FALSE)
        b_m <- matrix(rep(post_params$b, M), nrow = N, ncol = M,
            byrow = FALSE)
        life_m <- b_m * exp(a_m * S0)
        life_v <- post_params$b * exp(post_params$a * S0)
    }
    if (ls == "Arrhenius") {
        Ea_m <- matrix(rep(post_params$Ea, M), nrow = N, ncol = M,
            byrow = FALSE)
        b_m <- matrix(rep(post_params$b, M), nrow = N, ncol = M,
            byrow = FALSE)
        life_m <- b_m * exp(Ea_m/((8.617385e-05) * S0))
        life_v <- post_params$b * exp(post_params$Ea/((8.617385e-05) *
            S0))
    }
    if (ls == "Eyring") {
        a_m <- matrix(rep(post_params$a, M), nrow = N, ncol = M,
            byrow = FALSE)
        b_m <- matrix(rep(post_params$b, M), nrow = N, ncol = M,
            byrow = FALSE)
        life_m <- (b_m/S0) * exp(a_m/S0)
        life_v <- (post_params$b/S0) * exp(post_params$a/S0)
    }
    if (ls == "Eyring2") {
        a_m <- matrix(rep(post_params$a, M), nrow = N, ncol = M,
            byrow = FALSE)
        b_m <- matrix(rep(post_params$b, M), nrow = N, ncol = M,
            byrow = FALSE)
        life_m <- (1/S0) * exp(-(a_m - (b_m/S0)))
        life_v <- (1/S0) * exp(-(post_params$a - (post_params$b/S0)))
    }
    if (ls == "Power") {
        a_m <- matrix(rep(post_params$a, M), nrow = N, ncol = M,
            byrow = FALSE)
        b_m <- matrix(rep(post_params$b, M), nrow = N, ncol = M,
            byrow = FALSE)
        life_m <- b_m * (S0^a_m)
        life_v <- post_params$b * (S0^post_params$a)
    }
    if (ls == "InversePower") {
        a_m <- matrix(rep(post_params$a, M), nrow = N, ncol = M,
            byrow = FALSE)
        b_m <- matrix(rep(post_params$b, M), nrow = N, ncol = M,
            byrow = FALSE)
        life_m <- b_m * (S0^-a_m)
        life_v <- post_params$b * (S0^-post_params$a)
    }
    if (ls == "Logarithmic") {
        a_m <- matrix(rep(post_params$a, M), nrow = N, ncol = M,
            byrow = FALSE)
        b_m <- matrix(rep(post_params$b, M), nrow = N, ncol = M,
            byrow = FALSE)
        life_m <- a_m * log(S0) + b_m
        life_v <- post_params$a * log(S0) + post_params$b
    }
    if (ls == "MultiStress") {
        lsparams <- "vector a; real b;"
        lsparamsvec <- c("a", "b")
        pr1 <- paste(c("a ~ ", priors[ishift + 1], ";"), collapse = "")
        pr2 <- paste(c("b ~ ", priors[ishift + 2], ";"), collapse = "")
        lspriors <- paste(c(pr1, pr2), collapse = " ")
        lifeF <- exp(theta[ishift + 1:length(SF) + ishift + 1] \%*\%
            c(1, SF))
        function(theta) {
            exp(theta[ishift + 1:length(SF) + ishift + 1] \%*\%
                c(1, SF))
        }
        loglifeF <- function(theta) {
            theta[ishift + 1:length(SF) + ishift + 1] \%*\% c(1,
                SF)
        }
        if (missing(Tc) == FALSE) {
            lifeC <- function(theta) {
                exp(theta[ishift + 1:length(Sc) + ishift + 1] \%*\%
                  c(1, Sc))
            }
            loglifeC <- function(theta) {
                theta[ishift + 1:length(Sc) + ishift + 1] \%*\%
                  c(1, Sc)
            }
        }
    }
    if (ls == "TempHumidity") {
        A_m <- matrix(rep(post_params$A, M), nrow = N, ncol = M,
            byrow = FALSE)
        a_m <- matrix(rep(post_params$a, M), nrow = N, ncol = M,
            byrow = FALSE)
        b_m <- matrix(rep(post_params$b, M), nrow = N, ncol = M,
            byrow = FALSE)
        life_m <- A_m * exp((a_m/S0[, 1]) + (b_m/S0[, 2]))
    }
    if (ls == "TempNonthermal") {
        lsparams <- "real A; real a; real b;"
        lsparamsvec <- c("A", "a", "b")
        pr1 <- paste(c("a ~ ", priors[ishift + 1], ";"), collapse = "")
        pr2 <- paste(c("b ~ ", priors[ishift + 2], ";"), collapse = "")
        pr3 <- paste(c("c ~ ", priors[ishift + 3], ";"), collapse = "")
        lspriors <- paste(c(pr1, pr2, pr3), collapse = " ")
        lifeF <- "c/((Sf[,2]^b)*exp(-a/Sf[,1]))"
        loglifeF <- "log(c) - b*log(Sf[,2]) + (a/Sf[,1])"
        if (missing(Tc) == FALSE) {
            lifeC <- "c/((Sc[,2]^b)*exp(-a/Sc[,1]))"
            loglifeC <- "log(c) - b*log(Sc[,2]) + (a/Sc[,1])"
        }
        if (missing(S0) == FALSE) {
            lifepost <- "b + S0ref*a"
            loglifepost <- "log(b + S0ref*a)"
        }
    }
    if (ls == "Eyring3") {
        lsparams <- "real a; real b; real c; real d"
        lsparamsvec <- c("a", "b", "c", "d")
        pr1 <- paste(c("a ~ ", priors[ishift + 1], ";"), collapse = "")
        pr2 <- paste(c("b ~ ", priors[ishift + 2], ";"), collapse = "")
        pr3 <- paste(c("c ~ ", priors[ishift + 3], ";"), collapse = "")
        pr4 <- paste(c("d ~ ", priors[ishift + 4], ";"), collapse = "")
        lspriors <- paste(c(pr1, pr2, pr3, pr4), collapse = " ")
        lifeF <- "(1/Sf[,1])*exp((a + (b/Sf[,1])) + (c + (d/Sf[,1]))*Sf[,2])"
        loglifeF <- "-log(Sf[,1]) + a + (b/Sf[,1]) + (c + (d/Sf[,1]))*Sf[,2]"
        if (missing(Tc) == FALSE) {
            lifeC <- "(1/Sc[,1])*exp((a + (b/Sc[,1])) + (c + (d/Sc[,1]))*Sc[,2])"
            loglifeC <- "-log(Sc[,1]) + a + (b/Sc[,1]) + (c + (d/Sc[,1]))*Sc[,2]"
        }
        if (missing(S0) == FALSE) {
            lifepost <- "(1/S0ref[,1])*exp((a + (b/S0ref[,1])) + (c + (d/S0ref[,1]))*S0ref[,2])"
            loglifepost <- "-log(S0ref[,1]) + a + (b/S0ref[,1]) + (c + (d/S0ref[,1]))*S0ref[,2]"
        }
    }
    if (dist == "Weibull") {
        beta_m <- matrix(rep(post_params$beta, M), nrow = N,
            ncol = M, byrow = FALSE)
        R_m <- 1 - pweibull(X_m, life_m, beta_m)
        if (missing(t0) == FALSE) {
            R_t0 <- 1 - pweibull(t0, life_v, post_params$beta)
        }
        if (missing(R0) == FALSE) {
            t_R0 <- life_v * (-log(R0))^(1/post_params$beta)
        }
    }
    if (dist == "Lognormal") {
        sigma_t_m <- matrix(rep(post_params$sigma_t, M), nrow = N,
            ncol = M, byrow = FALSE)
        R_m <- 1 - plnorm(X_m, log(life_m), sigma_t_m)
        if (missing(t0) == FALSE) {
            R_t0 <- 1 - plnorm(t0, log(life_v), post_params$sigma_t)
        }
        if (missing(R0) == FALSE) {
            t_R0 <- exp(log(life_v) + post_params$sigma_t * qnorm(1 -
                R0))
        }
    }
    if (dist == "Normal") {
        sigma_m <- matrix(rep(post_params$sigma, M), nrow = N,
            ncol = M, byrow = FALSE)
        R_m <- 1 - pnorm(X_m, life_m, sigma_m)
        if (missing(t0) == FALSE) {
            R_t0 <- 1 - pnorm(t0, life_v, post_params$sigma)
        }
        if (missing(R0) == FALSE) {
            t_R0 <- life_v + post_params$sigma * qnorm(1 - R0)
        }
    }
    if (dist == "Exponential") {
        R_m <- 1 - pexp(X_m, 1/life_m)
        if (missing(t0) == FALSE) {
            R_t0 <- 1 - pexp(t0, 1/life_v)
        }
        if (missing(R0) == FALSE) {
            t_R0 <- -life_v * log(R0)
        }
    }
    if (dist == "2PExponential") {
        sigma_m <- matrix(rep(post_params$sigma, M), nrow = N,
            ncol = M, byrow = FALSE)
        R_m <- 1 - pexp(X_m - life_m, 1/sigma_m)
        if (missing(t0) == FALSE) {
            R_t0 <- 1 - pexp(t0 - life_v, 1/post_params$sigma)
        }
        if (missing(R0) == FALSE) {
            t_R0 <- life_v - post_params$sigma * log(R0)
        }
    }
    R_mean <- colMeans(R_m)
    R_low <- colQuantiles(R_m, (1 - confid)/2)
    R_high <- colQuantiles(R_m, 1 - (1 - confid)/2)
    if (missing(t0) == TRUE && missing(R0) == TRUE) {
        df <- data.frame(X = X_v, Rmean = R_mean, Rmin = R_low,
            Rmax = R_high)
    }
    if (missing(t0) == FALSE && missing(R0) == TRUE) {
        df <- data.frame(X = X_v, Rmean = R_mean, Rmin = R_low,
            Rmax = R_high, reliability_at_t0 = R_t0)
    }
    if (missing(t0) == TRUE && missing(R0) == FALSE) {
        df <- data.frame(X = X_v, Rmean = R_mean, Rmin = R_low,
            Rmax = R_high, life_at_R0 = t_R0)
    }
    if (missing(t0) == FALSE && missing(R0) == FALSE) {
        df <- data.frame(X = X_v, Rmean = R_mean, Rmin = R_low,
            Rmax = R_high, life_at_R0 = t_R0, reliability_at_t0 = R_t0)
    }
    plotout <- ggplot(data = df, aes(X, Rmean)) + geom_line(color = "red",
        size = 0.9) + geom_ribbon(aes(ymin = Rmin, ymax = Rmax,
        x = X), alpha = 0.5, fill = "red") + xlab(Xlab) + ylab("Reliability")
    if (missing(t0) == FALSE) {
        stats_R_t0 <- c(mean(R_t0), std(R_t0), quantile(R_t0,
            c((1 - confid)/2, 0.5, 1 - (1 - confid)/2)))
        plotout2 <- ggplot(data = df, aes(reliability_at_t0,
            fill = "blue", color = "blue")) + geom_density(alpha = 0.1)
    }
    if (missing(R0) == FALSE) {
        stats_t_R0 <- c(mean(t_R0), std(t_R0), quantile(t_R0,
            c((1 - confid)/2, 0.5, 1 - (1 - confid)/2)))
        plotout3 <- ggplot(data = df, aes(life_at_R0, fill = "blue",
            color = "blue")) + geom_density(alpha = 0.1)
    }
    if (missing(t0) == TRUE && missing(R0) == TRUE) {
        return(list(plotout, X_v, R_mean, R_low, R_high))
    }
    if (missing(t0) == FALSE && missing(R0) == TRUE) {
        return(list(plotout, plotout2, X_v, R_mean, R_low, R_high,
            stats_R_t0))
    }
    if (missing(t0) == TRUE && missing(R0) == FALSE) {
        return(list(plotout, plotout3, X_v, R_mean, R_low, R_high,
            stats_t_R0))
    }
    if (missing(t0) == FALSE && missing(R0) == FALSE) {
        return(list(plotout, plotout2, plotout3, X_v, R_mean,
            R_low, R_high, stats_R_t0, stats_t_R0))
    }
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
