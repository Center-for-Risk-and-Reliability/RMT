\name{lifestress.BAYESest}
\alias{lifestress.BAYESest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Life-Stress Estimator
}
\description{
Computes Bayesian posterior of Life-Stress model parameters.  The evaluation is performed using R's Stan interface library RStan.
}
\usage{
lifestress.BAYESest(pt_est, ls, dist, TTF, SF, Tc, Sc, confid, priors, nsamples, burnin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pt_est}{
The  initial estimate for the Bayesian parameter estimation process to begin with (the life distribution parameter leads followed by the life-stress parameters).  Make sure the order matches the model being entered.  That is, if the output paramters are the form c(\eqn{\beta, a, b}), then the input must also be some form of c(\eqn{\beta, a, b}).
}
  \item{ls}{
The life-stress model definition for analysis.  Choose from the list of life-models outlined in \code{\link{lifestress.select}}.
}
  \item{dist}{
  The distribution that represents time or life.
  \enumerate{
  \item \strong{Weibull} (Enter as "Weibull".  See more on \code{\link{probplot.wbl}}.)
  \item \strong{Three Parameter Weibull} (Enter as "3PWeibull".  See more on \code{\link{probplot.wbl3P}}.)
  \item \strong{Exponential} (Enter as "Exponential".  See more on \code{\link{probplot.exp}}.)
  \item \strong{Two Parameter Exponential} (Enter as "2PExponential".  See more on \code{\link{probplot.exp2P}}.)
  \item \strong{Normal} (Enter as "Normal".  See more on \code{\link{probplot.nor}}.)
  \item \strong{Lognormal} (Enter as "Lognormal".  See more on \code{\link{probplot.logn}}.)
  \item \strong{Gamma} (Enter as "Gamma".  See more on \code{\link{probplot.gam}}.)
  \item \strong{Logistic} (Enter as "Logistic".  See more on \code{\link{probplot.logist}}.)
  \item \strong{Log-logistic} (Enter as "Loglogistic".  See more on \code{\link{probplot.loglogist}}.)
  \item \strong{Gumbel} (Enter as "Gumbel".  See more on \code{\link{probplot.gumb}}.)
  }
}
  \item{TTF}{
A vector consisting of the failure or primary event data in a given set.  Output[[1]] of \code{sort.xircstressdata}.
}
  \item{SF}{
A vector consisting of the stress at the failure or primary event data \code{xi} in a given set.  In the event that there are multiple stress factors, the output will be a matrix.  Output[[3]] of \code{sort.xircstressdata}.
}
  \item{Tc}{
A vector consisting of the right censored data of a given set of data (NOTE: If no censored data exisits, this will return \code{NULL}).  Output[[2]] of \code{sort.xircstressdata}.
}
  \item{Sc}{
A vector consisting of the stress at the right censored data \code{rc} of a given set of data.  In the event that there are multiple stress factors, the output will be a matrix (NOTE: If no censored data exisits, this will return \code{NULL}).  Output[[4]] of \code{sort.xircstressdata}.
}
  \item{confid}{
The confidence in the maximum likelihood parameter estimates entered as a probability between \eqn{0} and \eqn{1}.  Default is \eqn{0.95} (or \eqn{95\%})
}
  \item{priors}{
A vector of prior distributions assigned to each of the parameters for the model (must be entered in the same order as \code{ptest}).  Because this procedure makes use of RStan, the priors shall adhere to the notation there-in.

\enumerate{
  \item \strong{Uniform (non-informative)} ~ \code{uniform(a low bound, b upper bound)}
  \item \strong{Normal} ~ {normal(mu mean, sigma standard deviation)}
  \item \strong{Lognormal} ~ \code{lognormal(mu logmean, sigma log standard deviation)}
  \item \strong{Exponential} ~ \code{exponential(lambda inverse scale)}
  \item \strong{Weibull} ~ \code{weibull(beta shape, alpha scale)}
  \item \strong{Gamma} ~ \code{gamma(alpha shape, beta scale)}
  \item \strong{Beta} ~ \code{beta(alpha shape, beta scale)}
}

Example: \code{priors<-c("normal(3,4)","normal(1,4)", "lognormal(-2,3)")}.  Any other distribution types that may not be on the list type in \code{lookup()} to see RStan's notation.
}
  \item{nsamples}{
The number of samples each Markov chain will produce.
}
  \item{burnin}{
The number of initial iterations to drop from the final Markov chain.
}
}
\details{
This script uses R's Stan library to formulate a stan file specific to the entered life-stress model and parameter priors therein.  The created Stan file runs a full Bayesian analysis based on four Markov chains per run (number of will be controlled in future versions of the library).
}
\value{
The output includes a summary of the completed Bayesian analysis in addition to the following list:
\item{fit}{
  The full summary of the Stan analysis as presented in the console and other output.  Output[[1]]
  }
\item{stats}{Description of 'comp2'Output[[2]]}
\item{
  dataout}{The mean and standard deviation of the posterior parameter markov chains as well as the log-likelihood.  Output[[3]]
  }
\item{quanout}{
  The 2.5, 5, 10, 25, 50, 75, 95, and 97.5\% quantile of the life-stress posterior parameters.  Output[[4]]
  }
\item{
  plot1_MCtrace}{Plot of the Markov chains for each life-stress posterior parameter.  Output[[5]]
  }
\item{
  plot2_hist}{Plot of the histograms for each life-stress posterior parameter.  Output[[6]]
  }
\item{
  plot3_density}{Plot of the density for each life-stress posterior parameter.  Output[[7]]
  }
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2021-2022)
}
\note{
As of version 1.2.6, this tool is operational.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example: A BAYESIAN ESTIMATION ALT DATA ANALYSIS EXAMPLE
## Takes prior data from a fatigue test (cycle to failure data vs. stress amplitude)
time_prior <- c(45000, 240000, 800000, 1500000, 2700000, 7800000, 10000000, 26000000, 12000000, 22000000)
cens_prior <- c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0)
stress_prior <- c(78.9, 74.02, 68.16, 63.27, 62.05, 59.61, 59.61, 58.63, 57.65, 57.41)
mat1 <- matrix(c(time_prior, cens_prior, stress_prior), nrow = length(time_prior), ncol = 3, byrow = FALSE,dimnames = list(NULL,c("cycles to failure","cens","stress amplitude")))
est1 <- lifestress.LSQest("InversePower","Lognormal",probplot.logn(mat1,"Blom","cycles to failure")[[1]])[[3]]
# Update data for Bayesian analysis
ttf1 <- c(1400000, 9000000)
ttrc1 <- c(10000000, 10000000, 10000000)
Sf1 <- c(58.7, 56.2)
Src1 <- c(57.2, 55.3, 45.1)
output <- lifestress.BAYESest(est1,ls="InversePower",dist="Lognormal",ttf1,Sf1,ttrc1,Src1,0.90,c("lognormal(-0.696369535,0.306574275)","lognormal(2.816114853,0.002729075)","lognormal(96.2360357,7.943598083)"),200000)



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
