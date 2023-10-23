\name{distribution.BAYESest}
\alias{distribution.BAYESest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Updater for Probability Distributions
}
\description{
Computes Bayesian posterior of basic distribution model parameters. The evaluation is performed using R's Stan interface library RStan.
}
\usage{
distribution.BAYESest(pt_est, dist, TTF, Tc = NULL, confid = 0.95, priors, nsamples = 20000, burnin, nchains = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pt_est}{
    The  initial estimate for the Bayesian parameter estimation process to begin with.  Make sure the order matches the model being entered.  That is, if the output paramters are the form c(\eqn{\alpha, \beta}), then the input must also be some form of c(\eqn{\alpha, \beta}).
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
  \item{Tc}{
A vector consisting of the right censored data of a given set of data (NOTE: If no censored data exisits, this will return \code{NULL}).  Output[[2]] of \code{sort.xircstressdata}.
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

Example: \code{priors<-c("normal(3,4)", "lognormal(-2,3)")}.  Any other distribution types that may not be on the list type in \code{lookup()} to see RStan's notation.
}
  \item{nsamples}{
The number of samples each Markov chain will produce.
}
  \item{burnin}{
The number of initial iterations to drop from the final Markov chain.
}
  \item{nchains}{
Number of Markov chains that the program runs.  (default is set to 4)
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
## Example: A BAYESIAN ESTIMATION DATA ANALYSIS EXAMPLE
## Take a set of test data where three failures occur at 135, 143, and 540 hours and three data are censored at 500, 600, and 600 hours.  Say that we have a prior distribution for Weibull parameters alpha ~ normal(106.11,6.61) and beta normal(5.16,1.311).

priorset <- c("normal(106.11,6.61)","normal(5.16,1.311)")
output <- distribution.BAYESest(pt_est = c(100,6),"Weibull",c(135,143,540),c(500,600,600),priors = priorset,nsamples = 100000,burnin = 10000)

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
