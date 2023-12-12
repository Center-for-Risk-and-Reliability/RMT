\name{stepstress.BAYESest}
\alias{stepstress.BAYESest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Step-Stress Estimator
}
\description{
Computes Bayesian posterior of life step-stress model parameters.
}
\usage{
stepstress.BAYESest(pt_est, data, stepstresstable, ls, dist, confid, priors, nsamples, burnin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pt_est}{
The  initial estimate for the Bayesian parameter estimation process to begin with (the life distribution parameter leads followed by the life step-stress parameters).  Make sure the order matches the model being entered.  That is, if the output paramters are the form c(\eqn{\beta, a, b}), then the input must also be some form of c(\eqn{\beta, a, b}).
}
  \item{data}{
A matrix that represents a data set series which may consist of the following:
\itemize{
  \item{
  Column 1 consists of the failure (or primary event) data and the right censored data in a given set.  Data must consist of real and finite numbers.
  }
  \item{
  Column 2 gives the censored status of the data in column 1.  Demark failed data as \eqn{1} and right-censored data as \eqn{0}.
  }
  \item{
  Column(s) 3 (and on for multi-type stress) states the stress (or step-stress) values for each data in column 1.  List one stress type per column (i.e. Column 3 - Temperature and Column 4 - Humidity).
  }
}
}
  \item{stepstresstable}{
An \eqn{n x m} matrix that represents the stress and test time at a set of \eqn{n} step levels where \eqn{n} is the number of rows.
\itemize{
  \item{
  Column(s) 1 (through \eqn{m-1} for multi-type stress) states the  step-stress values for each possible step.  List one stress type per column (i.e. Column 1 - Temperature and Column 2 - Humidity).
  }
  \item{
  Column \eqn{m} is the overall test time for each step level.  Data must consist of real and finite numbers.
  }
}
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
This script uses R's Stan library to formulate a stan file specific to the entered life step-stress model and parameter priors therein.  The created Stan file runs a full Bayesian analysis based on four Markov chains per run (number of will be controlled in future versions of the library).
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
  The 2.5, 5, 10, 25, 50, 75, 95, and 97.5\% quantile of the life step-stress posterior parameters.  Output[[4]]
  }
\item{
  plot1_MCtrace}{Plot of the Markov chains for each life step-stress posterior parameter.  Output[[5]]
  }
\item{
  plot2_hist}{Plot of the histograms for each life step-stress posterior parameter.  Output[[6]]
  }
\item{
  plot3_density}{Plot of the density for each life step-stress posterior parameter.  Output[[7]]
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

\section{Warning }{
%% post warnings
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1: Single Stress Step-Stress
## Load the CSVs for step-stress data and the step-stress table from the Github site by read.csv.  The
## first table is made up of the step-stress data made up of 11 units tested from voltages 2V to
## 7V and the failures that occur on each step.  While the second table is the step description
## based on voltage and test time for each step.  This is based on Table 4.8 in the text.
datStepStress1 <- read.csv("https://raw.githubusercontent.com/Center-for-Risk-and-Reliability/RMT/main/CSVExampleData/Step_Stress_Data_1_by_Voltage_V_step_based_time_Example_Table_4_8_and_Example_4_17.csv")
tableStepStress1 <- read.csv("https://raw.githubusercontent.com/Center-for-Risk-and-Reliability/RMT/main/CSVExampleData/Step_Stress_Table_1_by_Voltage_V_step_based_time_Example_Table_4_8_and_Example_4_17.csv")
est1<-stepstress.LSQest(datStepStress1,tableStepStress1, "InversePower", "Weibull", "Blom")[[2]]
MLEest1<-stepstress.MLEest(est1,datStepStress1,tableStepStress1,"InversePower","Weibull",0.95)[[1]]
stepstress.BAYESest(MLEest1,datStepStress1,tableStepStress1,"InversePower","Weibull",0.95,c("uniform(1,100)","uniform(1,100)","uniform(1,100000)"),20000,100)


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
