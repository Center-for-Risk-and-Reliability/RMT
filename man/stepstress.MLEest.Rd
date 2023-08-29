\name{stepstress.MLEest}
\alias{stepstress.MLEest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Likelihood Step-Stress Estimator
}
\description{
Computes maximum likelihood estimate of life step-stress model parameters.
}
\usage{
stepstress.MLEest(LSQest, data, stepstresstable, ls, dist, confid, sided)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LSQest}{
The least-squares or another initial estimate for the maximum likelihood process to begin with (the life distribution parameter leads followed by the life step-stress parameters).  Make sure the order matches the model being entered.  That is, if the output paramters are the form c(\eqn{\beta, a, b}), then the input must also be some form of c(\eqn{\beta, a, b}).
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
  \item{sided}{
Indicates whether the limit calculations are based on a two-sided confidence or one-sided confidence.  The three options for this heading are \code{"twosided"} for two-sided confidence (default), \code{"onesidedlow"} for one-sided (low) confidence, and \code{"onesidedhigh"} for one-sided (high) confidence.
}
}
\details{
The maximum likelihood estimation is performed on the negated log-likelihood of the entered step-stress accelerated life model using one of the built-in functions defined in the \code{\link{MLE.var.covar.select}} function.
}
\value{
Outputs a printed summary of the life step-stress distribution model, its maximum likelihood estimate of the parameters, and average life at each step-stress level.  The list output includes the following:
%%  If it is a LIST, use
\item{Output[[1]] - theta.hat}{
  Maximum likelihood mean estimate of the parameters (in order of \code{LSQest}).
}
\item{Output[[2]] - inv.fish}{
  The variance-covariance matrix as computed by non-linear minimization (in order of \code{LSQest}).
}
\item{Output[[3]] - conflim}{
  Vector-by-vector output of the confidence bounds for each parameter (in order of \code{LSQest}).
}
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2021-2022)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1: Single Stress Step-Stress
## Load the CSVs for step-stress data and the step-stress table from the Github site by read.csv. The
## first table is made up of the step-stress data made up of 11 units tested from voltages 2V to
## 7V and the failures that occur on each step.  While the second table is the step description
## based on voltage and test time for each step.  This is based on Table 4.8 in the text.
datStepStress1 <- read.csv("https://raw.githubusercontent.com/ReuS009/reliabilitytk/main/CSVExampleData/Step_Stress_Data_1_by_Voltage_V_step_based_time_Example_Table_4_8_and_Example_4_17.csv")
tableStepStress1 <- read.csv("https://raw.githubusercontent.com/ReuS009/reliabilitytk/main/CSVExampleData/Step_Stress_Table_1_by_Voltage_V_step_based_time_Example_Table_4_8_and_Example_4_17.csv")
est1<-stepstress.LSQest(datStepStress1,tableStepStress1, "InversePower", "Weibull", "Blom")[[2]]
stepstress.MLEest(est1,datStepStress1,tableStepStress1,"InversePower","Weibull",0.95)

## Example 2: Dual Stress Step-Stress
## Load the CSVs for step-stress data and the step-stress table from the Github site by read.csv.  The
## first table is made up of the step-stress data made up of 15 units tested from relative
## humidity 70\% to 90\% and temperature 525K to 555K.  While the second table is the step
## description based on relative humidity and temperature and test time for each step.  This is
## based on Table 4.18 and Example 4.16 in the text.
datStepStress2 <- read.csv("https://raw.githubusercontent.com/ReuS009/reliabilitytk/main/CSVExampleData/Step_Stress_Data_3_by_temperature_Kelvin_and_relative_humidity_percent_step_based_time_Example_4_16.csv")
tableStepStress2 <- read.csv("https://raw.githubusercontent.com/ReuS009/reliabilitytk/main/CSVExampleData/Step_Stress_Table_3_by_temperature_Kelvin_and_relative_humidity_percent_step_based_time_Example_4_16.csv")
est2<-stepstress.LSQest(datStepStress2,tableStepStress2, "TempHumidity", "Weibull", "Blom")[[2]]
stepstress.MLEest(est2,datStepStress2,tableStepStress2,"TempHumidity","Weibull",0.95)

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
