\name{stepstress.LSQest}
\alias{stepstress.LSQest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Least-Squares Step-Stress Estimator
}
\description{
Computes least-squares parameter estimate of a step-stress life distribution model.
}
\usage{
stepstress.LSQest(data, stepstresstable, ls, dist, pp, therm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
  }
}
  \item{pp}{
Name of the plotting position to select for computation.  See \code{\link{plotposit.select}} for plotting position definitions.
}
  \item{therm}{
Designation of the primary thermal stress in the original data.  If it is the first column after the censored column then the value is 1 and the alternate is 2 (this is set to default).  If it is the second column after the censored column then the value is 2 and the alternate is 1.  Only necessary to input for dual-stress life-stress models.
}
}
\details{
The least-squares estimate of an entered step-stress life distribution model is based primarily of the step-stress-based probability plot of the entered data (usually in the \code{probplot} subfunction).  The step-stress life parameters are estimated by the relationship linear analysis between stress and a parameter that is based on life (for Weibull or three-parameter Weibull that is \eqn{\alpha}, for normal that is \eqn{\mu}, for lognormal that is \eqn{\mu_t}, for exponential that is \eqn{\lambda}, and for two parameter exponential that is \eqn{\sigma}).  The final distribution paramter is obtained by averaging through the multiple step-stress lines.
}
\value{
The output of the least-squares analysis for a given step-stress life model include a substantial list of items needed for both basic estimates from a leeast-squares fit and for further evalutation for a refined maximum-likelihood estimate \code{\link{stepstress.MLEest}}.
%%  If it is a LIST, use
  \item{Output[[1]]}{
  The initial step-stress life model parameter estimate based on the retabluated data from \code{\link{stepstress.data}}.  The life distribution parameter leads followed by the life-stress parameters.
  }
  \item{Output[[2]]}{
  The post optimization step-stress life model parameter estimate based on the retabluated data from \code{\link{stepstress.data}} and optimization by adjusted acceleration factors.  The life distribution parameter leads followed by the life-stress parameters.
  }
  \item{Output[[3]]}{
  The step-stress based mean life based on the optimized parameter estimates.
  }
  \item{Output[[4]]}{
  An \eqn{n x n} matrix that represents the "from" step-stress state where stresses are listed by column and repeated.  Needed for both least-squares and maximum likelihood analysis.
  }
  \item{Output[[5]]}{
  An \eqn{n x n} matrix that represents the "to" step-stress state where stresses are listed by row and repeated.  Needed for both least-squares and maximum likelihood analysis.
  }
  \item{Output[[6]]}{
  An \eqn{n x n} matrix that represents the "from/to" acceleration factors based on the "from" and "to" step-stress matricies.  Needed for both least-squares and maximum likelihood analysis.
  }
  \item{Output[[7]]}{
  An matrix that represents the "from/to" adjusted acceleration factors based on the equivalent times and the cumulative damage.  Only relevant for steps where damage occurs.
  }
  \item{Output[[8]]}{
  The adjusted step-stress life parameters based on the adjusted acceleration factors for all \eqn{n} steps.  The average of which is used to update the next iteration of the solution.
  }
  \item{Output[[9]]}{
  An \eqn{n x n} matrix that represents the test times for each step-stress state where times are listed by column and repeated.  Needed for both least-squares and maximum likelihood analysis.
  }
  \item{Output[[10]]}{
  The step equivalent times for all \eqn{n} steps.
  }
  \item{Output[[11]]}{
  The difference between the step equivalent times and the oriinal test times for all \eqn{n} steps.
  }
  \item{Output[[12]]}{
  The adjusted failure data table for probability plotting based on the optimized step-stress life least-square analysis.
  }
%% ...
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
## Load the CSVs for step-stress data and the step-stress table from the Github site by read.csv.  The
## first table is made up of the step-stress data made up of 11 units tested from voltages 2V to
## 7V and the failures that occur on each step.  While the second table is the step description
## based on voltage and test time for each step.  This is based on Table 4.8 in the text.
datStepStress1 <- read.csv("https://raw.githubusercontent.com/ReuS009/reliabilitytk/main/CSVExampleData/Step_Stress_Data_1_by_Voltage_V_step_based_time_Example_Table_4_8_and_Example_4_17.csv")
tableStepStress1 <- read.csv("https://raw.githubusercontent.com/ReuS009/reliabilitytk/main/CSVExampleData/Step_Stress_Table_1_by_Voltage_V_step_based_time_Example_Table_4_8_and_Example_4_17.csv")
stepstress.LSQest(datStepStress1,tableStepStress1, "InversePower", "Weibull", "Blom")

## Example 2: Dual Stress Step-Stress
## Load the CSVs for step-stress data and the step-stress table from the Github site by read.csv.  The
## first table is made up of the step-stress data made up of 15 units tested from relative
## humidity 70\% to 90\% and temperature 525K to 555K.  While the second table is the step
## description based on relative humidity and temperature and test time for each step.  This is
## based on Table 4.18 and Example 4.16 in the text.
datStepStress2 <- read.csv("https://raw.githubusercontent.com/ReuS009/reliabilitytk/main/CSVExampleData/Step_Stress_Data_3_by_temperature_Kelvin_and_relative_humidity_percent_step_based_time_Example_4_16.csv")
tableStepStress2 <- read.csv("https://raw.githubusercontent.com/ReuS009/reliabilitytk/main/CSVExampleData/Step_Stress_Table_3_by_temperature_Kelvin_and_relative_humidity_percent_step_based_time_Example_4_16.csv")
stepstress.LSQest(datStepStress2,tableStepStress2, "TempHumidity", "Weibull", "Blom")

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
