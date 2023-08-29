\name{stepstress.data.cum}
\alias{stepstress.data.cum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Step-Stress Cumulative Data Retabulator
}
\description{
Prepares the step-stress data into a cumulative form based on given step.  Failure times are documented from reference of the first step.
}
\usage{
stepstress.data.cum(data, stepstresstable)
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
}
\details{
An expanded form of \code{stepstress.data} where the failure times are documented as cumulative from the reference of the first step.  The function includes multiple output for the purpose of step-stress evaluation using maximum likelihood estimation.
}
\value{
Returns the following list:
%%  If it is a LIST, use
  \item{Output[[1]]}{
    A matrix of the reprocessed data that includes the censored times from all of the steps.  This includes steps where no failures are observed.
  }
  \item{Output[[2]]}{
  Vector comprised of the initial time of a given step.
  }
  \item{Output[[3]]}{
  Vector comprised of the initial time of a given step where failed units are observed.
  }
  \item{Output[[4]]}{
  Vector comprised of the initial time of a given step where units are censored.
  }
  \item{Output[[5]]}{
  Restatement of the step-stress table \code{stepstresstable}
  }
  \item{Output[[6]]}{
  Vector comprised of the steps where failed units are observed.  Relates with \code{Output[[3]]}.
  }
  \item{Output[[7]]}{
  Vector comprised of the steps where units are censored.  Relates with \code{Output[[4]]}.
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
stepstress.data.cum(datStepStress1,tableStepStress1)

## Example 2: Dual Stress Step-Stress
## Load the CSVs for step-stress data and the step-stress table from the Github site by read.csv.  The
## first table is made up of the step-stress data made up of 15 units tested from relative
## humidity 70\% to 90\% and temperature 525K to 555K.  While the second table is the step
## description based on relative humidity and temperature and test time for each step.  This is
## based on Table 4.18 in the text.
datStepStress2 <- read.csv("https://raw.githubusercontent.com/ReuS009/reliabilitytk/main/CSVExampleData/Step_Stress_Data_3_by_temperature_Kelvin_and_relative_humidity_percent_step_based_time_Example_4_16.csv")
tableStepStress2 <- read.csv("https://raw.githubusercontent.com/ReuS009/reliabilitytk/main/CSVExampleData/Step_Stress_Table_3_by_temperature_Kelvin_and_relative_humidity_percent_step_based_time_Example_4_16.csv")
stepstress.data.cum(datStepStress2,tableStepStress2)


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
