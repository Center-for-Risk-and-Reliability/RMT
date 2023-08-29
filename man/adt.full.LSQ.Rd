\name{adt.full.LSQ}
\alias{adt.full.LSQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Least-Squares Accelerated Degradation Testing Estimator
}
\description{
Computes least squares estimate of degradation-life model parameters based on a linearized degradation fit.
}
\usage{
adt.full.LSQ(data, lifedam, D0, Tuse)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A matrix that represents a data set series which may consist of the following:
\itemize{
  \item{
  Column 1 consists of the time or life entries which correspond to a degradation value.  Data must consist of real and finite numbers.
  }
  \item{
  Column 2 is the degradation value that corresponds to a specific time or life value (Column 1).
  }
  \item{
  Column 3 defines the label for the specific unit number for a given set of degradation data.  For example, the first four entries of this column may be labeled 'Unit 1' which labels the degradation for that particular unit while the following five entries may be labeled 'Unit 2' which represents that unit's degradation.
  }
  \item{
  Column(s) 4 (and on for multi-type stress) states the stress (or step-stress) values for each data in columns 1 and 2.  List one stress type per column (i.e. Column 4 - Temperature and Column 5 - Humidity).
  }
  }
  }
  \item{lifedam}{
The degradation-life model definition for analysis.  Choose from the list of degradation-life models outlined in \code{\link{adt.rank}}.}
  \item{D0}{
The degradation limit (or endurance limit) where "pseudo" failure times are defined.
}
  \item{Tuse}{
(\emph{Optional}) The use-level temperature which is needed only if the \emph{Hamada degradation-life model} is selected in \code{lifedam}.  If left blank, room temperature (293.15 \eqn{K}) is assumed for this temperature.
}
}
\details{
Estimator provides the least-squares estimate for the unit-based degradation-life model fit and the pseudo-failure time based on the degradation limit given in \code{D0}.
}
\value{
The output consists of two tables made up of the following:
\item{Output[[1]]}{
  A row number-by-four matrix where the number of rows are the number of units, the first and second columns are made up of the unit-specific LSQ parameter estimates for \eqn{a} and \eqn{b} respectively, the third column are the pseudo failure times, and column four are the coefficient of determination values for each unit.
}
\item{Output[[2]]}{
A matrix that represents a data set series which may consist of the following:
\itemize{
  \item{
  Column 1 consists of the pseudo failure (or primary event) data and the right censored data in a given set.  Data must consist of real and finite numbers.
  }
  \item{
  Column 2 gives the censored status of the data in column 1.  All data in this column are demarked as \eqn{1} because they are pseudo and \emph{assumed} to have all failed.
  }
  \item{
  Column(s) 3 (and on for multi-type stress) states the stress values for each data in column 1.  List one stress type per column (i.e. Column 3 - Temperature and Column 4 - Humidity).  Column 3 is not necessary if all data is from the same stress level however.
  }
}
Data is set up such that it may be used directly in probablity plotting functions: : (1) \code{\link{probplot.wbl}}, (2) \code{\link{probplot.wbl3P}}, (3) \code{\link{probplot.exp}}, (4) \code{\link{probplot.exp2P}}, (5) \code{\link{probplot.nor}}, and (6) \code{\link{probplot.logn}}.
}
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2021-2022)
}
\note{
Suggest running the ADT Rank Calculator \code{\link{adt.rank}} before running this script.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1:
## Download the CSV containing the metallic alloy test data from Table 5.2
## Computes the Least Squared estimate of the
## pseudo failure times based on reaching a 50 microgram scar width.
datADT1 <- read.csv("https://raw.githubusercontent.com/ReuS009/reliabilitytk/main/CSVExampleData/Degradation_Data_1_Mass_Loss_by_Weight_gms_Example_5_2.csv")
OUT <- adt.full.LSQ(datADT1,"Power",50)


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
