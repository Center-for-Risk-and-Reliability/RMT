\name{adt.full.addME.MLE}
\alias{adt.full.addME.MLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Likelihood Estimation Accelerated Degradation Testing Estimator (with additive measurement error)
}
\description{
Computes maximum likelihood estimate of degradation model parameters based on additive measurement error and normal degradation fit.
}
\usage{
adt.full.addME.MLE(data, lifedam, D0, confid, sided, addME, Tuse)
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
The degradation-life model definition for analysis.  Choose from the list of degradation-life models outlined in \code{\link{adt.rank}}.
}
  \item{D0}{
The degradation limit (or endurance limit) where "pseudo" failure times are defined.
}
  \item{confid}{
The confidence in the maximum likelihood parameter estimates entered as a probability between \eqn{0} and \eqn{1}.  Default is \eqn{0.95} (or \eqn{95\%})
}
  \item{sided}{
Indicates whether the limit calculations are based on a two-sided confidence or one-sided confidence.  The three options for this heading are \code{"twosided"} for two-sided confidence (default), \code{"onesidedlow"} for one-sided (low) confidence, and \code{"onesidedhigh"} for one-sided (high) confidence.
}
  \item{addME}{
The additive measurement error that is given as a numeric offset.
}
  \item{Tuse}{
(\emph{Optional}) The use-level temperature which is needed only if the \emph{Hamada degradation-life model} is selected in \code{lifedam}.  If left blank, room temperature (293.15 \eqn{K}) is assumed for this temperature.
}
}
\details{
This is a variation on the \code{\link{adt.full.MLE}} function in that it computes the maximum likelihood estimation is performed on the negated log-likelihood of the entered degradation-life model with a known degrdation measurement (additive) error.
}
\value{
Outputs a printed summary of the degradation-life distribution model and its maximum likelihood estimate of the parameters.  The list output includes the following:
%%  If it is a LIST, use
\item{theta.hat}{
  Maximum likelihood mean estimate of the parameters (in order of \code{LSQest}).
}
\item{inv.fish}{
  The variance-covariance matrix as computed by non-linear minimization (in order of \code{LSQest}).
}
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2022)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1:
## Download the CSV 'DegradationData1.csv' from the CRR website which is made up of the
## metallic alloy test data from Table 5.2.  Computes the maximum likelihood estimate of the
## pseudo failure times based on reaching a 50 microgram scar width with the condition that
## there is a 0.01 additive measurement error.
datADT1 <- read.csv("c:/Download/DegradationData1.csv")
adt.full.addME.MLE(datADT1,"Power",50,0.9,"twosided",0.01)
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
