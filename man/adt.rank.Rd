\name{adt.rank}
\alias{adt.rank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Accelerated Degradation Testing Rank System
}
\description{
Ranks the fitness of a set of degradation data based on a list of degradation model.
}
\usage{
adt.rank(data)
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
  }}
}
}
\details{
This script identifies which of the following degradation-life models is a the best fit for a set of degradation data based on the coefficient of determination average obtained from the least-squares estimate (from \code{\link{adt.full.LSQ}}).  The following models are defined such that \eqn{D} is degradation, \eqn{L} is life, and \eqn{a} and \eqn{b} are model parameters:
  \enumerate{
  \item \strong{Linear} (Enter as "Linear"):
  \deqn{D = a + b L} \emph{Parameter order} - c(a,b)

  \item \strong{Exponential} (Enter as "Exponential")
  \deqn{D = b \exp(a L) } \emph{Parameter order} - c(a,b)

  \item \strong{Square-Root} (Enter as "SquareRoot")
  \deqn{D = (a + b L)^2} or \deqn{D^{1/2} = a + b L}
  \emph{Parameter order} - c(a,b)

  \item \strong{Square-Root v.2} (Enter as "SquareRoot2")
  \deqn{D = a + b (L^{1/2})}
  \emph{Parameter order} - c(a,b)

  \item \strong{Power} (Enter as "Power")
  \deqn{D = b L^a} \emph{Parameter order} - c(a,b)

  \item \strong{Logarithmic} (Enter as "Logarithmic")
  \deqn{D = a + b \ln(L)} \emph{Parameter order} - c(a,b)

  \item \strong{Lloyd-Lipow} (Enter as "LloydLipow")
  \deqn{D = a - \frac{b}{L}} \emph{Parameter order} - c(a,b)

  \item \strong{Mitsuom} (Enter as "Mitsuom")
  \deqn{D = \frac{1}{1 + b L^a}} \emph{Parameter order} - c(a,b) (\emph{NOTE}: The Mitsuom model only works for ratio based degradation that is degradation between 0 and 1.)
}
}
\value{
Outputs a printed summary of the degradation ranking by model.  The list output includes the following:
%%  If it is a LIST, use
\item{R2set}{
  A table made up of the model and unit based coefficient of determination.
}
\item{rankset}{
  A matrix made up of the model and unit-based ranking where 1 is the best rank.
}
\item{avgrankset}{
  The average model ranking for the degradation models.
}
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2021-2026)
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1:
## Download the CSV containing the metallic alloy test data from Table 5.2
datADT1 <- read.csv("https://raw.githubusercontent.com/ReuS009/reliabilitytk/main/CSVExampleData/Degradation_Data_1_Mass_Loss_by_Weight_gms_Example_5_2.csv")
adt.rank(datADT1)
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
