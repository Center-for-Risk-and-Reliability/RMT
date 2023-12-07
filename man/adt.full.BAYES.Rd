\name{adt.full.BAYES}
\alias{adt.full.BAYES}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Accelerated Degradation Testing Estimator
}
\description{
Computes Bayesian posterior of degradation model parameters based on either a lognormal or normal degradation fit.
}
\usage{
adt.full.BAYES(data, lifedam, dist, D0, Tuse, confid, priors, nsamples, burnin)
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
  \item{dist}{
  The distribution that represents time or life.
  \enumerate{
  \item \strong{Normal} (Enter as "Normal".  See more on \code{\link{probplot.nor}}.)
  \item \strong{Lognormal} (Enter as "Lognormal".  See more on \code{\link{probplot.logn}}.)
  }
}
  \item{D0}{
The degradation limit (or endurance limit) where "pseudo" failure times are defined.
}
  \item{Tuse}{
(\emph{Optional}) The use-level temperature which is needed only if the \emph{Hamada degradation-life model} is selected in \code{lifedam}.  If left blank, room temperature (293.15 \eqn{K}) is assumed for this temperature.
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
This script uses R's Stan library to formulate a stan file specific to the entered degradation-life model and parameter priors therein.  The created Stan file runs a full Bayesian analysis based on four Markov chains per run (number of will be controlled in future versions of the library).
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
  The 2.5, 5, 10, 25, 50, 75, 95, and 97.5\% quantile of the degradation-life posterior parameters.  Output[[4]]
  }
\item{
  plot1_MCtrace}{Plot of the Markov chains for each degradation-life posterior parameter.  Output[[5]]
  }
\item{
  plot2_hist}{Plot of the histograms for each degradation-life posterior parameter.  Output[[6]]
  }
\item{
  plot3_density}{Plot of the density for each degradation-life posterior parameter.  Output[[7]]
  }
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2022)
}
\note{
Suggest running the ADT Rank Calculator \code{\link{adt.rank}} before running this script.  As of version 1.2.6, this tool is operational.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1:
## Download the following CSV from GitHub which is made up of the
## full set of LED luminosity test data based on Table 5.10.  Per Example 5.5, failure occurs
## when luminosity drops to 0.5.  Normal use temperature is given as 20 deg C and the priors
## for
datADT4 <- read.csv("https://raw.githubusercontent.com/Center-for-Risk-and-Reliability/RMT/main/CSVExampleData/Degradation_Data_4_Luminosity_Drop_by_Temperature_Kelvin_Example_5_5.csv")
adt.full.BAYES(ptest=c(1,1,1,1),datADT4,"Hamada","Lognormal",0.5,293.15,0.9,c("uniform(0.000001,100)","uniform(0.00001,100)","uniform(-100,100)","uniform(0.00001,100)"), 5000,100)
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
