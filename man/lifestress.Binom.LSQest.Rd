\name{lifestress.Binom.LSQest}
\alias{lifestress.Binom.LSQest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Least-Squares Binomial Accelerated Life-Stress Estimator
}
\description{
Computes least-squares estimate of the Binomial acceleration Life-Stress model parameters.
}
\usage{
lifestress.Binom.LSQest(data, S0, interact_stress)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A matrix made up of necessary test data.  The first column is made up of the \eqn{n} data and the second column is made up of the \eqn{k} data.  The remaining columns are made up of the stress differential data (single or multiple).  If data is weighed however, that goes in the last column after the stress differential data columns.
}
  \item{interact_stress}{
A status which determines whether or not  the estimator will consider interactions between the various stresses or not.  If evaluation checks for interactions \code{interact_stress}, then set to 1, otherwise set to 0 (default).
}
  \item{weight0}{
A status which determines whether or not there is weight considered in the evaluation or not (for MLE step primarirly).  If evaluation accounts for weight, then set to 1, otherwise set to 0 (default).
}
}
\details{
This tool computes least-squares estimate of the Binomial Acceleration Life-Stress model parameters as defined.  The Binomial Acceleration distribution is discrete by definition and bound between 0 and 1.  \eqn{k} is a random variable indicating a value between 0 and \eqn{n} trials.  Parameters consist of a use-level probability of failure event \eqn{p_0}, accelerated stress differential levels \eqn{\DeltaS} = \eqn{S_acc - S_use}, and \eqn{\theta_j} coefficent parameter for each stress differential (and stress differential interaction if considered).

\deqn{Pr(k) = (n|k) [1 - (1 - p_0) exp[-\Sigma_(j=1)^m \theta_j \DeltaS_j]]^k [(1 - p_0) exp(\Sigma_(j=1)^m \theta_j \DeltaS_j)]}

Estimates are obtained through \code{\link{nls}}
}
\value{
Outputs a printed summary of the life-stress distribution model, its least-squares estimate, and average life at each stress level.  The list output includes the following:
\item{Output[[1]]}{
  The least squares estimate for the model parameters (\eqn{p_0} parameter leads followed by the stress-differential coefficient parameters)
}
\item{Output[[2]]}{
  The vector or matrix of accelerated stress differentials
}
\item{Output[[3]]}{
  The model parameters used in the analysis
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mohammad Modarres and Reuel Smith (2022)
}
\note{
Some data may produce a \eqn{p_0} estimate that is out of bounds.  That may indicate outlier(s) or incorrect input.  \code{lifestress.Binom.LSQest} will ask to reevaluate input data \code{data} in that event.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1:
## Download the CSV 'BinomialData2.csv' from the CRR website which is made up of
## binomial test data where the stress differential is temperature in Kelvin.  Weight
## is considered so it is listed in the last column.
datbinom2 <- read.csv("c:/Download/BinomialData2.csv")
OUT <- lifestress.Binom.LSQest(datbinom2,0,1)

## Example 2:
## Download the CSV 'BinomialData2b.csv' from the CRR website which is made up of
## binomial test data where the stress differentials are temperature in Kelvin
## (column 3) and voltage (column 4).  Weight is considered so it is listed in
## the last column.
datbinom2b <- read.csv("c:/Download/BinomialData2b.csv")
OUT <- lifestress.Binom.LSQest(datbinom2b,0,1)

## Example 3:
## Repeat example 2 but with consideration of stress differential interaction.
OUT <- lifestress.Binom.LSQest(datbinom2b,1,1)
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
