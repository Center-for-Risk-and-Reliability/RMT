\name{lifestress.Binom.MLEest}
\alias{lifestress.Binom.MLEest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Likelihood Binomial Accelerated Life-Stress Estimator
}
\description{
Computes maximum likelihood estimate of the Binomial accelerated Life-Stress model parameters.
}
\usage{
lifestress.Binom.MLEest(data, interact_stress, weight0, confid, sided)
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
  \item{confid}{
The confidence in the maximum likelihood parameter estimates entered as a probability between \eqn{0} and \eqn{1}.  Default is \eqn{0.95} (or \eqn{95\%})
}
  \item{sided}{
Indicates whether the limit calculations are based on a two-sided confidence or one-sided confidence.  The three options for this heading are \code{"twosided"} for two-sided confidence (default), \code{"onesidedlow"} for one-sided (low) confidence, and \code{"onesidedhigh"} for one-sided (high) confidence.
}
}
\details{
The maximum likelihood estimation is performed on the negated log-likelihood of the entered life-stress accelerated life model using the built-in \code{\link{nlm}} function, a non-linear minimization script that is based on a Newton-type algorithm.  Certain output will not respond well to a poor initial estimate or data \code{data}.

The Binomial Acceleration distribution is discrete by definition and bound between 0 and 1.  \eqn{k} is a random variable indicating a value between 0 and \eqn{n} trials.  Parameters consist of a use-level probability of failure event \eqn{p_0}, accelerated stress differential levels \eqn{\DeltaS} = \eqn{S_acc - S_use}, and \eqn{\theta_j} coefficent parameter for each stress differential (and stress differential interaction if considered).

\deqn{Pr(k) = (n|k) [1 - (1 - p_0) exp[-\Sigma_(j=1)^m \theta_j \DeltaS_j]]^k [(1 - p_0) exp(\Sigma_(j=1)^m \theta_j \DeltaS_j)]}

}
\value{
Outputs a printed summary of the life-stress distribution model, its maximum likelihood estimate of the parameters, and average life at each stress level.  The list output includes the following:
%%  If it is a LIST, use
\item{Output[[1]]}{
  Maximum likelihood mean estimate of the parameters (in the order listed in the console output).
}
\item{Output[[2]]}{
  The variance-covariance matrix as computed by non-linear minimization (in the order listed in the console output).  NOTE: The variance for \eqn{p_0} is currently given as \eqn{ln(-ln(p_0))}.
}
\item{Output[[3]]}{
  Vector-by-vector output of the confidence bounds for each parameter (in the order listed in the console output).
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mohammad Modarres and Reuel Smith (2022)
}
\note{
MLE operation performed on the \eqn{p_0} variable requires a pre-evaluation transformation to the pseudovariable \eqn{ln(-ln(p_0))} in order to maintain \eqn{p_0}'s bound between \eqn{0} and \eqn{1}.
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
OUT <- lifestress.Binom.MLEest(datbinom2,0,1,0.9)

## Example 2:
## Download the CSV 'BinomialData2b.csv' from the CRR website which is made up of
## binomial test data where the stress differentials are temperature in Kelvin
## (column 3) and voltage (column 4).  Weight is considered so it is listed in
## the last column.
datbinom2b <- read.csv("c:/Download/BinomialData2b.csv")
OUT <- lifestress.Binom.MLEest(datbinom2b,0,1,0.9)

## Example 3:
## Repeat example 2 but with consideration of stress differential interaction.
OUT <- lifestress.Binom.MLEest(datbinom2b,1,1,0.9)


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
