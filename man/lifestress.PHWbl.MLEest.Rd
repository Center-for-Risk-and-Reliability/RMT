\name{lifestress.PHWbl.MLEest}
\alias{lifestress.PHWbl.MLEest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Likelihood Proportional Hazard-Weibull Life-Stress Estimator
}
\description{
Computes maximum likelihood estimate of the Proportional-Hazard Weibull Life-Stress model parameters.
}
\usage{
lifestress.PHWbl.MLEest(LSQest, data, interact_stress, confid, sided)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LSQest}{
The least-squares or another initial estimate for the maximum likelihood process to begin with (the life distribution parameter leads followed by the life-stress parameters).  Make sure the order matches the model being entered.  That is, if the output paramters are the form c(\eqn{\beta, a, b}), then the input must also be some form of c(\eqn{\beta, a, b}).
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
  Column(s) 3 (and on for multi-type stress) states the stress values for each data in column 1.  List one stress type per column (i.e. Column 3 - Temperature and Column 4 - Humidity).  Column 3 is not necessary if all data is from the same stress level however.
  }
}
}
  \item{interact_stress}{
A status which determines whether or not  the estimator will consider interactions between the various stresses or not.  If evaluation checks for interactions \code{interact_stress}, then set to 1, otherwise set to 0 (default).
}
  \item{confid}{
The confidence in the maximum likelihood parameter estimates entered as a probability between \eqn{0} and \eqn{1}.  Default is \eqn{0.95} (or \eqn{95\%})
}
  \item{sided}{
Indicates whether the limit calculations are based on a two-sided confidence or one-sided confidence.  The three options for this heading are \code{"twosided"} for two-sided confidence (default), \code{"onesidedlow"} for one-sided (low) confidence, and \code{"onesidedhigh"} for one-sided (high) confidence.
}
}
\details{
The maximum likelihood estimation is performed on the negated log-likelihood of the entered life-stress accelerated life model using the built-in \code{\link{nlm}} function, a non-linear minimization script that is based on a Newton-type algorithm.  Certain output will not respond well to a poor initial estimate \code{LSQest}.

The proportional hazard Weibull distribution is positively bounded (\eqn{x \geq 0}) consisting of a Weibull shape parameter \eqn{\beta}, and \eqn{S} and \eqn{\theta} are stress and parameters (respectively) of a modifier funciton \eqn{g(S;\theta)}.

  PDF:
  \deqn{f(x) = \beta x^{\beta - 1} \exp[\Sigma_{j=0}^m \theta_j S_j - x_j^\beta \exp(\Sigma_{j=0}^m \theta_j S_j)]}

  CDF:

  \deqn{F(x) = 1 - \exp[- x_j^\beta \exp(\Sigma_{j=0}^m \theta_j S_j))]}

  Modifier Function:

  \deqn{g(S;\theta) = \exp(\Sigma_{j=0}^m \theta_j S_j)}

The proportional hazard function is formed by multiplying the baseline hazard function (for a typical Weibull distribution) by the modifier funciton, from which we get the relations for PDF and CDF.
}
\value{
Outputs a printed summary of the life-stress distribution model, its maximum likelihood estimate of the parameters, and average life at each stress level.  The list output includes the following:
%%  If it is a LIST, use
\item{theta.hat}{
  Maximum likelihood mean estimate of the parameters (in order of \code{LSQest}).
}
\item{inv.fish}{
  The variance-covariance matrix as computed by non-linear minimization (in order of \code{LSQest}).
}
\item{conflim}{
  Vector-by-vector output of the confidence bounds for each parameter (in order of \code{LSQest}).
}
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2022)
}
\note{
Some of the confidence intervals may produce an NaN because of a negative value in the computed variance of the variance-covariance matrix.  Full assessment of non-linear minimization functions in R are in process.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1: Parametric PH-Weibull model with two continuous covariates (inverse temperature K^-1 and inverse relative humidity)
time <- c(102, 115, 151, 196, 210, 184, 193, 244, 265, 281, 298, 311, 337, 347, 405)
cens <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
invtemp <- c(0.002544529, 0.002544529, 0.002544529, 0.002544529, 0.002544529, 0.002832861, 0.002832861, 0.002832861, 0.002832861, 0.002832861, 0.002832861, 0.002832861, 0.002832861, 0.002832861, 0.002832861)
invRH <- c(1.666666667, 1.666666667, 1.666666667, 1.666666667, 1.666666667, 1.176470588, 1.176470588, 1.176470588, 1.176470588, 1.176470588, 1.666666667, 1.666666667, 1.666666667, 1.666666667, 1.666666667)
mat1 <- matrix(c(time, cens, invtemp,invRH), nrow = length(time), ncol = 4, byrow = FALSE,dimnames = list(NULL,c("time (hrs)","cens","invtemp","invrelhumidity")))
est1 <- lifestress.PHWbl.LSQest(probplot.wbl(mat1,"Blom","hours")[[1]], 0)[[3]]
lifestress.PHWbl.MLEest(est1,mat1,0,0.9)

## Example 2: Parametric PH-Weibull model with two continuous covariates (inverse temperature K^-1 and inverse relative humidity), but in this instance, consider the interactions between the two stresses
est2 <- lifestress.PHWbl.LSQest(probplot.wbl(mat1,"Blom","hours")[[1]], 1)[[3]]
lifestress.PHWbl.MLEest(est2,mat1,1,0.9)


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
