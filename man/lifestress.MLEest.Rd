\name{lifestress.MLEest}
\alias{lifestress.MLEest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Likelihood Life-Stress Estimator
}
\description{
Computes maximum likelihood estimate of Life-Stress model parameters.
}
\usage{
lifestress.MLEest(LSQest,ls,dist,TTF,SF,Tc=NULL,Sc=NULL,Suse=NULL,confid=0.95,sided="twosided",pp="Blom",xlabel1="X")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LSQest}{
The least-squares or another initial estimate for the maximum likelihood process to begin with (the life distribution parameter leads followed by the life-stress parameters).  Make sure the order matches the model being entered.  That is, if the output paramters are the form c(\eqn{\beta, a, b}), then the input must also be some form of c(\eqn{\beta, a, b}).
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
  \item{TTF}{
A vector consisting of the failure or primary event data in a given set.  Output[[1]] of \code{sort.xircstressdata}.
}
  \item{SF}{
A vector consisting of the stress at the failure or primary event data \code{xi} in a given set.  In the event that there are multiple stress factors, the output will be a matrix.  Output[[3]] of \code{sort.xircstressdata}.
}
  \item{Tc}{
A vector consisting of the right censored data of a given set of data (NOTE: If no censored data exisits, this will return \code{NULL}).  Output[[2]] of \code{sort.xircstressdata}.
}
  \item{Sc}{
A vector consisting of the stress at the right censored data \code{rc} of a given set of data.  In the event that there are multiple stress factors, the output will be a matrix (NOTE: If no censored data exisits, this will return \code{NULL}).  Output[[4]] of \code{sort.xircstressdata}.
}
  \item{Suse}{
(optional) The use or nominal stress level for the accletated life test scenario.  If dual or multiple stress, enter as a vector starting with temperature based stress.  Default input is NULL.
}
  \item{confid}{
The confidence in the maximum likelihood parameter estimates entered as a probability between \eqn{0} and \eqn{1}.  Default is \eqn{0.95} (or \eqn{95\%})
}
  \item{sided}{
Indicates whether the limit calculations are based on a two-sided confidence or one-sided confidence.  The three options for this heading are \code{"twosided"} for two-sided confidence (default), \code{"onesidedlow"} for one-sided (low) confidence, and \code{"onesidedhigh"} for one-sided (high) confidence.
}
  \item{pp}{
Name of the plotting position to select for computation.  See \code{\link{plotposit.select}} for plotting position definitions.  The defaulr is "Blom".
}
  \item{xlabel1}{
Label for the x-axis for the probability plot.  If left blank, the default label is "X".
}
}
\details{
The maximum likelihood estimation is performed on the negated log-likelihood of the entered life-stress accelerated life model using one of the built-in functions defined in the \code{\link{MLE.var.covar.select}} function.
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
\item{uselife}{
  The mean use life based on the use/nominal stress level and the mean maximum likelihood estimate of the parameters.  This output will only be presented if \code{Suse} is set as an input.
}
\item{conflim}{
  Vector-by-vector output of the confidence bounds for each parameter (in order of \code{LSQest}).
}
\item{uselifelim}{
  Thee confidence bounds for the use life based on the use/nominal stress level and the mean maximum likelihood estimate of the parameters.  This output will only be presented if \code{Suse} is set as an input.
}
\item{AIC}{
  The Akaike Information Criterion value.
}
\item{BIC}{
  The Bayesian Information Criterion value.
}
\item{$plotoutput}{
  The probability plot of the data by stress level based on the maximum likelihood estimate of the parameters.
}
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2021-2024)
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1: Single Stress value (Temperature in Kelvin)
time <- c(52, 80, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 20, 28, 35, 50, 50, 50, 50, 50, 50, 50, 4, 11, 12, 18, 25, 25, 25)
cens <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0)
temp <- c(175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 250, 250, 250, 250, 250, 250, 250) + 273.15
mat1 <- cbind(time, cens, temp)
est1 <- lifestress.LSQest(mat1,"Arrhenius","Weibull","Blom")[[3]]
ttf1 <- sort.xircstressdata(mat1)[[1]]
ttrc1 <- sort.xircstressdata(mat1)[[2]]
Sf1 <- sort.xircstressdata(mat1)[[3]]
Src1 <- sort.xircstressdata(mat1)[[4]]
lifestress.MLEest(est1,"Arrhenius","Weibull",ttf1,Sf1,ttrc1,Src1,confid=0.90)

## Example 2: USe the failure data from Example 4.1 of the textbook to find the MLE for the parameters and the use life given the use temperature of 333 Kelvin.
time <- c(248,456,528,731,813,965,972,1528,164,176,289,319,386,459,92,105,155,184,219,235)
cens <- rep(1,20)
temp <- c(rep(406,8),rep(436,6),rep(466,6))
mat2 <- cbind(time,cens,temp)

# Estimate the LSQ parameters for MLE estimate
est2<-lifestress.LSQest(mat2,"Arrhenius","Weibull","Blom")[[3]]
lifestress.MLEest(est2,"Arrhenius","Weibull",time,temp,Suse = 333)

## Example 3: Dual Stress values (Temperature in Kelvin and Relative Humidity (RH))
time <- c(102, 115, 151, 196, 210, 184, 193, 244, 265, 281, 298, 311, 337, 347, 405)
cens <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
temp <- c(393, 393, 393, 393, 393, 353, 353, 353, 353, 353, 353, 353, 353, 353, 353)
RH <- c(0.6, 0.6, 0.6, 0.6, 0.6, 0.85, 0.85, 0.85, 0.85, 0.85, 0.6, 0.6, 0.6, 0.6, 0.6)
mat3 <- cbind(time, cens, temp, RH)
est3 <- lifestress.LSQest(mat3,"TempHumidity","Weibull","Blom")[[3]]
ttf3 <- sort.xircstressdata(mat3)[[1]]
ttrc3 <- sort.xircstressdata(mat3)[[2]]
Sf3 <- sort.xircstressdata(mat3)[[3]]
Src3 <- sort.xircstressdata(mat3)[[4]]
lifestress.MLEest(est3,"TempHumidity","Weibull",ttf3,Sf3,ttrc3,Src3,confid=0.90)

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
