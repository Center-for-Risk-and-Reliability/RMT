\name{MLE.var.covar.select}
\alias{MLE.var.covar.select}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Likelihood Variance-Covariance Matrix Select
}
\description{
A short code designed to perform maximum likelihood estimation by way of several R-based non-linear estimation tools.  As a result of the complexity of the log-likelihoods that are formed by the life-stress distribution models created by the reliabilitytk package, the possibility of an improper estimate based on certain procedures is possible (example: a negative variance is often a result of an improper non-linear estimation procedure).  Therefore this short code chooses the procedure that produces a proper estimate.
}
\usage{
MLE.var.covar.select(loglik, LSQest)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loglik}{
The log-likelihood function entered as a function of a parameter vector \code{theta}.
}
  \item{LSQest}{
The initial estimate for the log-likelihood parameter vector.
}
}
\details{
The \code{MLE.var.covar.select} function operates on a tiered system of acceptance with regard to the maximum likelihood (ML) estimate  and its variance-covariance matrix where the selected set is contingent on the latter output.  Tier (1) uses the \code{\link{nlm}} function which computes ML estimates based on a Newton-type algorithm (see references in  \code{\link{nlm}} for more).  Tier (2) uses the \code{\link[ucminf]{ucminf}} function which is based on a quasi-Newton type algorithm with a BFSG (Broyden-Fletcher-Goldfarb-Shanno) updating procedure (see references in  \code{\link[ucminf]{ucminf}} for more).  The final tier (3) employs the \code{\link{optim}} function which is based on Nelder–Mead, quasi-Newton and conjugate-gradient algorithms to perform a general optimization (see references in  \code{\link{optim}} for more).
}
\value{
Returns the following:
%%  If it is a LIST, use
\item{theta.hat}{
  The chosen maximum likelihood parameter estimate based on the diagonal of the variance-covariance matrix.
}
\item{inv.fish}{
  The chosen variance-covariance matrix based on its diagonal (non-negative).
}
}
\references{
%% ~put references to the literature/web site here ~
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
## Example 1: Single Stress value (Temperature in Kelvin)
## Estimate based on Arrhenius-Weibull model (see example 1 in lifestress.LSQest help file)
time <- c(52, 80, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 20, 28, 35, 50, 50, 50, 50, 50, 50, 50, 4, 11, 12, 18, 25, 25, 25)
cens <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0)
temp <- c(175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 250, 250, 250, 250, 250, 250, 250) + 273.15
est1<-c(1.8662677563,0.5399077595,0.0001537369)
K<-8.617385e-5
lifeF1 <- function(theta) {
  theta[3]*exp(theta[2]/(K*temp[which(cens==1)]))
}
loglifeF1 <- function(theta) {
  log(theta[3]) + theta[2]/(K*temp[which(cens==1)])
}
lifeC1 <- function(theta) {
  theta[3]*exp(theta[2]/(K*temp[which(cens==0)]))
}
loglik1 <- function(theta){
  -sum(log(theta[1]) + (theta[1]-1)*log(time[which(cens==1)]) - theta[1]*loglifeF1(theta) - ((time[which(cens==1)]/lifeF1(theta))^theta[1])) - sum(- ((time[which(cens==0)]/lifeC1(theta))^theta[1]))
}
MLE.var.covar.select(loglik1,est1)


## Example 2: Dual Stress values (Temperature in Kelvin and Relative Humidity (RH))
## Estimate based on Temperature-Humidity-Weibull model (see example 2 in lifestress.LSQest help file)
time <- c(102, 115, 151, 196, 210, 184, 193, 244, 265, 281, 298, 311, 337, 347, 405)
temp <- c(393, 393, 393, 393, 393, 353, 353, 353, 353, 353, 353, 353, 353, 353, 353)
RH <- c(0.6, 0.6, 0.6, 0.6, 0.6, 0.85, 0.85, 0.85, 0.85, 0.85, 0.6, 0.6, 0.6, 0.6, 0.6)
est2 <- c(5.784249e,8.681083e-02,2.518336e+03,7.142889e-01)
lifeF2 <- function(theta) {
  theta[2]*exp((theta[3]/temp) + (theta[4]/RH))
}
loglifeF2 <- function(theta) {
  log(theta[2]) + (theta[3]/temp) + (theta[4]/RH)
}
loglik2 <- function(theta){
  -sum(log(theta[1]) + (theta[1]-1)*log(time) - theta[1]*loglifeF2(theta) - ((time/lifeF2(theta))^theta[1]))
}
MLE.var.covar.select(loglik2,est2)

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
