\name{optimizer.BAYES.alt}
\alias{optimizer.BAYES.alt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Accelerated Life Testing Bayes Optimization
}
\description{
\strong{NOTE: Function and help file not complete yet as of version 0.1.1.}


An accelerated life testing optimization tool that computes the posterior Bayesian estimation of a predefined two-stress-level life-stress test.
}
\usage{
optimizer.BAYES.alt(Strset,n,tc,p,pt_est,ls,dist,confid,priors,nruns,nsamples,burnin)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Strset}{
  A vector of your design stress levels that are defined c(S0,S1,S2).  S0 is the given use or nominal stress, S1 is the first accelerated stress and S2 is the second accelerated stress at time \eqn{\tau} the stress changing point.  For dual or greater life-stress models enter vector as c(S0,H0,S1,H1,S2,H2) where (S,H) represents a dual stress type.
}
\item{n}{
  Number of items to be tested.
}
\item{tc}{
  Predetermined censoring time for the test.
}
\item{p}{
  The given percentile that is used to determine the time to failure at the use or nominal stress t_p.
}
\item{pt_est}{
  The initial point estimate for your model.
}
\item{ls}{
  The life-stress model definition for analysis.  Choose from the following life-models where \eqn{L} is life, \eqn{S} is stress, and \eqn{a,b,c,d} are model parameters:
  \enumerate{
  \item \strong{Linear} (Enter as "Linear"):
  \deqn{L = b + a S} \emph{Parameter order} - c(a,b)

  \item \strong{Exponential} (Enter as "Exponential")
  \deqn{L = b exp(a S)} \emph{Parameter order} - c(a,b)

  \item \strong{Arrhenius} (Enter as "Arrhenius")
  \deqn{L = b exp(E_a/(K S))} \emph{Parameter order} - c(E_a,b) (\emph{NOTE}: K = \eqn{8.617385 x 10^-5 eV/K})

  \item \strong{Exponential} (Enter as "Exponential")
  \deqn{L = b exp(a S)} \emph{Parameter order} - c(a,b)

  \item \strong{Eyring} (Enter as "Eyring")
  \deqn{L = (b/S) exp(a/S)} \emph{Parameter order} - c(a,b)

  \item \strong{Eyring \emph{(Alternative Notation)}} (Enter as "Eyring2")
  \deqn{L = (1/S) exp(-(a - (b/S)))} \emph{Parameter order} - c(a,b)

  \item \strong{Power} (Enter as "Power")
  \deqn{L = b S^a} \emph{Parameter order} - c(a,b)

  \item \strong{Inverse Power} (Enter as "InversePower")
  \deqn{L = b S^-a or L = b / S^a} \emph{Parameter order} - c(a,b)

  \item \strong{Logarithmic} (Enter as "Logarithmic")
  \deqn{L = b + a ln(S)} \emph{Parameter order} - c(a,b)

  \item \strong{Temperature-Humidity} (Enter as "TempHumidity")
  \deqn{L = A exp((a/S) + (b/H))} \emph{Parameter order} - c(A,a,b) (\emph{NOTE}: \eqn{S} is temperature and \eqn{H} is humidity)

  \item \strong{Temperature-Nonthermal (\emph{or Power-Exponential})} (Enter as "TempNonthermal")
  \deqn{L = c/((U^b) exp(-a/S))} \emph{Parameter order} - c(a,b,c) (\emph{NOTE}: \eqn{S} is temperature and \eqn{U} is the nonthermal stress)

  \item \strong{Generalized Eyring (\emph{two-stress})} (Enter as "Eyring3")
  \deqn{L = (1/S) exp((a + (b/S)) + (c + (d/S)) U)} \emph{Parameter order} - c(a,b,c,d) (\emph{NOTE}: \eqn{S} is temperature and \eqn{U} is the nonthermal stress)

  \item \strong{General Multi-Stress} (Enter as "MultiStress")
  \deqn{L = exp(a_0 + a_1*S_1 + a_2*S_2 + ...+ a_n*S_n)} \emph{Parameter order} - c(a_1,a_2,...,a_n)
  }
}
\item{dist}{
  The distribution that represents time or life.
  \enumerate{
  \item \strong{Weibull} (Enter as "Weibull")

  CDFs:

  \deqn{F(0 < t <= \tau) = 1 - exp[-(t/L_1)^\beta]}
  \deqn{F(\tau < t) = 1 - exp{-[\tau/L_1 + (t - \tau)/L_2]^\beta}}

  PDFs:
  \deqn{f(0 < t <= \tau) = (\beta/L_1) (t/L_1)^(\beta - 1) exp[-(t/L_1)^\beta]}
  \deqn{f(\tau < t) = (\beta/L_2) [\tau/L_1 + (t - \tau)/L_2]^(\beta - 1) exp{-[\tau/L_1 + (t - \tau)/L_2]^\beta}}

  \item \strong{Three Parameter Weibull} (Enter as "3PWeibull")
  \item \strong{Exponential} (Enter as "Exponential")
  \item \strong{Two Parameter Exponential} (Enter as "2PExponential")
  \item \strong{Normal} (Enter as "Normal")
  \item \strong{Lognormal} (Enter as "Lognormal")
  }
}
\item{confid}{
  Confidence interval \eqn{100(1 - \alpha)} used for the posterior parameters and percentile time (or time unit) to failure.
}
\item{priors}{
  The prior distribution for each of the life-stress model parameters entered as a vector.  Each prior is entered as a quoted statement and must adhere to RStan notation.

  Example: priors <- c("normal(3,4)","normal(1,4)", "lognormal(-2,3)")
}
\item{nruns}{
  Number of iterations that the analysis will run.
}
\item{nsamples}{
  Number of samples generated for each Markov chain.
}
\item{burnin}{
  The number of samples each iteration burns through and throws out initially.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
List of Life-Stress Models
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Dr. Reuel Smith
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
# Accelerated Life Testing Bayesian Optimizer
# Developed by Dr. Reuel Smith, 2022

optimizer.BAYES.alt <- function(Strset,n,tc,p,pt_est,ls,dist,confid,priors,nruns,nsamples,burnin){
  # NOTE: For right now I will stick to a few life distributions (mainly the Weibull,
  # Normal, Exponential, and Lognormal), and maintain the life-stress model variety.
  # - V1.0 will also only stick to a two-stress-level ALT optimization as outlined
  #   in the current rendition of the textbook and Yuan, Liu, & Kuo (2012).
  # - V2.o should implement three-stress-level analysis to update its usefulness.
  #   This too is outlined in Yuan, Liu, & Kuo (2012).
  # - V3.0 should have an n-stress-level general analysis.

  #Load pracma library for erf
  library(pracma)
  library(StanHeaders)
  library(rstan)
  library(ggplot2)
  library(shinystan)
  library(bayesplot)

  # Add input to this to include prior estimates for LS parameters.
  # Example: priors<-c("normal(3,4),normal(1,4), lognormal(-2,3)")
  # I will have to cite the Rstan text for distributions in the code.  Use lookup("") for the translation.
  # The code takes these and separates them so that they are written into the stan file.
  # Then the code will run the program and compute the Bayes estimation

  # Check to see if dist="Exponential" so you can exclude life
  # distribution parameters.
  if (dist=="Exponential") {
    ishift<-0
  } else {
    ishift<-1
  }
  # Check to see if confidence exists
  if(missing(confid)){
    conf.level <- 0.95
  } else {
    conf.level <- confid
  }

  # Check to see if estimate exists
  if(missing(pt_est)){
    pt_est <- 'random'
  } else {
    pt_est <- pt_est
  }

  # Check to see if burn-in exists
  if(missing(burnin)){
    burnin <- floor(nsamples/2)
  } else {
    burnin <- burnin
  }

  # Initialize life-stress parameter estimates for theta
  if (ls=="Linear") {
    # theta[1] - parameter a, theta[2] - parameter b
    lsparams <- "real a; real b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*((a*Si[i] + b)/(a*Sni[i] + b));"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*((a*Sj[j] + b)/(a*Snj[j] + b));"

    # Life functions
    life1F <- "Life1i[i] = b + S1i[i]*a"
    loglife1F <- "Life1i[i] = log(b + S1i[i]*a)"
    life2F <- "Life2j[j] = b + S2j[j]*a"
    loglife2F <- "Life2j[j] = log(b + S2j[j]*a)"
    lifeC <- "Lifek[k] = b + Sk[k]*a"
    loglifeC <- "Lifek[k] = log(b + Sk[k]*a)"
  }

  if (ls=="Exponential"){
    # theta[1] - parameter a, theta[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*exp(a*(Si[i] - Sni[i]));"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*exp(a*(Sj[j] - Snj[j]));"


    # Life functions
    life1F <- "Life1i[i] = b*exp(a*S1i[i])"
    loglife1F <- "Life1i[i] = log(b) + a*S1i[i]"
    life2F <- "Life2j[j] = b*exp(a*S2j[j])"
    loglife2F <- "Life2j[j] = log(b) + a*S2j[j]"
    lifeC <- "Lifek[k] = b*exp(a*Sk[k])"
    loglifeC <- "Lifek[k] = log(b) + a*Sk[k]"
  }

  if (ls=="Arrhenius") {
    # lsparams[1] - parameter Ea, lsparams[2] - parameter b
    # Temperature HAS to be in Kelvin for this to work
    lsparams <- "real<lower=0> Ea; real<lower=0> b;"
    lsparamsvec <- c("Ea","b")
    pr1<-paste(c("Ea ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    life1F <- "b*exp(Ea/((8.617385e-5)*Si))"
    loglife1F <- "log(b) + (Ea/((8.617385e-5)*Si))"
    life2F <- "b*exp(Ea/((8.617385e-5)*Sj))"
    loglife2F <- "log(b) + (Ea/((8.617385e-5)*Sj))"
    lifeC <- "b*exp(Ea/((8.617385e-5)*Sk))"
    loglifeC <- "log(b) + (Ea/((8.617385e-5)*Sk))"

    K<-8.617385e-5
    lifeF <- function(theta,S) {
      theta[ishift+2]*exp(theta[ishift+1]/(K*S))
    }
    loglifeF <- function(theta,S) {
      log(theta[ishift+2]) + theta[ishift+1]/(K*S)
    }
  }

  if (ls=="Eyring") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*exp((Ea/8.617385e-5)*((1/Si[i]) - (1/Sni[i])));"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*exp((Ea/8.617385e-5)*((1/Sj[j]) - (1/Snj[j])));"

    # Life functions
    lifeF <- "Lifei[i] = (b/Si[i])*exp(a/Si[i])"
    loglifeF <- "Lifei[i] = log(b) - log(Si[i]) + (a/Si[i])"
    lifeC <- "Lifej[j] = (b/Sj[j])*exp(a/Sj[j])"
    loglifeC <- "Lifej[j] = log(b) - log(Sj[j]) + (a/Sj[j])"
  }

  if (ls=="Eyring2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*exp((Ea/8.617385e-5)*((1/Si[i]) - (1/Sni[i])));"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*exp((Ea/8.617385e-5)*((1/Sj[j]) - (1/Snj[j])));"

    # Life functions
    lifeF <- "Lifei[i] = (1/Si[i])*exp(-(a - (b/Si[i])))"
    loglifeF <- "Lifei[i] = -log(Si[i]) - a + (b/Si[i])"
    lifeC <- "Lifej[j] = (1/Sj[j])*exp(-(a - (b/Sj[j])))"
    loglifeC <- "Lifej[j] = -log(Sj[j]) - a + (b/Sj[j])"
  }

  if (ls=="Power") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*((Si[i]/Sni[i])^a);"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*((Sj[j]/Snj[j])^a);"

    # Life functions
    lifeF <- "Lifei[i] = b*(Si[i]^a)"
    loglifeF <- "Lifei[i] = log(b) + a*log(Si[i])"
    lifeC <- "Lifej[j] = b*(Sj[j]^a)"
    loglifeC <- "Lifej[j] = log(b) + a*log(Sj[j])"
  }

  if (ls=="InversePower") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*((Sni[i]/Si[i])^a);"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*((Snj[j]/Sj[j])^a);"

    # Life functions
    lifeF <- "Lifei[i] = b *(Si[i]^-a);"
    loglifeF <- "Lifei[i] = log(b) - a*log(Si[i]);"
    lifeC <- "Lifej[j] = b*(Sj[j]^-a);"
    loglifeC <- "Lifej[j] = log(b) - a*log(Sj[j]);"
  }

  if (ls=="Logarithmic") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*((a*log(Si[i]) + b)/(a*log(Sni[i]) + b));"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*((a*log(Sj[j]) + b)/(a*log(Snj[j]) + b));"

    # Life functions
    lifeF <- "Lifei[i] = a*log(Si[i]) + b"
    loglifeF <- "Lifei[i] = log(a*log(Si[i]) + b)"
    lifeC <- "Lifej[j] = a*log(Sj[j]) + b"
    loglifeC <- "Lifej[j] = log(a*log(Sj[j]) + b)"
  }

  if (ls=="MultiStress") {
    # CHECK THIS LAST
    # lsparams - c(ao,a1,...,an), S - c(1,S1,...,Sn)
    lsparams <- "vector a; real b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- exp(theta[ishift+1:length(SF)+ishift+1]%*%c(1,SF))
    function(theta) {
      exp(theta[ishift+1:length(SF)+ishift+1]%*%c(1,SF))
    }

    loglifeF <- function(theta) {
      theta[ishift+1:length(SF)+ishift+1]%*%c(1,SF)
    }
    if(missing(Tc)==FALSE){
      lifeC <- function(theta) {
        exp(theta[ishift+1:length(Sc)+ishift+1]%*%c(1,Sc))
      }
      loglifeC <- function(theta) {
        theta[ishift+1:length(Sc)+ishift+1]%*%c(1,Sc)
      }
    }
  }

  if (ls=="TempHumidity") {
    # lsparams[1] - parameter A, lsparams[2] - parameter a, lsparams[3] - parameter b
    lsparams <- "real A; real a; real b;"
    lsparamsvec <- c("A","a","b")
    pr1<-paste(c("A ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("a ~ ",priors[ishift+2],";"),collapse = "")
    pr3<-paste(c("b ~ ",priors[ishift+3],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2,pr3),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*exp((Ea/8.617385e-5)*((1/Si[i]) - (1/Sni[i])));"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*exp((Ea/8.617385e-5)*((1/Sj[j]) - (1/Snj[j])));"

    # Life functions
    lifeF <- "Lifei[i] = A*exp((a/Sf[,1]) + (b/Sf[,2]))"
    loglifeF <- "Lifei[i] = log(A) + (a/Sf[,1]) + (b/Sf[,2])"
    lifeC <- "Lifej[j] = A*exp((a/Sc[,1]) + (b/Sc[,2]))"
    loglifeC <- "Lifej[j] = log(A) + (a/Sc[,1]) + (b/Sc[,2])"
  }

  if (ls=="TempNonthermal") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
    lsparams <- "real A; real a; real b;"
    lsparamsvec <- c("A","a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    pr3<-paste(c("c ~ ",priors[ishift+3],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2,pr3),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*exp((Ea/8.617385e-5)*((1/Si[i]) - (1/Sni[i])));"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*exp((Ea/8.617385e-5)*((1/Sj[j]) - (1/Snj[j])));"

    # Life functions
    lifeF <- "c/((Sf[,2]^b)*exp(-a/Sf[,1]))"
    loglifeF <- "log(c) - b*log(Sf[,2]) + (a/Sf[,1])"
    lifeC <- "c/((Sc[,2]^b)*exp(-a/Sc[,1]))"
    loglifeC <- "log(c) - b*log(Sc[,2]) + (a/Sc[,1])"
  }

  if (ls=="Eyring3") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # lsparams[3] - parameter c, lsparams[4] - parameter d
    lsparams <- "real a; real b; real c; real d"
    lsparamsvec <- c("a","b","c","d")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    pr3<-paste(c("c ~ ",priors[ishift+3],";"),collapse = "")
    pr4<-paste(c("d ~ ",priors[ishift+4],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2,pr3,pr4),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*exp((Ea/8.617385e-5)*((1/Si[i]) - (1/Sni[i])));"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*exp((Ea/8.617385e-5)*((1/Sj[j]) - (1/Snj[j])));"

    # Life functions
    lifeF <- "(1/Sf[,1])*exp((a + (b/Sf[,1])) + (c + (d/Sf[,1]))*Sf[,2])"
    loglifeF <- "-log(Sf[,1]) + a + (b/Sf[,1]) + (c + (d/Sf[,1]))*Sf[,2]"
    lifeC <- "(1/Sc[,1])*exp((a + (b/Sc[,1])) + (c + (d/Sc[,1]))*Sc[,2])"
    loglifeC <- "-log(Sc[,1]) + a + (b/Sc[,1]) + (c + (d/Sc[,1]))*Sc[,2]"
  }

  # Fit to log-likelihood distributions
  if (dist=="Weibull") {
    distparam <-"real<lower=0> beta;"
    distpriors<-paste(c("beta ~ ",priors[ishift],";"),collapse = "")

    stress2_lpdf <- "real weibull2alp_lpdf(vector y, real beta, real alpha1, real alpha2, real tau){  vector[num_elements(y)] lprob_i;  real lprob;  for (i in 1:num_elements(y)){  lprob_i[i] <- log(beta) - log(alpha1) + (beta - 1)*log((tau/alpha1) + ((y[i] - tau)/alpha2)) - (((tau/alpha1) + ((y[i] - tau)/alpha2))^beta);  }  lprob <- sum(lprob_i);  return lprob;}"
    stressC_lccdf <- "real weibull2alp_lccdf(vector y, real beta, real alpha1, real alpha2, real tau){  vector[num_elements(y)] lprob_i;  real lprob;  for (i in 1:num_elements(y)){    lprob_i[i] <- - ((tau/alpha1) + ((y[i] - tau)/alpha2))^beta;  }  lprob <- sum(lprob_i);  return lprob;}"

    transformedparams <- paste(c(" real<lower=0> alpha1; real<lower=0> alpha2; alpha1 = ",life1F,"; alpha2 = ",life2F,";"),collapse = "")

    timestress1 <- function(theta,S,p) {
      exp((1/theta[1])*log(-log(1-p)) + loglifeF(theta,S))
    }


    loglik <- paste(c("target += weibull_lpdf(Ti | beta, alpha1) + weibull2alp_lpdf(Tj | beta, alpha1, alpha2, tau) + weibull2alp_lccdf(Tk | beta, alpha1, alpha2, tau);"),collapse = "")
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("beta",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")

  }
  if (dist=="Lognormal") {
    distparam <-"real<lower=0> sigma_t;"
    distpriors<-paste(c("sigma_t ~ ",priors[ishift],";"),collapse = "")

    loglik <- paste(c("target += lognormal_lpdf(tiadj | Lifei, sigma_t) + lognormal_lccdf(tjadj | Lifej, sigma_t);"),collapse = "")
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("sigma_t",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")
  }
  if (dist=="Normal") {
    distparam <-"real<lower=0> sigma;"
    distpriors<-paste(c("sigma ~ ",priors[ishift],";"),collapse = "")

    loglik <- paste(c("target += normal_lpdf(tiadj | Lifei, sigma) + normal_lccdf(tjadj | Lifej, sigma);"),collapse = "")
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("sigma",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")
  }
  if (dist=="Exponential") {
    loglik <- paste(c("target += exponential_lpdf(tiadj | 1/Lifei) + exponential_lccdf(tjadj 1/Lifej);"),collapse = "")
    params <- lsparams
    paramsvec <- sparamsvec
    priors <- lspriors
  }
  if (dist=="2PExponential") {
    distparam <-"real<lower=0> sigma;"
    distpriors<-paste(c("sigma ~ ",priors[ishift],";"),collapse = "")

    loglik <- paste(c("target += double_exponential_lpdf(Ti | Lifei, sigma) + double_exponential_lccdf(tjadj | Lifej, sigma);"),collapse = "")
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("sigma",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")
  }

  # Define stancode here
  # Data Block
  block0 <- paste(c("functions {",stress2_lpdf," ",stressC_lccdf,"}"),collapse = " ")
  block1 <- paste(c("data {int<lower=0> n1; int<lower=0> n2; int<lower=0> n3; vector[n1] Ti; vector[n2] Tj; vector[n3] Tk; real tau; real Si; real Sj; real Sk; }"),collapse = " ")
  # Parameter Block
  # block2 <- paste(c("parameters {",params,"}"),collapse = " ")
  block2 <- paste(c("parameters {",params,"} transformed parameters {",transformedparams,"}"),collapse = " ")
  for(j in 1:nruns){
    # Model Block
    # Randomize data for run #j here

    priorno<-gsub(".*\\(", "", gsub("\\)", "", set1[3]))

    n1<-sample.int(n,1)
    n2<-sample.int(n-n1,1)
    n3<-n-n1-n2
    tau<-sample.int(tc,1)
    ttf1<-sample(0:tau,n1)
    ttf2<-sample(tau:tc,n2)
    tC<-rep(tc,n3)

    # if (dist=="Lognormal"){
    #   block3 <- paste(c("model {vector[n] tiadj; vector[m] tjadj; vector[n] Lifei; vector[m] Lifej; ",priors," for(i in 1:n){",Ti_adj,loglifeF,"} for(j in 1:m){",Tj_adj,loglifeC,"}",loglik,"}"),collapse = " ")
    # } else{
    #   block3 <- paste(c("model {vector[n] tiadj; vector[m] tjadj; vector[n] Lifei; vector[m] Lifej; ",priors," for(i in 1:n){",Ti_adj,lifeF,"} for(j in 1:m){",Tj_adj,lifeC,"}",loglik,"}"),collapse = " ")
    # }
    block3 <- paste(c("model {",priors,loglik,"}"),collapse = " ")
    # block3 <- paste(c("model {vector T_1[n1]; vector[n2] T_2;  vector[n3] T_C;",priors," for(i in 1:n1){T_1[i] = Ti[i];} for(j in 1:n2){T_2[j] = Tj[j]; for(k in 1:n3){T_C[k] = Tk[k];}}",loglik,"}"),collapse = " ")
    # Data Input Block
    datablock <- list(n1 = n1, n2 = n2, n3 = n3, Ti = ttf1, Tj = ttf2, Tk = tC, tau = tau, Si = Strset[2], Sj = Strset[3], Sk = Strset[3])
    # NOT RUN {
    stanlscode <- paste(c(block0,block1,block2,block3),collapse=" ")
    # Build or compile Stan code to C++
    lsmod <- stan_model(model_code = stanlscode, verbose = TRUE)
    fit <- sampling(lsmod, data = datablock, iter = nsamples, warmup = burnin, init = pt_est)
    # }
    # Print results.  I need to get this as an output
    stats <- print(fit, pars = paramsvec, probs=c((1-confid)/2,.5,1-(1-confid)/2))
    dataout <- fit@.MISC[["summary"]][["msd"]]
    quanout <- fit@.MISC[["summary"]][["quan"]]
  }

  # Trace the Markov Chains for each parameter
  # plot1_MCtrace <- traceplot(fit, pars = paramsvec, inc_warmup = TRUE, nrow = 3)
  plot1_MCtrace <- mcmc_trace(as.matrix(fit),pars=paramsvec, facet_args = list(nrow = length(paramsvec), labeller = label_parsed))
  plot2_hist <- stan_hist(fit)
  plot3_density <- stan_dens(fit)

  return(list(fit,stats,dataout,quanout,plot1_MCtrace,plot2_hist,plot3_density))
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
}
