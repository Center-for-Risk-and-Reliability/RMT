# Bayesian Life-Stress Estimator
# Developed by Dr. Reuel Smith, 2021-2022

lifestress.BAYESest <- function(pt_est,ls,dist,TTF,SF,Tc,Sc,confid,priors,nsamples,burnin){
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

    lifeF <- "b + Sf*a"
    loglifeF <- "log(b + Sf*a)"
    if(missing(Tc)==FALSE){
      lifeC <- "b + Sc*a"
      loglifeC <- "log(b + Sc*a)"
    }
  }

  if (ls=="Exponential"){
    # theta[1] - parameter a, theta[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "b*exp(a*Sf)"
    loglifeF <- "log(b) + a*Sf"
    if(missing(Tc)==FALSE){
      lifeC <- "b*exp(a*Sc)"
      loglifeC <- "log(b) + a*Sc"
    }
  }

  if (ls=="Arrhenius") {
    # lsparams[1] - parameter Ea, lsparams[2] - parameter b
    # Temperature HAS to be in Kelvin for this to work
    lsparams <- "real<lower=0> Ea; real<lower=0> b;"
    lsparamsvec <- c("Ea","b")
    pr1<-paste(c("Ea ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "b*exp(Ea/((8.617385e-5)*Sf))"
    loglifeF <- "log(b) + (Ea/((8.617385e-5)*Sf))"
    if(missing(Tc)==FALSE){
      lifeC <- "b*exp(Ea/((8.617385e-5)*Sc))"
      loglifeC <- "log(b) + (Ea/((8.617385e-5)*Sc))"
    }
  }

  if (ls=="Eyring") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "(b/Sf)*exp(a/Sf)"
    loglifeF <- "log(b) - log(Sf) + (a/Sf)"
    if(missing(Tc)==FALSE){
      lifeC <- "(b/Sc)*exp(a/Sc)"
      loglifeC <- "log(b) - log(Sc) + (a/Sc)"
    }
  }

  if (ls=="Eyring2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "(1/Sf)*exp(-(a - (b/Sf)))"
    loglifeF <- "-log(Sf) - a + (b/Sf)"

    if(missing(Tc)==FALSE){
      lifeC <- "(1/Sc)*exp(-(a - (b/Sc)))"
      loglifeC <- "-log(Sc) - a + (b/Sc)"
    }
  }

  if (ls=="Power") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "b*(Sf^a)"
    loglifeF <- "log(b) + a*log(Sf)"
    if(missing(Tc)==FALSE){
      lifeC <- "b*(Sc^a)"
      loglifeC <- "log(b) + a*log(Sc)"
    }
  }

  if (ls=="InversePower") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "b*(Sf^-a)"
    loglifeF <- "log(b) - a*log(Sf)"
    if(missing(Tc)==FALSE){
      lifeC <- "b*(Sc^-a)"
      loglifeC <- "log(b) - a*log(Sc)"
    }
  }

  if (ls=="Logarithmic") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "a*log(Sf) + b"
    loglifeF <- "log(a*log(Sf) + b)"
    if(missing(Tc)==FALSE){
      lifeC <- "a*log(Sc) + b"
      loglifeC <- "log(a*log(Sc) + b)"
    }
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

    lifeF <- "A*exp((a/Sf[,1]) + (b/Sf[,2]))"
    loglifeF <- "log(A) + (a/Sf[,1]) + (b/Sf[,2])"
    if(missing(Tc)==FALSE){
      lifeC <- "A*exp((a/Sc[,1]) + (b/Sc[,2]))"
      loglifeC <- "log(A) + (a/Sc[,1]) + (b/Sc[,2])"
    }
  }

  if (ls=="TempNonthermal") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
    lsparams <- "real A; real a; real b;"
    lsparamsvec <- c("A","a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    pr3<-paste(c("c ~ ",priors[ishift+3],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2,pr3),collapse = " ")

    lifeF <- "c/((Sf[,2]^b)*exp(-a/Sf[,1]))"
    loglifeF <- "log(c) - b*log(Sf[,2]) + (a/Sf[,1])"
    if(missing(Tc)==FALSE){
      lifeC <- "c/((Sc[,2]^b)*exp(-a/Sc[,1]))"
      loglifeC <- "log(c) - b*log(Sc[,2]) + (a/Sc[,1])"
    }
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

    lifeF <- "(1/Sf[,1])*exp((a + (b/Sf[,1])) + (c + (d/Sf[,1]))*Sf[,2])"
    loglifeF <- "-log(Sf[,1]) + a + (b/Sf[,1]) + (c + (d/Sf[,1]))*Sf[,2]"
    if(missing(Tc)==FALSE){
      lifeC <- "(1/Sc[,1])*exp((a + (b/Sc[,1])) + (c + (d/Sc[,1]))*Sc[,2])"
      loglifeC <- "-log(Sc[,1]) + a + (b/Sc[,1]) + (c + (d/Sc[,1]))*Sc[,2]"
    }
  }

  # Fit to log-likelihood distributions
  if (dist=="Weibull") {
    distparam <-"real<lower=0> beta;"
    distpriors<-paste(c("beta ~ ",priors[ishift],";"),collapse = "")

    if(missing(Tc)){
      loglik <- paste(c("target += weibull_lpdf(TTF | beta,",lifeF,");"),collapse = "")
    } else{
      loglik <- paste(c("target += weibull_lpdf(TTF | beta",lifeF,") + weibull_lccdf(TTS | beta",lifeC,");"),collapse = "")
    }
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("beta",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")
  }
  if (dist=="Lognormal") {
    distparam <-"real<lower=0> sigma_t;"
    distpriors<-paste(c("sigma_t ~ ",priors[ishift],";"),collapse = "")

    if(missing(Tc)){
      loglik <- paste(c("target += lognormal_lpdf(TTF |",loglifeF,", sigma_t);"),collapse = "")
    } else{
      loglik <- paste(c("target += lognormal_lpdf(TTF |",loglifeF,", sigma_t) + lognormal_lccdf(TTS |",loglifeC,", sigma_t);"),collapse = "")
    }
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("sigma_t",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")
  }
  if (dist=="Normal") {
    distparam <-"real<lower=0> sigma;"
    distpriors<-paste(c("sigma ~ ",priors[ishift],";"),collapse = "")
    if(missing(Tc)){
      loglik <- paste(c("target += normal_lpdf(TTF |",lifeF,", sigma);"),collapse = "")
    } else{
      loglik <- paste(c("target += normal_lpdf(TTF |",lifeF,", sigma) + normal_lccdf(TTS |",lifeC,", sigma);"),collapse = "")
    }
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("sigma",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")
  }
  if (dist=="Exponential") {
    if(missing(Tc)){
      loglik <- paste(c("target += exponential_lpdf(TTF | 1/(",lifeF,"));"),collapse = "")
    } else{
      loglik <- paste(c("target += exponential_lpdf(TTF | 1/(",lifeF,")) + exponential_lccdf(TTS 1/(",lifeC,"));"),collapse = "")
    }
    params <- lsparams
    paramsvec <- sparamsvec
    priors <- lspriors
  }
  if (dist=="2PExponential") {
    distparam <-"real<lower=0> sigma;"
    distpriors<-paste(c("sigma ~ ",priors[ishift],";"),collapse = "")
    if(missing(Tc)){
      loglik <- paste(c("target += double_exponential_lpdf(TTF |",lifeF,", sigma);"),collapse = "")
    } else{
      loglik <- paste(c("target += double_exponential_lpdf(TTF |",lifeF,", sigma) + double_exponential_lccdf(TTS |",lifeC,", sigma);"),collapse = "")
    }
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("sigma",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")
  }

  # Define stancode here
  if(missing(Tc)){
    block1 <- "data {int<lower=0> n; vector[n] TTF; vector[n] Sf;}"
    datablock <- list(n = length(TTF), TTF = TTF, Sf = SF)
  } else{
    block1 <- "data {int<lower=0> n; int<lower=0> m; vector[n] TTF; vector[m] TTS; vector[n] Sf; vector[m] Sc;}"
    datablock <- list(n = length(TTF), m = length(Tc), TTF = TTF, Sf = SF, TTS = Tc, Sc = Sc)
  }
  block2 <- paste(c("parameters {",params,"}"),collapse = " ")
  block3 <- paste(c("model {",priors,loglik,"}"),collapse = " ")
  # NOT RUN {
  stanlscode <- paste(c(block1,block2,block3),collapse=" ")
  # Build or compile Stan code to C++
  lsmod <- stan_model(model_code = stanlscode, verbose = TRUE)
  fit <- sampling(lsmod, data = datablock, iter = nsamples, warmup = burnin, init = pt_est)
  # }
  # Print results.  I need to get this as an output
  stats <- print(fit, pars = paramsvec, probs=c((1-confid)/2,.5,1-(1-confid)/2))
  dataout <- fit@.MISC[["summary"]][["msd"]]
  quanout <- fit@.MISC[["summary"]][["quan"]]

  # Trace the Markov Chains for each parameter
  # plot1_MCtrace <- traceplot(fit, pars = paramsvec, inc_warmup = TRUE, nrow = 3)
  plot1_MCtrace <- mcmc_trace(as.matrix(fit),pars=paramsvec, facet_args = list(nrow = length(paramsvec), labeller = label_parsed))
  plot2_hist <- stan_hist(fit)
  plot3_density <- stan_dens(fit)

  return(list(fit,stats,dataout,quanout,plot1_MCtrace,plot2_hist,plot3_density))
}
