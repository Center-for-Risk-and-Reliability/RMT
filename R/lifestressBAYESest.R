# Bayesian Life-Stress Estimator
# Developed by Dr. Reuel Smith, 2021-2023

lifestress.BAYESest <- function(pt_est,ls,dist,TTF,SF,Tc=NULL,Sc=NULL,confid=0.95,priors,nsamples=20000,burnin=10000,nchains=4){
  #Load pracma library for erf
  library(pracma)
  library(StanHeaders)
  library(rstan)
  library(ggplot2)
  library(shinystan)
  library(cmdstanr)
  library(bayesplot)

  # (UPDATE 11/14/2023): Adding an IPL form of l = 1/(b x S^a) and Exponential form b x exp(a/S)

  # Add input to this to include prior estimates for LS parameters.
  # Example: priors<-c("normal(3,4)","normal(1,4)", "lognormal(-2,3)")
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
  conf.level <- confid
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
    if(is.null(Tc)==FALSE){
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
    if(is.null(Tc)==FALSE){
      lifeC <- "b*exp(a*Sc)"
      loglifeC <- "log(b) + a*Sc"
    }
  }

  if (ls=="Exponential2"){
    # theta[1] - parameter a, theta[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "b*exp(a/Sf)"
    loglifeF <- "log(b) + a/Sf"
    if(is.null(Tc)==FALSE){
      lifeC <- "b*exp(a/Sc)"
      loglifeC <- "log(b) + a/Sc"
    }
  }

  if (ls=="Arrhenius") {
    # lsparams[1] - parameter Ea, lsparams[2] - parameter b
    # Temperature HAS to be in Kelvin for this to work
    lsparams <- "real Ea; real<lower=0> b;"
    lsparamsvec <- c("Ea","b")
    pr1<-paste(c("Ea ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "b*exp(Ea/((8.617385e-5)*Sf))"
    loglifeF <- "log(b) + (Ea/((8.617385e-5)*Sf))"
    if(is.null(Tc)==FALSE){
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

    # NOTE (11/15/2023): Eyring life-stress model needs to be entered as a vector because it multiplies two Sf or Sc vectors and
    # this is how Stan typically treats vectors multiplied or divided from other vectors.
    lifeF <- "Lifei[i] = (b/Sf[i])*exp(a/Sf[i]);"
    loglifeF <- "Lifei[i] = log(b) - log(Sf[i]) + (a/Sf[i]);"
    if(is.null(Tc)==FALSE){
      lifeC <- "Lifej[j] = (b/Sc[j])*exp(a/Sc[j]);"
      loglifeC <- "Lifej[j] = log(b) - log(Sc[j]) + (a/Sc[j]);"
    }
  }

  if (ls=="Eyring2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "Lifei[i] = (1/Sf[i])*exp(-(a - (b/Sf[i])));"
    loglifeF <- "Lifei[i] = -log(Sf[i]) - a + (b/Sf[i]);"

    if(is.null(Tc)==FALSE){
      lifeC <- "Lifej[j] = (1/Sc[j])*exp(-(a - (b/Sc[j])));"
      loglifeC <- "Lifej[j] = -log(Sc[j]) - a + (b/Sc[j]);"
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
    if(is.null(Tc)==FALSE){
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
    if(is.null(Tc)==FALSE){
      lifeC <- "b*(Sc^-a)"
      loglifeC <- "log(b) - a*log(Sc)"
    }
  }

  if (ls=="InversePower2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "1/(b*(Sf^a))"
    loglifeF <- "-log(b) - a*log(Sf)"
    if(is.null(Tc)==FALSE){
      lifeC <- "1/(b*(Sc^a))"
      loglifeC <- "-log(b) - a*log(Sc)"
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
    if(is.null(Tc)==FALSE){
      lifeC <- "a*log(Sc) + b"
      loglifeC <- "log(a*log(Sc) + b)"
    }
  }

  if (ls=="MultiStress") {
    # CHECK THIS LAST
    # lsparams - c(ao,a1,...,an), S - c(1,S1,...,Sn)
    if((length(priors)-ishift)==2){
      lsparams <- "real a0; real a1; "
      lsparamsvec <- c("a0","a1")
      pr1<-paste(c("a0 ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("a1 ~ ",priors[ishift+2],";"),collapse = "")
      lspriors <- paste(c(pr1,pr2),collapse = " ")

      lifeF <- "Lifei[i] = exp(a0 + a1*Sf[i]);"
      loglifeF <- "Lifei[i] = a0 + a1*Sf[i];"
      if(is.null(Tc)==FALSE){
        lifeC <- "Lifej[j] = exp(a0 + a1*Sc[j]);"
        loglifeC <- "Lifej[j] = a0 + a1*Sc[j];"
      }
    }
    if((length(priors)-ishift)==3){
      lsparams <- "real a0; real a1; real a2; "
      lsparamsvec <- c("a0","a1","a2")
      pr1<-paste(c("a0 ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("a1 ~ ",priors[ishift+2],";"),collapse = "")
      pr3<-paste(c("a2 ~ ",priors[ishift+3],";"),collapse = "")
      lspriors <- paste(c(pr1,pr2,pr3),collapse = " ")

      lifeF <- "Lifei[i] = exp(a0 + a1*Sf[i,1] + a2*Sf[i,2]);"
      loglifeF <- "Lifei[i] = a0 + a1*Sf[i,1] + a2*Sf[i,2];"
      if(is.null(Tc)==FALSE){
        lifeC <- "Lifej[j] = exp(a0 + a1*Sc[j,1] + a2*Sc[j,2]);"
        loglifeC <- "Lifej[j] = a0 + a1*Sc[j,1] + a2*Sc[j,2];"
      }
    }
    if((length(priors)-ishift)==4){
      lsparams <- "real a0; real a1; real a2; real a3;"
      lsparamsvec <- c("a0","a1","a2","a3")
      pr1<-paste(c("a0 ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("a1 ~ ",priors[ishift+2],";"),collapse = "")
      pr3<-paste(c("a2 ~ ",priors[ishift+3],";"),collapse = "")
      pr4<-paste(c("a3 ~ ",priors[ishift+4],";"),collapse = "")
      lspriors <- paste(c(pr1,pr2,pr3,pr4),collapse = " ")

      lifeF <- "Lifei[i] = exp(a0 + a1*Sf[i,1] + a2*Sf[i,2] + a3*Sf[i,3]);"
      loglifeF <- "Lifei[i] = a0 + a1*Sf[i,1] + a2*Sf[i,2] + a3*Sf[i,3];"
      if(is.null(Tc)==FALSE){
        lifeC <- "Lifej[j] = exp(a0 + a1*Sc[j,1] + a2*Sc[j,2] + a3*Sc[j,3]);"
        loglifeC <- "Lifej[j] = a0 + a1*Sc[j,1] + a2*Sc[j,2] + a3*Sc[j,3];"
      }
    }
    if((length(priors)-ishift)==5){
      lsparams <- "real a0; real a1; real a2; real a3; real a4;"
      lsparamsvec <- c("a0","a1","a2","a3")
      pr1<-paste(c("a0 ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("a1 ~ ",priors[ishift+2],";"),collapse = "")
      pr3<-paste(c("a2 ~ ",priors[ishift+3],";"),collapse = "")
      pr4<-paste(c("a3 ~ ",priors[ishift+4],";"),collapse = "")
      pr5<-paste(c("a4 ~ ",priors[ishift+5],";"),collapse = "")
      lspriors <- paste(c(pr1,pr2,pr3,pr4,pr5),collapse = " ")

      lifeF <- "Lifei[i] = exp(a0 + a1*Sf[i,1] + a2*Sf[i,2] + a3*Sf[i,3] + a4*Sf[i,4]);"
      loglifeF <- "Lifei[i] = a0 + a1*Sf[i,1] + a2*Sf[i,2] + a3*Sf[i,3] + a4*Sf[i,4];"
      if(is.null(Tc)==FALSE){
        lifeC <- "Lifej[j] = exp(a0 + a1*Sc[j,1] + a2*Sc[j,2] + a3*Sc[j,3] + a4*Sc[j,4]);"
        loglifeC <- "Lifej[j] = a0 + a1*Sc[j,1] + a2*Sc[j,2] + a3*Sc[j,3] + a4*Sc[j,4];"
      }
    }
  }

  if (ls=="TempHumidity") {
    # lsparams[1] - parameter A, lsparams[2] - parameter a, lsparams[3] - parameter b
    lsparams <- "real<lower=0> A; real a; real b;"
    lsparamsvec <- c("A","a","b")
    pr1<-paste(c("A ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("a ~ ",priors[ishift+2],";"),collapse = "")
    pr3<-paste(c("b ~ ",priors[ishift+3],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2,pr3),collapse = " ")

    lifeF <- "Lifei[i] = A*exp((a/Sf[i,1]) + (b/Sf[i,2]));"
    loglifeF <- "Lifei[i] = log(A) + (a/Sf[i,1]) + (b/Sf[i,2]);"
    if(is.null(Tc)==FALSE){
      lifeC <- "Lifej[j] = A*exp((a/Sc[j,1]) + (b/Sc[j,2]));"
      loglifeC <- "Lifej[j] = log(A) + (a/Sc[j,1]) + (b/Sc[j,2]);"
    }
  }

  if (ls=="TempNonthermal") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
    lsparams <- "real a; real b; real<lower=0> c;"
    lsparamsvec <- c("a","b","c")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    pr3<-paste(c("c ~ ",priors[ishift+3],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2,pr3),collapse = " ")

    lifeF <- "Lifei[i] = c/((Sf[i,2]^b)*exp(-a/Sf[i,1]));"
    loglifeF <- "Lifei[i] = log(c) - b*log(Sf[i,2]) + (a/Sf[i,1]);"
    if(is.null(Tc)==FALSE){
      lifeC <- "Lifej[j] = c/((Sc[j,2]^b)*exp(-a/Sc[j,1]));"
      loglifeC <- "Lifej[j] = log(c) - b*log(Sc[j,2]) + (a/Sc[j,1]);"
    }
  }

  if (ls=="Eyring3") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # lsparams[3] - parameter c, lsparams[4] - parameter d
    lsparams <- "real a; real b; real c; real d;"
    lsparamsvec <- c("a","b","c","d")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    pr3<-paste(c("c ~ ",priors[ishift+3],";"),collapse = "")
    pr4<-paste(c("d ~ ",priors[ishift+4],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2,pr3,pr4),collapse = " ")

    lifeF <- "Lifei[i] = (1/Sf[i,1])*exp((a + (b/Sf[i,1])) + (c + (d/Sf[i,1]))*Sf[i,2]);"
    loglifeF <- "Lifei[i] = -log(Sf[i,1]) + a + (b/Sf[i,1]) + (c + (d/Sf[i,1]))*Sf[i,2];"
    if(is.null(Tc)==FALSE){
      lifeC <- "Lifej[j] = (1/Sc[j,1])*exp((a + (b/Sc[j,1])) + (c + (d/Sc[j,1]))*Sc[j,2]);"
      loglifeC <- "Lifej[j] = -log(Sc[j,1]) + a + (b/Sc[j,1]) + (c + (d/Sc[j,1]))*Sc[j,2];"
    }
  }

  # Fit to log-likelihood distributions
  if (dist=="Weibull") {
    distparam <-"real<lower=0> beta;"
    distpriors<-paste(c("beta ~ ",priors[ishift],";"),collapse = "")

    if(ls=="Eyring" || ls=="Eyring2" || ls=="Eyring3" || ls=="TempHumidity" || ls=="TempNonthermal" || ls=="MultiStress"){
      if(missing(Tc)){
        loglik <- paste(c("target += weibull_lpdf(TTF | beta, Lifei);"),collapse = "")
      } else{
        loglik <- paste(c("target += weibull_lpdf(TTF | beta, Lifei) + weibull_lccdf(TTS | beta, Lifej);"),collapse = "")
      }
    } else{
      if(missing(Tc)){
        loglik <- paste(c("target += weibull_lpdf(TTF | beta,",lifeF,");"),collapse = "")
      } else{
        loglik <- paste(c("target += weibull_lpdf(TTF | beta,",lifeF,") + weibull_lccdf(TTS | beta,",lifeC,");"),collapse = "")
      }
    }
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("beta",lsparamsvec)
    outputparamset <- c("\U03B2",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")
  }
  if (dist=="Lognormal") {
    distparam <-"real<lower=0> sigma_t;"
    distpriors<-paste(c("sigma_t ~ ",priors[ishift],";"),collapse = "")

    if(ls=="Eyring" || ls=="Eyring2" || ls=="Eyring3" || ls=="TempHumidity" || ls=="TempNonthermal" || ls=="MultiStress"){
      if(missing(Tc)){
        loglik <- paste(c("target += lognormal_lpdf(TTF |Lifei, sigma_t);"),collapse = "")
      } else{
        loglik <- paste(c("target += lognormal_lpdf(TTF |Lifei, sigma_t) + lognormal_lccdf(TTS |Lifej, sigma_t);"),collapse = "")
      }
    } else{
      if(missing(Tc)){
        loglik <- paste(c("target += lognormal_lpdf(TTF |",loglifeF,", sigma_t);"),collapse = "")
      } else{
        loglik <- paste(c("target += lognormal_lpdf(TTF |",loglifeF,", sigma_t) + lognormal_lccdf(TTS |",loglifeC,", sigma_t);"),collapse = "")
      }
    }

    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("sigma_t",lsparamsvec)
    outputparamset <- c("\U03C3_t",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")
  }
  if (dist=="Normal") {
    distparam <-"real<lower=0> sigma;"
    distpriors<-paste(c("sigma ~ ",priors[ishift],";"),collapse = "")
    if(ls=="Eyring" || ls=="Eyring2" || ls=="Eyring3" || ls=="TempHumidity" || ls=="TempNonthermal" || ls=="MultiStress"){
      if(missing(Tc)){
        loglik <- paste(c("target += normal_lpdf(TTF | Lifei, sigma);"),collapse = "")
      } else{
        loglik <- paste(c("target += normal_lpdf(TTF | Lifei, sigma) + normal_lccdf(TTS | Lifej, sigma);"),collapse = "")
      }
    } else{
      if(missing(Tc)){
        loglik <- paste(c("target += normal_lpdf(TTF |",lifeF,", sigma);"),collapse = "")
      } else{
        loglik <- paste(c("target += normal_lpdf(TTF |",lifeF,", sigma) + normal_lccdf(TTS |",lifeC,", sigma);"),collapse = "")
      }
    }
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("sigma",lsparamsvec)
    outputparamset <- c("\U03C3",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")
  }
  if (dist=="Exponential") {
    if(ls=="Eyring" || ls=="Eyring2" || ls=="Eyring3" || ls=="TempHumidity" || ls=="TempNonthermal" || ls=="MultiStress"){
      if(missing(Tc)){
        loglik <- paste(c("target += exponential_lpdf(TTF | 1/(Lifei));"),collapse = "")
      } else{
        loglik <- paste(c("target += exponential_lpdf(TTF | 1/(Lifei)) + exponential_lccdf(TTS 1/(Lifej));"),collapse = "")
      }
    } else{
      if(missing(Tc)){
        loglik <- paste(c("target += exponential_lpdf(TTF | 1/(",lifeF,"));"),collapse = "")
      } else{
        loglik <- paste(c("target += exponential_lpdf(TTF | 1/(",lifeF,")) + exponential_lccdf(TTS 1/(",lifeC,"));"),collapse = "")
      }
    }
    params <- lsparams
    paramsvec <- lsparamsvec
    outputparamset <- lsparamsvec
    priors <- lspriors
  }
  if (dist=="2PExponential") {
    distparam <-"real<lower=0> sigma;"
    distpriors<-paste(c("sigma ~ ",priors[ishift],";"),collapse = "")
    if(ls=="Eyring" || ls=="Eyring2" || ls=="Eyring3" || ls=="TempHumidity" || ls=="TempNonthermal" || ls=="MultiStress"){
      if(missing(Tc)){
        loglik <- paste(c("target += double_exponential_lpdf(TTF | Lifei, sigma);"),collapse = "")
      } else{
        loglik <- paste(c("target += double_exponential_lpdf(TTF | Lifei, sigma) + double_exponential_lccdf(TTS | Lifej, sigma);"),collapse = "")
      }
    } else{
      if(missing(Tc)){
        loglik <- paste(c("target += double_exponential_lpdf(TTF |",lifeF,", sigma);"),collapse = "")
      } else{
        loglik <- paste(c("target += double_exponential_lpdf(TTF |",lifeF,", sigma) + double_exponential_lccdf(TTS |",lifeC,", sigma);"),collapse = "")
      }
    }
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("sigma",lsparamsvec)
    outputparamset <- c("\U03C3",lsparamsvec)
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
  if ((ls=="Eyring" || ls=="Eyring2"  || ls=="Eyring3" || ls=="TempHumidity" || ls=="TempNonthermal" || ls=="MultiStress") && dist == "Lognormal"){
    if(missing(Tc)==TRUE){
      block3 <- paste(c("model { vector[n] Lifei; vector[m] Lifej; ",priors," for(i in 1:n){",loglifeF,"}",loglik,"}"),collapse = " ")
    }
    if(missing(Tc)==FALSE){
      block3 <- paste(c("model { vector[n] Lifei; vector[m] Lifej; ",priors," for(i in 1:n){",loglifeF,"} for(j in 1:m){",loglifeC,"}",loglik,"}"),collapse = " ")
    }
  }
  if ((ls=="Eyring" || ls=="Eyring2"  || ls=="Eyring3" || ls=="TempHumidity" || ls=="TempNonthermal" || ls=="MultiStress") && (dist == "Normal" || dist=="Weibull" || dist=="Exponential")){
    if(missing(Tc)==TRUE){
      block3 <- paste(c("model { vector[n] Lifei; ",priors," for(i in 1:n){",lifeF,"}",loglik,"}"),collapse = " ")
    }
    if(missing(Tc)==FALSE){
      block3 <- paste(c("model { vector[n] Lifei; vector[m] Lifej; ",priors," for(i in 1:n){",lifeF,"} for(j in 1:m){",lifeC,"}",loglik,"}"),collapse = " ")
    }
  }
  # NOT RUN {
  stanlscode <- paste(c(block1,block2,block3),collapse=" ")
  stanlsfile <- write_stan_file(stanlscode)
  print(stanlsfile)
  # Generate initial list (one list per chain)
  names(pt_est) <- paramsvec
  pt_estlist <- as.list(pt_est)
  init_pt_est <- vector("list",nchains)
  for(i in 1:nchains){
    init_pt_est[[i]] <- pt_estlist
  }
  # Build or compile Stan code to C++
  # return(list(stanlscode,stanlsfile))

  # lsmod <- stan_model(model_code = stanlscode, verbose = TRUE)
  lsmod <- cmdstan_model(stanlsfile)
  # return(lsmod)
  # fit <- sampling(lsmod, data = datablock, iter = nsamples, warmup = burnin, init = pt_est)
  fit <- lsmod$sample(data = datablock, init = init_pt_est, chains = nchains, iter_warmup = burnin, iter_sampling = nsamples)
  # }
  # return(fit)
  # Print results.  I need to get this as an output
  # stats <- print(fit, pars = paramsvec, probs=c((1-confid)/2,.5,1-(1-confid)/2))
  # dataout <- fit@.MISC[["summary"]][["msd"]]
  conflim_txt<-c(paste(c("Lower ",100*conf.level,"%"),collapse = ""),paste(c("Upper ",100*conf.level,"%"),collapse = ""))
  stats <- fit$summary(variables = paramsvec)
  dataout <- fit$draws(format = "df")
  confidbounds <- mcmc_intervals_data(fit$draws(variables = paramsvec),prob_outer = confid)
  outputtable <- matrix(c(stats[[2]],stats[[4]],confidbounds[[5]],stats[[3]],confidbounds[[9]],stats[[8]]), nrow = length(outputparamset), ncol = 6, byrow = FALSE,dimnames = list(outputparamset,c("Mean","Standard Deviation",conflim_txt[1],"Median",conflim_txt[2],"R\U005E")))


  # Trace the Markov Chains for each parameter
  # plot1_MCtrace <- traceplot(fit, pars = paramsvec, inc_warmup = TRUE, nrow = 3)
  # plot1_MCtrace <- mcmc_trace(as.matrix(fit),pars=paramsvec, facet_args = list(nrow = length(paramsvec), labeller = label_parsed))
  # plot2_hist <- stan_hist(fit)
  # plot3_density <- stan_dens(fit)
  plot1_MCtrace <- mcmc_trace(fit$draws(paramsvec))
  plot2_hist <- mcmc_hist(fit$draws(paramsvec))
  plot3_density <- mcmc_dens(fit$draws(paramsvec))
  plot4_densityoverlay <- mcmc_dens_overlay(fit$draws(paramsvec))

  # Produce some output text that summarizes the results
  cat(c("Posterior estimates for Bayesian Analysis.\n\n"),sep = "")
  print(outputtable)
  cat(c("\n"),sep = "")


  return(list(fit,stats,dataout,plot1_MCtrace,plot2_hist,plot3_density,plot4_densityoverlay))
}
