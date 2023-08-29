# Bayesian Accelerated Degradation Testing Estimator
# Developed by Dr. Reuel Smith, 2022

adt.full.BAYES <- function(data,lifedam,dist,D0,Tuse,confid,priors,nsamples,burnin){
  # (pt_est,ls,dist,TTF,SF,Tc,Sc,confid,priors,nsamples,burnin)
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

  # Check to see if confidence exists
  if(missing(confid)){
    conf.level <- 0.95
  } else {
    conf.level <- confid
  }

  # Check to see if burn-in exists
  if(missing(burnin)){
    burnin <- floor(nsamples/2)
  } else {
    burnin <- burnin
  }

  # Check the condition if the Hamada damage model is chosen and whether or not
  # a use temperature has been identified or not.  If not then room temperature
  # 293.15 K will be used by default.
  if(lifedam=="Hamada" & missing(Tuse)){
    Tuse <- 293.15
  }

  # Start with the pre-processing of the data in which you obtain an initial parameter
  # estimate based on the curve fit of the data
  if(lifedam=="Hamada"){
    TempDam <- data[,4]
    adtLSQ<-adt.full.LSQ(data,lifedam,D0,Tuse)[[1]]
  } else {
    adtLSQ<-adt.full.LSQ(data,lifedam,D0)[[1]]
  }

  # Check to see if estimate exists
  pt_est <- c(colMeans(adtLSQ)[1:(dim(adtLSQ)[2]-2)],1)

  # Pull necessary data damage time tDam and damage level Dam
  tDam <- data[,1]
  Dam <- data[,2]

  # Initialize life-stress parameter estimates for theta
  if (lifedam == "Linear") {
    # D = a + b*t
    # lifedamparams[1] - parameter a, lifedamparams[2] - parameter b
    lifedamparams <- "real a; real b;"
    lifedamparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[2],";"),collapse = "")
    lifedampriors <- paste(c(pr1,pr2),collapse = " ")

    DamT <- "Damagei[i] = b + t[i]*a"
    logDamT <- "Damagei[i] = log(b + t[i]*a)"
    sigparamno<-3
  }

  if (lifedam == "Exponential"){
    # D = b*exp(a*t)
    # lifedamparams[1] - parameter a, lifedamparams[2] - parameter b
    lifedamparams <- "real a; real<lower=0> b;"
    lifedamparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[2],";"),collapse = "")
    lifedampriors <- paste(c(pr1,pr2),collapse = " ")

    DamT <- "Damagei[i] = b*exp(a*t[i])"
    logDamT <- "Damagei[i] = log(b) + a*t[i]"
    sigparamno<-3
  }

  if(lifedam=="SquareRoot"){
    # D^(1/2) = a + b*t
    # theta[1] ~ parameter a, theta[2] ~ parameter b
    lifedamparams <- "real a; real b;"
    lifedamparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[2],";"),collapse = "")
    lifedampriors <- paste(c(pr1,pr2),collapse = " ")

    DamT <- "Damagei[i] = (a + b*t[i])^2"
    logDamT <- "Damagei[i] = 2*log(a + b*t[i])"
    sigparamno<-3
  }

  if(lifedam=="Power"){
    # D = b*(t^a)
    # lifedamparams[1] - parameter a, lifedamparams[2] - parameter b
    lifedamparams <- "real a; real<lower=0> b;"
    lifedamparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[2],";"),collapse = "")
    lifedampriors <- paste(c(pr1,pr2),collapse = " ")

    DamT <- "Damagei[i] = b*(t[i]^a);"
    logDamT <- "Damagei[i] = log(b) + a*log(t[i]);"
    sigparamno<-3
  }

  if (lifedam=="Logarithmic") {
    # D = a + b*ln(t)
    # lifedamparams[1] - parameter a, lifedamparams[2] - parameter b
    lifedamparams <- "real a; real b;"
    lifedamparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[2],";"),collapse = "")
    lifedampriors <- paste(c(pr1,pr2),collapse = " ")

    DamT <- "Damagei[i] = a*log(t[i]) + b;"
    logDamT <- "Damagei[i] = log(a*log(t[i]) + b);"
    sigparamno<-3
  }

  if (lifedam=="Gompertz") {
    # D = a + b^(c*t)
    # lifedamparams[1] ~ parameter a, lifedamparams[2] ~ parameter b, lifedamparams[3] ~ parameter c
    lifedamparams <- "real a; real b; real c"
    lifedamparamsvec <- c("a","b","c")
    pr1<-paste(c("a ~ ",priors[1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[2],";"),collapse = "")
    pr3<-paste(c("c ~ ",priors[3],";"),collapse = "")
    lifedampriors <- paste(c(pr1,pr2,pr3),collapse = " ")

    DamT <- "Damagei[i] = a + b^(c*t[i]);"
    logDamT <- "Damagei[i] = log(a + b^(c*t[i]));"
    sigparamno<-4
  }

  if (lifedam=="LloydLipow") {
    # D = a - b/t
    # lifedamparams[1] - parameter a, lifedamparams[2] - parameter b
    lifedamparams <- "real a; real b;"
    lifedamparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[2],";"),collapse = "")
    lifedampriors <- paste(c(pr1,pr2),collapse = " ")

    DamT <- "Damagei[i] = a - b/t[i];"
    logDamT <- "Damagei[i] = log(a - b/t[i]);"
    sigparamno<-3
  }

  if (lifedam=="Mitsuom") {
    # D = 1/(1 + b*(t^a))
    # lifedamparams[1] - parameter a, lifedamparams[2] - parameter b
    lifedamparams <- "real a; real b;"
    lifedamparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[2],";"),collapse = "")
    lifedampriors <- paste(c(pr1,pr2),collapse = " ")

    DamT <- "Damagei[i] = 1/(1 + b*(t[i]^a));"
    logDamT <- "Damagei[i] = -log(1 + b*(t[i]^a));"
    sigparamno<-3
  }

  if (lifedam=="Hamada") {
    # D = 1/(1 + beta1*(t*exp(beta3*11605*(1/Tu - 1/Ti)))^beta2)
    # lifedamparams[1] ~ parameter beta1, lifedamparams[2] ~ parameter beta2, lifedamparams[3] ~ parameter beta3
    lifedamparams <- "real<lower=0> beta1; real<lower=0> beta2; real beta3; "
    lifedamparamsvec <- c("beta1","beta2","beta3")
    pr1<-paste(c("beta1 ~ ",priors[1],";"),collapse = "")
    pr2<-paste(c("beta2 ~ ",priors[2],";"),collapse = "")
    pr3<-paste(c("beta3 ~ ",priors[3],";"),collapse = "")
    lifedampriors <- paste(c(pr1,pr2,pr3),collapse = " ")

    DamT <- "Damagei[i] = 1/(1 + beta1*(t[i]*exp(beta3*11605*(1/Tu - 1/Ti[i])))^beta2); "
    logDamT <- "Damagei[i] = -log(1 + beta1*(t[i]*exp(beta3*11605*(1/Tu - 1/Ti[i])))^beta2); "

    Tempf <- data[,4]
    sigparamno<-4
  }

  # Fit to log-likelihood distributions
  if (dist=="Lognormal") {
    distparam <-"real<lower=0> sigma_t;"
    distpriors<-paste(c("sigma_t ~ ",priors[sigparamno],";"),collapse = "")
    loglik <- paste(c("target += lognormal_lpdf(Dt | Damagei, sigma_t);"),collapse = "")
    params <- paste(c(distparam,lifedamparams),collapse = " ")
    paramsvec <- c("sigma_t",lifedamparamsvec)
    priors <- paste(c(distpriors,lifedampriors),collapse = " ")
  }

  if (dist=="Normal") {
    distparam <-"real<lower=0> sigma;"
    distpriors<-paste(c("sigma ~ ",priors[sigparamno],";"),collapse = "")
    loglik <- paste(c("target += normal_lpdf(Dt | Damagei, sigma);"),collapse = "")
    params <- paste(c(distparam,lifedamparams),collapse = " ")
    paramsvec <- c("sigma",lifedamparamsvec)
    priors <- paste(c(distpriors,lifedampriors),collapse = " ")
  }

  # Define stancode here

  if(lifedam=="Hamada"){
    block1 <- "data {int<lower=0> n; vector[n] Dt; vector[n] t; vector[n] Ti; real<lower=0> Tu; }"
    datablock <- list(n = length(Dam), Dt = Dam, t = tDam, Ti = Tempf, Tu = Tuse)
  } else{
    block1 <- "data {int<lower=0> n; vector[n] Dt; vector[n] t;}"
    datablock <- list(n = length(Dam), Dt = Dam, t = tDam)
  }
  block2 <- paste(c("parameters {",params,"}"),collapse = " ")
  if (dist=="Lognormal"){
    block3 <- paste(c("model { vector[n] Damagei;",priors," for(i in 1:n){",logDamT,"}",loglik,"}"),collapse = " ")
  } else {
    block3 <- paste(c("model { vector[n] Damagei;",priors," for(i in 1:n){",DamT,"}",loglik,"}"),collapse = " ")
  }
  # NOT RUN {
  stanlscode <- paste(c(block1,block2,block3),collapse=" ")
  # Build or compile Stan code to C++
  lifedammod <- stan_model(model_code = stanlscode, verbose = TRUE)
  fit <- sampling(lifedammod, data = datablock, iter = nsamples, warmup = burnin, init = pt_est)
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
