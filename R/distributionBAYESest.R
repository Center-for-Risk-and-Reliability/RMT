# Bayesian Updater for Probability Distributions
# Developed by Dr. Reuel Smith, 2021-2023

distribution.BAYESest <- function(pt_est,dist,TTF,Tc=NULL,confid=0.95,priors,nsamples=20000,burnin,nchains=4){
  #Load pracma library for erf
  library(pracma)
  library(StanHeaders)
  library(rstan)
  library(ggplot2)
  library(shinystan)
  library(cmdstanr)
  library(bayesplot)

  # Add input to this to include prior estimates for LS parameters.
  # Example: priors<-c("normal(3,4)","normal(1,4)", "lognormal(-2,3)")
  # I will have to cite the Rstan text for distributions in the code.  Use lookup("") for the translation.
  # The code takes these and separates them so that they are written into the stan file.
  # Then the code will run the program and compute the Bayes estimation

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

  # Fit to log-likelihood distributions
  if (dist=="Weibull") {
    distparam <-"real<lower=0> alpha; real<lower=0> beta;"
    distpriors<-paste(c("alpha ~ ",priors[1],";","beta ~ ",priors[2],";"),collapse = "")

    if(missing(Tc)){
      loglik <- paste(c("target += weibull_lpdf(TTF | beta, alpha);"),collapse = "")
    } else{
      loglik <- paste(c("target += weibull_lpdf(TTF | beta, alpha) + weibull_lccdf(TTS | beta, alpha);"),collapse = "")
    }
    params <- paste(c(distparam),collapse = " ")
    paramsvec <- c("alpha","beta")
    priors <- paste(c(distpriors),collapse = " ")
  }
  if (dist=="Lognormal") {
    distparam <-"real mu_t; real<lower=0> sigma_t;"
    distpriors<-paste(c("mu_t ~ ",priors[1],";","sigma_t ~ ",priors[2],";"),collapse = "")

    if(missing(Tc)){
      loglik <- paste(c("target += lognormal_lpdf(TTF |mu_t, sigma_t);"),collapse = "")
    } else{
      loglik <- paste(c("target += lognormal_lpdf(TTF |mu_t, sigma_t) + lognormal_lccdf(TTS |mu_t, sigma_t);"),collapse = "")
    }
    params <- paste(c(distparam),collapse = " ")
    paramsvec <- c("mu_t","sigma_t")
    priors <- paste(c(distpriors),collapse = " ")
  }
  if (dist=="Normal") {
    distparam <-"real mu; real<lower=0> sigma;"
    distpriors<-paste(c("mu ~ ",priors[1],";","sigma ~ ",priors[2],";"),collapse = "")
    if(missing(Tc)){
      loglik <- paste(c("target += normal_lpdf(TTF |mu, sigma);"),collapse = "")
    } else{
      loglik <- paste(c("target += normal_lpdf(TTF |mu, sigma) + normal_lccdf(TTS |mu, sigma);"),collapse = "")
    }
    params <- paste(c(distparam),collapse = " ")
    paramsvec <- c("mu","sigma")
    priors <- paste(c(distpriors),collapse = " ")
  }
  if (dist=="Exponential") {
    distparam <-"real<lower=0> lambda;"
    distpriors<-paste(c("lambda ~ ",priors[1],";"),collapse = "")
    if(missing(Tc)){
      loglik <- paste(c("target += exponential_lpdf(TTF | 1/(lambda));"),collapse = "")
    } else{
      loglik <- paste(c("target += exponential_lpdf(TTF | 1/(lambda)) + exponential_lccdf(TTS 1/(lambda));"),collapse = "")
    }
    params <- paste(c(distparam),collapse = " ")
    paramsvec <- c("lambda")
    priors <- paste(c(distpriors),collapse = " ")
  }
  if (dist=="2PExponential") {
    distparam <-"real<lower=0> theta; real<lower=0> sigma;"
    distpriors<-paste(c("sigma ~ ",priors[ishift],";"),collapse = "")
    if(missing(Tc)){
      loglik <- paste(c("target += double_exponential_lpdf(TTF |theta, sigma);"),collapse = "")
    } else{
      loglik <- paste(c("target += double_exponential_lpdf(TTF |theta, sigma) + double_exponential_lccdf(TTS |theta, sigma);"),collapse = "")
    }
    params <- paste(c(distparam),collapse = " ")
    paramsvec <- c("theta","sigma")
    priors <- paste(c(distpriors),collapse = " ")
  }

  # Define stancode here
  if(missing(Tc)){
    block1 <- "data {int<lower=0> n; vector[n] TTF;}"
    datablock <- list(n = length(TTF), TTF = TTF)
  } else{
    block1 <- "data {int<lower=0> n; int<lower=0> m; vector[n] TTF; vector[m] TTS;}"
    datablock <- list(n = length(TTF), m = length(Tc), TTF = TTF, TTS = Tc)
  }
  block2 <- paste(c("parameters {",params,"}"),collapse = " ")
  block3 <- paste(c("model {",priors,loglik,"}"),collapse = " ")
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
  stats <- fit$summary(variables = paramsvec)
  dataout <- fit$draws(format = "df")

  # Trace the Markov Chains for each parameter
  # plot1_MCtrace <- traceplot(fit, pars = paramsvec, inc_warmup = TRUE, nrow = 3)
  # plot1_MCtrace <- mcmc_trace(as.matrix(fit),pars=paramsvec, facet_args = list(nrow = length(paramsvec), labeller = label_parsed))
  # plot2_hist <- stan_hist(fit)
  # plot3_density <- stan_dens(fit)
  plot1_MCtrace <- mcmc_trace(fit$draws(paramsvec))
  plot2_hist <- mcmc_hist(fit$draws(paramsvec))
  plot3_density <- mcmc_dens(fit$draws(paramsvec))
  plot4_densityoverlay <- mcmc_dens_overlay(fit$draws(paramsvec))


  return(list(fit,stats,dataout,plot1_MCtrace,plot2_hist,plot3_density,plot4_densityoverlay))
}
