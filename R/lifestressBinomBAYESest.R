# Bayesian Binomial Accelerated Life-Stress Estimator
# Developed by Mohammad Modarres and Reuel Smith, 2022

lifestress.Binom.BAYESest <- function(data,interact_stress,weight0,confid,priors,nsamples,burnin){
  # (pt_est,ls,dist,TTF,SF,Tc,Sc,confid,priors,nsamples,burnin)
  #Load pracma library for erf
  library(pracma)
  library(StanHeaders)
  library(rstan)
  library(ggplot2)
  library(shinystan)
  library(bayesplot)
  library(stringr)

  # Add input to this to include prior estimates for LS parameters.
  # Example: priors<-c("normal(3,4)","normal(1,4)", "lognormal(-2,3)")
  # I will have to cite the Rstan text for distributions in the code.  Use lookup("") for the translation.
  # The code takes these and separates them so that they are written into the stan file.
  # Then the code will run the program and compute the Bayes estimation
  # NOTE: The prior for p_0 must always be uniform(0,1) so only the theta coefficients need definition.

  # Sets the default for weight use to zero if input is not used
  if(missing(weight0)){
    weight0<-0
  }

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

  # Check and see that there are multiple stress levels
  if(dim(data)[2]<3 && weight0==0) {
    stop('Need one or more stress levels to generate estimates')
  }
  if(dim(data)[2]<4 && weight0==1) {
    stop('Need one or more stress levels to generate estimates')
  }

  # Check to see if there are more than one set of data
  if(dim(data)[1]==1) {
    stop('Need more than one data entry to generate estimates')
  }

  # Pull specific data from data matrix
  ni <- data[,1]
  ki <- data[,2]
  if(weight0==1){
    wi <- data[,dim(data)[2]]
  }

  # Start with the pre-processing of the data in which you obtain an initial parameter
  # estimate based on the curve fit of the data
  LSQout <- lifestress.Binom.LSQest(data,interact_stress,weight0)
  pt_est <- LSQout[[1]]
  Dstress <- LSQout[[2]]
  lifedamparamsvec <- LSQout[[3]]
  Dstress_txt <- str_replace_all(lifedamparamsvec[2:length(lifedamparamsvec)],"\U03B8","DS")
  theta_txt <- str_replace_all(lifedamparamsvec[2:length(lifedamparamsvec)],"\U03B8","theta")

  # Initialize life-stress parameter estimates for theta
  for(i in 1:length(pt_est)){
    if(i==1){
      lifedamparams <- "real<lower=0> p_0;"
      lifedampriors <- "p_0 ~ uniform(0,1); "
    } else{
      # lifedamparams <- paste(c(lifedamparams," real ",lifedamparamsvec[i],";"),collapse = "")
      lifedamparams <- paste(c(lifedamparams," real ",theta_txt[i-1],";"),collapse = "")
      # lifedampriors <- paste(c(lifedampriors,lifedamparamsvec[i]," ~ ",priors[i-1],";"),collapse = "")
      lifedampriors <- paste(c(lifedampriors,theta_txt[i-1]," ~ ",priors[i-1],";"),collapse = "")
      if(i==2){
        # stresssum <- paste(c(lifedamparamsvec[i],"*DSij[i,",i-1,"]"),collapse = "")
        stresssum <- paste(c(theta_txt[i-1],"*DSij[i,",i-1,"]"),collapse = "")
      } else{
        # stresssum <- paste(c(stresssum," + ",lifedamparamsvec[i],"*DSij[i,",i-1,"]"),collapse = "")
        stresssum <- paste(c(stresssum," + ",theta_txt[i-1],"*DSij[i,",i-1,"]"),collapse = "")
      }
    }
  }

  P_acc <- paste(c("P_i[i] = (1 - (1 - p_0)*exp(-(",stresssum,")));"),collapse = "")
  w_binom <- paste(c("wbinomial_lpmf[i] = Wi[i]*binomial_lpmf(Ki[i] | Ni[i], P_i[i]);"),collapse = "")

  # return(list(lifedamparamsvec,lifedamparams,lifedampriors,Dstress_txt,theta_txt,P_acc))

  # Fit to log-likelihood distributions
  if(weight0==0){
    loglik <- paste(c("target += binomial_lpmf(Ki | Ni, P_i);"),collapse = "")
  } else{
    # loglik <- paste(c("target += Wi*binomial_lpmf(Ki | Ni, P_i);"),collapse = "")
    loglik <- paste(c("target += wbinomial_lpmf;"),collapse = "")
  }

  params <- lifedamparams
  paramsvec <- lifedamparamsvec
  priors <- lifedampriors

  # Define stancode here

  if(weight0==0){
    block1 <- "data {int<lower=0> N; int<lower=0> M; vector[N] Ni; vector[N] Ki; matrix[N,M] DSij; vector[N] Wi;}"
    if(is.null(dim(Dstress))==TRUE){
      datablock <- list(N = length(ni), M = 1, Ni = ni, Ki = ki, DSij = Dstress)
    } else{
      datablock <- list(N = length(ni), M = dim(Dstress)[2], Ni = ni, Ki = ki, DSij = Dstress)
    }

  } else{
    block1 <- "data {int<lower=0> N; int<lower=0> M; vector[N] Ni; vector[N] Ki; matrix[N,M] DSij; vector[N] Wi;}"
    datablock <- list(N = length(ni), M = dim(Dstress)[2], Ni = ni, Ki = ki, DSij = Dstress, Wi = wi)
    # if(is.null(dim(Dstress))==TRUE){
    #   datablock <- list(N = length(ni), M = 1, Ni = ni, Ki = ki, DSij = Dstress)
    # } else{
    #   datablock <- list(N = length(ni), M = dim(Dstress)[2], Ni = ni, Ki = ki, DSij = Dstress)
    # }
  }
  block2 <- paste(c("parameters {",params,"}"),collapse = " ")
  block3 <- paste(c("model { vector[N] P_i; vector[N] wbinomial_lpmf;",priors," for(i in 1:N){",P_acc,w_binom,"}",loglik,"}"),collapse = " ")

  # NOT RUN {
  stanlscode <- paste(c(block1,block2,block3),collapse=" ")
  # return(stanlscode)
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
