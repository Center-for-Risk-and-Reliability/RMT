# Bayesian Updater for Probability Distributions
# Developed by Dr. Reuel Smith, 2021-2023

distribution.BAYESest <- function(pt_est,dist,TTF,Tc=NULL,Tlc=NULL,confid=0.95,priors,nsamples=20000,burnin,nchains=4,setbeta=NULL,N=NULL,Kf=NULL,TTFf=NULL,Tcf=NULL){
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

  # RCS04012025 - Adding TTFf and Tcf for failure time and censored time frequency respectively

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
  if (dist=="Weibull" && is.null(setbeta) == TRUE) {
    distparam <-"real<lower=0> alpha; real<lower=0> beta;"
    distpriors<-paste(c("alpha ~ ",priors[1],";","beta ~ ",priors[2],";"),collapse = "")

    # Adjust IF statement to include right-censored only data cases
    if(is.null(TTF) == FALSE && is.null(Tc) == TRUE){ # Only failure data
      loglik <- paste(c("target += weibull_lpdf(TTF | beta, alpha);"),collapse = "")
    }
    if(is.null(TTF) == TRUE && is.null(Tc) == FALSE){ # only right-censored data
      loglik <- paste(c("target += weibull_lccdf(TTS | beta, alpha);"),collapse = "")
    }
    if(is.null(TTF) == FALSE && is.null(Tc) == FALSE){ # both failure data and right-censored data
      loglik <- paste(c("target += weibull_lpdf(TTF | beta, alpha) + weibull_lccdf(TTS | beta, alpha);"),collapse = "")
    }
    params <- paste(c(distparam),collapse = " ")
    paramsvec <- c("alpha","beta")
    outputparamset <- c("\U03B1 (Scale)","\U03B2 (Shape)")
    priors <- paste(c(distpriors),collapse = " ")
  }
  # RCS03312025 - Add case where beta is known and given.  We treat this as an Exponential likelihood
  if (dist=="Weibull" && is.null(setbeta) == FALSE) {
    distparam <-"real<lower=0> alpha;"
    distpriors<-paste(c("alpha ~ ",priors[1],";"),collapse = "")

    # Adjust IF statement to include right-censored only data cases
    if(is.null(TTF) == FALSE && is.null(Tc) == TRUE && is.null(Tlc) == TRUE){ # Only failure data
      loglik <- paste(c("target += weibull_lpdf(TTF | beta, alpha);"),collapse = "")
    }
    if(is.null(TTF) == TRUE && is.null(Tc) == FALSE && is.null(Tlc) == TRUE){ # only right-censored data
      loglik <- paste(c("target += weibull_lccdf(TTS | beta, alpha);"),collapse = "")
    }
    if(is.null(TTF) == TRUE && is.null(Tc) == TRUE && is.null(Tlc) == FALSE){ # only left-censored data
      loglik <- paste(c("target += weibull_lccdf(TTS | beta, alpha);"),collapse = "")
    }
    if(is.null(TTF) == FALSE && is.null(Tc) == FALSE && is.null(Tlc) == TRUE){ # both failure data and right-censored data
      loglik <- paste(c("target += weibull_lpdf(TTF | beta, alpha) + weibull_lccdf(TTS | beta, alpha);"),collapse = "")
    }
    if(is.null(TTF) == FALSE && is.null(Tc) == FALSE && is.null(Tlc) == FALSE){ # failure data, right-censored data, and left-censored data
      loglik <- paste(c("target += weibull_lpdf(TTF | beta, alpha) + weibull_lccdf(TTS | beta, alpha) + weibull_lcdf(TTLC | beta, alpha);"),collapse = "")
    }

    params <- paste(c(distparam),collapse = " ")
    paramsvec <- c("alpha")
    # paramsvec0 <- c(paramsvec,"lambda")

    outputparamset <- c("\U03B1 (Scale)")
    priors <- paste(c(distpriors),collapse = " ")
  }
  if (dist=="3PWeibull") {
    distparam <-"real<lower=0> alpha; real<lower=0> beta; real gamma;"
    distpriors<-paste(c("alpha ~ ",priors[1],";","beta ~ ",priors[2],";","gamma ~ ",priors[3],";"),collapse = "")

    if(missing(Tc)){
      loglik <- paste(c("for(i in 1:n){	target += log(beta) + (beta-1)*log(TTF[i] - gamma) - beta*log(alpha) - (((TTF[i] -gamma)/alpha)^beta);}"),collapse = "")
    } else{
      loglik <- paste(c("for(i in 1:n){	target += log(beta) + (beta-1)*log(TTF[i] - gamma) - beta*log(alpha) - (((TTF[i] -gamma)/alpha)^beta);} for(j in 1:m){	target += - (((TTS[j] - gamma)/alpha)^beta);}"),collapse = "")
    }
    params <- paste(c(distparam),collapse = " ")
    paramsvec <- c("alpha","beta","gamma")
    outputparamset <- c("\U03B1 (Scale)","\U03B2 (Shape)","\U03B3 (Location)")
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
    outputparamset <- c("\U03BC_x (Scale)","\U03C3_X (Shape)")
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
    outputparamset <- c("\U03BC (Location)","\U03C3 (Scale)")
    priors <- paste(c(distpriors),collapse = " ")
  }
  if (dist=="Exponential") {
    distparam <-"real<lower=0> lambda;"
    distpriors<-paste(c("lambda ~ ",priors[1],";"),collapse = "")
    if(missing(Tc)){
      loglik <- paste(c("target += exponential_lpdf(TTF | lambda);"),collapse = "")
    } else{
      loglik <- paste(c("target += exponential_lpdf(TTF | lambda) + exponential_lccdf(TTS | lambda);"),collapse = "")
    }
    params <- paste(c(distparam),collapse = " ")
    paramsvec <- c("lambda")
    outputparamset <- c("\U03BB (Rate)")
    priors <- paste(c(distpriors),collapse = " ")
  }
  if (dist=="2PExponential") {
    distparam <-"real<lower=0> theta; real<lower=0> sigma;"
    distpriors<-paste(c("sigma ~ ",priors[ishift],";"),collapse = "")
    if(missing(Tc)){
      loglik <- paste(c("target += double_exponential_lpdf(TTF |theta, sigma);"),collapse = "")
    } else{
      loglik <- paste(c("target += double_exponential_lpdf(TTF |theta, sigma) + double_exponential_lccdf(TTS | theta, sigma);"),collapse = "")
    }
    params <- paste(c(distparam),collapse = " ")
    paramsvec <- c("theta","sigma")
    outputparamset <- c("\U03B8 (Location)","\U03C3 (Scale)")
    priors <- paste(c(distpriors),collapse = " ")
  }
  # NEW 11/26/2024 - Binomial Likelihood for updating probability p
  # If Binomial likelihood then TTF = NULL and so will Tc.  Set a new input for N and Kf (set to NULL normally)
  # and they will be single input or vectors depending on the scenario.
  if (dist=="Binomial") {
    distparam <-"real<lower=0> p;"
    distpriors<-paste(c("p ~ ",priors[1],";"),collapse = "")

    loglik <- paste(c("target += binomial_lpmf(nf | N, p);"),collapse = "")

    # if(missing(Tc)){
    #   loglik <- paste(c("target += weibull_lpdf(TTF | beta, alpha);"),collapse = "")
    # } else{
    #   loglik <- paste(c("target += weibull_lpdf(TTF | beta, alpha) + weibull_lccdf(TTS | beta, alpha);"),collapse = "")
    # }
    params <- paste(c(distparam),collapse = " ")
    paramsvec <- c("p")
    outputparamset <- c("p")
    priors <- paste(c(distpriors),collapse = " ")
  }
  if (dist=="Poisson") {
    distparam <-"real<lower=0> lambda;"
    distpriors<-paste(c("lambda ~ ",priors[1],";"),collapse = "")

    loglik <- paste(c("target += poisson_lpmf(nf | mu);"),collapse = "")

    params <- paste(c(distparam),collapse = " ")
    paramsvec <- c("lambda")
    outputparamset <- c("\U03BB (rate)")
    priors <- paste(c(distpriors),collapse = " ")
  }
  if (dist=="Jetliner") {
    distparam <-"real<lower=0> a; real<lower=0> b;"
    distpriors<-paste(c("a ~ ",priors[1],";","b ~ ",priors[2],";"),collapse = "")

    if(missing(Tc)){
      loglik <- paste(c("target += log(a) + log(b) + (b - 1)*log(TTF) - 2*log(1 + a*TTF^b);"),collapse = "")
    } else{
      loglik <- paste(c("target += log(a) + log(b) + (b - 1)*log(TTF) - 2*log(1 + a*TTF^b);"),collapse = "")
    }
    params <- paste(c(distparam),collapse = " ")
    paramsvec <- c("a","b")
    outputparamset <- c("a","b")
    priors <- paste(c(distpriors),collapse = " ")
  }

  # Define stancode here
  if((is.null(TTF)==FALSE && is.null(Tc)==TRUE && is.null(Tlc) == TRUE) && (dist != "Binomial" || dist != "Poisson")){
    block1 <- "data {int<lower=0> n; vector[n] TTF;}"
    datablock <- list(n = length(TTF), TTF = TTF)
    if (is.null(TTFf) == TRUE && dist=="Weibull" && is.null(setbeta) == FALSE){
      block1 <- "data {int<lower=0> n; vector[n] TTF; int<lower=0> beta;}"
      datablock <- list(n = length(TTF), TTF = TTF, beta = setbeta)
    }
    if (is.null(TTFf) == FALSE && dist=="Weibull" && is.null(setbeta) == FALSE){
      block1 <- "data {int<lower=0> n; vector[n] TTF; vector[n] Nf; int<lower=0> beta;}"
      datablock <- list(n = length(TTF), TTF = TTF, Nf = TTFf, beta = setbeta)
    }
  }
  if((is.null(TTF)==TRUE && is.null(Tc)==FALSE && is.null(Tlc) == TRUE) && (dist != "Binomial" || dist != "Poisson")){
    block1 <- "data {int<lower=0> m; vector[m] TTS;}"
    datablock <- list(m = length(Tc), TTS = Tc)
    if (is.null(Tcf) == TRUE&& dist=="Weibull" && is.null(setbeta) == FALSE){
      block1 <- "data {int<lower=0> m; vector[m] TTS; int<lower=0> beta;}"
      datablock <- list(m = length(Tc), TTS = Tc, beta = setbeta)
    }
    if(is.null(Tcf) == FALSE && dist=="Weibull" && is.null(setbeta) == FALSE){
      block1 <- "data {int<lower=0> m; vector[m] TTS; vector[m] Nc; int<lower=0> beta;}"
      datablock <- list(m = length(Tc), TTS = Tc, Nc = Tcf, beta = setbeta)
    }
  }
  if((is.null(TTF)==FALSE && is.null(Tc)==FALSE && is.null(Tlc) == TRUE) && (dist != "Binomial" || dist != "Poisson")){
    block1 <- "data {int<lower=0> n; int<lower=0> m; vector[n] TTF; vector[m] TTS;}"
    datablock <- list(n = length(TTF), m = length(Tc), TTF = TTF, TTS = Tc)
    if ((is.null(TTFf) == TRUE && is.null(Tcf) == TRUE) && dist=="Weibull" && is.null(setbeta) == FALSE){
      block1 <- "data {int<lower=0> n; int<lower=0> m; vector[n] TTF; vector[m] TTS; real<lower=0> beta;}"
      datablock <- list(n = length(TTF), m = length(Tc), TTF = TTF, TTS = Tc, beta = setbeta)
    }
    if ((is.null(TTFf) == FALSE && is.null(Tcf) == FALSE) && dist=="Weibull" && is.null(setbeta) == FALSE){
      block1 <- "data {int<lower=0> n; int<lower=0> m; vector[n] TTF; vector[m] TTS; vector[n] Nf; vector[m] Nc; real<lower=0> beta;}"
      datablock <- list(n = length(TTF), m = length(Tc), TTF = TTF, TTS = Tc, Nf = TTFf, Nc = Tcf, beta = setbeta)
    }
  }
  if((is.null(TTF)==FALSE && is.null(Tc)==FALSE && is.null(Tlc) == FALSE) && (dist != "Binomial" || dist != "Poisson")){
    block1 <- "data {int<lower=0> n; int<lower=0> m; int<lower=0> q; vector[n] TTF; vector[m] TTS; vector[q] TTLC;}"
    datablock <- list(n = length(TTF), m = length(Tc), q = length(Tlc), TTF = TTF, TTS = Tc, TTLC = Tlc)
    if ((is.null(TTFf) == TRUE && is.null(Tcf) == TRUE) && dist=="Weibull" && is.null(setbeta) == FALSE){
      block1 <- "data {int<lower=0> n; int<lower=0> m; int<lower=0> q; vector[n] TTF; vector[m] TTS; vector[q] TTLC; real<lower=0> beta;}"
      datablock <- list(n = length(TTF), m = length(Tc), q = length(Tlc), TTF = TTF, TTS = Tc, TTLC = Tlc, beta = setbeta)
    }
    if ((is.null(TTFf) == FALSE && is.null(Tcf) == FALSE) && dist=="Weibull" && is.null(setbeta) == FALSE){
      block1 <- "data {int<lower=0> n; int<lower=0> m; int<lower=0> q; vector[n] TTF; vector[m] TTS; vector[q] TTLC; vector[n] Nf; vector[m] Nc; real<lower=0> beta;}"
      datablock <- list(n = length(TTF), m = length(Tc), q = length(Tlc), TTF = TTF, TTS = Tc, TTLC = Tlc, Nf = TTFf, Nc = Tcf, beta = setbeta)
    }
  }
  # return(datablock)
  # if (dist=="Weibull" && is.null(setbeta) == FALSE){
  #   block2b <- paste(c("transformed parameters { real<lower=0> alphanew; alphanew = (1/alpha)^(1/beta);}"),collapse = " ")
  #   paramsvec <- "alphanew"
  # }
  if(dist == "Binomial"){
    block1 <- "data {int N; int nf;}"
    datablock <- list(N = as.integer(sum(N)), nf = as.integer(sum(Kf)))
  }
  if(dist == "Poisson"){
    block1 <- "data {real T; int nf;}"
    datablock <- list(T = TTF, nf = as.integer(N))
  }
  block2 <- paste(c("parameters {",params,"}"),collapse = " ")
  if(dist=="Poisson"){
    block2b <- paste(c("transformed parameters { real<lower=0> mu; mu = lambda*T;}"),collapse = " ")
    paramsvec0 <- c(paramsvec,"mu")
  }
  block3 <- paste(c("model {",priors,loglik,"}"),collapse = " ")
  # NOT RUN {
  stanlscode <- paste(c(block1,block2,block3),collapse=" ")
  if(dist=="Poisson"){
    stanlscode <- paste(c(block1,block2,block2b,block3),collapse=" ")
  }
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
  # Set up confidence limit text for output table
  conflim_txt<-c(paste(c("Lower ",100*conf.level,"%"),collapse = ""),paste(c("Upper ",100*conf.level,"%"),collapse = ""))
  # ==================================================================================================
  # NOTE RCS01102026 - The following block works under cmdstanr which is not functioning at the moment
  # ==================================================================================================
  # lsmod <- cmdstan_model(stanlsfile)
  # fit <- lsmod$sample(data = datablock, init = init_pt_est, chains = nchains, iter_warmup = burnin, iter_sampling = nsamples)
  # ==================================================================================================
  # PATCH RCS01102026 - The following block works under rstan which IS functioning at the moment
  # ==================================================================================================
  lsmod <- stan_model(model_code = stanlscode, verbose = TRUE)
  fit <- sampling(lsmod, data = datablock, iter = nsamples, warmup = burnin, init = pt_est)

  # return(lsmod)
  # }
  # return(fit)
  # Print results.  I need to get this as an output
  # ==================================================================================================
  # NOTE RCS01102026 - The following block works under cmdstanr which is not functioning at the moment
  # ==================================================================================================
  # stats <- fit$summary(variables = paramsvec)
  # confidbounds <- mcmc_intervals_data(fit$draws(variables = paramsvec),prob_outer = confid)
  # outputtable <- matrix(c(stats[[2]],stats[[4]],confidbounds[[5]],stats[[3]],confidbounds[[9]],stats[[8]]), nrow = length(outputparamset), ncol = 6, byrow = FALSE,dimnames = list(outputparamset,c("Mean","Standard Deviation",conflim_txt[1],"Median",conflim_txt[2],"R\U005E")))
  # ==================================================================================================
  # PATCH RCS01102026 - The following block works under rstan which IS functioning at the moment
  # ==================================================================================================
  stats.mean.sd <- summary(fit)$summary[,c(1,3)]
  stats.Rhat <- rhat(fit)
  confidbounds <- mcmc_intervals_data(data.frame(extract(fit, paramsvec)),prob_outer = confid)
  outputtable <- matrix(c(unname(stats.mean.sd)[1:length(paramsvec),1],unname(stats.mean.sd)[1:length(paramsvec),2],confidbounds[[5]],confidbounds[[7]],confidbounds[[9]],unname(stats.Rhat)[1:length(paramsvec)]), nrow = length(paramsvec), ncol = 6, byrow = FALSE,dimnames = list(paramsvec,c("Mean","Standard Deviation",conflim_txt[1],"Median",conflim_txt[2],"R\U005E")))

  # stats <- print(fit, pars = paramsvec, probs=c((1-confid)/2,.5,1-(1-confid)/2))
  # dataout <- fit@.MISC[["summary"]][["msd"]]
  # dataout <- fit$draws(format = "df")

  # NEW Post Bayes analysis plotting of posterior (12/12/25)
  if (dist=="Weibull" && is.null(setbeta) == TRUE){
    posterior_beta <- density(extract(fit,c("beta"))$beta)
    posterior_alpha <- density(extract(fit,c("alpha"))$alpha)
    df_posterior <- data.frame(x = c(posterior_alpha$x,posterior_beta$x),
                               ymin = rep(0,(length(posterior_alpha$x)+length(posterior_beta$x))),
                               ymax = c(posterior_alpha$y,posterior_beta$y),
                               distlabel = c(rep("α",length(posterior_alpha$x)),rep("β",length(posterior_beta$x))))

    # Density plot for Exponential rate parameter posterior
    plot3_density <- ggplot() + geom_ribbon(data = df_posterior, aes(x=x, ymin = ymin, ymax = ymax), fill = "red" ,alpha = 0.5) +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      facet_wrap(~distlabel, dir="v", scales = "free") +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0)) +
      xlab(" ")
      ylab("density")
  }
  if (dist=="Weibull" && is.null(setbeta) == FALSE){
    posterior_alpha <- density(extract(fit,c("alpha"))$alpha)
    df_posterior_alpha <- data.frame(x = c(posterior_alpha$x), ymin = rep(0,(length(posterior_alpha$x))), ymax = c(posterior_alpha$y), distlabel = c(rep("posterior",length(posterior_alpha$x))))
    #
    # Density plot for Exponential rate parameter posterior
    plot3_density <- ggplot() + geom_ribbon(data = df_posterior_alpha, aes(x=x, ymin = ymin, ymax = ymax), fill = "red" ,alpha = 0.5) +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0)) +
      xlab("α") +
      ylab("density")
  }
  if (dist=="3PWeibull") {
    posterior_beta <- density(extract(fit,c("beta"))$beta)
    posterior_alpha <- density(extract(fit,c("alpha"))$alpha)
    posterior_gamma <- density(extract(fit,c("gamma"))$gamma)
    df_posterior <- data.frame(x = c(posterior_alpha$x,posterior_beta$x,posterior_gamma$x),
                               ymin = rep(0,(length(posterior_alpha$x)+length(posterior_beta$x)+length(posterior_gamma$x))),
                               ymax = c(posterior_alpha$y,posterior_beta$y,posterior_gamma$y),
                               distlabel = c(rep("α",length(posterior_alpha$x)),rep("β",length(posterior_beta$x)),rep("γ",length(posterior_gamma$x))))

    # Density plot for Exponential rate parameter posterior
    plot3_density <- ggplot() + geom_ribbon(data = df_posterior, aes(x=x, ymin = ymin, ymax = ymax), fill = "red" ,alpha = 0.5) +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      facet_wrap(~distlabel, dir="v", scales = "free") +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0)) +
      xlab(" ")
    ylab("density")
  }
  if (dist=="Lognormal") {
    posterior_mu_t <- density(extract(fit,c("mu_t"))$mu_t)
    posterior_sigma_t <- density(extract(fit,c("sigma_t"))$sigma_t)
    df_posterior <- data.frame(x = c(posterior_mu_t$x,posterior_sigma_t$x),
                               ymin = rep(0,(length(posterior_mu_t$x)+length(posterior_sigma_t$x))),
                               ymax = c(posterior_mu_t$y,posterior_sigma_t$y),
                               distlabel = c(rep("μ_t",length(posterior_mu_t$x)),rep("σ_t",length(posterior_sigma_t$x))))

    # Density plot for Exponential rate parameter posterior
    plot3_density <- ggplot() + geom_ribbon(data = df_posterior, aes(x=x, ymin = ymin, ymax = ymax), fill = "red" ,alpha = 0.5) +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      facet_wrap(~distlabel, dir="v", scales = "free") +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0)) +
      xlab(" ")
    ylab("density")
  }
  if (dist=="Normal") {
    posterior_mu <- density(extract(fit,c("mu"))$mu)
    posterior_sigma <- density(extract(fit,c("sigma"))$sigma)
    df_posterior <- data.frame(x = c(posterior_mu$x,posterior_sigma$x),
                               ymin = rep(0,(length(posterior_mu$x)+length(posterior_sigma$x))),
                               ymax = c(posterior_mu$y,posterior_sigma$y),
                               distlabel = c(rep("μ",length(posterior_mu$x)),rep("σ",length(posterior_sigma$x))))

    # Density plot for Exponential rate parameter posterior
    plot3_density <- ggplot() + geom_ribbon(data = df_posterior, aes(x=x, ymin = ymin, ymax = ymax), fill = "red" ,alpha = 0.5) +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      facet_wrap(~distlabel, dir="v", scales = "free") +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0)) +
      xlab(" ")
    ylab("density")
  }
  if (dist=="Exponential") {
    posterior_lambda <- density(extract(fit,c("lambda"))$lambda)
    df_posterior_lambda <- data.frame(x = c(posterior_lambda$x), ymin = rep(0,(length(posterior_lambda$x))), ymax = c(posterior_lambda$y), distlabel = c(rep("posterior",length(posterior_lambda$x))))

    # Density plot for Exponential rate parameter posterior
    plot3_density <- ggplot() + geom_ribbon(data = df_posterior_lambda, aes(x=x, ymin = ymin, ymax = ymax), fill = "red" ,alpha = 0.5) +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0)) +
      xlab("λ") +
      ylab("density")
  }



  # Trace the Markov Chains for each parameter
  # plot1_MCtrace <- traceplot(fit, pars = paramsvec, inc_warmup = TRUE, nrow = 3)
  # plot1_MCtrace <- mcmc_trace(as.matrix(fit),pars=paramsvec, facet_args = list(nrow = length(paramsvec), labeller = label_parsed))
  # plot2_hist <- stan_hist(fit)
  # plot3_density <- stan_dens(fit)
  plot1_MCtrace <- mcmc_trace(fit,paramsvec) +
    theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4))
  plot2_hist <- mcmc_hist(fit,paramsvec) +
    theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4))

  # plot1_MCtrace <- mcmc_trace(fit$draws(paramsvec)) +
  #   theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4))
  # plot2_hist <- mcmc_hist(fit$draws(paramsvec)) +
  #   theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4))
 # Produce some output text that summarizes the results
  cat(c("Posterior estimates for Bayesian Analysis.\n\n"),sep = "")
  print(outputtable)
  cat(c("\n"),sep = "")


  return(list(posterior.fit=fit,post.stats=stats.mean.sd,MC.trace=plot1_MCtrace,post.histogram=plot2_hist,post.density=plot3_density,stanlscode))
}
