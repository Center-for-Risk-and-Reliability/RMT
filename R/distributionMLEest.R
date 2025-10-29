# Maximum Likelihood Estimator for Probability Distributions
# Developed by Dr. Reuel Smith, 2021-2023

distribution.MLEest <- function(LSQest,dist,TTF,Tc=NULL,confid=0.95,sided="twosided"){
  #Load pracma library for erf
  library(pracma)
  library(matrixcalc)
  library(zipfR)

  # Check to see if confidence exists
  conf.level <- confid

  # Setup positivity check vector for parameters
  positivity_v<-rep(0,length(LSQest))

  # Fit to log-likelihood distributions
  if (dist=="Weibull") {
    # positivity_v[1]<-1
    # positivity_v[2]<-1
    # shift alpha to log alpha
    LSQest[1] <- log(LSQest[1])
    # shift beta to log beta
    LSQest[2] <- log(LSQest[2])

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(theta[2] + (exp(theta[2])-1)*log(TTF) - exp(theta[2])*theta[1] - ((TTF/exp(theta[1]))^exp(theta[2])))
      }
    } else{
      loglik <- function(theta){
        -sum(theta[2] + (exp(theta[2])-1)*log(TTF) - exp(theta[2])*theta[1] - ((TTF/exp(theta[1]))^exp(theta[2]))) - sum(- ((Tc/exp(theta[1]))^exp(theta[2])))
      }
    }
  }

  if (dist=="3PWeibull") {
    # positivity_v[1]<-1
    # positivity_v[2]<-1
    # Check Gamma first and if it is greater than any of the times, reset
    if(LSQest[3] > min(TTF) || (is.null(Tc) == FALSE && LSQest[3] > min(Tc))){
      if(is.null(Tc)==TRUE){
        LSQest[3] <- 0.5*min(TTF)
      } else{
        LSQest[3] <- 0.5*min(Tc)
      }
      xiRFblock_list <- plotposit.select((TTF-LSQest[3]),(Tc-LSQest[3]),pp="Blom")
      LSQest[1:2] <- c(probplotparam.wbl(xiRFblock_list[,1],xiRFblock_list[,3])[[3]])
    }
    # return(LSQest)
    # shift alpha to log alpha
    LSQest[1] <- log(LSQest[1])
    # shift beta to log beta
    LSQest[2] <- log(LSQest[2])

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(theta[2] + (exp(theta[2])-1)*log(TTF-theta[3]) - exp(theta[2])*theta[1] - (((TTF-theta[3])/exp(theta[1]))^exp(theta[2])))
      }
    } else{
      loglik <- function(theta){
        -sum(theta[2] + (exp(theta[2])-1)*log(TTF - theta[3]) - exp(theta[2])*theta[1] - (((TTF - theta[3])/exp(theta[1]))^exp(theta[2]))) - sum(- (((Tc - theta[3])/exp(theta[1]))^exp(theta[2])))
      }
    }
  }

  if (dist=="Lognormal") {
    # positivity_v[2]<-1
    # shift sigma to log sigma
    LSQest[2] <- log(LSQest[2])

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(-theta[2] - log(TTF) - 0.5*log(2*pi) - 0.5*(exp(theta[2])^-2)*((log(TTF) - theta[1])^2))
      }
    } else{
      loglik <- function(theta){
        -sum(-theta[2] - 0.5*log(2*pi) - 0.5*(exp(theta[2])^-2)*((log(TTF) - theta[1])^2)) - sum(log(0.5 - 0.5*erf((2^-0.5)*(exp(theta[2])^-1)*(log(Tc) - theta[1]))))
      }
    }
  }
  if (dist=="Normal") {
    # positivity_v[2]<-1
    # shift sigma to log sigma
    LSQest[2] <- log(LSQest[2])

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(-theta[2] - 0.5*log(2*pi) - 0.5*(exp(theta[2])^-2)*((TTF - theta[1])^2))
      }
    } else{
      loglik <- function(theta){
        -sum(-theta[2] - 0.5*log(2*pi) - 0.5*(exp(theta[2])^-2)*((TTF - theta[1])^2)) - sum(log(0.5 - 0.5*erf((2^-0.5)*(exp(theta[2])^-1)*(Tc - theta[1]))))
      }
    }
  }
  if (dist=="Exponential") {
    # shift lambda to log lambda
    LSQest[1] <- log(LSQest[1])

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(theta - TTF*exp(theta))
      }
    } else{
      loglik <- function(theta){
        -sum(theta - TTF*exp(theta)) + sum( Tc*exp(theta))
      }
    }
  }
  if (dist=="2PExponential") {
    positivity_v<-rep(1,length(LSQest))

    if(is.null(Tc)){
      # Estimate directly theta
      LSQest[1] <- (length(TTF)*min(TTF) - mean(TTF))/(length(TTF) - 1)
      # and sigma
      LSQest[2] <- mean(TTF) - LSQest[1]
      # and variance
      SIGMA_est <- cbind(c((LSQest[2]^2)/(length(TTF)*(length(TTF)-1)),0),c(0,(LSQest[2]^2)/(length(TTF)-1)))
      # loglikelihood
      loglik <- function(theta){
        -sum(log(min(c(TTF)) - theta[1]) - log(theta[2]) - (theta[2]^-1)*(TTF - theta[1]))
      }
    } else{
      # Estimate directly theta
      LSQest[1] <- ((length(TTF)+length(Tc))*min(TTF) - (1/length(TTF))*sum(c(TTF,Tc)))/(length(TTF) + length(Tc) - 1 - (length(Tc)/length(TTF)))
      # and sigma
      LSQest[2] <- (1/length(TTF))*(sum(c(TTF,Tc))) - LSQest[1]*(length(Tc)/length(TTF))
      # loglikelihood
      SIGMA_est <- cbind(c((length(TTF)*(LSQest[2]^2))/(((length(TTF)+length(Tc))^2)*(length(TTF)-1)),0),c(0,(LSQest[2]^2)/(length(TTF)-1)))
      # loglikelihood
      loglik <- function(theta){
        -sum(log(min(c(TTF,Tc)) - theta[1]) - log(theta[2]) - (theta[2]^-1)*(TTF - theta[1])) - sum(-(theta[2]^-1)*(Tc - theta[1]))
      }
    }
  }
  if (dist=="Gamma") {
    # Re-parameterize alpha and beta to mu and lambda
    LSQest0 <- c(0,0)
    # mu = ln beta + ln alpha (-Inf to Inf)
    LSQest0[1] <- log(LSQest[2]) + log(LSQest[1])
    # lambda = 1/sqrt(alpha) (positive so make loglambda)
    LSQest0[2] <- log(1/sqrt(LSQest[1]))
    LSQest <- LSQest0

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(-log(gamma(1/(exp(theta[2]))^2)) - (1/(exp(theta[2]))^2)*(theta[1] + log((exp(theta[2]))^2)) + ((1/(exp(theta[2]))^2) - 1)*log(TTF) - (1/(exp(theta[2]))^2)*exp(log(TTF) - theta[1]))
      }
    } else{
      loglik <- function(theta){
        -sum(-log(gamma(1/(exp(theta[2]))^2)) - (1/(exp(theta[2]))^2)*(theta[1] + log((exp(theta[2]))^2)) + ((1/(exp(theta[2]))^2) - 1)*log(TTF) - (1/(exp(theta[2]))^2)*exp(log(TTF) - theta[1])) - sum(log(1 - Rgamma(1/(exp(theta[2]))^2,exp(log(Tc) - theta[1])/(exp(theta[2]))^2),lower=TRUE))
      }
    }
  }
  if (dist=="3PGamma") {
    # Re-parameterize alpha, beta, and gamma to mu, lambda, and sigma
    LSQest0 <- c(0,0,0)
    # mu = ln beta + (1/gamma)*ln alpha (-Inf to Inf)
    LSQest0[1] <- log(LSQest[2]) + (1/LSQest[3])*log(LSQest[1])
    # lambda = 1/sqrt(alpha) (positive so make loglambda)
    LSQest0[2] <- log(1/sqrt(LSQest[1]))
    # sigma = 1/gamma*sqrt(alpha) (positive so make logsigma)
    LSQest0[3] <- log(1/(LSQest[3]*sqrt(LSQest[1])))
    LSQest <- LSQest0

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(log(abs(exp(theta[2]))) - log(exp(theta[3])) - log(gamma(1/(exp(theta[2]))^2)) - (1/(exp(theta[2]))^2)*(((theta[1]*exp(theta[2]))/exp(theta[3])) + log((exp(theta[2]))^2)) + ((1/(exp(theta[2])*exp(theta[3]))) - 1)*log(TTF) - (1/(exp(theta[2]))^2)*exp((exp(theta[2])/exp(theta[3]))*(log(TTF) - theta[1])))
      }
    } else{
      loglik <- function(theta){
        -sum(log(abs(exp(theta[2]))) - log(exp(theta[3])) - log(gamma(1/(exp(theta[2]))^2)) - (1/(exp(theta[2]))^2)*(((theta[1]*exp(theta[2]))/exp(theta[3])) + log((exp(theta[2]))^2)) + ((1/(exp(theta[2])*exp(theta[3]))) - 1)*log(TTF) - (1/(exp(theta[2]))^2)*exp((exp(theta[2])/exp(theta[3]))*(log(TTF) - theta[1]))) - sum(log(1 - Rgamma(1/(exp(theta[2]))^2,exp((exp(theta[2])/exp(theta[3]))*(log(Tc) - theta[1]))/(exp(theta[2]))^2),lower=TRUE))
      }
    }
  }
  if (dist=="Logistic") {
    # positivity_v[2]<-1
    LSQest[2] <- log(LSQest[2])

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(((TTF - theta[1])/exp(theta[2])) - theta[2] - 2*log(1 + exp((TTF - theta[1])/exp(theta[2]))))
      }
    } else{
      loglik <- function(theta){
        -sum(((TTF - theta[1])/exp(theta[2])) - theta[2] - 2*log(1 + exp((TTF - theta[1])/exp(theta[2])))) + sum(log(1 + exp((Tc - theta[1])/exp(theta[2]))))
      }
    }
  }
  if (dist=="Loglogistic") {
    # positivity_v[2]<-1
    LSQest[2] <- log(LSQest[2])

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(((log(TTF) - theta[1])/exp(theta[2])) - theta[2] - log(TTF*((1 + exp((log(TTF) - theta[1])/exp(theta[2])))^2)))
      }
    } else{
      loglik <- function(theta){
        -sum(((log(TTF) - theta[1])/exp(theta[2])) - theta[2] - log(TTF*((1 + exp((log(TTF) - theta[1])/exp(theta[2])))^2))) + sum(log(1 + exp((log(Tc) - theta[1])/exp(theta[2]))))
      }
    }
  }
  if (dist=="Gumbel") {
    # positivity_v[2]<-1
    # shift sigma to log sigma
    LSQest[2] <- log(LSQest[2])

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(-theta[2] + ((TTF - theta[1])/exp(theta[2])) - exp((TTF - theta[1])/exp(theta[2])))
      }
    } else{
      loglik <- function(theta){
        -sum(-theta[2] + ((TTF - theta[1])/exp(theta[2])) - exp((TTF - theta[1])/exp(theta[2]))) + sum(exp((Tc - theta[1])/exp(theta[2])))
      }
    }
  }
  # return(list(loglik,LSQest))

  if (dist=="3PWeibull"){
    # Debugger check to see what initial conditions are under 3PWeibull
    # print(cbind(c(LSQest,0.99*min(c(TTF,Tc)))))
    MLEandvar <- MLE.var.covar.select(loglik,LSQest,0.99*min(c(TTF,Tc)))
  }
  if (dist=="2PExponential"){
    MLEandvar <- list(LSQest,SIGMA_est)
  }
  if (dist=="Weibull" || dist=="Normal" || dist=="Lognormal" || dist=="Exponential" || dist=="Gumbel" || dist=="Logistic" || dist=="Loglogistic" || dist=="Gamma" || dist=="3PGamma"){
    MLEandvar <- MLE.var.covar.select(loglik,LSQest)
  }
  # return(MLEandvar)

  theta.hat <- MLEandvar[[1]]
  inv.fish  <- MLEandvar[[2]]
  loglik.hat <- -loglik(theta.hat)
  likeli.hat <- exp(loglik.hat)

  # return(list(loglik,MLEandvar))
  crit <- qnorm((1 + conf.level)/2)
  crit2 <- qnorm(conf.level)
  conflim<-vector(mode = "list", length = length(LSQest))
  fulllimset<-vector(mode = "list", length = length(LSQest))

  for(i in 1:length(theta.hat)){
    if(sided == "twosided"){
      if(positivity_v[i]==0){
        conflim[[i]]<-theta.hat[i] + c(-1, 1) * crit * sqrt(inv.fish[i, i])
      } else if(positivity_v[i]==1){
        conflim[[i]]<-theta.hat[i]*exp(c(-1, 1) * crit * (sqrt(inv.fish[i, i])/theta.hat[i]))
      }
      if((dist=="Weibull" || dist=="3PWeibull") && (i == 1 || i == 2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(dist=="Exponential" && i == 1){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if ((dist=="Normal" || dist=="Lognormal" || dist=="Logistic" || dist=="Loglogistic" || dist=="Gamma" || dist=="3PGamma" || dist=="Gumbel") && i == 2){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(dist=="3PGamma" && i == 3){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      conflim_txt<-c(paste(c("Lower ",100*conf.level,"%"),collapse = ""),paste(c("Upper ",100*conf.level,"%"),collapse = ""))
    }
    if(sided == "onesidedlow"){
      if(positivity_v[i]==0){
        conflim[[i]]<-theta.hat[i] - crit2 * sqrt(inv.fish[i, i])
      } else if(positivity_v[i]==1){
        conflim[[i]]<-theta.hat[i]*exp(-crit2 * (sqrt(inv.fish[i, i])/theta.hat[i]))
      }
      if((dist=="Weibull" || dist=="3PWeibull") && (i == 1 || i == 2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(dist=="Exponential" && i == 1){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if ((dist=="Normal" || dist=="Lognormal" || dist=="Logistic" || dist=="Loglogistic" || dist=="Gamma" || dist=="3PGamma" || dist=="Gumbel") && i == 2){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(dist=="3PGamma" && i == 3){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      conflim_txt<-paste(c("One-Sided Low ",100*conf.level,"%"),collapse = "")
    }
    if(sided == "onesidedhigh"){
      if(positivity_v[i]==0){
        conflim[[i]]<-theta.hat[i] + crit2 * sqrt(inv.fish[i, i])
      } else if(positivity_v[i]==1){
        conflim[[i]]<-theta.hat[i]*exp(crit2 * (sqrt(inv.fish[i, i])/theta.hat[i]))
      }
      if((dist=="Weibull" || dist=="3PWeibull") && (i == 1 || i == 2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(dist=="Exponential" && i == 1){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if ((dist=="Normal" || dist=="Lognormal" || dist=="Logistic" || dist=="Loglogistic" || dist=="Gamma" || dist=="3PGamma" || dist=="Gumbel") && i == 2){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(dist=="3PGamma" && i == 3){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      conflim_txt<-paste(c("One-Sided High ",100*conf.level,"%"),collapse = "")
    }
  }

  if(dist=="Gamma"){
    conflim0 <- conflim
    conflim[[1]] <- 1/(conflim0[[2]])^2
    conflim[[2]] <- exp(conflim0[[1]] + log((conflim0[[2]])^2))
  }
  if(dist=="3PGamma"){
    conflim0 <- conflim
    conflim[[1]] <- 1/(conflim0[[2]])^2
    conflim[[2]] <- exp(conflim0[[1]] + (conflim0[[3]]/conflim0[[2]])*log((conflim0[[2]])^2))
    conflim[[3]] <- conflim0[[2]]/conflim0[[3]]
  }

  AIC = 2*length(theta.hat) + 2*loglik(theta.hat)
  BIC = 2*log(length(TTF)+length(Tc)) + 2*loglik(theta.hat)

  # Recompute necessary output
  if(dist=="Weibull" || dist=="3PWeibull"){
    theta.hat[1] <- exp(theta.hat[1])
    theta.hat[2] <- exp(theta.hat[2])
  }
  if(dist=="Gamma"){
    theta.hat0 <- theta.hat
    theta.hat[1] <- 1/exp(theta.hat0[2])^2
    theta.hat[2] <- exp(theta.hat0[1] + log(exp(theta.hat0[2])^2))
  }
  if(dist=="3PGamma"){
    theta.hat0 <- theta.hat
    theta.hat[1] <- 1/(theta.hat0[2])^2
    theta.hat[2] <- exp(theta.hat0[1] + (theta.hat0[3]/theta.hat0[2])*log((theta.hat0[2])^2))
    theta.hat[3] <- theta.hat0[2]/theta.hat0[3]
  }
  if(dist=="Exponential"){
    theta.hat[1] <- exp(theta.hat[1])
  }
  if(dist=="Normal" || dist=="Lognormal" || dist=="Logistic" || dist=="Loglogistic" || dist=="Gumbel"){
    theta.hat[2] <- exp(theta.hat[2])
  }

  return(list(theta.hat,inv.fish,conflim,loglik.hat,likeli.hat,AIC,BIC))
}
