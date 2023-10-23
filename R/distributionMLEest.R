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
    positivity_v[1]<-1
    positivity_v[2]<-1

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(log(theta[2]) + (theta[2]-1)*log(TTF) - theta[2]*log(theta[1]) - ((TTF/theta[1])^theta[2]))
      }
    } else{
      loglik <- function(theta){
        -sum(log(theta[2]) + (theta[2]-1)*log(TTF) - theta[2]*log(theta[1]) - ((TTF/theta[1])^theta[2])) - sum(- ((Tc/theta[1])^theta[2]))
      }
    }
  }

  if (dist=="3PWeibull" && LSQest[2] > 2) {
    positivity_v[1]<-1
    positivity_v[2]<-1

    if(is.null(Tc)){
      loglik <- function(theta){
        ll <- -sum(log(theta[2]) + (theta[2]-1)*log(TTF-theta[3]) - theta[2]*log(theta[1]) - (((TTF-theta[3])/theta[1])^theta[2]))
        attr(ll, "gradient") <- -c(sum(-(theta[2]/theta[1]) + (theta[2]/(theta[1]^(1+theta[2])))*((TTF-theta[3])^theta[2])),
                                  sum((1/theta[2]) - log(theta[1]) + log(TTF-theta[3]) - log((TTF-theta[3])/theta[1])*(((TTF-theta[3])/theta[1])^theta[2])),
                                  sum(-((theta[2]-1)/(TTF-theta[3])) + (theta[2]/theta[1])*(((TTF-theta[3])/theta[1])^(theta[2]-1))))
        attr(ll,"hessian") <- -rbind(c(sum(theta[2]/theta[1]^2 + (2*theta[2]*(theta[3] - TTF)*(-(theta[3] - TTF)/theta[1])^(theta[2] - 1))/theta[1]^3 - (theta[2]*(theta[3] - TTF)^2*(-(theta[3] - TTF)/theta[1])^(theta[2] - 2)*(theta[2] - 1))/theta[1]^4),
                                       sum(- 1/theta[1] - ((theta[3] - TTF)*(-(theta[3] - TTF)/theta[1])^(theta[2] - 1))/theta[1]^2 - (theta[2]*log(-(theta[3] - TTF)/theta[1])*(theta[3] - TTF)*(-(theta[3] - TTF)/theta[1])^(theta[2] - 1))/theta[1]^2),
                                       sum((theta[2]*(theta[3] - TTF)*(-(theta[3] - TTF)/theta[1])^(theta[2] - 2)*(theta[2] - 1))/theta[1]^3 - (theta[2]*(-(theta[3] - TTF)/theta[1])^(theta[2] - 1))/theta[1]^2)),
                                     c(sum(- 1/theta[1] - ((theta[3] - TTF)*(-(theta[3] - TTF)/theta[1])^(theta[2] - 1))/theta[1]^2 - (theta[2]*log(-(theta[3] - TTF)/theta[1])*(theta[3] - TTF)*(-(theta[3] - TTF)/theta[1])^(theta[2] - 1))/theta[1]^2),
                                       sum(- log(-(theta[3] - TTF)/theta[1])^2*(-(theta[3] - TTF)/theta[1])^theta[2] - 1/theta[2]^2),
                                       sum(1/(theta[3] - TTF) - (-(theta[3] - TTF)/theta[1])^theta[2]/(theta[3] - TTF) + (theta[2]*log(-(theta[3] - TTF)/theta[1])*(-(theta[3] - TTF)/theta[1])^(theta[2] - 1))/theta[1])),
                                     c(sum((theta[2]*(theta[3] - TTF)*(-(theta[3] - TTF)/theta[1])^(theta[2] - 2)*(theta[2] - 1))/theta[1]^3 - (theta[2]*(-(theta[3] - TTF)/theta[1])^(theta[2] - 1))/theta[1]^2),
                                       sum(1/(theta[3] - TTF) - (-(theta[3] - TTF)/theta[1])^theta[2]/(theta[3] - TTF) + (theta[2]*log(-(theta[3] - TTF)/theta[1])*(-(theta[3] - TTF)/theta[1])^(theta[2] - 1))/theta[1]),
                                       sum(- (theta[2] - 1)/(theta[3] - TTF)^2 - (theta[2]*(-(theta[3] - TTF)/theta[1])^(theta[2] - 2)*(theta[2] - 1))/theta[1]^2)))
        ll
      }
    } else{
      loglik <- function(theta){
        ll <- -sum(log(theta[2]) + (theta[2]-1)*log(TTF-theta[3]) - theta[2]*log(theta[1]) - (((TTF-theta[3])/theta[1])^theta[2])) - sum(- (((Tc-theta[3])/theta[1])^theta[2]))
        ll
        # -sum(log(theta[2]) + (theta[2]-1)*log(TTF-exp(theta[3])) - theta[2]*log(theta[1]) - (((TTF-exp(theta[3]))/theta[1])^theta[2])) - sum(- (((Tc-exp(theta[3]))/theta[1])^theta[2]))
      }
    }
  }

  if (dist=="3PWeibull" && LSQest[2] < 2) {
    positivity_v[1]<-1
    positivity_v[2]<-1
    gammaset <- LSQest[3]
    LSQest1 <- LSQest[1:2]

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(log(theta[2]) + (theta[2]-1)*log(TTF-gammaset) - theta[2]*log(theta[1]) - (((TTF-gammaset)/theta[1])^theta[2]))
      }
    } else{
      loglik <- function(theta){
        -sum(log(theta[2]) + (theta[2]-1)*log(TTF-gammaset) - theta[2]*log(theta[1]) - (((TTF-gammaset)/theta[1])^theta[2])) - sum(- (((Tc-gammaset)/theta[1])^theta[2]))
      }
    }
  }

  if (dist=="Lognormal") {
    positivity_v[2]<-1

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(-log(theta[2]) - 0.5*log(2*pi) - 0.5*(theta[2]^-2)*((log(TTF) - theta[1])^2))
      }
    } else{
      loglik <- function(theta){
        -sum(-log(theta[2]) - 0.5*log(2*pi) - 0.5*(theta[2]^-2)*((log(TTF) - theta[1])^2)) - sum(log(0.5 - 0.5*erf((2^-0.5)*(theta[2]^-1)*(log(Tc) - theta[1]))))
      }
    }
  }
  if (dist=="Normal") {
    positivity_v[2]<-1

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(-log(theta[2]) - 0.5*log(2*pi) - 0.5*(theta[2]^-2)*((TTF - theta[1])^2))
      }
    } else{
      loglik <- function(theta){
        -sum(-log(theta[2]) - 0.5*log(2*pi) - 0.5*(theta[2]^-2)*((TTF - theta[1])^2)) - sum(log(0.5 - 0.5*erf((2^-0.5)*(theta[2]^-1)*(Tc - theta[1]))))
      }
    }
  }
  if (dist=="Exponential") {
    positivity_v[1]<-1
    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(log(theta) - TTF*theta)
      }
    } else{
      loglik <- function(theta){
        -sum(log(theta) - TTF*theta) + sum( Tc*theta)
      }
    }
  }
  if (dist=="2PExponential") {
    positivity_v[2]<-1

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(-log(theta[2]) - (theta[2]^-1)*(TTF - theta[1]) - 1)
      }
    } else{
      loglik <- function(theta){
        -sum(-log(theta[2]) - (theta[2]^-1)*(TTF - theta[1]) - 1) - sum(-(theta[2])*(Tc - theta[1]) - 1)
      }
    }
  }
  if (dist=="Gamma") {
    positivity_v[1]<-1
    positivity_v[2]<-1

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(-log(gamma(theta[1])) + theta[1]*log(theta[2]) + (theta[1] - 1)*log(TTF) - theta[2]*TTF)
      }
    } else{
      loglik <- function(theta){
        -sum(-log(gamma(theta[1])) + theta[1]*log(theta[2]) + (theta[1] - 1)*log(TTF) - theta[2]*TTF) - sum(1 - Rgamma(theta[1],Tc*theta[2]))
      }
    }
  }
  if (dist=="Logistic") {
    positivity_v[2]<-1

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(((TTF - theta[1])/theta[2]) - log(theta[2]) - 2*log(1 + exp((TTF - theta[1])/theta[2])))
      }
    } else{
      loglik <- function(theta){
        -sum(((TTF - theta[1])/theta[2]) - log(theta[2]) - 2*log(1 + exp((TTF - theta[1])/theta[2]))) + sum(log(1 + exp((Tc - theta[1])/theta[2])))
      }
    }
  }
  if (dist=="Loglogistic") {
    positivity_v[2]<-1

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(((log(TTF) - theta[1])/theta[2]) - log(theta[2]) - log(TTF*((1 + exp((log(TTF) - theta[1])/theta[2]))^2)))
      }
    } else{
      loglik <- function(theta){
        -sum(((log(TTF) - theta[1])/theta[2]) - log(theta[2]) - log(TTF*((1 + exp((log(TTF) - theta[1])/theta[2]))^2))) + sum(log(1 + exp((log(Tc) - theta[1])/theta[2])))
      }
    }
  }
  if (dist=="Gumbel") {
    positivity_v[2]<-1

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(-log(theta[2]) + ((TTF - theta[1])/theta[2]) - exp((TTF - theta[1])/theta[2]))
      }
    } else{
      loglik <- function(theta){
        -sum(-log(theta[2]) + ((TTF - theta[1])/theta[2]) - exp((TTF - theta[1])/theta[2])) + sum(exp((Tc - theta[1])/theta[2]))
      }
    }
  }


  if (dist=="3PWeibull" && LSQest[2] < 2){
    MLEandvar <- MLE.var.covar.select(loglik,LSQest1)
  } else {
    MLEandvar <- MLE.var.covar.select(loglik,LSQest)
  }

  theta.hat <- MLEandvar[[1]]
  inv.fish  <- MLEandvar[[2]]
  loglik.hat <- -loglik(theta.hat)
  likeli.hat <- exp(loglik.hat)


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
      conflim_txt<-c(paste(c("Lower ",100*conf.level,"%"),collapse = ""),paste(c("Upper ",100*conf.level,"%"),collapse = ""))
    }
    if(sided == "onesidedlow"){
      if(positivity_v[i]==0){
        conflim[[i]]<-theta.hat[i] - crit2 * sqrt(inv.fish[i, i])
      } else if(positivity_v[i]==1){
        conflim[[i]]<-theta.hat[i]*exp(-crit2 * (sqrt(inv.fish[i, i])/theta.hat[i]))
      }
      conflim_txt<-paste(c("One-Sided Low ",100*conf.level,"%"),collapse = "")
    }
    if(sided == "onesidedhigh"){
      if(positivity_v[i]==0){
        conflim[[i]]<-theta.hat[i] + crit2 * sqrt(inv.fish[i, i])
      } else if(positivity_v[i]==1){
        conflim[[i]]<-theta.hat[i]*exp(crit2 * (sqrt(inv.fish[i, i])/theta.hat[i]))
      }
      conflim_txt<-paste(c("One-Sided High ",100*conf.level,"%"),collapse = "")
    }
    # fulllimset[[i]]<-c(theta.hat[i],conflim[[i]])
  }

  if (dist=="3PWeibull" && LSQest[2] < 2){
    theta.hat <- c(theta.hat,gammaset)
    conflim[[3]]<-c(NA,NA)
  }

  return(list(theta.hat,inv.fish,conflim,loglik.hat,likeli.hat))
}
