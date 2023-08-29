# MLE Life-Model Estimator
# Developed by Dr. Reuel Smith, 2021-2022

life.MLEest <- function(LSQest,dist,TTF,Tc,confid,sided){
  #Load pracma library for erf
  library(pracma)
  library(matrixcalc)

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

  # Check to see if two or one sided entry exists
  if(missing(sided)){
    sided <- "twosided"
  } else {
    sided <- sided
  }

  # Setup positivity check vector for parameters
  positivity_v<-rep(0,length(LSQest))

  # Fit to log-likelihood distributions
  if (dist=="Weibull") {
    positivity_v[1]<-1
    positivity_v[2]<-1

    if(missing(Tc)){
      loglik <- function(theta){
        -sum(log(theta[2]) + (theta[2]-1)*log(TTF) - theta[2]*log(theta[1]) - ((TTF/theta[1])^theta[2]))
      }
    } else{
      loglik <- function(theta){
        -sum(log(theta[2]) + (theta[2]-1)*log(TTF) - theta[2]*log(theta[1]) - ((TTF/theta[1])^theta[2])) - sum(- ((Tc/theta[1])^theta[2]))
      }
    }
  }
  if (dist=="Lognormal") {
    positivity_v[2]<-1

    if(missing(Tc)){
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

    if(missing(Tc)){
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
    if(missing(Tc)){
      loglik <- function(theta){
        -sum(-loglifeF(theta) - TTF*theta)
      }
    } else{
      loglik <- function(theta){
        -sum(-loglifeF(theta) - TTF*theta) - sum(-Tc*theta)
      }
    }
  }
  if (dist=="2PExponential") {
    positivity_v[2]<-1

    if(missing(Tc)){
      loglik <- function(theta){
        -sum(-log(theta[2]) - (theta[2]^-1)*(TTF - theta[1]) - 1)
      }
    } else{
      loglik <- function(theta){
        -sum(-log(theta[2]) - (theta[2]^-1)*(TTF - theta[1]) - 1) - sum(-(theta[2])*(Tc - theta[1]) - 1)
      }
    }
  }

  out <- nlm(loglik, theta <- LSQest, hessian=TRUE)
  fish <- out$hessian

  # Check singularity.  If it is singular, then use pseudoinverse
  if (is.singular.matrix(fish)==TRUE){
    inv.fish <- pinv(fish)
  } else {
    inv.fish <- solve(fish)
  }

  theta.hat <- out$estimate
  crit <- qnorm((1 + conf.level)/2)
  crit2 <- qnorm(conf.level)
  conflim<-vector(mode = "list", length = length(LSQest))

  for(i in 1:length(LSQest)){
    if(sided == "twosided"){
      if(positivity_v[i]==0){
        conflim[[i]]<-theta.hat[i] + c(-1, 1) * crit * sqrt(inv.fish[i, i])
      } else if(positivity_v[i]==1){
        conflim[[i]]<-theta.hat[i]*exp(c(-1, 1) * crit * (sqrt(inv.fish[i, i])/theta.hat[i]))
      }
    }
    if(sided == "onesidedlow"){
      if(positivity_v[i]==0){
        conflim[[i]]<-theta.hat[i] - crit2 * sqrt(inv.fish[i, i])
      } else if(positivity_v[i]==1){
        conflim[[i]]<-theta.hat[i]*exp(-crit2 * (sqrt(inv.fish[i, i])/theta.hat[i]))
      }
    }
    if(sided == "onesidedhigh"){
      if(positivity_v[i]==0){
        conflim[[i]]<-theta.hat[i] + crit2 * sqrt(inv.fish[i, i])
      } else if(positivity_v[i]==1){
        conflim[[i]]<-theta.hat[i]*exp(crit2 * (sqrt(inv.fish[i, i])/theta.hat[i]))
      }
    }
  }
  return(list(theta.hat,conflim))
}
