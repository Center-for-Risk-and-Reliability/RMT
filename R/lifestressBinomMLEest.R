# Maximum Likelihood Binomial Accelerated Life-Stress Estimator
# Developed by Mohammad Modarres and Reuel Smith, 2022

lifestress.Binom.MLEest <- function(data,interact_stress,weight0,confid,sided){
  #Load pracma library for pseudoinverse
  library(pracma)

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

  # Check to see if two or one sided entry exists
  if(missing(sided)){
    sided <- "twosided"
  } else {
    sided <- sided
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

  # Isolate the stress differentials
  # Dstressi <- data[,3:dim(data)[2]-weight0]
  # names(Dstressi) <- NULL

  # Pull specific data from data matrix
  ni <- data[,1]
  ki <- data[,2]
  if(weight0==1){
    wi <- data[,dim(data)[2]]
  }

  Dstressi <- data[,3:(dim(data)[2]-weight0)]
  names(Dstressi) <- NULL


  Dstress_txt<-paste("\U03B8", 1:(dim(data)[2]-2-weight0), sep = "_")

  # Preprocess stress matrices (will add to them if interaction is evaluated later)
  if(is.null(dim(Dstressi))==TRUE){
    # One Stress Differential case
    Dstressi_v<-Dstressi
    Dstress_m<-matrix(Dstressi_v,nrow = dim(data)[1],ncol = 1,byrow = FALSE)

  } else{
    # Multi Stress Differential
    for(i in 1:(dim(data)[2]-2-weight0)){
      if(i==1){
        Dstressi_v<-Dstressi[,i]
      } else {
        Dstressi_v<-c(Dstressi_v,Dstressi[,i])
      }
    }

    if(interact_stress==0){
      # No Stress dependency check
      Dstress_m<-matrix(Dstressi_v,nrow = dim(data)[1],ncol = dim(Dstressi)[2],byrow = FALSE)
    }
    if(interact_stress==1){
      # Stress dependency check
      # Permutations of stress pairs
      Dstresspairs <- combn(dim(Dstressi)[2],2)
      for(i in 1:dim(Dstresspairs)[2]){
        if(i==1){
          Dstresspair_v <- Dstressi[, Dstresspairs[1,i]]*Dstressi[, Dstresspairs[2,i]]
          Dstressjoint_txt<-paste("\U03B8", paste(Dstresspairs[,i],collapse=""), sep = "_")
        } else{
          Dstresspair_v <- c(Dstresspair_v,Dstressi[, Dstresspairs[1,i]]*Dstressi[, Dstresspairs[2,i]])
          Dstressjoint_txt<-c(Dstressjoint_txt,paste("\U03B8", paste(Dstresspairs[,i],collapse=""), sep = "_"))
        }
      }
      Dstress_m<-matrix(c(Dstressi_v,Dstresspair_v),nrow = dim(data)[1],ncol = dim(Dstressi)[2]+dim(Dstresspairs)[2],byrow = FALSE)
      Dstress_txt<-c(Dstress_txt,Dstressjoint_txt)
    }
  }

  # Initiate LSQ estimate from raw data
  LSQest <- lifestress.Binom.LSQest(data,interact_stress,weight0)[[1]]

  # If p0 initial is negative, reset to a workable initial
  # Note however that this p_0,0 has to be between 0 and 1 and it also as to produce a non
  # non-NaN number for log(1 - (1 - p_0)*exp(-sum(c(Dstress_m%*%theta[2:length(LSQest)]))))
  # or rather p_0 > 1 - exp(sum(c(Dstress_m%*%theta[2:length(LSQest)])))
  # if(LSQest[1]<0){
  #   p_0_low <- max(1 - exp((c(Dstress_m%*%LSQest[2:length(LSQest)]))))
  #   LSQest[1]<-log(-log(.5*(1-p_0_low)))
  # }

  # Set p_0 bounds between 0 and 1
  LSQest[1]<-log(-log(LSQest[1]))
  for(i in 2:length(LSQest)){
    if(LSQest[i]<10e-10){
      LSQest[i] <- max(LSQest[2:length(LSQest)])
    }
  }
  # Set theta bounds to positive
  LSQest[2:length(LSQest)] <- log(LSQest[2:length(LSQest)])


  # Setup positivity check vector for parameters
  positivity_v<-rep(0,length(LSQest))
  # positivity_v[2:length(LSQest)] <- 1

  # If p0 initial is negative, reset to a workable initial

  # # Initialize life-stress parameter estimates for theta
  ls_txt<-"Binomial Acceleration"
  pdf_txt<-"Pr(k) = (n|k) \U2219 [1 - (1 - p_0) \U2219 exp(-\U03A3_(j=1)^m \U03B8_j \U0394S_j)]^k \U2219 [(1 - p_0) \U2219 exp(-\U03A3_(j=1)^m \U03B8_j \U0394S_j)]^(n - k)"
  distparam_txt<-"p_0"
  params_txt<-c(distparam_txt,Dstress_txt)

  # Fit to log-likelihood distributions
  # positivity_v[2:length(LSQest)]<-1

  if (weight0==0){
    loglik <- function(theta){
      # -sum(lfactorial(ni)) + sum(lfactorial(ki)) + sum(lfactorial(ni-ki)) - sum(ki*log(1 - (1 - theta[1])*exp(-c(Dstress_m%*%theta[2:length(LSQest)])))) - log(1 - theta[1])*sum(ni - ki) + sum((ni - ki)*c(Dstress_m%*%theta[2:length(LSQest)]))
      # -sum(lfactorial(ni)) + sum(lfactorial(ki)) + sum(lfactorial(ni-ki)) - sum(ki*log(1 - (1 - exp(-exp(theta[1])))*exp(-c(Dstress_m%*%theta[2:length(LSQest)])))) - log(1 - exp(-exp(theta[1])))*sum(ni - ki) + sum((ni - ki)*c(Dstress_m%*%theta[2:length(LSQest)]))
      -sum(lfactorial(ni)) + sum(lfactorial(ki)) + sum(lfactorial(ni-ki)) - sum(ki*log(1 - (1 - exp(-exp(theta[1])))*exp(-c(Dstress_m%*%exp(theta[2:length(LSQest)]))))) - log(1 - exp(-exp(theta[1])))*sum(ni - ki) + sum((ni - ki)*c(Dstress_m%*%exp(theta[2:length(LSQest)])))
    }
  } else{
    loglik <- function(theta){
      # -sum(wi*lfactorial(ni)) + sum(wi*lfactorial(ki)) + sum(wi*lfactorial(ni-ki)) - sum(wi*ki*log(1 - (1 - theta[1])*exp(-c(Dstress_m%*%theta[2:length(LSQest)])))) - log(1 - theta[1])*sum(wi*(ni - ki)) + sum(wi*(ni - ki)*c(Dstress_m%*%theta[2:length(LSQest)]))
      # -sum(wi*lfactorial(ni)) + sum(wi*lfactorial(ki)) + sum(wi*lfactorial(ni-ki)) - sum(wi*ki*log(1 - (1 - exp(-exp(theta[1])))*exp(-c(Dstress_m%*%theta[2:length(LSQest)])))) - log(1 - exp(-exp(theta[1])))*sum(wi*(ni - ki)) + sum(wi*(ni - ki)*c(Dstress_m%*%theta[2:length(LSQest)]))
      -sum(wi*lfactorial(ni)) + sum(wi*lfactorial(ki)) + sum(wi*lfactorial(ni-ki)) - sum(wi*ki*log(1 - (1 - exp(-exp(theta[1])))*exp(-c(Dstress_m%*%exp(theta[2:length(LSQest)]))))) - log(1 - exp(-exp(theta[1])))*sum(wi*(ni - ki)) + sum(wi*(ni - ki)*c(Dstress_m%*%exp(theta[2:length(LSQest)])))
    }
  }
  # return(list(loglik,LSQest))

  MLEandvar <- MLE.var.covar.select(loglik,LSQest)
  theta.hat <- MLEandvar[[1]]
  inv.fish  <- MLEandvar[[2]]

  # theta.hat <- out$estimate
  crit <- qnorm((1 + conf.level)/2)
  crit2 <- qnorm(conf.level)
  conflim<-vector(mode = "list", length = length(LSQest))
  fulllimset<-vector(mode = "list", length = length(LSQest))

  for(i in 1:length(LSQest)){
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
    fulllimset[[i]]<-c(theta.hat[i],conflim[[i]])
  }

  theta.hat.res<-theta.hat
  theta.hat.res[1]<-exp(-exp(theta.hat[1]))
  theta.hat.res[2:length(LSQest)]<-exp(theta.hat[2:length(LSQest)])
  conflim.res<-conflim
  conflim.res[[1]][1]<-exp(-exp(conflim[[1]][2]))
  conflim.res[[1]][2]<-exp(-exp(conflim[[1]][1]))
  fulllimset.res<-fulllimset
  fulllimset.res[[2]]<-exp(fulllimset[[2]])
  fulllimset.res[[1]][1]<-exp(-exp(fulllimset[[1]][1]))
  fulllimset.res[[1]][2]<-exp(-exp(fulllimset[[1]][3]))
  fulllimset.res[[1]][3]<-exp(-exp(fulllimset[[1]][2]))



  # Produce some output text that summariZes the results
  cat(c("Maximum-Likelihood estimates for the ",ls_txt," Life-Stress model.\n\n"),sep = "")
  print(matrix(unlist(fulllimset.res), nrow = length(unlist(fulllimset.res))/length(LSQest), ncol = length(LSQest), byrow = FALSE,dimnames = list(c("Life-Stress Parameters Mean",conflim_txt),params_txt)))
  # print(matrix(unlist(fulllimset), nrow = length(unlist(fulllimset))/length(LSQest), ncol = length(LSQest), byrow = FALSE,dimnames = list(c("Life-Stress Parameters Mean",conflim_txt),params_txt)))
  # return(list(theta.hat,theta.hat.res,inv.fish,conflim,conflim.res))
  return(list(theta.hat.res,inv.fish,conflim.res,fulllimset.res))

  # return(list(theta.hat,inv.fish,conflim))
  # return(list(Dstressi,Dstress_m,Dstress_txt, LSQest,ni,ki,loglik))
}
