# MLE PH-Weibull Life-Stress Estimator
# Developed by Dr. Reuel Smith, 2022

lifestress.PHWbl.MLEest <- function(LSQest,data,interact_stress,confid,sided){
  #Load pracma library for erf
  library(pracma)
  library(matrixcalc)
  library(ucminf)

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

  # Sets the default computation for stresses to simply single evaluation
  # and not check for interactions
  if(missing(interact_stress)){
    interact_stress<-0
  }

  # Pull specific data from data matrix
  TTF<-data[which(data[,2]==1),1]
  TTS<-data[which(data[,2]==0),1]
  SF<-data[which(data[,2]==1),3:dim(data)[2]]
  Sc<-data[which(data[,2]==0),3:dim(data)[2]]

  stress_txt<-paste("\U03B8", 0:dim(SF)[2], sep = "_")

  # Preprocess stress matrices (will add to them if interaction is evaluated later)
  for(i in 1:(1+dim(SF)[2])){
    if(i==1){
      SF_v<-rep(1,length(TTF))
      if(length(TTS)>0){
        Sc_v<-rep(1,length(TTS))
      }
    } else {
      SF_v<-c(SF_v,SF[,i-1])
      if(length(TTS)>0){
        Sc_v<-c(Sc_v,Sc[,i-1])
      }
    }
  }
  if(interact_stress==0){
    # No Stress dependency check
    SF_m<-matrix(SF_v,nrow = length(TTF),ncol = 1+dim(SF)[2],byrow = FALSE)
    if(length(TTS)>0){
      Sc_m<-matrix(Sc_v,nrow = length(TTS),ncol = 1+dim(Sc)[2],byrow = FALSE)
    }
  }
  if(interact_stress==1){
    # Stress dependency check
    # Permutations of stress pairs
    Spairs <- combn(dim(SF)[2],2)
    for(i in 1:dim(Spairs)[2]){
      if(i==1){
        SFpair_v <- SF[, Spairs[1,i]]*SF[, Spairs[2,i]]
        if(length(TTS)>0){
          Scpair_v <- Sc[, Spairs[1,i]]*Sc[, Spairs[2,i]]
        }
        stressjoint_txt<-paste("\U03B8", paste(Spairs[,i],collapse=""), sep = "_")
      } else{
        SFpair_v <- c(SFpair_v,S[, Spairs[1,i]]*S[, Spairs[2,i]])
        if(length(TTS)>0){
          Scpair_v <- c(Scpair_v,Sc[, Spairs[1,i]]*Sc[, Spairs[2,i]])
        }
        stressjoint_txt<-c(stressjoint_txt,paste("\U03B8", paste(Spairs[,i],collapse=""), sep = "_"))
      }
    }
    SF_m<-matrix(c(SF_v,SFpair_v),nrow = length(TTF),ncol = 1+dim(SF)[2]+dim(Spairs)[2],byrow = FALSE)
    if(length(TTS)>0){
      Sc_m<-matrix(c(Sc_v,Scpair_v),nrow = length(TTS),ncol = 1+dim(Sc)[2]+dim(Spairs)[2],byrow = FALSE)
    }
    stress_txt<-c(stress_txt,stressjoint_txt)
  }

  # Setup positivity check vector for parameters
  positivity_v<-rep(0,length(LSQest))

  # Initialize life-stress parameter estimates for theta
  ls_txt<-"Parametric Proportional Hazard"
  pdf_txt<-"\U03B2\U2219 t^(\U03B2-1)\U2219 t\U2219 exp[\U03A3_(j=0)^m \U03B8_j S_j - t^\U03B2 exp(\U03A3_(j=0)^m \U03B8_j S_j)]"
  dist_txt<-"Weibull"
  distparam_txt<-"\U03B2"
  params_txt<-c(distparam_txt,stress_txt)

  # Fit to log-likelihood distributions
  positivity_v[1]<-1

  if(length(TTS)==0){
    loglik <- function(theta){
      -sum(log(theta[1]) + (theta[1]-1)*log(TTF) + c(SF_m%*%theta[2:length(LSQest)]) - (TTF^theta[1])*exp(c(SF_m%*%theta[2:length(LSQest)])))
    }
  } else{
    loglik <- function(theta){
      -sum(log(theta[1]) + (theta[1]-1)*log(TTF) + c(SF_m%*%theta[2:length(LSQest)]) - (TTF^theta[1])*exp(c(SF_m%*%theta[2:length(LSQest)]))) - sum((TTS^theta[1])*exp(c(Sc_m%*%theta[2:length(LSQest)])))
    }
  }


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

  # Produce some output text that summariZes the results
  cat(c("Maximum-Likelihood estimates for the ",ls_txt,"-",dist_txt," Life-Stress model.\n\n"),sep = "")
  print(matrix(unlist(fulllimset), nrow = length(unlist(fulllimset))/length(LSQest), ncol = length(LSQest), byrow = FALSE,dimnames = list(c("Life-Stress Parameters Mean",conflim_txt),params_txt)))
  return(list(theta.hat,inv.fish,conflim))
}
