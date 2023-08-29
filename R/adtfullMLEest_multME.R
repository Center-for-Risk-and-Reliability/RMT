# MLE Accelerated Degradation Testing Estimator (with multiplicative measurement error)
# Developed by Dr. Reuel Smith, 2022

adt.full.multME.MLE <- function(data,lifedam,D0,confid=0.95,sided="twosided",multME,Tuse=293.15){
  # Load pracma library for pseudo-inverse
  library(pracma)
  library(dplyr)
  library(plotly)
  library(matrixcalc)

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

  # Pulls the unit designations from column 3 of the input data
  unitnames <- unique(data[,3])

  # Pull the time and damage data from the input
  TimeDam <- data[,1]
  Dam <- data[,2]

  # Setup positivity check vector for parameters
  # For when measurement error is a single value throughout
  if(length(multME)==1){
    positivity_v<-rep(0,dim(adtLSQ)[2]-1)
  }
  # For when measurement error is a distribution ME~ logn(b_e,sig_e)
  # Enter as a vector c(b_e,sig_e)
  if(length(multME)==2){
    positivity_v<-rep(0,dim(adtLSQ)[2])
  }

  # Will need to consider a check to see whether the data all has the same time stamp or not
  # for each set of degradation data

  # Life-Degradation (damage) models
  if(lifedam=="Linear"){
    # D = a + b*t
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"D = a + b*t"
    damage <- function(theta) {
      theta[1] + TimeDam*theta[2]
    }
    damage_mod <- function(theta,time_unit) {
      theta[1] + time_unit*theta[2]
    }
    logdamage <- function(theta) {
      log(theta[1] + TimeDam*theta[2])
    }
    lifedam_txt<-lifedam
    params_txt<-c("a","b")
    LSQest<-c(mean(adtLSQ[,1]), mean(adtLSQ[,2]), 1)
    sigparamno<-3
  }

  if(lifedam=="Exponential"){
    # D = b*exp(a*t)
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"b*exp(a*t)"
    damage <- function(theta) {
      theta[2]*exp(TimeDam*theta[1])
    }
    damage_mod <- function(theta,time_unit) {
      theta[2]*exp(time_unit*theta[1])
    }
    logdamage <- function(theta) {
      log(theta[2]) + TimeDam*theta[1]
    }
    lifedam_txt<-lifedam
    params_txt<-c("a","b")
    LSQest<-c(mean(adtLSQ[,1]), mean(adtLSQ[,2]), 1)
    positivity_v[2]<-1
    sigparamno<-3
  }

  if(lifedam=="SquareRoot"){
    # D^(1/2) = a + b*t
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"(a + b*t)^2"
    damage <- function(theta) {
      (theta[1] + TimeDam*theta[2])^2
    }
    damage_mod <- function(theta,time_unit) {
      (theta[1] + time_unit*theta[2])^2
    }
    logdamage <- function(theta) {
      2*log(theta[1] + TimeDam*theta[2])
    }
    lifedam_txt<-lifedam
    params_txt<-c("a","b")
    LSQest<-c(mean(adtLSQ[,1]), mean(adtLSQ[,2]), 1)
    sigparamno<-3
  }

  if(lifedam=="Power"){
    # D = b*(t^a)
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"b*(t^a)"
    damage <- function(theta) {
      theta[2]*(TimeDam^theta[1])
    }
    damage_mod <- function(theta,time_unit) {
      theta[2]*(time_unit^theta[1])
    }
    logdamage <- function(theta) {
      log(theta[2]) + theta[1]*log(TimeDam)
    }
    lifedam_txt<-lifedam
    params_txt<-c("a","b")
    LSQest<-c(mean(adtLSQ[,1]), mean(adtLSQ[,2]), 0.1)
    positivity_v[2]<-1
    sigparamno<-3
  }

  if(lifedam=="Logarithmic"){
    # D = a + b*ln(t)
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"a + b*ln(t)"
    damage <- function(theta) {
      theta[1] + theta[2]*log(TimeDam)
    }
    damage_mod <- function(theta,time_unit) {
      theta[1] + theta[2]*log(time_unit)
    }
    logdamage <- function(theta) {
      log(theta[1] + theta[2]*log(TimeDam))
    }
    lifedam_txt<-lifedam
    params_txt<-c("a","b")
    LSQest<-c(mean(adtLSQ[,1]), mean(adtLSQ[,2]), 1)
    sigparamno<-3
  }

  # if(lifedam=="Gompertz"){
  #   # D = a + b^(c*t)
  #   # theta[1] ~ a, theta[2] ~ b, theta[3] ~ c
  #   Dt_txt<-"a + b^(c*t)"
  #   damage <- function(theta) {
  #     theta[1] + theta[2]*log(TimeDam)
  #   }
  #   logdamage <- function(theta) {
  #     log(theta[1] + theta[2]*log(TimeDam))
  #   }
  #   lifedam_txt<-lifedam
  #   params_txt<-c("a","b")
  #   LSQest<-c(mean(adtLSQ[,1]), mean(adtLSQ[,2]), 1)
  #   sigparamno<-4
  # }

  if(lifedam=="LloydLipow"){
    # D = a - b/t
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"a + b/t"
    damage <- function(theta) {
      theta[1] - (theta[2]/TimeDam)
    }
    damage_mod <- function(theta,time_unit) {
      theta[1] - (theta[2]/time_unit)
    }
    logdamage <- function(theta) {
      log(theta[1] - (theta[2]/TimeDam))
    }
    lifedam_txt<-lifedam
    params_txt<-c("a","b")
    LSQest<-c(mean(adtLSQ[,1]), mean(adtLSQ[,2]), 1)
    sigparamno<-3
  }

  if(lifedam=="Mitsuom"){
    # D = 1/(1 + b*(t^a))
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"1/(1 + b*(t^a))"
    damage <- function(theta) {
      1/(1 + theta[2]*(TimeDam^theta[1]))
    }
    damage_mod <- function(theta,time_unit) {
      1/(1 + theta[2]*(time_unit^theta[1]))
    }
    logdamage <- function(theta) {
      -log(1 + theta[2]*(TimeDam^theta[1]))
    }
    lifedam_txt<-lifedam
    params_txt<-c("a","b")
    LSQest<-c(mean(adtLSQ[,1]), mean(adtLSQ[,2]), 1)
    sigparamno<-3
  }

  if(lifedam=="Hamada"){
    # D = 1/(1 + beta1*(t*exp(beta3*11605*(1/Tu - 1/Ti)))^beta2)
    # theta[1] ~ beta1, theta[2] ~ beta2, theta[3] ~ beta3
    TimeDam
    Dt_txt<-"1/(1 + \U03B2_1*(t*exp(\U03B2_3*11605*(1/Tuse - 1/Temp)))^\U03B2_2)"
    damage <- function(theta) {
      1/(1 + theta[1]*(TimeDam*exp(theta[3]*11605*(1/Tuse - 1/TempDam)))^theta[2])
    }
    logdamage <- function(theta) {
      -log(1 + theta[1]*(TimeDam*exp(theta[3]*11605*(1/Tuse - 1/TempDam)))^theta[2])
    }
    lifedam_txt<-lifedam
    params_txt<-c("\U03B2_1","\U03B2_2","\U03B2_3")
    LSQest<-c(mean(adtLSQ[,1]), mean(adtLSQ[,2]), mean(adtLSQ[,3]), 1)
    positivity_v[1]<-1
    positivity_v[2]<-1
    sigparamno<-4
  }

  # Obtain the model values for the input data based on the LSQ estimates for each unit
  # if we are dealing with a measurement error distribution


  # Then compute the MLE estimate based on the lognormal distribution (default for multiplicative measurement error)
  # for the damage model
  dist_txt<-"Lognormal"

  if(length(multME)==1){
    positivity_v[sigparamno]<-1
    loglik <- function(theta){
      -sum(-log(theta[sigparamno]) - 0.5*log(2*pi) - log(multME) - log(Dam) - 0.5*(theta[sigparamno]^-2)*((log(multME) + log(Dam) - logdamage(theta))^2))
    }
    distparam_txt<-"\U03C3_t"
  }
  if(length(multME)==2){
    # Obtain the model values for the input data based on the LSQ estimates for each unit
    # if we are dealing with a measurement error distribution
    for(i in 1:length(unitnames)){
      if(i==1){
        D_model<-damage_mod(adtLSQ[i,],data[,1][which(data[,3]==unitnames[i])])
      } else{
        D_model<-c(D_model,damage_mod(adtLSQ[i,],data[,1][which(data[,3]==unitnames[i])]))
      }
    }
    # Compute model error initials based on LSQ
    LSQmodelerror<-sort((data[,2]/D_model)*rlnorm(length(data[,2]),multME[1],multME[2]))
    F1<-plotposit.blom(rankcalc(LSQmodelerror),LSQmodelerror)[,2]
    LSQmodelerrorparams<-c(probplotparam.logn(LSQmodelerror,F1)[[3]])
    positivity_v[sigparamno + 1]<-1
    loglik <- function(theta){
      -sum(-log(theta[sigparamno+1]^2 + multME[2]^2) - 0.5*log(2*pi) - log(Dam) + logdamage(theta) - 0.5*(1/(theta[sigparamno+1]^2 + multME[2]^2))*((log(Dam) - logdamage(theta) - (theta[sigparamno] - multME[1]))^2))
    }
    distparam_txt<-c("b_m","\U03C3_m")
    # adt.full.MLE(data,lifedam,"Lognormal",D0,confid,sided)[[1]]
    # LSQest<-c(LSQest[1:length(LSQest)-1],LSQmodelerrorparams)
    LSQest<-c(adt.full.MLE(data,lifedam,"Lognormal",D0,confid,sided)[[1]][1:2],LSQmodelerrorparams)
  }

  MLEandvar <- MLE.var.covar.select(loglik,LSQest)
  theta.hat <- MLEandvar[[1]]
  inv.fish  <- MLEandvar[[2]]

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

  params_txt<-c(params_txt,distparam_txt)

  # Produce some output text that summariZes the results
  cat(c("Maximum-Likelihood estimates for the ",lifedam_txt,"-",dist_txt," Degradation model.\n\nD(t) = ",Dt_txt,"\n\n"),sep = "")
  print(matrix(unlist(fulllimset), nrow = length(unlist(fulllimset))/length(LSQest), ncol = length(LSQest), byrow = FALSE,dimnames = list(c("Degradation Model Parameters Mean",conflim_txt),params_txt)))
  return(list(theta.hat,inv.fish))
}
