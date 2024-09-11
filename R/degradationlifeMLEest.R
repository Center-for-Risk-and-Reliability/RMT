# Maximum Likelihood Accelerated Degradation Testing Estimator
# Developed by Dr. Reuel Smith, 2021-2022

degradationlife.MLEest <- function(data,dl,dist="Lognormal",pp="Blom",D0,correl=0,ME=0,modelstress=NULL,Tuse=NULL,confid=0.95,sided="twosided",xlabel=NULL,ylabel=NULL){
  # Load pracma library for pseudo-inverse
  library(pracma)
  library(dplyr)
  library(plotly)
  library(matrixcalc)
  library(ucminf)
  library(MASS)
  library(mclust)
  library(mvtnorm)
  library(ggplot2)

  # Set confidence
  conf.level <- confid

  # Check the condition if the Hamada damage model is chosen and whether or not
  # a use temperature has been identified or not.  If not then room temperature
  # 293.15 K will be used by default.
  if(dl=="Hamada" & is.null(Tuse)==TRUE){
    Tuse <- 293.15
  }

  # Start with the pre-processing of the data in which you obtain an initial parameter
  # estimate based on the curve fit of the data.  Also used to obtain initial Sigma and
  # bounds for calculating the integral based likelihood
  # NEW - Now consider model stress for parameters and include that as output
  if(is.null(modelstress) == TRUE){
    adtLSQ<-degradationlife.LSQest(data,dl,dist,pp,D0,modelstress,xlabel,ylabel,Tuse)[[1]]
  }
  if(is.null(modelstress) == FALSE){
    adtLSQ<-degradationlife.LSQest(data,dl,dist,pp,D0,modelstress,xlabel,ylabel,Tuse)[[1]]
    adtLSQ2<-degradationlife.LSQest(data,dl,dist,pp,D0,modelstress,xlabel,ylabel,Tuse)[[2]]
  }

  # UPDATE (2/12/2024) Now check to see if any time values (column 1) are zero with respect to
  # Power, Logarithmic, Lloyd-Lipow, and Mitsuom degradation life models (Add more where ln(L) is involved)
  data0 <- data
  if((dl=="Power" || dl=="Logarithmic" || dl=="LloydLipow" || dl=="Mitsuom") && min(data[,1])==0){
    # Nix time data at zero
    data <- data0[which(data0[,1]!=0),]
  }

  # if(dl=="Hamada"){
  #   TempDam <- data[,4]
  #   adtLSQ<-adt.full.LSQ(data,dl,D0,Tuse)[[1]]
  # } else {
  #   adtLSQ<-adt.full.LSQ(data,dl,D0)[[1]]
  # }

  # Pulls the unit designations from column 3 of the input data
  unitnames <- unique(data[,3])

  # # Pulls the stresses from column 4 of the input data
  # unitstress <- unique(data[,4])
  # Pulls stress from columns 4 and up
  stresscount<-dim(data)[2]-3
  if(stresscount>1){
    stressvals<-vector(mode = "list", length = stresscount)
  }

  for(i in 1:length(unitnames)){
    stressgroup<-which(data[,3]==unitnames[i])
    stressvals0<-data[stressgroup[1],4:dim(data)[2]]
    names(stressvals0)<-NULL
    if(stresscount==1){
      if(i==1){
        stressvals<-stressvals0
      } else if(i>1){
        stressvals<-c(stressvals,stressvals0)
      }
    } else {
      for(j in 1:stresscount){
        if(i==1){
          stressvals[[j]]<-stressvals0[j]
        } else if(i>1){
          stressvals[[j]]<-c(stressvals[[j]],stressvals0[j])
        }
      }
    }
  }

  # Pull the time and damage data from the input
  TimeFULL <- data[,1]
  DegradationFULL <- data[,2]
  StresSULL <- data[,4:dim(data)[2]]
  N <- length(data[,1])


  # NEW PLAN 2/26/2024 (Happy 80th Daddy ♥♥♥)
  # 1. We will want to do the MLE for each unit based on dl regardless of the modelstress
  # 2. If additive error is selected (=1) then we do the full likelihood with parameter covariance
  # with normal distribution fit. Use mvn fit to initialize that part and numerical analysis for the rest.
  # No error (=0) means we compute this without correlation.
  # 3. If multiplicative error is selected (=2) then we do this with a lognormal distribution and the multiplicative
  # model error definition.

  # Initialize parameter-stress relation if modelstress is not NULL (default)
  # First tick off which modelstress type it is based on the name (parameter-based - "1" or AF-based - "2")
  if(is.null(modelstress) == FALSE){
    if(modelstress == "Linear" || modelstress == "Exponential" || modelstress == "Exponential2" ||
       modelstress == "Arrhenius" || modelstress == "Eyring" || modelstress == "Eyring2" || modelstress == "Eyring3" || modelstress == "Eyring4" ||
       modelstress=="Power" || modelstress=="InversePower" || modelstress=="InversePower2" || modelstress=="Logarithmic" ||
       modelstress=="MultiStress" || modelstress=="TempHumidity" || modelstress=="TempNonthermal"){
      modelstresstype <- 1
    }
    if(modelstress == "LinearAF" || modelstress == "ExponentialAF" || modelstress == "ExponentialAF2" ||
       modelstress == "ArrheniusAF" || modelstress == "EyringAF" || modelstress == "EyringAF2" || modelstress == "EyringAF3" || modelstress == "EyringAF4"  ||
       modelstress=="PowerAF" || modelstress=="InversePowerAF" || modelstress=="InversePowerAF2" || modelstress=="LogarithmicAF" ||
       modelstress=="MultiStressAF" || modelstress=="TempHumidityAF" || modelstress=="TempNonthermalAF"){
      modelstresstype <- 2
    }
  } else{
    modelstresstype <- 0
  }

  # Establish the model stress relations for parameters


  if (is.null(modelstress) == FALSE && modelstress=="Linear"){
    # lsparams[1] - parameter a_0, lsparams[2] - parameter b_0
    LSQest3 <- adtLSQ2[1:2]
    paramstress <- function(theta,S) {
      theta[2] + S*theta[1]
    }
    logparamstress <- function(theta,S) {
      log(theta[2] + S*theta[1])
    }
    # Writeup for the output text
    ps_txt<-modelstress
    params_txt3<-c("a_0","b_0")
    params_txt2<-"(b_0 + S*a_0)"
    logparam_txt<-"ln(b_0 + S*a_0)"
  }

  if (is.null(modelstress) == FALSE && modelstress=="Exponential"){
    # lsparams[1] - parameter a_0, lsparams[2] - parameter b_0, lsparams[3] - R^2
    # shift b to log b
    # positivity_v[2]<-1
    LSQest3 <- c(adtLSQ2[1],log(adtLSQ2[2]))

    paramstress <- function(theta,S) {
      exp(theta[2])*exp(S*theta[1])
    }
    logparamstress <- function(theta,S) {
      theta[2] + S*theta[1]
    }
    # Writeup for the output text
    ps_txt<-modelstress
    params_txt3<-c("a_0","b_0")
    param_txt2<-"b_0*exp(a_0*S)"
    logparam_txt<-"(log(b_0) + a_0*S)"
  }
  if (is.null(modelstress) == FALSE && modelstress=="Exponential2"){
    # lsparams[1] - parameter a_0, lsparams[2] - parameter b_0
    # shift b to log b
    # positivity_v[2]<-1
    LSQest3 <- c(adtLSQ2[1],log(adtLSQ2[2]))

    paramstress <- function(theta,S) {
      exp(theta[2])*exp(theta[1]/S)
    }
    logparamstress <- function(theta,S) {
      theta[2] + theta[1]/S
    }
    # Writeup for the output text
    ps_txt<-modelstress
    params_txt3<-c("a_0","b_0")
    param_txt2<-"b_0*exp(a_0/S)"
    logparam_txt<-"(log(b_0) + a_0/S)"
  }
  if (is.null(modelstress) == FALSE && modelstress=="Arrhenius"){
    # psparams[1] - parameter Ea, psparams[2] - parameter b
    # Temperature HAS to be in Kelvin for this to work
    # shift b to log b
    # positivity_v[2]<-1
    LSQest3 <- c(adtLSQ2[1],log(adtLSQ2[2]))

    K<-8.617385e-5
    paramstress <- function(theta,S) {
      exp(theta[2])*exp(theta[1]/(K*S))
    }
    logparamstress <- function(theta,S) {
      theta[2] + theta[1]/(K*S)
    }
    # Writeup for the output text
    ps_txt<-modelstress
    params_txt3<-c("E_a_0","b_0")
    param_txt2<-"b_0*exp(E_a_0)/(K*S))"
    logparam_txt<-"(log(b_0) + (E_a_0/(K*S)))"
  }
  if (is.null(modelstress) == FALSE && modelstress=="Eyring"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # First redefine b parameter as logb
    LSQest3 <-  c(adtLSQ2[1],log(adtLSQ2[2]))
    paramstress <- function(theta,S) {
      (exp(theta[2])/S)*exp(theta[1]/S)
    }
    logparamstress <- function(theta,S) {
      theta[2] - log(S) + theta[1]/S
    }
    # Writeup for the output text
    ps_txt<-modelstress
    params_txt3<-c("a_0","b_0")
    param_txt2<-"(b_0/S)*exp(a_0/S)"
    logparam_txt<-"(log(b_0) - log(S) + (a_0/S))"
  }
  if (is.null(modelstress) == FALSE && modelstress=="Eyring2"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    LSQest3 <- adtLSQ2[1:2]
    paramstress <- function(theta,S) {
      (1/S)*exp(-(theta[1] - (theta[2]/S)))
    }
    logparamstress <- function(theta,S) {
      -log(S) - theta[1] + theta[2]/S
    }
    # Writeup for the output text
    ps_txt<-"Eyring (Type-2)"
    params_txt3<-c("a_0","b_0")
    param_txt2<-"(1/S)*exp(-(a_0 - (b_0/S)))"
    logparam_txt<-"(-log(S) - a_0 + (b_0/S))"
  }
  if (is.null(modelstress) == FALSE && modelstress=="Power"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # First redefine b parameter as logb
    LSQest3 <- c(adtLSQ2[1],log(adtLSQ2[2]))
    # positivity_v[2]<-1

    paramstress <- function(theta,S) {
      exp(theta[2])*(S^theta[1])
    }
    logparamstress <- function(theta,S) {
      theta[2] + theta[1]*log(S)
    }
    # Writeup for the output text
    ps_txt<-modelstress
    params_txt3<-c("a_0","b_0")
    param_txt2<-"b_0*(S^a_0)"
    logparam_txt<-"(ln(b_0) + a_0ln(S))"
  }
  if (is.null(modelstress) == FALSE && modelstress=="InversePower"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # First redefine b parameter as logb
    LSQest3 <- c(adtLSQ2[1],log(adtLSQ2[2]))
    # positivity_v[2]<-1

    paramstress <- function(theta,S) {
      exp(theta[2])*(S^-theta[1])
    }
    logparamstress <- function(theta,S) {
      theta[2] - theta[1]*log(S)
    }
    # Writeup for the output text
    ps_txt<-"Inverse Power"
    params_txt3<-c("a_0","b_0")
    param_txt2<-"b_0*(S^-a_0)"
    logparam_txt<-"(ln(b_0) - a_0ln(S))"
  }
  if (is.null(modelstress) == FALSE && modelstress=="InversePower2"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # First redefine b parameter as logb
    LSQest3 <- c(adtLSQ2[1],log(adtLSQ2[2]))
    # positivity_v[2]<-1

    paramstress <- function(theta,S) {
      1/(exp(theta[2])*(S^theta[1]))
    }
    logparamstress <- function(theta,S) {
      -theta[2] - theta[1]*log(S)
    }
    # Writeup for the output text
    ps_txt<-"Inverse Power"
    params_txt3<-c("a_0","b_0")
    param_txt2<-"1/[b_0*(S^a_0)]"
    logparam_txt<-"(-ln(b_0) - a_0ln(S))"
  }

  if (is.null(modelstress) == FALSE && modelstress=="Logarithmic"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    LSQest3 <- adtLSQ2[1:2]

    paramstress <- function(theta,S) {
      theta[1]*log(S) + theta[2]
    }
    logparamstress <- function(theta,S) {
      log(theta[1]*log(S) + theta[2])
    }

    # Writeup for the output text
    ps_txt<-modelstress
    params_txt3<-c("a_0","b_0")
    param_txt2<-"(b_0 + a_0*ln(S))"
    logparam_txt<-"ln(b_0 + a_0*ln(S))"
  }

  if (is.null(modelstress) == FALSE && modelstress=="MultiStress"){
    # lsparams - c(ao,a1,...,an), S - c(1,S1,...,Sn)
    LSQest3 <- adtLSQ2[1:dim(S)[2]]

    paramstress <- function(theta,S) {
      if(dim(S)[2]==2){
        eqn1<-exp(theta[1]+theta[2]*S[,1]+theta[3]*S[,2])
      }
      if(dim(S)[2]==3){
        eqn1<-exp(theta[1]+theta[2]*S[,1]+theta[3]*S[,2]+theta[4]*S[,3])
      }
      if(dim(S)[2]==4){
        eqn1<-exp(theta[1]+theta[2]*S[,1]+theta[3]*S[,2]+theta[4]*S[,3]+theta[5]*S[,4])
      }
      if(dim(S)[2]==5){
        eqn1<-exp(theta[1]+theta[2]*S[,1]+theta[3]*S[,2]+theta[4]*S[,3]+theta[5]*S[,4]+theta[4]*S[,5])
      }
      if(dim(S)[2]==6){
        eqn1<-exp(theta[1]+theta[2]*S[,1]+theta[3]*S[,2]+theta[4]*S[,3]+theta[5]*S[,4]+theta[4]*S[,5]+theta[5]*S[,6])
      }
      if(dim(S)[2]==7){
        eqn1<-exp(theta[1]+theta[2]*S[,1]+theta[3]*S[,2]+theta[4]*S[,3]+theta[5]*S[,4]+theta[4]*S[,5]+theta[5]*S[,6]+theta[6]*S[,7])
      }
      if(dim(S)[2]==8){
        eqn1<-exp(theta[1]+theta[2]*S[,1]+theta[3]*S[,2]+theta[4]*S[,3]+theta[5]*S[,4]+theta[4]*S[,5]+theta[5]*S[,6]+theta[6]*S[,7]+theta[7]*S[,8])
      }
      if(dim(S)[2]==9){
        eqn1<-exp(theta[1]+theta[2]*S[,1]+theta[3]*S[,2]+theta[4]*S[,3]+theta[5]*S[,4]+theta[4]*S[,5]+theta[5]*S[,6]+theta[6]*S[,7]+theta[7]*S[,8]+theta[8]*S[,9])
      }
      if(dim(S)[2]==10){
        eqn1<-exp(theta[1]+theta[2]*S[,1]+theta[3]*S[,2]+theta[4]*S[,3]+theta[5]*S[,4]+theta[4]*S[,5]+theta[5]*S[,6]+theta[6]*S[,7]+theta[7]*S[,8]+theta[8]*S[,9]+theta[9]*S[,10])
      }
      return(eqn1)
    }
    logparamstress <- function(theta,S) {
      if(dim(S)[2]==2){
        eqn2<-theta[1]*rep(1,length(S[,1]))+theta[2]*S[,1]+theta[3]*S[,2]
      }
      if(dim(S)[2]==3){
        eqn2<-theta[1]*rep(1,length(S[,1]))+theta[2]*S[,1]+theta[3]*S[,2]+theta[4]*S[,3]
      }
      if(dim(S)[2]==4){
        eqn2<-theta[1]+theta[2]*S[,1]+theta[3]*S[,2]+theta[4]*S[,3]+theta[5]*S[,4]
      }
      if(dim(S)[2]==5){
        eqn2<-theta[1]+theta[2]*S[,1]+theta[3]*S[,2]+theta[4]*S[,3]+theta[5]*S[,4]+theta[4]*S[,5]
      }
      if(dim(S)[2]==6){
        eqn2<-theta[1]+theta[2]*S[,1]+theta[3]*S[,2]+theta[4]*S[,3]+theta[5]*S[,4]+theta[4]*S[,5]+theta[5]*S[,6]
      }
      if(dim(S)[2]==7){
        eqn2<-theta[1]+theta[2]*S[,1]+theta[3]*S[,2]+theta[4]*S[,3]+theta[5]*S[,4]+theta[4]*S[,5]+theta[5]*S[,6]+theta[6]*S[,7]
      }
      if(dim(S)[2]==8){
        eqn2<-theta[1]+theta[2]*S[,1]+theta[3]*S[,2]+theta[4]*S[,3]+theta[5]*S[,4]+theta[4]*S[,5]+theta[5]*S[,6]+theta[6]*S[,7]+theta[7]*S[,8]
      }
      if(dim(S)[2]==9){
        eqn2<-theta[1]+theta[2]*S[,1]+theta[3]*S[,2]+theta[4]*S[,3]+theta[5]*S[,4]+theta[4]*S[,5]+theta[5]*S[,6]+theta[6]*S[,7]+theta[7]*S[,8]+theta[8]*S[,9]
      }
      if(dim(S)[2]==10){
        eqn2<-theta[1]+theta[2]*S[,1]+theta[3]*S[,2]+theta[4]*S[,3]+theta[5]*S[,4]+theta[4]*S[,5]+theta[5]*S[,6]+theta[6]*S[,7]+theta[7]*S[,8]+theta[8]*S[,9]+theta[9]*S[,10]
      }
      return(eqn2)
    }
    # Writeup for the output text
    ps_txt<-"Multi-Stress"
    params_txt3<-paste("a_",c(0:length(S[1,])),sep="")
    param_txt2<-"exp(a_0 + a_1*S_1 + a_2*S_2 + ...+ a_n*S_n)"
    logparam_txt<-"a_0 + a_1*S_1 + a_2*S_2 + ...+ a_n*S_n"
  }

  if (is.null(modelstress) == FALSE && modelstress=="TempHumidity"){
    # lsparams[1] - parameter A, lsparams[2] - parameter a, lsparams[3] - parameter b
    # First redefine A parameter as logA
    LSQest3 <- c(log(adtLSQ2[1],adtLSQ2[2:3]))
    # positivity_v[1]<-1

    paramstress <- function(theta,S) {
      exp(theta[1])*exp((theta[2]/S[,1]) + (theta[3]/S[,2]))
    }
    logparamstress <- function(theta,S) {
      theta[1] + (theta[2]/S[,1]) + (theta[3]/S[,2])
    }
    # Writeup for the output text
    ps_txt<-"Temperature-Humidity"
    params_txt3<-c("A_0","a_0","b_0")
    param_txt2<-"A_0 exp(a_0/S + b_0/H)"
    logparam_txt<-"ln(A_0) + a_0/S + b_0/H"
  }

  if (is.null(modelstress) == FALSE && modelstress=="TempNonthermal"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
    # First redefine c parameter as logc
    LSQest3 <- c(adtLSQ2[1:2],log(adtLSQ2[3]))
    # positivity_v[3]<-1

    paramstress <- function(theta,S) {
      exp(theta[3])/((S[,2]^theta[2])*exp(-theta[1]/S[,1]))
    }
    logparamstress <- function(theta,S) {
      theta[3] - theta[2]*log(S[,2]) + (theta[1]/S[,1])
    }
    # Writeup for the output text
    ps_txt<-"Temperature-Non-thermal"
    params_txt3<-c("a_0","b_0","c_0")
    param_txt2<-"c_0/(U^b_0 * exp(-a_0/S))"
    logparam_txt<-"a_0(1/S) - b_0*ln(U) + ln(c_0)"
  }

  if (is.null(modelstress) == FALSE && modelstress=="Eyring3"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # lsparams[3] - parameter c, lsparams[4] - parameter d
    LSQest3 <- adtLSQ2

    paramstress <- function(theta,S) {
      (1/S[,1])*exp((theta[1] + (theta[2]/S[,1])) + (theta[3] + (theta[4]/S[,1]))*S[,2])
    }
    logparamstress <- function(theta,S) {
      -log(S[,1]) + theta[1] + (theta[2]/S[,1]) + (theta[3] + (theta[4]/S[,1]))*S[,2]
    }
    # Writeup for the output text
    ps_txt<-"Eyring (Type 3)"
    params_txt3<-c("a_0","b_0","c_0","d_0")
    param_txt2<-"(1/S) exp((a_0 + (b_0/S)) + (c_0 + (d_0/S)) U)"
    logparam_txt<-"-ln(S) + (a_0 + (b_0/S)) + (c_0 + (d_0/S)) U"
  }

  if (is.null(modelstress) == FALSE && modelstress=="Eyring4"){
    # lsparams[1] - parameter A, lsparams[2] - parameter b
    # lsparams[3] - parameter Ea
    # Temperature HAS to be in Kelvin for this to work
    LSQest3 <- adtLSQ2
    K<-8.617385e-5
    paramstress <- function(theta,S) {
      theta[1]*exp(theta[3]/(K*S[,1]))*(S[,2]^-theta[2])
    }
    logparamstress <- function(theta,S) {
      log(theta[1]) + theta[3]/(K*S[,1]) - theta[2]*log(S[,2])
    }

    ps_txt<-"Eyring (Type 3)"
    params_txt3<-c("A_0","b_0","E_a_0")
    param_txt2<-"A_0 exp(E_a_0/(K*S)) U^-b_0"
    logparam_txt<-"ln(A_0) + (E_a_0/(K*S)) - b_0 ln(U)"
  }

  # ADJUSTMENT
  # Adding nlme into tool to account for mixed effects so I will need to make some
  # code changes
  # Make a copy of the data block to unify naming structure used by nlme
  # Time column = "LIFE", damage column = "DAM", unit column = "UNIT", stress columns
  # stay the same

  # ===================================
  # COMMENT OUT FOR NOW 2/26/2024
  # ===================================
  # dataCOPY <- data
  # names(dataCOPY)[1:3] <- c("LIFE","DAM","UNIT")
  # dataGROUP <- split(cbind(dataCOPY$LIFE,dataCOPY$DAM),dataCOPY$UNIT)
  # # Resort the data groups into a list of pairs for each unit for readability
  # dat_list <- function(x1){
  #   mapply(c, as.list(x1[1:(0.5*length(x1))]),as.list(x1[(1+(0.5*length(x1))):length(x1)]), SIMPLIFY=FALSE)
  # }
  # dataGROUP2 <- lapply(dataGROUP,dat_list)
  # ===================================

  # Test functions for likelihoods
  # x is grouped data, a and b are damage life model parameters that will be subjected
  # to numerical integration

  # x input will be the list of grouped data (life/damage) but each must now be processed
  # as its own list to allow for matrix a and b input needed to compute the integral.
  # Define primary function where input is grouped data and set of a and b matrices (same nxn size)
  # Use this for Power model and alternate ones for other models

  # ===================================
  # COMMENT OUT FOR NOW 2/26/2024
  # ===================================
  # f_ab1 <- function(x,a,b){
  #   # Function for summation part of f(a,b) equation
  #   damage1 <- function(x1){
  #     subdamage <- function(x2){
  #       # (D - b*l^a)^2 Power damage-life model
  #       (x2[2] - b*(x2[1]^a))^2
  #     }
  #     lapply(x1,subdamage)
  #   }
  #   parts <- lapply(x,damage1)
  #   # sum matrices per unit
  #   sum_mat <- function(x3){
  #     Reduce("+",x3)
  #   }
  #   # SUM (D - b*l^a)^2 Power damage-life model
  #   parts2 <- lapply(parts,sum_mat)
  #   # In one equation
  #   parts <- lapply(lapply(x,damage1),sum_mat)
  #   return(parts)
  # }
  #
  # # Second part of the exponential expression where the multivariate normal part is defined.  May place as general function
  # # after damage-life model section is stated.
  #
  # f_ab2 <- function(x,a,b,theta){
  #   # Number of elements per unit
  #   mi <- lapply(x,length)
  #   # Covariance of a and b as a function of variance of a (theta[4]),
  #   # variance of b (theta[5]), and correlation between a and b (theta[6])
  #   cov_ab <- theta[6]*sqrt(theta[4]*theta[5])
  #   # MVN matrix
  #   MVN_mat <- ((((a - theta[2])^2)*theta[5]) + (((b - theta[3])^2)*theta[4]) - 2*(a - theta[2])*(b - theta[3])*theta[6]*sqrt(theta[4]*theta[5]))/(theta[4]*theta[5]*(1 - (theta[6]^2)))
  #   MVN_mat<-lapply(lapply(mi,function(x){x/2}),function(x2){x2*MVN_mat})
  #   return(MVN_mat)
  # }
  #
  # f_ab <- function(x,a,b,theta){
  #   # Number of elements per unit
  #   fab <- lapply(mapply(function(x2,y2){x2+y2}, lapply(f_ab1(x,a,b),function(x1){0.5*(1/theta[1]^2)*x1}), f_ab2(x,a,b,theta), SIMPLIFY=FALSE),
  #                 function(x3){exp(-x3)})
  #   return(fab)
  # }
  #
  # # Number of integral trapezoids per parameter
  # n <- 1000
  # # Setup the h-step matrix
  # h_mat <- ones(n)
  # h_mat[,1] <- rep(0.5,n)
  # h_mat[,n] <- rep(0.5,n)
  # h_mat[1,] <- rep(0.5,n)
  # h_mat[n,] <- rep(0.5,n)
  # h_mat[1,1] <- 0.25
  # h_mat[1,n] <- 0.25
  # h_mat[n,1] <- 0.25
  # h_mat[n,n] <- 0.25
  #
  # # New log-likelihood equation
  # loglik <- function(theta){
  #   # neg sum of LN of integral (double or triple) of lognpdfxmvnpdf
  #   N*log(theta[1]) + 1.5*N*log(2*pi) + N*log(theta[4]*theta[5]*(1 - (theta[6]^2)))
  # }
  # ===================================



  # Setup positivity check vector for parameters
  if(dl=="Hamada"){
    positivity_v<-rep(0,dim(adtLSQ)[2]-1+6)
  } else {
    positivity_v<-rep(0,dim(adtLSQ)[2]-1+3)
  }

  # Life-Degradation (damage) models
  # Modify to accept different time values for unit based MLE
  if(dl=="Linear"){
    # D = a + b*t
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"D = a + b*t"
    damage <- function(theta,Time) {
      theta[1] + Time*theta[2]
    }
    logdamage <- function(theta,Time) {
      log(theta[1] + Time*theta[2])
    }
    damfit <- function(TimeDamfit,params){
      params[1] + TimeDamfit*params[2]
    }
    dl_txt<-dl
    params_txt<-c("a","b")
    LSQest1 <- cbind(rep(0,length(unitnames)),adtLSQ[,1],adtLSQ[,2])
    LSQest<-c(0,mean(adtLSQ[,1]), mean(adtLSQ[,2]))
    sigparamno<-3

    if(modelstresstype == 1 || modelstresstype == 2){
      damage1 <- function(theta,Time,S) {
        (paramstress(theta[2:length(theta)],S)) + Time*theta[1]
      }
      logdamage1 <- function(theta,Time,S) {
        log((paramstress(theta[2:length(theta)],S)) + Time*theta[1])
      }
      LSQest2<-c(0,mean(adtLSQ[,2]),LSQest3)
      params_txt<-c("b",params_txt3)
    }
  }

  if(dl=="Exponential"){
    # D = b*exp(a*t)
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"b*exp(a*t)"
    damage <- function(theta,Time) {
      exp(theta[2])*exp(Time*theta[1])
    }
    logdamage <- function(theta,Time) {
      theta[2] + Time*theta[1]
    }
    damfit <- function(TimeDamfit,params){
      exp(params[2])*exp(TimeDamfit*params[1])
    }
    dl_txt<-dl
    params_txt<-c("a","b")
    LSQest1 <- cbind(rep(0,length(unitnames)),adtLSQ[,1],log(adtLSQ[,2]))
    LSQest<-c(0,mean(adtLSQ[,1]), mean(log(adtLSQ[,2])))
    sigparamno<-3

    if(modelstresstype == 1 || modelstresstype == 2){
      damage1 <- function(theta,Time,S) {
        (paramstress(theta[2:length(theta)],S))*exp(Time*theta[1])
      }
      logdamage1 <- function(theta,Time,S) {
        logparamstress(theta[2:length(theta)],S) + Time*theta[1]
      }
      LSQest2<-c(0,mean(adtLSQ[,1]),LSQest3)
      params_txt<-c("a",params_txt3)
    }
  }

  if(dl=="SquareRoot"){
    # D^(1/2) = a + b*t
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"(a + b*t)^2"
    damage <- function(theta,Time) {
      (theta[1] + Time*theta[2])^2
    }
    logdamage <- function(theta,Time) {
      2*log(theta[1] + Time*theta[2])
    }
    damfit <- function(TimeDamfit,params){
      (params[1] + TimeDamfit*params[2])^2
    }
    dl_txt<-dl
    params_txt<-c("a","b")
    LSQest1 <- cbind(rep(0,length(unitnames)),adtLSQ[,1],adtLSQ[,2])
    LSQest<-c(0,mean(adtLSQ[,1]), mean(adtLSQ[,2]))
    sigparamno<-3

    if(modelstresstype == 1 || modelstresstype == 2){
      damage1 <- function(theta,Time,S) {
        (paramstress(theta[2:length(theta)],S) + Time*theta[1])^2
      }
      logdamage1 <- function(theta,Time,S) {
        (paramstress(theta[2:length(theta)],S) + TimeDamfit*params[1])^2
      }
      LSQest2<-c(0,mean(adtLSQ[,2]),LSQest3)
      params_txt<-c("b",params_txt3)
    }
  }

  if(dl=="Power"){
    # D = b*(t^a)
    # theta[1] ~ a, theta[2] ~ b
    # shift b to log b
    Dt_txt<-"b*(t^a)"
    damage <- function(theta,Time) {
      exp(theta[2])*(Time^theta[1])
    }
    logdamage <- function(theta,Time) {
      theta[2] + theta[1]*log(Time)
    }
    damfit <- function(TimeDamfit,params){
      params[2]*(TimeDamfit^params[1])
    }
    dl_txt<-dl
    params_txt<-c("a","b")
    LSQest1 <- cbind(rep(0,length(unitnames)),adtLSQ[,1],log(adtLSQ[,2]))
    LSQest<-c(0,mean(adtLSQ[,1]), mean(log(adtLSQ[,2])))
    sigparamno<-3

    if(modelstresstype == 1 || modelstresstype == 2){
      damage1 <- function(theta,Time,S) {
        (paramstress(theta[2:length(theta)],S))*(Time^theta[1])
      }
      logdamage1 <- function(theta,Time,S) {
        logparamstress(theta[2:length(theta)],S) + theta[1]*log(Time)
      }
      LSQest2<-c(0,mean(adtLSQ[,1]),LSQest3)
      params_txt<-c("a",params_txt3)
    }
  }
  # return(list(loglik,LSQest,loglik(LSQest),Dij,Lij))

  if(dl=="Logarithmic"){
    # D = a + b*ln(t)
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"a + b*ln(t)"
    damage <- function(theta,Time) {
      theta[1] + theta[2]*log(Time)
    }
    logdamage <- function(theta,Time) {
      log(theta[1] + theta[2]*log(Time))
    }
    damfit <- function(TimeDamfit,params){
      params[1] + params[2]*log(TimeDamfit)
    }

    dl_txt<-dl
    params_txt<-c("a","b")
    LSQest1 <- cbind(rep(0,length(unitnames)),adtLSQ[,1],adtLSQ[,2])
    LSQest<-c(0,mean(adtLSQ[,1]), mean(adtLSQ[,2]))
    sigparamno<-3

    if(modelstresstype == 1 || modelstresstype == 2){
      damage1 <- function(theta,Time,S) {
        (paramstress(theta[2:length(theta)],S)) + theta[1]*log(Time)
      }
      logdamage1 <- function(theta,Time,S) {
        log((paramstress(theta[2:length(theta)],S)) + theta[1]*log(Time))
      }
      LSQest2<-c(0,mean(adtLSQ[,1]),LSQest3)
      params_txt<-c("b",params_txt3)
    }
  }

  # if(dl=="Gompertz"){
  #   # D = a + b^(c*t)
  #   # theta[1] ~ a, theta[2] ~ b, theta[3] ~ c
  #   Dt_txt<-"a + b^(c*t)"
  #   damage <- function(theta) {
  #     theta[1] + theta[2]*log(TimeDam)
  #   }
  #   logdamage <- function(theta) {
  #     log(theta[1] + theta[2]*log(TimeDam))
  #   }
  #   dl_txt<-dl
  #   params_txt<-c("a","b")
  #   LSQest<-c(mean(adtLSQ[,1]), mean(adtLSQ[,2]), 1)
  #   sigparamno<-4
  # }

  if(dl=="LloydLipow"){
    # D = a - b/t
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"a + b/t"
    damage <- function(theta,Time) {
      theta[1] - (theta[2]/Time)
    }
    logdamage <- function(theta,Time) {
      log(theta[1] - (theta[2]/Time))
    }
    damfit <- function(TimeDamfit,params){
      params[1] - (params[2]/TimeDamfit)
    }
    dl_txt<-dl
    params_txt<-c("a","b")
    LSQest1 <- cbind(rep(0,length(unitnames)),adtLSQ[,1],adtLSQ[,2])
    LSQest<-c(1,mean(adtLSQ[,1]), mean(adtLSQ[,2]))
    sigparamno<-3

    if(modelstresstype == 1 || modelstresstype == 2){
      damage1 <- function(theta,Time,S) {
        (paramstress(theta[2:length(theta)],S)) - (theta[1]/Time)
      }
      logdamage1 <- function(theta,Time,S) {
        log((paramstress(theta[2:length(theta)],S)) - (theta[1]/Time))
      }
      LSQest2<-c(0,mean(adtLSQ[,2]),LSQest3)
      params_txt<-c("b",params_txt3)
    }
  }

  if(dl=="Mitsuom"){
    # D = 1/(1 + b*(t^a))
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"1/(1 + b*(t^a))"
    damage <- function(theta,Time) {
      1/(1 + theta[2]*(Time^theta[1]))
    }
    logdamage <- function(theta,Time) {
      -log(1 + theta[2]*(Time^theta[1]))
    }
    damfit <- function(TimeDamfit,params){
      1/(1 + params[2]*(TimeDamfit^params[1]))
    }
    dl_txt<-dl
    params_txt<-c("a","b")
    LSQest1 <- cbind(rep(0,length(unitnames)),adtLSQ[,1],adtLSQ[,2])
    LSQest<-c(1,mean(adtLSQ[,1]), mean(adtLSQ[,2]))
    sigparamno<-3

    if(modelstresstype == 1 || modelstresstype == 2){
      damage1 <- function(theta,Time,S) {
        1/(1 + (paramstress(theta[2:length(theta)],S))*(Time^theta[1]))
      }
      logdamage1 <- function(theta,Time,S) {
        -log(1 + (paramstress(theta[2:length(theta)],S))*(Time^theta[1]))
      }
      LSQest2<-c(0,mean(adtLSQ[,1]),LSQest3)
      params_txt<-c("a",params_txt3)
    }
  }

  if(dl=="Hamada"){
    # D = 1/(1 + beta1*(t*exp(beta3*11605*(1/Tu - 1/Ti)))^beta2)
    # theta[1] ~ beta1, theta[2] ~ beta2, theta[3] ~ beta3
    TimeDam
    Dt_txt<-"1/(1 + \U03B2_1*(t*exp(\U03B2_3*11605*(1/Tuse - 1/Temp)))^\U03B2_2)"
    damage <- function(theta,Time) {
      1/(1 + theta[1]*(Time*exp(theta[3]*11605*(1/Tuse - 1/TempDam)))^theta[2])
    }
    logdamage <- function(theta,Time) {
      -log(1 + theta[1]*(Time*exp(theta[3]*11605*(1/Tuse - 1/TempDam)))^theta[2])
    }
    damfit <- function(TimeDamfit,TempDamfit,params){
      1/(1 + params[1]*(TimeDamfit*exp(params[3]*11605*(1/Tuse - 1/TempDamfit)))^params[2])
    }
    dl_txt<-dl
    params_txt<-c("\U03B2_1","\U03B2_2","\U03B2_3")
    LSQest1 <- cbind(rep(0,length(unitnames)),adtLSQ[,1],adtLSQ[,2])
    LSQest<-c(1,mean(adtLSQ[,1]), mean(adtLSQ[,2]), mean(adtLSQ[,3]))
    positivity_v[1]<-1
    positivity_v[2]<-1
    sigparamno<-4
  }
  # life damage model "CrackProp1" is an assumed effective zero initial crack length
  # model adjusted for an initial crack length of 1 mm or 0.001 m.  It was used in the
  # ENRE 641 2022 Final Exam with an assumed DS of 200 MPa.
  if(dl=="CrackProp1"){
    # D = a_0(=0.001 m) + exp((2/(2-m))*LN(N) + (2/(2-m))*(LN(C*sqrt-pi) + LN(1 - 0.5*m) + LN(DS)))
    # theta[1] ~ logC, theta[2] ~ m (C is too small to control so take the natural log and process that)
    Dt_txt<-"a_0(=0.001 m) + exp((2/(2-m))*LN(N) + (2/(2-m))*(LN(C*sqrt-pi) + LN(1 - 0.5*m) + LN(DS)))"
    damage <- function(theta) {
      0.001 + exp((2/(2 - theta[2]))*log(TimeDam) + (2/(2 - theta[2]))*(log(exp(theta[1])*sqrt(pi)) + log(1 - 0.5*theta[2]) + log(200)))
    }
    logdamage <- function(theta) {
      log(0.001 + exp((2/(2 - theta[2]))*log(TimeDam) + (2/(2 - theta[2]))*(log(exp(theta[1])*sqrt(pi)) + log(1 - 0.5*theta[2]) + log(200))))
    }
    damfit <- function(TimeDamfit,params){
      0.001 + exp((2/(2 - params[2]))*log(TimeDamfit) + (2/(2 - params[2]))*(log(params[1]*sqrt(pi)) + log(1 - 0.5*params[2]) + log(200)))
    }
    dl_txt<-dl
    params_txt<-c("C","m")
    LSQest<-c(0.1,mean(log(adtLSQ[,1])), (mean(adtLSQ[,2])))
    sigparamno<-3
  }

  if(dl=="CrackProp2"){
    # D = a_0(=0.001 m)*exp(pi*C*DS^2*N_a)
    # theta[1] ~ logC (C is too small to control so take the natural log and process that)
    Dt_txt<-"a_0(=0.001 m)*exp(pi*C*DS^2*N_a)"
    damage <- function(theta) {
      0.001*exp(pi*exp(theta[1])*(200^2)*TimeDam)
    }
    logdamage <- function(theta) {
      log(0.001) + (pi*exp(theta[1])*(200^2)*TimeDam)
    }
    damfit <- function(TimeDamfit,params){
      0.001*exp(pi*params[1]*(200^2)*TimeDamfit)
    }
    dl_txt<-dl
    params_txt<-c("C")
    LSQest<-c(0.1, mean(log(adtLSQ[,1])))
    sigparamno<-2
  }

  ## NEW (3/4/2024) Initialize variance-covariance matrix in the case of parameter correlation
  ## effects.  Consider getting the initial variance-covariance matrix from a non-correlated run
  ## executed within the script.

  # # Establish cases where degradation-life model is only 2 parameters
  # if((correl == 1 || correl == 2) && modelstresstype == 0){
  #   VARCOV <- function(theta){
  #     # VAR(1) - theta[1]
  #     # COV(1,2), COV(2,1) - theta[2]
  #     # VAR(2) - theta[3]
  #     cbind(c(theta[1],theta[2]),
  #           c(theta[2],theta[3]))
  #   }
  #   # Set the order for the LSQ elements of the Variance-Covariance
  #   LSQVARCOVSET <- function(VARCOV){
  #     c(VARCOV[1,1],VARCOV[1,2],
  #       VARCOV[2,2])
  #   }
  #   # Number of parameters in correlation
  #   d <- 2
  #   # Establish function for numeric integration
  #   Funct2 <- function(theta,x,y){
  #     exp(-0.5*(exp(theta[1])^-2)*((DegradationFULL - damage(c(x,y),TimeFULL))^2) - 0.5*c(t(c(x,y) - theta[2:3])%*%inv(VARCOV(theta[4:6]))%*%(c(x,y) - theta[2:3])))
  #   }
  #   # Establish the numeric integration function
  #   NUMINT <- function(theta){
  #     alp <- sqrt(2/3)
  #     bet <- 1/sqrt(5)
  #     # First bounds
  #     h1 <- 0.5*(b-a)
  #     m1 <- 0.5*(b+a)
  #     # Second bounds
  #     h2 <- 0.5*(d-c)
  #     m2 <- 0.5*(d+c)
  #
  #     est <- ((h1*h2)/(1470^2))*((77^2)*(Func2(theta, a, c) + Func2(theta, a, d) +
  #                                          Func2(theta, b, c) + Func2(theta, b, d)) +
  #                                  (432^2)*(Func2(theta, m1 + alp*h1, m2 + alp*h2) + Func2(theta, m1 + alp*h1, m2 - alp*h2) +
  #                                             Func2(theta, m1 - alp*h1, m2 + alp*h2) + Func2(theta, m1 - alp*h1, m2 - alp*h2)) +
  #                                  (625^2)*(Func2(theta, m1 + bet*h1, m2 + bet*h2) + Func2(theta, m1 + bet*h1, m2 - bet*h2) +
  #                                             Func2(theta, m1 - bet*h1, m2 + bet*h2) + Func2(theta, m1 - bet*h1, m2 - bet*h2)) +
  #                                  (672^2)*Func2(theta, m1, m2) +
  #                                  (77*432)*(Func2(theta, a, m2 + alp*h2) + Func2(theta, a, m2 - alp*h2) +
  #                                              Func2(theta, b, m2 + alp*h2) + Func2(theta, b, m2 - alp*h2) +
  #                                              Func2(theta, m1 + alp*h1, c) + Func2(theta, m1 - alp*h1, c) +
  #                                              Func2(theta, m1 + alp*h1, d) + Func2(theta, m1 - alp*h1, d)) +
  #                                  (77*625)*(Func2(theta, a, m2 + bet*h2) + Func2(theta, a, m2 - bet*h2) +
  #                                              Func2(theta, b, m2 + bet*h2) + Func2(theta, b, m2 - bet*h2) +
  #                                              Func2(theta, m1 + bet*h1, c) + Func2(theta, m1 - bet*h1, c) +
  #                                              Func2(theta, m1 + bet*h1, d) + Func2(theta, m1 - bet*h1, d)) +
  #                                  (77*672)*(Func2(theta, a, m2) + Func2(theta, b, m2) +
  #                                              Func2(theta, m1, c) + Func2(theta, m1, d)) +
  #                                  (432*625)*(Func2(theta, m1 + alp*h1, m2 + bet*h2) + Func2(theta, m1 + bet*h1, m2 + alp*h2) +
  #                                               Func2(theta, m1 + alp*h1, m2 - bet*h2) + Func2(theta, m1 - alp*h1, m2 + bet*h2) +
  #                                               Func2(theta, m1 + bet*h1, m2 - alp*h2) + Func2(theta, m1 - bet*h1, m2 + alp*h2) +
  #                                               Func2(theta, m1 - alp*h1, m2 - bet*h2) + Func2(theta, m1 - bet*h1, m2 - alp*h2)) +
  #                                  (432*672)*(Func2(theta, m1 + alp*h1, m2) + Func2(theta, m1, m2 + alp*h2) +
  #                                               Func2(theta, m1 - alp*h1, m2) + Func2(theta, m1, m2 - alp*h2)) +
  #                                  (625*672)*(Func2(theta, m1 + bet*h1, m2) + Func2(theta, m1, m2 + bet*h2) +
  #                                               Func2(theta, m1 - bet*h1, m2) + Func2(theta, m1, m2 - bet*h2)))
  #     return(est)
  #   }
  # }
  # Establish cases where degradation-life model is only 3 parameters
  if((correl == 1 || correl == 2) && modelstresstype == 1 &&
     (modelstress=="Linear" || modelstress=="Exponential" || modelstress=="Exponential2" ||
      modelstress=="Arrhenius" || modelstress=="Eyring" || modelstress=="Eyring2" ||
      modelstress=="Power" || modelstress=="InversePower" || modelstress=="InversePower2" ||
      modelstress=="Logarithmic")){
    VARCOV <- function(theta){
      # VAR(1) - theta[1]
      # COV(1,2), COV(2,1) - theta[2]
      # COV(1,3), COV(3,1) - theta[3]
      # VAR(2) - theta[4]
      # COV(2,3), COV(3,2) - theta[5]
      # VAR(3) - theta[6]
      cbind(c(theta[1],theta[2],theta[3]),
            c(theta[2],theta[4],theta[5]),
            c(theta[3],theta[5],theta[6]))
    }
    # Set the order for the LSQ elements of the Variance-Covariance
    LSQVARCOVSET <- function(VARCOV){
      c(VARCOV[1,1],VARCOV[1,2],VARCOV[1,3],
        VARCOV[2,2],VARCOV[2,3],
        VARCOV[3,3])
    }
    # Number of parameters in correlation
    d <- 3
    # Establish the numeric integration function
    NUMINT <- function(x,y,z){
      est <- ((h1*h2*h3)/(1470^3))*((77^3)*(Func3(a, c, e) + Func3(a, c, f) + Func3(a, d, e) + Func3(b, c, e) +
                                              Func3(a, d, f) + Func3(b, c, f) + Func3(b, d, e) + Func3(b, d, f)) +
                                      ((77^2)*672)*(Func3(a, c, m3) + Func3(a, d, m3) + Func3(b, c, m3) + Func3(a, m2, e) +
                                                 Func3(b, d, m3) + Func3(a, m2, f) + Func3(b, m2, e) + Func3(m1, c, e) +
                                                 Func3(b, m2, f) + Func3(m1, c, f) + Func3(m1, d, e) + Func3(m1, d, f)) +
                                      (77*(672^2))*(Func3(a, m2, m3) + Func3(b, m2, m3) + Func3(m1, c, m3) +
                                                  Func3(m1, d, m3) + Func3(m1, m2, e) + Func3(m1, m2, f)) +
                                      (672^3)*Func3(m1, m2, m3) +
                                      3705625*(Func3(a, c, m3 + bet*h3) + Func3(a, c, m3 - bet*h3) + Func3(a, d, m3 + bet*h3) + Func3(b, c, m3 + bet*h3) +
                                                 Func3(a, d, m3 - bet*h3) + Func3(b, c, m3 - bet*h3) + Func3(a, m2 + bet*h2, e) + Func3(a, m2 - bet*h2, e) +
                                                 Func3(b, d, m3 + bet*h3) + Func3(b, d, m3 - bet*h3) + Func3(a, m2 + bet*h2, f) + Func3(b, m2 + bet*h2, e) +
                                                 Func3(a, m2 - bet*h2, f) + Func3(b, m2 - bet*h2, e) + Func3(b, m2 + bet*h2, f) + Func3(b, m2 - bet*h2, f) +
                                                 Func3(m1 + bet*h1, c, f) + Func3(m1 + bet*h1, d, e) + Func3(m1 - bet*h1, c, f) + Func3(m1 - bet*h1, d, e) +
                                                 Func3(m1 + bet*h1, d, f) + Func3(m1 - bet*h1, d, f)) +
                                      2561328*(Func3(a, d, m3 + alp*h3) + Func3(b, c, m3 + alp*h3) + Func3(a, d, m3 - alp*h3) + Func3(b, c, m3 - alp*h3) +
                                                 Func3(a, m2 + alp*h2, e) + Func3(a, m2 - alp*h2, e) + Func3(b, d, m3 + alp*h3) + Func3(b, d, m3 - alp*h3) +
                                                 Func3(a, m2 + alp*h2, f) + Func3(b, m2 + alp*h2, e) + Func3(a, m2 - alp*h2, f) + Func3(b, m2 - alp*h2, e) +
                                                 Func3(m1 + alp*h1, c, e) + Func3(m1 - alp*h1, c, e) + Func3(b, m2 + alp*h2, f) + Func3(b, m2 - alp*h2, f) +
                                                 Func3(m1 + alp*h1, c, f) + Func3(m1 + alp*h1, d, e) + Func3(m1 - alp*h1, c, f) + Func3(m1 - alp*h1, d, e) +
                                                 Func3(m1 + alp*h1, d, f) + Func3(m1 - alp*h1, d, f) + Func3(a, c, m3 + alp*h3) + Func3(a, c, m3 - alp*h3)) +
                                      22353408*(Func3(a, m2 + alp*h2, m3) + Func3(a, m2, m3 + alp*h3) + Func3(a, m2 - alp*h2, m3) + Func3(a, m2, m3 - alp*h3) +
                                                  Func3(b, m2 + alp*h2, m3) + Func3(b, m2, m3 + alp*h3) + Func3(b, m2 - alp*h2, m3) + Func3(b, m2, m3 - alp*h3) +
                                                  Func3(m1 + alp*h1, c, m3) + Func3(m1 - alp*h1, c, m3) + Func3(m1, c, m3 + alp*h3) + Func3(m1, c, m3 - alp*h3) +
                                                  Func3(m1 + alp*h1, d, m3) + Func3(m1 - alp*h1, d, m3) + Func3(m1, d, m3 + alp*h3) + Func3(m1, d, m3 - alp*h3) +
                                                  Func3(m1 + alp*h1, m2, e) + Func3(m1, m2 + alp*h2, e) + Func3(m1 - alp*h1, m2, e) + Func3(m1, m2 - alp*h2, e) +
                                                  Func3(m1 + alp*h1, m2, f) + Func3(m1, m2 + alp*h2, f) + Func3(m1 - alp*h1, m2, f) + Func3(m1, m2 - alp*h2, f)) +
                                      32340000*(Func3(a, m2 + bet*h2, m3) + Func3(a, m2, m3 + bet*h3) +
                                                  Func3(a, m2 - bet*h2, m3) + Func3(a, m2, m3 - bet*h3) +
                                                  Func3(b, m2 + bet*h2, m3) + Func3(b, m2, m3 + bet*h3) +
                                                  Func3(b, m2 - bet*h2, m3) + Func3(b, m2, m3 - bet*h3) +
                                                  Func3(m1 + bet*h1, c, m3) + Func3(m1 - bet*h1, c, m3) +
                                                  Func3(m1, c, m3 + bet*h3) + Func3(m1, c, m3 - bet*h3) +
                                                  Func3(m1 + bet*h1, d, m3) + Func3(m1 - bet*h1, d, m3) +
                                                  Func3(m1, d, m3 + bet*h3) + Func3(m1, d, m3 - bet*h3) +
                                                  Func3(m1 + bet*h1, m2, e) + Func3(m1, m2 + bet*h2, e) +
                                                  Func3(m1 - bet*h1, m2, e) + Func3(m1, m2 - bet*h2, e) +
                                                  Func3(m1 + bet*h1, m2, f) + Func3(m1, m2 + bet*h2, f) +
                                                  Func3(m1 - bet*h1, m2, f) + Func3(m1, m2 - bet*h2, f)) +
                   195084288*Func3(m1 + alp*h1, m2, m3) + 195084288*Func3(m1, m2 + alp*h2, m3) +
                   195084288*Func3(m1 - alp*h1, m2, m3) + 195084288*Func3(m1, m2, m3 + alp*h3) +
                   195084288*Func3(m1, m2 - alp*h2, m3) + 195084288*Func3(m1, m2, m3 - alp*h3) +
                   282240000*Func3(m1 + bet*h1, m2, m3) + 282240000*Func3(m1, m2 + bet*h2, m3) +
                   282240000*Func3(m1 - bet*h1, m2, m3) + 282240000*Func3(m1, m2, m3 + bet*h3) +
                   282240000*Func3(m1, m2 - bet*h2, m3) + 282240000*Func3(m1, m2, m3 - bet*h3) +
                   80621568*Func3(m1 + alp*h1, m2 + alp*h2, m3 + alp*h3) + 80621568*Func3(m1 + alp*h1, m2 + alp*h2, m3 - alp*h3) +
                   80621568*Func3(m1 + alp*h1, m2 - alp*h2, m3 + alp*h3) + 80621568*Func3(m1 - alp*h1, m2 + alp*h2, m3 + alp*h3) +
                   80621568*Func3(m1 + alp*h1, m2 - alp*h2, m3 - alp*h3) + 80621568*Func3(m1 - alp*h1, m2 + alp*h2, m3 - alp*h3) +
                   80621568*Func3(m1 - alp*h1, m2 - alp*h2, m3 + alp*h3) + 80621568*Func3(m1 - alp*h1, m2 - alp*h2, m3 - alp*h3) +
                   116640000*Func3(m1 + alp*h1, m2 + alp*h2, m3 + bet*h3) + 116640000*Func3(m1 + alp*h1, m2 + bet*h2, m3 + alp*h3) +
                   116640000*Func3(m1 + bet*h1, m2 + alp*h2, m3 + alp*h3) + 116640000*Func3(m1 + alp*h1, m2 + alp*h2, m3 - bet*h3) +
                   116640000*Func3(m1 + alp*h1, m2 - alp*h2, m3 + bet*h3) + 116640000*Func3(m1 + alp*h1, m2 + bet*h2, m3 - alp*h3) +
                   116640000*Func3(m1 + alp*h1, m2 - bet*h2, m3 + alp*h3) + 116640000*Func3(m1 - alp*h1, m2 + alp*h2, m3 + bet*h3) +
                   116640000*Func3(m1 - alp*h1, m2 + bet*h2, m3 + alp*h3) + 116640000*Func3(m1 + bet*h1, m2 + alp*h2, m3 - alp*h3) +
                   116640000*Func3(m1 + bet*h1, m2 - alp*h2, m3 + alp*h3) + 116640000*Func3(m1 - bet*h1, m2 + alp*h2, m3 + alp*h3) +
                   116640000*Func3(m1 + alp*h1, m2 - alp*h2, m3 - bet*h3) + 116640000*Func3(m1 + alp*h1, m2 - bet*h2, m3 - alp*h3) +
                   116640000*Func3(m1 - alp*h1, m2 + alp*h2, m3 - bet*h3) + 116640000*Func3(m1 - alp*h1, m2 - alp*h2, m3 + bet*h3) +
                   116640000*Func3(m1 - alp*h1, m2 + bet*h2, m3 - alp*h3) + 116640000*Func3(m1 - alp*h1, m2 - bet*h2, m3 + alp*h3) +
                   116640000*Func3(m1 + bet*h1, m2 - alp*h2, m3 - alp*h3) + 116640000*Func3(m1 - bet*h1, m2 + alp*h2, m3 - alp*h3) +
                   116640000*Func3(m1 - bet*h1, m2 - alp*h2, m3 + alp*h3) + 116640000*Func3(m1 - alp*h1, m2 - alp*h2, m3 - bet*h3) +
                   116640000*Func3(m1 - alp*h1, m2 - bet*h2, m3 - alp*h3) + 116640000*Func3(m1 - bet*h1, m2 - alp*h2, m3 - alp*h3) +
                   168750000*Func3(m1 + alp*h1, m2 + bet*h2, m3 + bet*h3) + 168750000*Func3(m1 + bet*h1, m2 + alp*h2, m3 + bet*h3) +
                   168750000*Func3(m1 + bet*h1, m2 + bet*h2, m3 + alp*h3) + 168750000*Func3(m1 + alp*h1, m2 + bet*h2, m3 - bet*h3) +
                   168750000*Func3(m1 + alp*h1, m2 - bet*h2, m3 + bet*h3) + 168750000*Func3(m1 - alp*h1, m2 + bet*h2, m3 + bet*h3) +
                   168750000*Func3(m1 + bet*h1, m2 + alp*h2, m3 - bet*h3) + 168750000*Func3(m1 + bet*h1, m2 - alp*h2, m3 + bet*h3) +
                   168750000*Func3(m1 + bet*h1, m2 + bet*h2, m3 - alp*h3) + 168750000*Func3(m1 + bet*h1, m2 - bet*h2, m3 + alp*h3) +
                   168750000*Func3(m1 - bet*h1, m2 + alp*h2, m3 + bet*h3) + 168750000*Func3(m1 - bet*h1, m2 + bet*h2, m3 + alp*h3) +
                   168750000*Func3(m1 + alp*h1, m2 - bet*h2, m3 - bet*h3) + 168750000*Func3(m1 - alp*h1, m2 + bet*h2, m3 - bet*h3) +
                   168750000*Func3(m1 - alp*h1, m2 - bet*h2, m3 + bet*h3) + 168750000*Func3(m1 + bet*h1, m2 - alp*h2, m3 - bet*h3) +
                   168750000*Func3(m1 + bet*h1, m2 - bet*h2, m3 - alp*h3) + 168750000*Func3(m1 - bet*h1, m2 + alp*h2, m3 - bet*h3) +
                   168750000*Func3(m1 - bet*h1, m2 - alp*h2, m3 + bet*h3) + 168750000*Func3(m1 - bet*h1, m2 + bet*h2, m3 - alp*h3) +
                   168750000*Func3(m1 - bet*h1, m2 - bet*h2, m3 + alp*h3) + 168750000*Func3(m1 - alp*h1, m2 - bet*h2, m3 - bet*h3) +
                   168750000*Func3(m1 - bet*h1, m2 - alp*h2, m3 - bet*h3) + 168750000*Func3(m1 - bet*h1, m2 - bet*h2, m3 - alp*h3) +
                   244140625*Func3(m1 + bet*h1, m2 + bet*h2, m3 + bet*h3) + 244140625*Func3(m1 + bet*h1, m2 + bet*h2, m3 - bet*h3) +
                   244140625*Func3(m1 + bet*h1, m2 - bet*h2, m3 + bet*h3) + 244140625*Func3(m1 - bet*h1, m2 + bet*h2, m3 + bet*h3) +
                   244140625*Func3(m1 + bet*h1, m2 - bet*h2, m3 - bet*h3) + 244140625*Func3(m1 - bet*h1, m2 + bet*h2, m3 - bet*h3) +
                   244140625*Func3(m1 - bet*h1, m2 - bet*h2, m3 + bet*h3) + 244140625*Func3(m1 - bet*h1, m2 - bet*h2, m3 - bet*h3) +
                     14370048*Func3(a, m2 + alp*h2, m3 + alp*h3) + 14370048*Func3(a, m2 + alp*h2, m3 - alp*h3) + 14370048*Func3(a, m2 - alp*h2, m3 + alp*h3) +
                   14370048*Func3(a, m2 - alp*h2, m3 - alp*h3) + 20790000*Func3(a, m2 + alp*h2, m3 + bet*h3) + 20790000*Func3(a, m2 + bet*h2, m3 + alp*h3) + 20790000*Func3(a, m2 + alp*h2, m3 - bet*h3) +
                   20790000*Func3(a, m2 - alp*h2, m3 + bet*h3) + 20790000*Func3(a, m2 + bet*h2, m3 - alp*h3) + 20790000*Func3(a, m2 - bet*h2, m3 + alp*h3) + 20790000*Func3(a, m2 - alp*h2, m3 - bet*h3) +
                   20790000*Func3(a, m2 - bet*h2, m3 - alp*h3) + 14370048*Func3(b, m2 + alp*h2, m3 + alp*h3) + 14370048*Func3(b, m2 + alp*h2, m3 - alp*h3) + 14370048*Func3(b, m2 - alp*h2, m3 + alp*h3) +
                   14370048*Func3(b, m2 - alp*h2, m3 - alp*h3) + 30078125*Func3(a, m2 + bet*h2, m3 + bet*h3) + 30078125*Func3(a, m2 + bet*h2, m3 - bet*h3) + 30078125*Func3(a, m2 - bet*h2, m3 + bet*h3) +
                   30078125*Func3(a, m2 - bet*h2, m3 - bet*h3) + 20790000*Func3(b, m2 + alp*h2, m3 + bet*h3) + 20790000*Func3(b, m2 + bet*h2, m3 + alp*h3) + 20790000*Func3(b, m2 + alp*h2, m3 - bet*h3) +
                   20790000*Func3(b, m2 - alp*h2, m3 + bet*h3) + 20790000*Func3(b, m2 + bet*h2, m3 - alp*h3) + 20790000*Func3(b, m2 - bet*h2, m3 + alp*h3) + 20790000*Func3(b, m2 - alp*h2, m3 - bet*h3) +
                   20790000*Func3(b, m2 - bet*h2, m3 - alp*h3) + 14370048*Func3(m1 + alp*h1, c, m3 + alp*h3) + 14370048*Func3(m1 + alp*h1, c, m3 - alp*h3) + 14370048*Func3(m1 - alp*h1, c, m3 + alp*h3) +
                   14370048*Func3(m1 - alp*h1, c, m3 - alp*h3) + 30078125*Func3(b, m2 + bet*h2, m3 + bet*h3) + 30078125*Func3(b, m2 + bet*h2, m3 - bet*h3) + 30078125*Func3(b, m2 - bet*h2, m3 + bet*h3) +
                   30078125*Func3(b, m2 - bet*h2, m3 - bet*h3) + 20790000*Func3(m1 + alp*h1, c, m3 + bet*h3) + 20790000*Func3(m1 + bet*h1, c, m3 + alp*h3) + 20790000*Func3(m1 + alp*h1, c, m3 - bet*h3) +
                   20790000*Func3(m1 - alp*h1, c, m3 + bet*h3) + 20790000*Func3(m1 + bet*h1, c, m3 - alp*h3) + 20790000*Func3(m1 - bet*h1, c, m3 + alp*h3) + 20790000*Func3(m1 - alp*h1, c, m3 - bet*h3) +
                   20790000*Func3(m1 - bet*h1, c, m3 - alp*h3) + 14370048*Func3(m1 + alp*h1, d, m3 + alp*h3) + 14370048*Func3(m1 + alp*h1, d, m3 - alp*h3) + 14370048*Func3(m1 - alp*h1, d, m3 + alp*h3) +
                   14370048*Func3(m1 - alp*h1, d, m3 - alp*h3) + 30078125*Func3(m1 + bet*h1, c, m3 + bet*h3) + 30078125*Func3(m1 + bet*h1, c, m3 - bet*h3) + 30078125*Func3(m1 - bet*h1, c, m3 + bet*h3) +
                   30078125*Func3(m1 - bet*h1, c, m3 - bet*h3) + 20790000*Func3(m1 + alp*h1, d, m3 + bet*h3) + 20790000*Func3(m1 + bet*h1, d, m3 + alp*h3) + 20790000*Func3(m1 + alp*h1, d, m3 - bet*h3) +
                   20790000*Func3(m1 - alp*h1, d, m3 + bet*h3) + 20790000*Func3(m1 + bet*h1, d, m3 - alp*h3) + 20790000*Func3(m1 - bet*h1, d, m3 + alp*h3) + 20790000*Func3(m1 - alp*h1, d, m3 - bet*h3) +
                   20790000*Func3(m1 - bet*h1, d, m3 - alp*h3) + 14370048*Func3(m1 + alp*h1, m2 + alp*h2, e) + 14370048*Func3(m1 + alp*h1, m2 - alp*h2, e) + 14370048*Func3(m1 - alp*h1, m2 + alp*h2, e) +
                   14370048*Func3(m1 - alp*h1, m2 - alp*h2, e) + 30078125*Func3(m1 + bet*h1, d, m3 + bet*h3) + 30078125*Func3(m1 + bet*h1, d, m3 - bet*h3) + 30078125*Func3(m1 - bet*h1, d, m3 + bet*h3) +
                   30078125*Func3(m1 - bet*h1, d, m3 - bet*h3) + 20790000*Func3(m1 + alp*h1, m2 + bet*h2, e) + 20790000*Func3(m1 + bet*h1, m2 + alp*h2, e) + 20790000*Func3(m1 + alp*h1, m2 - bet*h2, e) +
                   20790000*Func3(m1 - alp*h1, m2 + bet*h2, e) + 20790000*Func3(m1 + bet*h1, m2 - alp*h2, e) + 20790000*Func3(m1 - bet*h1, m2 + alp*h2, e) + 20790000*Func3(m1 - alp*h1, m2 - bet*h2, e) +
                   20790000*Func3(m1 - bet*h1, m2 - alp*h2, e) + 14370048*Func3(m1 + alp*h1, m2 + alp*h2, f) + 14370048*Func3(m1 + alp*h1, m2 - alp*h2, f) + 14370048*Func3(m1 - alp*h1, m2 + alp*h2, f) +
                   14370048*Func3(m1 - alp*h1, m2 - alp*h2, f) + 30078125*Func3(m1 + bet*h1, m2 + bet*h2, e) + 30078125*Func3(m1 + bet*h1, m2 - bet*h2, e) + 30078125*Func3(m1 - bet*h1, m2 + bet*h2, e) +
                   30078125*Func3(m1 - bet*h1, m2 - bet*h2, e) + 20790000*Func3(m1 + alp*h1, m2 + bet*h2, f) + 20790000*Func3(m1 + bet*h1, m2 + alp*h2, f) + 20790000*Func3(m1 + alp*h1, m2 - bet*h2, f) +
                   20790000*Func3(m1 - alp*h1, m2 + bet*h2, f) + 20790000*Func3(m1 + bet*h1, m2 - alp*h2, f) + 20790000*Func3(m1 - bet*h1, m2 + alp*h2, f) + 20790000*Func3(m1 - alp*h1, m2 - bet*h2, f) +
                   20790000*Func3(m1 - bet*h1, m2 - alp*h2, f) + 30078125*Func3(m1 + bet*h1, m2 + bet*h2, f) + 30078125*Func3(m1 + bet*h1, m2 - bet*h2, f) + 30078125*Func3(m1 - bet*h1, m2 + bet*h2, f) +
                   30078125*Func3(m1 - bet*h1, m2 - bet*h2, f) + 125411328*Func3(m1 + alp*h1, m2 + alp*h2, m3) + 125411328*Func3(m1 + alp*h1, m2, m3 + alp*h3) + 125411328*Func3(m1 + alp*h1, m2 - alp*h2, m3) +
                   125411328*Func3(m1 - alp*h1, m2 + alp*h2, m3) + 125411328*Func3(m1, m2 + alp*h2, m3 + alp*h3) + 125411328*Func3(m1 + alp*h1, m2, m3 - alp*h3) + 125411328*Func3(m1 - alp*h1, m2, m3 + alp*h3) +
                   125411328*Func3(m1 - alp*h1, m2 - alp*h2, m3) + 125411328*Func3(m1, m2 + alp*h2, m3 - alp*h3) + 125411328*Func3(m1, m2 - alp*h2, m3 + alp*h3) + 125411328*Func3(m1 - alp*h1, m2, m3 - alp*h3) +
                   125411328*Func3(m1, m2 - alp*h2, m3 - alp*h3) + 181440000*Func3(m1 + alp*h1, m2 + bet*h2, m3) + 181440000*Func3(m1 + bet*h1, m2 + alp*h2, m3) + 181440000*Func3(m1 + alp*h1, m2, m3 + bet*h3) +
                   181440000*Func3(m1 + alp*h1, m2 - bet*h2, m3) + 181440000*Func3(m1 - alp*h1, m2 + bet*h2, m3) + 181440000*Func3(m1 + bet*h1, m2, m3 + alp*h3) + 181440000*Func3(m1 + bet*h1, m2 - alp*h2, m3) +
                   181440000*Func3(m1 - bet*h1, m2 + alp*h2, m3) + 181440000*Func3(m1, m2 + alp*h2, m3 + bet*h3) + 181440000*Func3(m1, m2 + bet*h2, m3 + alp*h3) + 181440000*Func3(m1 + alp*h1, m2, m3 - bet*h3) +
                   181440000*Func3(m1 - alp*h1, m2, m3 + bet*h3) + 181440000*Func3(m1 - alp*h1, m2 - bet*h2, m3) + 181440000*Func3(m1 + bet*h1, m2, m3 - alp*h3) + 181440000*Func3(m1 - bet*h1, m2, m3 + alp*h3) +
                   181440000*Func3(m1 - bet*h1, m2 - alp*h2, m3) + 181440000*Func3(m1, m2 + alp*h2, m3 - bet*h3) + 181440000*Func3(m1, m2 - alp*h2, m3 + bet*h3) + 181440000*Func3(m1, m2 + bet*h2, m3 - alp*h3) +
                   181440000*Func3(m1, m2 - bet*h2, m3 + alp*h3) + 181440000*Func3(m1 - alp*h1, m2, m3 - bet*h3) + 181440000*Func3(m1 - bet*h1, m2, m3 - alp*h3) + 181440000*Func3(m1, m2 - alp*h2, m3 - bet*h3) +
                   181440000*Func3(m1, m2 - bet*h2, m3 - alp*h3) + 262500000*Func3(m1 + bet*h1, m2 + bet*h2, m3) + 262500000*Func3(m1 + bet*h1, m2, m3 + bet*h3) + 262500000*Func3(m1 + bet*h1, m2 - bet*h2, m3) +
                   262500000*Func3(m1 - bet*h1, m2 + bet*h2, m3) + 262500000*Func3(m1, m2 + bet*h2, m3 + bet*h3) + 262500000*Func3(m1 + bet*h1, m2, m3 - bet*h3) + 262500000*Func3(m1 - bet*h1, m2, m3 + bet*h3) +
                   262500000*Func3(m1 - bet*h1, m2 - bet*h2, m3) + 262500000*Func3(m1, m2 + bet*h2, m3 - bet*h3) + 262500000*Func3(m1, m2 - bet*h2, m3 + bet*h3) + 262500000*Func3(m1 - bet*h1, m2, m3 - bet*h3) +
                   262500000*Func3(m1, m2 - bet*h2, m3 - bet*h3))
      return(est)
    }
  }
  # Establish cases where degradation-life model is only 4 parameters
  if((correl == 1 || correl == 2) && modelstresstype == 1 &&
     (modelstress=="TempHumidity" || modelstress=="TempNonthermal" ||
      modelstress=="Eyring4")){
    VARCOV <- function(theta){
      # VAR(1) - theta[1]
      # COV(1,2), COV(2,1) - theta[2]
      # COV(1,3), COV(3,1) - theta[3]
      # COV(1,4), COV(4,1) - theta[4]
      # VAR(2) - theta[5]
      # COV(2,3), COV(3,2) - theta[6]
      # COV(2,4), COV(4,2) - theta[7]
      # VAR(3) - theta[8]
      # COV(3,4), COV(4,3) - theta[9]
      # VAR(4) - theta[10]
      cbind(c(theta[1],theta[2],theta[3],theta[4]),
            c(theta[2],theta[5],theta[6],theta[7]),
            c(theta[3],theta[6],theta[8],theta[9]),
            c(theta[4],theta[7],theta[9],theta[10]))
    }
    # Set the order for the LSQ elements of the Variance-Covariance
    LSQVARCOVSET <- function(VARCOV){
      c(VARCOV[1,1],VARCOV[1,2],VARCOV[1,3],VARCOV[1,4],
        VARCOV[2,2],VARCOV[2,3],VARCOV[2,4],
        VARCOV[3,3],VARCOV[3,4],
        VARCOV[4,4])
    }
    # Number of parameters in correlation
    d <- 4
    # Establish the numeric integration function
    NUMINT <- function(x,y,z){
      return(est)
    }
  }
  # Establish cases where degradation-life model is only 5 parameters
  if((correl == 1 || correl == 2) && modelstresstype == 1 && modelstress=="Eyring3"){
    VARCOV <- function(theta){
      # VAR(1) - theta[1]
      # COV(1,2), COV(2,1) - theta[2]
      # COV(1,3), COV(3,1) - theta[3]
      # COV(1,4), COV(4,1) - theta[4]
      # COV(1,5), COV(5,1) - theta[5]
      # VAR(2) - theta[6]
      # COV(2,3), COV(3,2) - theta[7]
      # COV(2,4), COV(4,2) - theta[8]
      # COV(2,5), COV(5,2) - theta[9]
      # VAR(3) - theta[10]
      # COV(3,4), COV(4,3) - theta[11]
      # COV(3,5), COV(5,3) - theta[12]
      # VAR(4) - theta[13]
      # COV(4,5), COV(5,4) - theta[14]
      # VAR(5) - theta[15]
      cbind(c(theta[1],theta[2],theta[3],theta[4],theta[5]),
            c(theta[2],theta[6],theta[7],theta[8],theta[9]),
            c(theta[3],theta[7],theta[10],theta[11],theta[12]),
            c(theta[4],theta[8],theta[11],theta[13],theta[14]),
            c(theta[5],theta[9],theta[12],theta[14],theta[15]))
    }
    # Set the order for the LSQ elements of the Variance-Covariance
    LSQVARCOVSET <- function(VARCOV){
      c(VARCOV[1,1],VARCOV[1,2],VARCOV[1,3],VARCOV[1,4],VARCOV[1,5],
        VARCOV[2,2],VARCOV[2,3],VARCOV[2,4],VARCOV[2,5],
        VARCOV[3,3],VARCOV[3,4],VARCOV[3,5],
        VARCOV[4,4],VARCOV[4,5],
        VARCOV[5,5])
    }
    # Number of parameters in correlation
    d <- 5
    # Establish the numeric integration function
    NUMINT <- function(x,y,z){
    }

  }

  ## NEW (2/27/2024) Set up log-likelihoods for different cases and scenarios.  The only constant
  ## will be unit by unit MLE evaluation
  # 1. Run MLE for each unit based on dl regardless of the modelstress
  # List of uncorrected MLE and variance-covariance output
  MLE.theta.hat.perunit<-vector(mode = "list", length = length(unitnames))
  MLE.inv.fish.perunit<-vector(mode = "list", length = length(unitnames))
  for(i in 1:length(unitnames)){
    # Since each damage is dependent on a separate time stamp, a different loglikelihood will be
    # formed for each part of the loop
    # RCS 09/05/2024 - ADD WEIBULL
    if (dist=="Weibull") {
      loglik <- function(theta){
        -sum(theta[1] + (exp(theta[1])-1)*log(data[which(data[,3]==unitnames[i]),2]) - exp(theta[1])*logdamage(theta[2:3],data[which(data[,3]==unitnames[i]),1]) - ((data[which(data[,3]==unitnames[i]),2]/damage(theta[2:3],data[which(data[,3]==unitnames[i]),1]))^exp(theta[1])))
      }
      dist_txt<-dist
      distparam_txt<-"\U03B2"
    }
    if (dist=="Normal") {
      loglik <- function(theta){
        -sum(-log(exp(theta[1])) - 0.5*log(2*pi) - 0.5*(exp(theta[1])^-2)*((data[which(data[,3]==unitnames[i]),2] - damage(theta[2:3],data[which(data[,3]==unitnames[i]),1]))^2))
      }
      dist_txt<-dist
      distparam_txt<-"\U03C3"
    }
    if (dist=="Lognormal") {
      loglik <- function(theta){
        -sum(-log(exp(theta[1])) - 0.5*log(2*pi) - log(data[which(data[,3]==unitnames[i]),2]) - 0.5*(exp(theta[1])^-2)*((log(data[which(data[,3]==unitnames[i]),2]) - logdamage(theta[2:3],data[which(data[,3]==unitnames[i]),1]))^2))
      }
      dist_txt<-dist
      distparam_txt<-"\U03C3_t"
    }


    MLEandvar <- MLE.var.covar.select(loglik,LSQest1[i,])
    MLE.theta.hat.perunit[[i]] <- MLEandvar[[1]]
    MLE.inv.fish.perunit[[i]]  <- MLEandvar[[2]]
  }
  # 2. If correlation is selected (=1) then we do the full likelihood with parameter covariance
  # with normal distribution fit. Use mvn fit to initialize that part and numerical analysis for the rest.
  # No correlation (=0) means we compute this without correlation.  However we still need the uncorrelated
  # output as an initial LSQ for the MU and VARCOV Matrix.
  # No correlation Case
  if(correl == 0 || correl == 1 || correl == 2){
    if(modelstresstype == 0){
      if (dist=="Weibull") {
        loglik <- function(theta){
          -sum(theta[1] + (exp(theta[1])-1)*log(DegradationFULL) - exp(theta[1])*logdamage(theta[2:3],TimeFULL) - ((DegradationFULL/damage(theta[2:3],TimeFULL))^exp(theta[1])))
        }
      }
      if (dist=="Normal") {
        loglik <- function(theta){
          -sum(-theta[1] - 0.5*log(2*pi) - 0.5*(exp(theta[1])^-2)*((DegradationFULL - damage(theta[2:3],TimeFULL))^2))
        }
      }
      if (dist=="Lognormal") {
        loglik <- function(theta){
          -sum(-theta[1] - 0.5*log(2*pi) - log(DegradationFULL) - 0.5*(exp(theta[1])^-2)*((log(DegradationFULL) - logdamage(theta[2:3],TimeFULL))^2))
        }
      }
      # if(is.null(Tc)){
      #   loglik <- function(theta){
      #     -sum(-theta[1] - 0.5*log(2*pi) - log(TTF) - 0.5*(exp(theta[1])^-2)*((log(TTF) - loglifeF(theta))^2))
      #   }
      # }

      # if(is.null(Tc)){
      #   loglik <- function(theta){
      #     -sum(theta[1] + (exp(theta[1])-1)*log(TTF) - exp(theta[1])*loglifeF(theta) - ((TTF/lifeF(theta))^exp(theta[1])))
      #   }
      # }
      MLEandvar <- MLE.var.covar.select(loglik,LSQest)
    }
    if(modelstresstype == 1 || modelstresstype == 2){
      if (dist=="Weibull") {
        loglik <- function(theta){
          -sum(theta[1] + (exp(theta[1])-1)*log(DegradationFULL) - exp(theta[1])*logdamage1(theta[2:3],TimeFULL) - ((DegradationFULL/damage1(theta[2:3],TimeFULL))^exp(theta[1])))
        }
      }
      if (dist=="Normal") {
        loglik <- function(theta){
          -sum(-theta[1] - 0.5*log(2*pi) - 0.5*(exp(theta[1])^-2)*((DegradationFULL - damage1(theta[2:length(theta)],TimeFULL,StresSULL))^2))
        }
      }
      if (dist=="Lognormal") {
        loglik <- function(theta){
          -sum(-theta[1] - 0.5*log(2*pi) - log(DegradationFULL) - 0.5*(exp(theta[1])^-2)*((log(DegradationFULL) - logdamage1(theta[2:length(theta)],TimeFULL,StresSULL))^2))
        }
      }

      MLEandvar <- MLE.var.covar.select(loglik,LSQest2)
    }
    theta.hat <- MLEandvar[[1]]
    inv.fish  <- MLEandvar[[2]]
  }
  # Additive Measurement Error
  # Establish cases where degradation-life model is only 2 parameters
  if((correl == 1 || correl == 2) && modelstresstype == 0){
    VARCOV <- function(theta){
      # VAR(1) - theta[1]
      # COV(1,2), COV(2,1) - theta[2]
      # rho - theta[2]
      # VAR(2) - theta[3]
      cbind(c(exp(theta[1]),0.999*erf(theta[2])*sqrt(exp(theta[1])*exp(theta[3]))),
            c(erf(theta[2])*sqrt(exp(theta[1])*exp(theta[3])),exp(theta[3])))
    }
    # Set the order for the LSQ elements of the Variance-Covariance
    LSQVARCOVSET <- function(VARCOV){
      c(log(VARCOV[1,1]),erfinv(VARCOV[1,2]/sqrt(VARCOV[1,1]*VARCOV[2,2]))/0.999,
        log(VARCOV[2,2]))
    }
    # Number of parameters in correlation
    d <- 2
    # Establish function for numeric integration

    # Funct2 <- function(theta,x,y){
    #   exp(-0.5*(exp(theta[1])^-2)*((DegradationFULL - damage(c(x,y),TimeFULL))^2) - 0.5*c(t(c(x,y) - theta[2:3])%*%inv(VARCOV(theta[4:6]))%*%(c(x,y) - theta[2:3])))
    # }

    Funct2a <- function(theta){
      alp <- sqrt(2/3)
      bet <- 1/sqrt(5)
      # First bounds
      a <- qnorm(.00001,theta.hat[2],sqrt(inv.fish[2,2]))
      b <- qnorm(.99999,theta.hat[2],sqrt(inv.fish[2,2]))
      h1 <- 0.5*(b-a)
      m1 <- 0.5*(b+a)
      # Second bounds
      c <- qnorm(.00001,theta.hat[3],sqrt(inv.fish[3,3]))
      d <- qnorm(.99999,theta.hat[3],sqrt(inv.fish[3,3]))
      h2 <- 0.5*(d-c)
      m2 <- 0.5*(d+c)
      # Bound Pairs
      pair1 <- list(c(a,c),c(a,d),c(b,c),c(b,d))
      pair2 <- c(m1,m2)
      pair3 <- list(c(m1 + alp*h1, m2 + alp*h2),c(m1 + alp*h1, m2 - alp*h2),
                    c(m1 - alp*h1, m2 + alp*h2),c(m1 - alp*h1, m2 - alp*h2))
      pair4 <- list(c(a, m2),c(b, m2),c(m1, c),c(m1, d))
      pair5 <- list(c(m1 + bet*h1, m2 + bet*h2),c(m1 + bet*h1, m2 - bet*h2),
                    c(m1 - bet*h1, m2 + bet*h2),c(m1 - bet*h1, m2 - bet*h2))
      pair6 <- list(c(a, m2 + alp*h2),c(a, m2 - alp*h2),
                    c(b, m2 + alp*h2),c(b, m2 - alp*h2),
                    c(m1 + alp*h1, c),c(m1 - alp*h1, c),
                    c(m1 + alp*h1, d),c(m1 - alp*h1, d))
      pair7 <- list(c(a, m2 + bet*h2),c(a, m2 - bet*h2),
                    c(b, m2 + bet*h2),c(b, m2 - bet*h2),
                    c(m1 + bet*h1, c),c(m1 - bet*h1, c),
                    c(m1 + bet*h1, d),c(m1 - bet*h1, d))
      pair8 <- list(c(m1 + alp*h1, m2 + bet*h2),c(m1 + bet*h1, m2 + alp*h2),
                    c(m1 + alp*h1, m2 - bet*h2),c(m1 - alp*h1, m2 + bet*h2),
                    c(m1 + bet*h1, m2 - alp*h2),c(m1 - bet*h1, m2 + alp*h2),
                    c(m1 - bet*h1, m2 - alp*h2),c(m1 - bet*h1, m2 - alp*h2))
      pair9 <- list(c(m1 + alp*h1, m2),c(m1, m2 + alp*h2),
                    c(m1 - alp*h1, m2),c(m1, m2 - alp*h2))
      pair10 <- list(c(m1 + bet*h1, m2),c(m1, m2 + bet*h2),
                     c(m1 - bet*h1, m2),c(m1, m2 - bet*h2))
      Funct2 <- function(X){
        exp(-0.5*(exp(theta[1])^-2)*((DegradationFULL - damage(X,TimeFULL))^2) - 0.5*c(t(X - theta[2:3])%*%inv(VARCOV(theta[4:6]))%*%(X - theta[2:3])))
      }
      est <- ((h1*h2)/(1470^2))*((77^2)*(Reduce('+',lapply(pair1,Funct2))) +
                                   (672^2)*Funct2(pair2) +
                                   (432^2)*(Reduce('+',lapply(pair3,Funct2))) +
                                   (77*672)*(Reduce('+',lapply(pair4,Funct2))) +
                                   (625^2)*(Reduce('+',lapply(pair5,Funct2))) +
                                   (77*432)*(Reduce('+',lapply(pair6,Funct2))) +
                                   (77*625)*(Reduce('+',lapply(pair7,Funct2))) +
                                   (432*625)*(Reduce('+',lapply(pair8,Funct2))) +
                                   (432*672)*(Reduce('+',lapply(pair9,Funct2))) +
                                   (625*672)*(Reduce('+',lapply(pair10,Funct2))))
      # return(list(est,a,b,c,d,h1,m1,h2,m2))
      return(est)
    }
    # Establish the numeric integration function
  }
  if(correl == 1){
    LSQest4 <- c(theta.hat,LSQVARCOVSET(inv.fish[2:(d+1),2:(d+1)]))
    loglik <- function(theta){
      -sum(-log(exp(theta[1])) - 0.5*(1+d)*log(2*pi) - 0.5*log(det(VARCOV(theta[(1+length(theta.hat)):length(LSQest4)]))) + log(Funct2a(theta)))
    }
    # MLEandvar1 <- MLE.var.covar.select(loglik,LSQest4)
    return(list(loglik,LSQest4,VARCOV,DegradationFULL,TimeFULL))
    if(modelstresstype == 0){}
    if(modelstresstype == 1 || modelstresstype == 2){}
  }
  # Multiplicative Measurement Error
  if(correl == 2){
    LSQest4 <- c(theta.hat,LSQVARCOVSET(inv.fish[2:(d+1),2:(d+1)]))
    loglik <- function(theta){
      -sum(-log(exp(theta[1])) - 0.5*(1+d)*log(2*pi) - log(DegradationFULL) - 0.5*log(det(VARCOV(theta[(1+length(theta.hat)):length(LSQest4)]))))
    }
    return(list(loglik,LSQest4,DegradationFULL,TimeFULL))
    if(modelstresstype == 0){}
    if(modelstresstype == 1 || modelstresstype == 2){}
  }

  # return(list(MLE.theta.hat.perunit,MLE.inv.fish.perunit,theta.hat,inv.fish,params_txt))

  # # Then compute the MLE estimate based on the distribution (normal or lognormal)
  # # for the damage model
  # # NOTE: For Version 2.0 consider other distributions that might work like Weibull
  # if (dist=="Lognormal") {
  #   # positivity_v[sigparamno]<-1
  #   -sum(-log(theta[1]) - 0.5*log(2*pi) - log(TTF) - 0.5*(theta[1]^-2)*((log(TTF) - logparamstress(theta))^2))
  #
  #
  #   # loglik <- function(theta){
  #   #   # -sum(-log(theta[sigparamno]) - 0.5*log(2*pi) - log(Dam) - 0.5*(theta[sigparamno]^-2)*((log(Dam) - logdamage(theta))^2))
  #   #   -sum(-log(exp(theta[ishift])) - 0.5*log(2*pi) - log(Dam) - 0.5*(exp(theta[ishift])^-2)*((log(Dam) - logdamage(theta))^2))
  #   # }
  #   # IF case for Hamada since it is the only one with three parameters
  #   if(dl=="Hamada"){
  #     loglik <- function(theta){
  #       # neg sum of LN of integral (double or triple) of lognpdfxmvnpdf
  #       -sum(-log(exp(theta[ishift])) - 0.5*log(2*pi) - log(Dam) - 0.5*(exp(theta[ishift])^-2)*((log(Dam) - logdamage(theta))^2))
  #     }
  #   } else {
  #     # mainfunct <- function(theta){
  #     #   integral2(prod(dlnorm(Dam,logdamage_integral(a,b),theta[1])*dmvnorm(c(a,b),c(theta[2],theta[3]),cbind(c(theta[4],theta[5]),C(theta[5],theta[6])))),aL,aH,bL,bH)
  #     # }
  #     loglik <- function(theta){
  #       # neg sum of LN of integral (double or triple) of lognpdfxmvnpdf
  #       -sum(log(integral2(function(a,b) prod(dlnorm(Dam,logdamage_integral(a,b),theta[1])*dmvnorm(c(a,b),c(theta[2],theta[3]),cbind(c(theta[4],theta[5]),C(theta[5],theta[6])))),a_L,a_H,b_L,b_H)$Q))
  #       # -sum(-log(exp(theta[ishift])) - 0.5*log(2*pi) - log(Dam) - 0.5*(exp(theta[ishift])^-2)*((log(Dam) - logdamage(theta))^2))
  #     }
  #   }
  #
  #   dist_txt<-dist
  #   distparam_txt<-"\U03C3_t"
  # }
  # if (dist=="Normal") {
  #   # positivity_v[sigparamno]<-1
  #
  #   loglik <- function(theta){
  #     -sum(-log(exp(theta[1])) - 0.5*log(2*pi) - 0.5*(exp(theta[1])^-2)*((Dam - damage(theta))^2))
  #   }
  #   dist_txt<-dist
  #   distparam_txt<-"\U03C3"
  # }

  # NEW 8/11/2023
  # Second part of the exponential expression where the multivariate normal part is defined.  May place as general function
  # after damage-life model section is stated.
  # NOTE: The Hamada model will have a modified MVN_mat
  # f_ab2 <- function(x,a,b,theta){
  #   # Number of elements per unit
  #   mi <- lapply(x,length)
  #   # Covariance of a and b as a function of variance of a (theta[4]),
  #   # variance of b (theta[5]), and correlation between a and b (theta[6])
  #   cov_ab <- theta[6]*sqrt(theta[4]*theta[5])
  #   # MVN matrix
  #   MVN_mat <- ((((a - theta[2])^2)*theta[5]) + (((b - theta[3])^2)*theta[4]) - 2*(a - theta[2])*(b - theta[3])*theta[6]*sqrt(theta[4]*theta[5]))/(theta[4]*theta[5]*(1 - (theta[6]^2)))
  #   MVN_mat<-lapply(lapply(mi,function(x){x/2}),function(x2){x2*MVN_mat})
  #   return(MVN_mat)
  # }
  #
  # f_ab <- function(x,a,b,theta){
  #   # Number of elements per unit
  #   fab <- lapply(mapply(function(x2,y2){x2+y2}, lapply(f_ab1(x,a,b),function(x1){0.5*(1/theta[1]^2)*x1}), f_ab2(x,a,b,theta), SIMPLIFY=FALSE),
  #                 function(x3){exp(-x3)})
  #   return(fab)
  # }
  #
  # # New log-likelihood equation
  # loglik <- function(theta){
  #   # neg sum of LN of integral (double or triple) of lognpdfxmvnpdf
  #   N*log(theta[1]) + 1.5*N*log(2*pi) + N*log(theta[4]*theta[5]*(1 - (theta[6]^2))) - Reduce("+",lapply(f_ab(dataGROUP2,a_mat,b_mat,theta),function(x){sum(log(x*h_mat))}))
  # }



  # return(list(loglik,LSQest,a_mat,b_mat,f_ab))

  # Compute the boundaries for output

  crit <- qnorm((1 + conf.level)/2)
  crit2 <- qnorm(conf.level)
  conflim<-vector(mode = "list", length = length(theta.hat))
  fulllimset<-vector(mode = "list", length = length(theta.hat))

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
    fulllimset[[i]]<-c(theta.hat[i],conflim[[i]])
  }

  fulllimset[[1]] <- c(exp(fulllimset[[1]][1]),sort(exp(fulllimset[[1]][2:3])))

  return(list(MLE.theta.hat.perunit,MLE.inv.fish.perunit,theta.hat,inv.fish,params_txt,fulllimset))


  params_txt<-c(distparam_txt,params_txt)

  # Plot MLE distribution band
  # Processing Time
  if(dl=="Logarithmic" || dl=="LloydLipow" || dl=="CrackProp1"){
    timefit <- linspace(0.1,max(data[,1]),1000)
  } else{
    timefit <- linspace(0,max(data[,1]),1000)
  }
  # Processing parameter samples
  samples1 <- mvrnorm(1000,theta.hat,inv.fish)

  if (dist=="Lognormal"){

  }

  # samples1[,1]<- exp(samples1[,1])
  # samples1[,2]<- exp(samples1[,2])

  # Produce some output text that summariZes the results
  cat(c("Maximum-Likelihood estimates for the ",dl_txt,"-",dist_txt," Degradation model.\n\nD(t) = ",Dt_txt,"\n\n"),sep = "")
  print(matrix(unlist(fulllimset), nrow = length(unlist(fulllimset))/length(LSQest), ncol = length(LSQest), byrow = FALSE,dimnames = list(c("Degradation Model Parameters Mean",conflim_txt),params_txt)))
  return(list(theta.hat,inv.fish,fulllimset))
}
