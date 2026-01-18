# Maximum Likelihood Degradation Life Estimator
# Developed by Dr. Reuel Smith, 2021-2025

degradationlife.MLEest <- function(data,dl,dist="Normal",pp="Blom",D0,modelstress=NULL,Suse=NULL,
                                   confid=0.95,sided="twosided",
                                   xlabel="Time",ylabel="Degradation",Q=20,
                                   stressunit1 = NULL, stressunit2 = NULL){
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

  if(is.null(modelstress)==TRUE){ # RCS 01152026: Error message until we reactivate single stress ADT analysis
    stop("RMT 1.5.0.1 currently requires a modelstress entry for degradationlife.MLEest.")
  }

  # NOTE (11/10/2025): Beginning a massive effort to fix the degradationlife.MLEest tool for RAMS.
  # PHASE 1: Redefine the input/output parameters for a simple single stress-level MLE case.
  # Input:
  #          data
  #          dl = the degradation life model
  #          dist = "Normal" or "Lognormal"
  #          pp
  #          D0 = endurance limit
  #          conf = confidence
  #          sided
  # Output:
  #          parameter estimates and confidence of model and life
  #          loglik and likelihood
  #          AIC and BIC
  #          A plot of degradation and life distribution
  # PHASE 2: Add param-stress functions for multi-stress cases.  Stick to direct parameter-stress models only for now
  # with AF based models to come later.  Define,
  # New Input:
  #          modelstress = parameter-stress model
  #          Q = the number of weights being used (default between 20 and 40)
  # PHASE 3: Add multiplication error code
  # New Input:
  #          ME
  # Comment everything carefully and block out sections of the code so that it is easier to modify if the users wish

  # Legend colors
  col_legend <- c("red","blue","darkgreen","violet","aquamarine","orange","pink","darkblue","lightgreen","yellow","green")
  # Legend shapes
  shape_legend <- c(0:25)

  # UPDATE (11/19/2025) Check to see if any time values (column 1) are zero with respect to
  # Power, Logarithmic, Lloyd-Lipow, and Mitsuom degradation life models (Add more where ln(L) is involved)
  data0 <- data
  if((dl=="Power" || dl=="Logarithmic" || dl=="LloydLipow" || dl=="Mitsuom") && min(data[,1])==0){
    # Nix time data at zero
    data <- data0[which(data0[,1]!=0),]
  }

  # Set confidence
  conf.level <- confid

  # Start with the pre-processing of the data in which you obtain an initial parameter
  # estimate based on the curve fit of the data.  Also used to obtain initial Sigma and
  # bounds for calculating the integral based likelihood
  # NEW - Now consider model stress for parameters and include that as output
  if(is.null(modelstress) == TRUE){ # no model stress consideration.  Mainly for single stress level degradation
    adtLSQ<-degradationlife.LSQest(data=data,dl=dl,dist=dist,pp=pp,D0=D0,modelstress=modelstress,xlabel=xlabel,ylabel=ylabel,Suse=Suse)[[1]]
  }
  if(is.null(modelstress) == FALSE){ # model stress consideration for multi-accelerated stress degradation
    adtLSQ<-degradationlife.LSQest(data=data,dl=dl,dist=dist,pp=pp,D0=D0,modelstress=modelstress,xlabel=xlabel,ylabel=ylabel,Suse=Suse)[[1]]
    adtLSQ2<-degradationlife.LSQest(data=data,dl=dl,dist=dist,pp=pp,D0=D0,modelstress=modelstress,xlabel=xlabel,ylabel=ylabel,Suse=Suse)[[2]]
  }
  # return(adtLSQ)

  # UPDATE (2/12/2024) Now check to see if any time values (column 1) are zero with respect to
  # Power, Logarithmic, Lloyd-Lipow, and Mitsuom degradation life models (Add more where ln(L) is involved)
  data0 <- data
  if((dl=="Power" || dl=="Logarithmic" || dl=="LloydLipow" || dl=="Mitsuom") && min(data[,1])==0){
    data <- data0[which(data0[,1]!=0),] # Nix time data at zero
  }

  # Pulls the unit designations from column 3 of the input data
  unitnames <- unique(data[,3])

  # # Pulls the stresses from column 4 of the input data
  # unitstress <- unique(data[,4])
  # Pulls stress from columns 4 and up
  stresscount<-dim(data)[2]-3 # count the number of stress types (temperature, RH, etc.)
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
  # return(stressvals)

  # Pull the time and degradation data from the input
  TimeFULL <- data[,1]
  DegradationFULL <- data[,2]
  StresSFULL <- data[,4:dim(data)[2]]
  N <- length(data[,1])
  TimeLIST <- vector(mode = "list", length = length(unitnames))
  DegradationLIST <- vector(mode = "list", length = length(unitnames))
  StresSLIST <- vector(mode = "list", length = length(unitnames))

  # Form Matrix of Time, Degradation, and Stress with cbind
  # for(i in 2:length(unitnames)){
  #   if(i == 2){
  #     TimeMAT <- cbind(TimeFULL[which(data[,3]==unitnames[1])],TimeFULL[which(data[,3]==unitnames[2])])
  #     DegradationMAT <- cbind(DegradationFULL[which(data[,3]==unitnames[1])],DegradationFULL[which(data[,3]==unitnames[2])])
  #     StessMAT <- cbind(StresSFULL[which(data[,3]==unitnames[1])],StresSFULL[which(data[,3]==unitnames[2])])
  #   } else{
  #     TimeMAT <- cbind(TimeMAT,TimeFULL[which(data[,3]==unitnames[i])])
  #     DegradationMAT <- cbind(DegradationMAT,DegradationFULL[which(data[,3]==unitnames[i])])
  #     StessMAT <- cbind(StessMAT,StresSFULL[which(data[,3]==unitnames[i])])
  #   }
  # }
  #
  # for(i in 1:length(unitnames)){
  #   TimeLIST[[i]] <- TimeFULL[which(data[,3]==unitnames[i])]
  #   DegradationLIST[[i]] <- DegradationFULL[which(data[,3]==unitnames[i])]
  #   StresSLIST[[i]] <- StresSFULL[which(data[,3]==unitnames[i])]
  # }

  if(is.null(modelstress) == FALSE){ # If modelstress is active, we need to find the z values and weights per Gauss-Hermite integration
    # Suggested by Justin Ryan, this provides a closed form function to compute the MLE and Bayesian calculations
    # Weights need to be calculated though for Q instances (usually 20-40 but I will make 40 the default)

    # subfunction to generate weight vector given any size Q
    GH_weight <- function(Q){ # Gauss–Hermite quadrature weight
      poly_Hermite <- rep(0,(Q+1)) # Set up the Hermite polynomial https://en.wikipedia.org/wiki/Hermite_polynomials
      poly_Hermite[1] <- 2^Q       # Set up the first Hermite polynomial coefficient https://en.wikipedia.org/wiki/Hermite_polynomials

      for(i in 1:floor(Q/2)){ # Form the remaining Hermite polynomial coefficients
        poly_Hermite[(2*i+1)] <- ((((-1)^i)*fact(Q))/(fact(i)*fact(Q - 2*i)))*(2^(Q - 2*i))
      }

      z_i <- sort(roots(poly_Hermite)) # Roots obtain the z values
      w_i <- exp(-z_i^2)               # Weights w_i = exp(-Z^2)

      return(list(z_i,w_i)) # Return z-vector and weight vector
    }

    Z_set <- GH_weight(Q)[[1]] # Pull z vector per element
    W_set <- GH_weight(Q)[[2]] # pull w vector per element

    # Form Q x N matrices for Time, Degradation, Stress, Z, and W
    Time_MAT <- matrix(data = rep(TimeFULL,Q), nrow = Q, ncol = N,byrow = TRUE)
    Degradation_MAT <- matrix(data = rep(DegradationFULL,Q), nrow = Q, ncol = N,byrow = TRUE)
    Stress_MAT <- matrix(data = rep(StresSFULL,Q), nrow = Q, ncol = N,byrow = TRUE)
    Z_MAT <- matrix(data = rep(Z_set,N), nrow = Q, ncol = N,byrow = FALSE)
    W_MAT <- matrix(data = rep(W_set,N), nrow = Q, ncol = N,byrow = FALSE)
  }
  # return(list(Time_MAT,Degradation_MAT,Stress_MAT,Z_MAT,W_MAT))

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
    if(modelstress == "Linear" || modelstress == "LinearD0" || modelstress == "Exponential" || modelstress == "Exponential2" ||
       modelstress == "Arrhenius" || modelstress == "Eyring" || modelstress == "Eyring2" || modelstress == "Eyring3" || modelstress == "Eyring4" ||
       modelstress=="Power" || modelstress=="InversePower" || modelstress=="InversePower2" || modelstress=="Logarithmic" ||
       modelstress=="MultiStress" || modelstress=="TempHumidity" || modelstress=="TempNonthermal"){
      modelstresstype <- 1 # Intercept replacement
    }
    # if(modelstress == "LinearA" || modelstress == "ExponentialA" || modelstress == "Exponential2A" ||
    #    modelstress == "ArrheniusA" || modelstress == "EyringA" || modelstress == "Eyring2A" || modelstress == "Eyring3A" || modelstress == "Eyring4A" ||
    #    modelstress=="PowerA" || modelstress=="InversePowerA" || modelstress=="InversePower2A" || modelstress=="LogarithmicA" ||
    #    modelstress=="MultiStressA" || modelstress=="TempHumidityA" || modelstress=="TempNonthermalA"){
    #   modelstresstype <- 2
    # }
    # if(modelstress == "LinearB" || modelstress == "ExponentialB" || modelstress == "Exponential2B" ||
    #    modelstress == "ArrheniusB" || modelstress == "EyringB" || modelstress == "Eyring2B" || modelstress == "Eyring3B" || modelstress == "Eyring4B" ||
    #    modelstress=="PowerB" || modelstress=="InversePowerB" || modelstress=="InversePower2B" || modelstress=="LogarithmicB" ||
    #    modelstress=="MultiStressB" || modelstress=="TempHumidityB" || modelstress=="TempNonthermalB"){
    #   modelstresstype <- 3
    # }
    # if(modelstress == "LinearAF" || modelstress == "ExponentialAF" || modelstress == "ExponentialAF2" ||
    #    modelstress == "ArrheniusAF" || modelstress == "EyringAF" || modelstress == "EyringAF2" || modelstress == "EyringAF3" || modelstress == "EyringAF4"  ||
    #    modelstress=="PowerAF" || modelstress=="InversePowerAF" || modelstress=="InversePowerAF2" || modelstress=="LogarithmicAF" ||
    #    modelstress=="MultiStressAF" || modelstress=="TempHumidityAF" || modelstress=="TempNonthermalAF"){
    #   modelstresstype <- 5
    # }
  } else{
    modelstresstype <- 0
  }

  # ===================================================================
  # Establish the model-stress relations for parameters (usually based on the slope of linearized degradation life model)
  # ===================================================================
  if (is.null(modelstress) == FALSE && (modelstress=="Linear")){
    # psparams[1] - parameter a_0, psparams[2] - parameter b_0
    ishift2<-2                            # Second ishift for number of parameters
    LSQest3 <- adtLSQ2[1:2]               # Initial LSQ estimate for parameters a_0 and b_0
    paramstress <- function(theta,S) {    # Parameter-stress relation
      theta[2] + S*theta[1]
    }
    logparamstress <- function(theta,S) { # Log-parameter-stress relation
      log(theta[2] + S*theta[1])
    }
    # Writeup for the output text
    ps_txt<-modelstress
    params_txt<-c("a_0","b_0")
    params_txt3<-c("a_0","b_0")
    params_txt2<-"(b_0 + S*a_0)"
    logparam_txt<-"ln(b_0 + S*a_0)"
  }

  if (is.null(modelstress) == FALSE && (modelstress=="Exponential")){
    # psparams[1] - parameter a_0, psparams[2] - parameter b_0, lsparams[3] - R^2
    ishift2<-2                               # Second ishift for number of parameters
    LSQest3 <- c(adtLSQ2[1],log(adtLSQ2[2])) # Initial LSQ estimate for parameters a_0 and logb_0
    paramstress <- function(theta,S) {       # Parameter-stress relation
      exp(theta[2])*exp(S*theta[1])
    }
    logparamstress <- function(theta,S) {    # Log-parameter-stress relation
      theta[2] + S*theta[1]
    }
    # Writeup for the output text
    ps_txt<-modelstress
    params_txt<-c("a_0","b_0")
    params_txt3<-c("a_0","b_0")
    param_txt2<-"b_0*exp(a_0*S)"
    logparam_txt<-"(log(b_0) + a_0*S)"
  }
  if (is.null(modelstress) == FALSE && (modelstress=="Exponential2")){
    # psparams[1] - parameter a_0, psparams[2] - parameter b_0
    ishift2<-2                               # Second ishift for number of parameters
    LSQest3 <- c(adtLSQ2[1],log(adtLSQ2[2])) # Initial LSQ estimate for parameters a_0 and logb_0
    paramstress <- function(theta,S) {       # Parameter-stress relation
      exp(theta[2])*exp(theta[1]/S)
    }
    logparamstress <- function(theta,S) {    # Log-parameter-stress relation
      theta[2] + theta[1]/S
    }
    # Writeup for the output text
    ps_txt<-modelstress
    params_txt<-c("a_0","b_0")
    params_txt3<-c("a_0","b_0")
    param_txt2<-"b_0*exp(a_0/S)"
    logparam_txt<-"(log(b_0) + a_0/S)"
  }
  if (is.null(modelstress) == FALSE && (modelstress=="Arrhenius")){
    # psparams[1] - parameter E_a_0, psparams[2] - parameter b_0
    # Temperature HAS to be in Kelvin for this to work
    ishift2<-2                               # Second ishift for number of parameters
    LSQest3 <- c(adtLSQ2[1],log(adtLSQ2[2])) # Initial LSQ estimate for parameters E_a_0 and logb_0
    K<-8.617385e-5                           # Boltzmann Constant
    paramstress <- function(theta,S) {       # Parameter-stress relation
      exp(theta[2])*exp(theta[1]/(K*S))
    }
    logparamstress <- function(theta,S) {    # Log-parameter-stress relation
      theta[2] + theta[1]/(K*S)
    }
    # Writeup for the output text
    ps_txt<-modelstress
    params_txt<-c("E_a_0","b_0")
    params_txt3<-c("E_a_0","b_0")
    param_txt2<-"b_0*exp(E_a_0)/(K*S))"
    logparam_txt<-"(log(b_0) + (E_a_0/(K*S)))"
  }
  if (is.null(modelstress) == FALSE && (modelstress=="Eyring")){
    # lsparams[1] - parameter a_0, lsparams[2] - parameter b_0
    # First redefine b parameter as logb
    ishift2<-2                               # Second ishift for number of parameters
    LSQest3 <- c(adtLSQ2[1],log(adtLSQ2[2])) # Initial LSQ estimate for parameters a_0 and logb_0
    paramstress <- function(theta,S) {       # Parameter-stress relation
      (exp(theta[2])/S)*exp(theta[1]/S)
    }
    logparamstress <- function(theta,S) {    # Log-parameter-stress relation
      theta[2] - log(S) + theta[1]/S
    }
    # Writeup for the output text
    ps_txt<-modelstress
    params_txt<-c("a_0","b_0")
    params_txt3<-c("a_0","b_0")
    param_txt2<-"(b_0/S)*exp(a_0/S)"
    logparam_txt<-"(log(b_0) - log(S) + (a_0/S))"
  }
  if (is.null(modelstress) == FALSE && (modelstress=="Eyring2")){
    # lsparams[1] - parameter a_0, lsparams[2] - parameter b_0, lsparams[3] - R^2
    ishift2<-2                               # Second ishift for number of parameters
    LSQest3 <- adtLSQ2[1:2]                  # Initial LSQ estimate for parameters a_0 and logb_0
    paramstress <- function(theta,S) {       # Parameter-stress relation
      (1/S)*exp(-(theta[1] - (theta[2]/S)))
    }
    logparamstress <- function(theta,S) {    # Log-parameter-stress relation
      -log(S) - theta[1] + theta[2]/S
    }
    # Writeup for the output text
    ps_txt<-"Eyring (Type-2)"
    params_txt<-c("a_0","b_0")
    params_txt3<-c("a_0","b_0")
    param_txt2<-"(1/S)*exp(-(a_0 - (b_0/S)))"
    logparam_txt<-"(-log(S) - a_0 + (b_0/S))"
  }
  if (is.null(modelstress) == FALSE && (modelstress=="Power")){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # First redefine b parameter as logb
    ishift2<-2                               # Second ishift for number of parameters
    LSQest3 <- c(adtLSQ2[1],log(adtLSQ2[2])) # Initial LSQ estimate for parameters a_0 and logb_0
    paramstress <- function(theta,S) {       # Parameter-stress relation
      exp(theta[2])*(S^theta[1])
    }
    logparamstress <- function(theta,S) {    # Log-parameter-stress relation
      theta[2] + theta[1]*log(S)
    }
    # Writeup for the output text
    ps_txt<-modelstress
    params_txt<-c("a_0","b_0")
    params_txt3<-c("a_0","b_0")
    param_txt2<-"b_0*(S^a_0)"
    logparam_txt<-"(ln(b_0) + a_0ln(S))"
  }
  if (is.null(modelstress) == FALSE && (modelstress=="InversePower")){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    ishift2<-2                               # Second ishift for number of parameters
    LSQest3 <- c(adtLSQ2[1],log(adtLSQ2[2])) # Initial LSQ estimate for parameters a_0 and logb_0
    paramstress <- function(theta,S) {       # Parameter-stress relation
      exp(theta[2])*(S^-theta[1])
    }
    logparamstress <- function(theta,S) {    # Log-parameter-stress relation
      theta[2] - theta[1]*log(S)
    }
    # Writeup for the output text
    ps_txt<-"Inverse Power"
    params_txt<-c("a_0","b_0")
    params_txt3<-c("a_0","b_0")
    param_txt2<-"b_0*(S^-a_0)"
    logparam_txt<-"(ln(b_0) - a_0ln(S))"
  }
  if (is.null(modelstress) == FALSE && (modelstress=="InversePower2")){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    ishift2<-2                               # Second ishift for number of parameters
    LSQest3 <- c(adtLSQ2[1],log(adtLSQ2[2])) # Initial LSQ estimate for parameters a_0 and logb_0
    paramstress <- function(theta,S) {       # Parameter-stress relation
      1/(exp(theta[2])*(S^theta[1]))
    }
    logparamstress <- function(theta,S) {    # Log-parameter-stress relation
      -theta[2] - theta[1]*log(S)
    }
    # Writeup for the output text
    ps_txt<-"Inverse Power"
    params_txt<-c("a_0","b_0")
    params_txt3<-c("a_0","b_0")
    param_txt2<-"1/[b_0*(S^a_0)]"
    logparam_txt<-"(-ln(b_0) - a_0*ln(S))"
  }

  if (is.null(modelstress) == FALSE && (modelstress=="Logarithmic")){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    ishift2<-2                               # Second ishift for number of parameters
    LSQest3 <- adtLSQ2[1:2]                  # Initial LSQ estimate for parameters a_0 and b_0
    paramstress <- function(theta,S) {       # Parameter-stress relation
      theta[1]*log(S) + theta[2]
    }
    logparamstress <- function(theta,S) {    # Log-parameter-stress relation
      log(theta[1]*log(S) + theta[2])
    }
    # Writeup for the output text
    ps_txt<-modelstress
    params_txt<-c("a_0","b_0")
    params_txt3<-c("a_0","b_0")
    param_txt2<-"(b_0 + a_0*ln(S))"
    logparam_txt<-"ln(b_0 + a_0*ln(S))"
  }

  if (is.null(modelstress) == FALSE && (modelstress=="MultiStress")){
    # lsparams - c(ao,a1,...,an), S - c(1,S1,...,Sn)
    ishift2<-dim(S)[2]
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

  if (is.null(modelstress) == FALSE && (modelstress=="TempHumidity")){
    # psparams[1] - parameter A_0, psparams[2] - parameter a_0, psparams[3] - parameter b_0
    ishift2<-3                                   # Second ishift for number of parameters
    LSQest3 <- c(log(adtLSQ2[1]),adtLSQ2[2:3])   # Initial LSQ estimate for parameters logA_0, a_0, and b_0
    paramstress <- function(theta,S) {           # Parameter-stress relation
      exp(theta[1])*exp((theta[2]/S[,1]) + (theta[3]/S[,2]))
    }
    logparamstress <- function(theta,S) {        # Log-parameter-stress relation
      theta[1] + (theta[2]/S[,1]) + (theta[3]/S[,2])
    }
    # Writeup for the output text
    ps_txt<-"Temperature-Humidity"
    params_txt3<-c("A_0","a_0","b_0")
    param_txt2<-"A_0 exp(a_0/S + b_0/H)"
    logparam_txt<-"ln(A_0) + a_0/S + b_0/H"
  }

  if (is.null(modelstress) == FALSE && (modelstress=="TempNonthermal")){
    # psparams[1] - parameter a_0, psparams[2] - parameter b_0, psparams[3] - parameter c_0
    ishift2<-3                                   # Second ishift for number of parameters
    LSQest3 <- c(adtLSQ2[1:2],log(adtLSQ2[3]))   # Initial LSQ estimate for parameters a_0, b_0, and logc_0
    paramstress <- function(theta,S) {           # Parameter-stress relation
      exp(theta[3])/((S[,2]^theta[2])*exp(-theta[1]/S[,1]))
    }
    logparamstress <- function(theta,S) {        # Log-parameter-stress relation
      theta[3] - theta[2]*log(S[,2]) + (theta[1]/S[,1])
    }
    # Writeup for the output text
    ps_txt<-"Temperature-Non-thermal"
    params_txt3<-c("a_0","b_0","c_0")
    param_txt2<-"c_0/(U^b_0 * exp(-a_0/S))"
    logparam_txt<-"a_0(1/S) - b_0*ln(U) + ln(c_0)"
  }

  if (is.null(modelstress) == FALSE && (modelstress=="Eyring3")){
    # psparams[1] - parameter a_0, psparams[2] - parameter b_0
    # psparams[3] - parameter c_0, psparams[4] - parameter d_0
    ishift2<-4                                   # Second ishift for number of parameters
    LSQest3 <- adtLSQ2                           # Initial LSQ estimate for parameters a_0, b_0, c_0, and d_0
    paramstress <- function(theta,S) {           # Parameter-stress relation
      (1/S[,1])*exp((theta[1] + (theta[2]/S[,1])) + (theta[3] + (theta[4]/S[,1]))*S[,2])
    }
    logparamstress <- function(theta,S) {        # Log-parameter-stress relation
      -log(S[,1]) + theta[1] + (theta[2]/S[,1]) + (theta[3] + (theta[4]/S[,1]))*S[,2]
    }
    # Writeup for the output text
    ps_txt<-"Eyring (Type 3)"
    params_txt3<-c("a_0","b_0","c_0","d_0")
    param_txt2<-"(1/S) exp((a_0 + (b_0/S)) + (c_0 + (d_0/S)) U)"
    logparam_txt<-"-ln(S) + (a_0 + (b_0/S)) + (c_0 + (d_0/S)) U"
  }

  if (is.null(modelstress) == FALSE && (modelstress=="Eyring4")){
    # psparams[1] - parameter A_0, psparams[2] - parameter b_0
    # psparams[3] - parameter E_a_0
    # Temperature HAS to be in Kelvin for this to work
    ishift2<-3                                   # Second ishift for number of parameters
    LSQest3 <- c(log(adtLSQ2[1]),adtLSQ2[2:3])   # Initial LSQ estimate for parameters logA_0, b_0, and E_a_0
    K<-8.617385e-5                               # Boltzmann Constant
    paramstress <- function(theta,S) {           # Parameter-stress relation
      theta[1]*exp(theta[3]/(K*S[,1]))*(S[,2]^-theta[2])
    }
    logparamstress <- function(theta,S) {        # Log-parameter-stress relation
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
  # positivity_v<-rep(0,length(LSQest))
  # if(dl=="Hamada"){
  #   positivity_v<-rep(0,dim(adtLSQ)[2]-1+6)
  # } else {
  #   positivity_v<-rep(0,dim(adtLSQ)[2]-1+3)
  # }

  # ===================================================================
  # Degradation-Life (degradationlife) models
  # ===================================================================
  if(dl=="Linear"){
    # D = a + b*t
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"D = a + b*t"
    degradationlife <- function(theta,Time) {
      theta[1] + Time*theta[2]
    }
    logdegradationlife <- function(theta,Time) {
      log(theta[1] + Time*theta[2])
    }
    degfit <- function(TimeDegfit,params){
      params[1] + TimeDegfit*params[2]
    }
    dl_txt<-dl
    dlmodel_txt0<-"a + bl"
    params_txt<-c("a","b")
    LSQest1 <- cbind(rep(0,length(unitnames)),adtLSQ[,1],adtLSQ[,2]) # Unit specific initial parameter estimates
    LSQest<-c(0,mean(adtLSQ[,1]), mean(adtLSQ[,2]))

    if(modelstresstype == 1){ # Replace parameter(S) = a or a(S)
      dlparams_txt<-c("a","b")
      # dlparams_txt<-c(params_txt,"b")
      paramref_txt<-"a(S) = "
      paramref_txt2<-"b ~ NOR(\U03BC_b,\U03C3_b)"
      params_txt2<-c("\U03BC_b","\U03C3_b")
      degradationlife1 <- function(theta,Time,S,Z) {
        (paramstress(theta[2:length(theta)],S)) + Time*(theta[ishift2+2] + exp(theta[ishift2+3])*Z*sqrt(2))
      }
      logdegradationlife1 <- function(theta,Time,S,Z) {
        log((paramstress(theta[2:length(theta)],S)) + Time*(theta[ishift2+2] + exp(theta[ishift2+3])*Z*sqrt(2)))
      }
      LSQest2<-c(0,LSQest3,adtLSQ2[3],log(adtLSQ2[4]))
      params_txt<-c("b",params_txt3)
    }
    # <==================================================================================================================
    # DEVELOPER NOTE (RCS 11/19/2025): When assigning stress to your degradation-life model replace the parameter of choice with
    # "paramstress(theta[2:length(theta)],S)".  The other parameter is to be replaced by
    # "(theta[ishift2+2] + exp(theta[ishift2+3])*Z*sqrt(2))".  The stress-based parameter is usually the intercept in
    # a linearized model and the variable parameter is usually the slope.
    # <==================================================================================================================
  }

  if(dl=="Exponential"){
    # D = b*exp(a*t)
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"b*exp(a*t)"
    degradationlife <- function(theta,Time) {
      exp(theta[2])*exp(Time*theta[1])
    }
    logdegradationlife <- function(theta,Time) {
      theta[2] + Time*theta[1]
    }
    degfit <- function(TimeDegfit,params){
      exp(params[2])*exp(TimeDegfit*params[1])
    }
    dl_txt<-dl
    dlmodel_txt0<-"b exp(al)"
    params_txt<-c("a","b")
    LSQest1 <- cbind(rep(0,length(unitnames)),adtLSQ[,1],log(adtLSQ[,2])) # Unit specific initial parameter estimates
    LSQest<-c(0,mean(adtLSQ[,1]), mean(log(adtLSQ[,2])))
    sigparamno<-3

    if(modelstresstype == 1){ # Replace parameter(S) = b or b(S)
      dlparams_txt<-c("a","b")
      paramref_txt<-"b(S) = "
      paramref_txt2<-"a ~ NOR(\U03BC_a,\U03C3_a)"
      params_txt2<-c("\U03BC_a","\U03C3_a")
      degradationlife1 <- function(theta,Time,S,Z) {
        (paramstress(theta[2:length(theta)],S))*exp(Time*(theta[ishift2+2] + exp(theta[ishift2+3])*Z*sqrt(2)))
      }
      logdegradationlife1 <- function(theta,Time,S,Z) {
        logparamstress(theta[2:length(theta)],S) + Time*(theta[ishift2+2] + exp(theta[ishift2+3])*Z*sqrt(2))
      }
      LSQest2<-c(0,LSQest3,adtLSQ2[3],log(adtLSQ2[4]))
      params_txt<-c(params_txt3,"\U03BC_b","\U03C3_b")
    }
  }

  if(dl=="SquareRoot"){
    # D^(1/2) = a + b*t
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"(a + b*t)^2"

    # NOTE 11/25/2024 - Check for general slope (positive or negative)
    # NOTE 11/26/2024 - Also check for whether intercept is generally positive or if it can include negative.
    if (mean(adtLSQ[,2]) > 0 && min(adtLSQ[,1]) >= 0){ # Positive slope and positive intercept
      degradationlife <- function(theta,Time) { # degradation-life
        (exp(theta[1]) + Time*exp(theta[2]))^2
      }
      logdegradationlife <- function(theta,Time) { # Log degradation-life
        2*log(exp(theta[1]) + Time*exp(theta[2]))
      }
      degfit <- function(TimeDegfit,params){ # degradation-life fit
        (exp(params[1]) + TimeDegfit*exp(params[2]))^2
      }
      LSQest1 <- cbind(rep(0,length(unitnames)),log(adtLSQ[,1]),log(adtLSQ[,2])) # Unit specific initial parameter estimates
      LSQest<-c(0,log(mean(adtLSQ[,1])), log(mean(adtLSQ[,2])))
    }
    if (mean(adtLSQ[,2]) > 0  && min(adtLSQ[,1]) < 0){ # Positive slope and positive/negative intercept
      degradationlife <- function(theta,Time) { # degradation-life
        (theta[1] + Time*exp(theta[2]))^2
      }
      logdegradationlife <- function(theta,Time) { # Log degradation-life
        2*log(theta[1] + Time*exp(theta[2]))
      }
      degfit <- function(TimeDegfit,params){ # degradation-life Fit
        (params[1] + TimeDegfit*exp(params[2]))^2
      }
      LSQest1 <- cbind(rep(0,length(unitnames)),adtLSQ[,1],log(adtLSQ[,2])) # Unit specific initial parameter estimates
      LSQest<-c(0,mean(adtLSQ[,1]), log(mean(adtLSQ[,2])))
    }
    if (min(adtLSQ[,2]) < 0){ # Negative slope
      degradationlife <- function(theta,Time) { # degradation-life
        (theta[1] + Time*theta[2])^2
      }
      logdegradationlife <- function(theta,Time) { # Log degradation-life
        2*log(theta[1] + Time*theta[2])
      }
      degfit <- function(TimeDegfit,params){ # degradation-life Fit
        (params[1] + TimeDegfit*params[2])^2
      }
      LSQest1 <- cbind(rep(0,length(unitnames)),adtLSQ[,1],adtLSQ[,2]) # Unit specific initial parameter estimates
      LSQest<-c(0,mean(adtLSQ[,1]), mean(adtLSQ[,2]))
    }
    dl_txt<-"Square-Root"
    dlmodel_txt0<-"(a + bl)\U00B2"
    params_txt<-c("a","b")

    if(modelstresstype == 1){ # Replace parameter(S) = a or a(S)
      dlparams_txt<-c("a","b")
      paramref_txt<-"a(S) = "
      paramref_txt2<-"b ~ NOR(\U03BC_b,\U03C3_b)"
      params_txt2<-c("\U03BC_b","\U03C3_b")
      if (mean(adtLSQ[,2]) > 0 && min(adtLSQ[,1]) >= 0){ # Positive slope
        degradationlife1 <- function(theta,Time,S,Z) { # degradationlife
          (paramstress(theta,S) + Time*exp(A))^2
        }
        logdegradationlife1 <- function(theta,Time,S,Z) { # Log degradationlife
          2*log(paramstress(theta,S) + Time*(theta[ishift2+2] + exp(theta[ishift2+3])*Z*sqrt(2)))
        }
        LSQest2<-c(0,LSQest3,adtLSQ2[3],log(adtLSQ2[4]))
        # LSQest2<-c(0,LSQest3,adtLSQ2[length(adtLSQ2)-1],log(adtLSQ2[length(adtLSQ2)]))
      }
      if (min(adtLSQ[,2]) < 0){ # Negative slope
        degradationlife1 <- function(theta,Time,S,Z) { # degradationlife
          (paramstress(theta,S) + Time*(theta[ishift2+2] + exp(theta[ishift2+3])*Z*sqrt(2)))^2
        }
        logdegradationlife1 <- function(theta,Time,S,Z) { # Log degradationlife
          2*log(paramstress(theta,S) + Time*(theta[ishift2+2] + exp(theta[ishift2+3])*Z*sqrt(2)))
        }
        LSQest2<-c(0,LSQest3,adtLSQ2[3],log(adtLSQ2[4]))
        # LSQest2<-c(0,LSQest3,adtLSQ2[length(adtLSQ2)-1],log(adtLSQ2[length(adtLSQ2)]))
      }
      params_txt<-c(params_txt3,"\U03BC_b","\U03C3_b")
    }

  }

  if(dl=="Power"){
    # D = b*(t^a)
    # theta[1] ~ a, theta[2] ~ b
    # shift b to log b
    Dt_txt<-"b*(t^a)"
    degradationlife <- function(theta,Time) {    # Degradation-life model
      exp(theta[2])*(Time^theta[1])
    }
    logdegradationlife <- function(theta,Time) { # log-Degradation-life model
      theta[2] + theta[1]*log(Time)
    }
    degfit <- function(TimeDegfit,params){
      params[2]*(TimeDegfit^params[1])
    }
    dl_txt<-dl
    dlmodel_txt0<-"bl\U00AA"
    params_txt<-c("a","b")
    LSQest1 <- cbind(rep(0,length(unitnames)),adtLSQ[,1],log(adtLSQ[,2])) # Unit specific initial parameter estimates
    LSQest<-c(0,mean(adtLSQ[,1]), mean(log(adtLSQ[,2])))                  # General initial parameter estimate
    sigparamno<-3

    if(modelstresstype == 1){ # Replace parameter(S) = b or b(S)
      dlparams_txt<-c("a","b")
      paramref_txt<-"b(S) = "
      paramref_txt2<-"a ~ NOR(\U03BC_a,\U03C3_a)"
      params_txt2<-c("\U03BC_a","\U03C3_a")
      degradationlife1 <- function(theta,Time,S,Z) {
        (paramstress(theta[2:length(theta)],S))*(Time^(theta[ishift2+2] + exp(theta[ishift2+3])*Z*sqrt(2)))
      }
      logdegradationlife1 <- function(theta,Time,S,Z) {
        logparamstress(theta[2:length(theta)],S) + (theta[ishift2+2] + exp(theta[ishift2+3])*Z*sqrt(2))*log(Time)
      }
      lifedegradation1 <- function(theta,Deg,S,Z) {
        (Deg/(paramstress(theta[2:length(theta)],S)))^(1/(theta[ishift2+2] + exp(theta[ishift2+3])*Z*sqrt(2)))
      }
      LSQest2<-c(0,LSQest3,adtLSQ2[3],log(adtLSQ2[4]))
      params_txt<-c(params_txt3,"\U03BC_a","\U03C3_a")
    }
  }
  # return(list(LSQest2,params_txt))
  # return(list(loglik,LSQest,loglik(LSQest),Dij,Lij))

  if(dl=="Logarithmic"){
    # D = a + b*ln(t)
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"a + b*ln(t)"
    degradationlife <- function(theta,Time) {
      theta[1] + theta[2]*log(Time)
    }
    logdegradationlife <- function(theta,Time) {
      log(theta[1] + theta[2]*log(Time))
    }
    degfit <- function(TimeDegfit,params){
      params[1] + params[2]*log(TimeDegfit)
    }

    dl_txt<-dl
    dlmodel_txt0<-"a + b ln(l)"
    params_txt<-c("a","b")
    LSQest1 <- cbind(rep(0,length(unitnames)),adtLSQ[,1],adtLSQ[,2]) # Unit specific initial parameter estimates
    LSQest<-c(0,mean(adtLSQ[,1]), mean(adtLSQ[,2]))
    sigparamno<-3

    if(modelstresstype == 1){ # Replace parameter(S) = a or a(S)
      dlparams_txt<-c("a","b")
      paramref_txt<-"a(S) = "
      paramref_txt2<-"b ~ NOR(\U03BC_b,\U03C3_b)"
      params_txt2<-c("\U03BC_b","\U03C3_b")
      degradationlife1 <- function(theta,Time,S) {
        (paramstress(theta[2:length(theta)],S,Z)) + (theta[ishift2+2] + exp(theta[ishift2+3])*Z*sqrt(2))*log(Time)
      }
      logdegradationlife1 <- function(theta,Time,S,Z) {
        log((paramstress(theta[2:length(theta)],S)) + (theta[ishift2+2] + exp(theta[ishift2+3])*Z*sqrt(2))*log(Time))
      }
      LSQest2<-c(0,LSQest3,adtLSQ2[3],log(adtLSQ2[4]))
      params_txt<-c("b",params_txt3)
    }
  }

  # if(dl=="Gompertz"){
  #   # D = a + b^(c*t)
  #   # theta[1] ~ a, theta[2] ~ b, theta[3] ~ c
  #   Dt_txt<-"a + b^(c*t)"
  #   degradationlife <- function(theta) {
  #     theta[1] + theta[2]*log(TimeDam)
  #   }
  #   logdegradationlife <- function(theta) {
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
    degradationlife <- function(theta,Time) {
      theta[1] - (theta[2]/Time)
    }
    logdegradationlife <- function(theta,Time) {
      log(theta[1] - (theta[2]/Time))
    }
    degfit <- function(TimeDegfit,params){
      params[1] - (params[2]/TimeDegfit)
    }
    dl_txt<-dl
    dlmodel_txt0<-"a - b/l"
    params_txt<-c("a","b")
    LSQest1 <- cbind(rep(0,length(unitnames)),adtLSQ[,1],adtLSQ[,2]) # Unit specific initial parameter estimates
    LSQest<-c(1,mean(adtLSQ[,1]), mean(adtLSQ[,2]))
    sigparamno<-3

    if(modelstresstype == 1){ # Replace parameter(S) = a or a(S)
      dlparams_txt<-c("a","b")
      # dlparams_txt<-c(params_txt,"b")
      paramref_txt<-"a(S) = "
      paramref_txt2<-"b ~ NOR(\U03BC_b,\U03C3_b)"
      params_txt2<-c("\U03BC_b","\U03C3_b")
      degradationlife1 <- function(theta,Time,S,Z) {
        (paramstress(theta[2:length(theta)],S)) - ((theta[ishift2+2] + exp(theta[ishift2+3])*Z*sqrt(2))/Time)
      }
      logdegradationlife1 <- function(theta,Time,S,Z) {
        log((paramstress(theta[2:length(theta)],S)) - ((theta[ishift2+2] + exp(theta[ishift2+3])*Z*sqrt(2))/Time))
      }
      LSQest2<-c(0,LSQest3,adtLSQ2[3],log(adtLSQ2[4]))
      params_txt<-c("b",params_txt3)
    }
  }

  if(dl=="Mitsuom"){
    # D = 1/(1 + b*(t^a))
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"1/(1 + b*(t^a))"
    degradationlife <- function(theta,Time) {
      1/(1 + theta[2]*(Time^theta[1]))
    }
    logdegradationlife <- function(theta,Time) {
      -log(1 + theta[2]*(Time^theta[1]))
    }
    degfit <- function(TimeDegfit,params){
      1/(1 + params[2]*(TimeDegfit^params[1]))
    }
    dl_txt<-dl
    dlmodel_txt0<-"(1 + bl\U00AA)\U207B\U00B9"
    params_txt<-c("a","b")
    LSQest1 <- cbind(rep(0,length(unitnames)),adtLSQ[,1],adtLSQ[,2]) # Unit specific initial parameter estimates
    LSQest<-c(1,mean(adtLSQ[,1]), mean(adtLSQ[,2]))
    sigparamno<-3

    if(modelstresstype == 1){ # Replace parameter(S) = b or b(S)
      dlparams_txt<-c("a","b")
      # dlparams_txt<-c("a",params_txt)
      paramref_txt<-"b(S) = "
      paramref_txt2<-"a ~ NOR(\U03BC_a,\U03C3_a)"
      params_txt2<-c("\U03BC_a","\U03C3_a")
      degradationlife1 <- function(theta,Time,S) {
        1/(1 + (paramstress(theta[2:length(theta)],S))*(Time^(theta[ishift2+2] + exp(theta[ishift2+3])*Z*sqrt(2))))
      }
      logdegradationlife1 <- function(theta,Time,S) {
        -log(1 + (paramstress(theta[2:length(theta)],S))*(Time^(theta[ishift2+2] + exp(theta[ishift2+3])*Z*sqrt(2))))
      }
      LSQest2<-c(0,LSQest3,adtLSQ2[3],log(adtLSQ2[4]))
      params_txt<-c(params_txt3,"\U03BC_a","\U03C3_a")
    }
  }

  if(dl=="Hamada"){
    # D = 1/(1 + beta1*(t*exp(beta3*11605*(1/Tu - 1/Ti)))^beta2)
    # theta[1] ~ beta1, theta[2] ~ beta2, theta[3] ~ beta3
    TimeDam
    Dt_txt<-"1/(1 + \U03B2_1*(t*exp(\U03B2_3*11605*(1/Tuse - 1/Temp)))^\U03B2_2)"
    degradationlife <- function(theta,Time) {
      1/(1 + theta[1]*(Time*exp(theta[3]*11605*(1/Tuse - 1/TempDam)))^theta[2])
    }
    logdegradationlife <- function(theta,Time) {
      -log(1 + theta[1]*(Time*exp(theta[3]*11605*(1/Tuse - 1/TempDam)))^theta[2])
    }
    degfit <- function(TimeDegfit,TempDegfit,params){
      1/(1 + params[1]*(TimeDegfit*exp(params[3]*11605*(1/Tuse - 1/TempDegfit)))^params[2])
    }
    dl_txt<-dl
    params_txt<-c("\U03B2_1","\U03B2_2","\U03B2_3")
    LSQest1 <- cbind(rep(0,length(unitnames)),adtLSQ[,1],adtLSQ[,2]) # Unit specific initial parameter estimates
    LSQest<-c(1,mean(adtLSQ[,1]), mean(adtLSQ[,2]), mean(adtLSQ[,3]))
    positivity_v[1]<-1
    positivity_v[2]<-1
    sigparamno<-4
  }
  # life degradationlife model "CrackProp1" is an assumed effective zero initial crack length
  # model adjusted for an initial crack length of 1 mm or 0.001 m.  It was used in the
  # ENRE 641 2022 Final Exam with an assumed DS of 200 MPa.
  if(dl=="CrackProp1"){
    # D = a_0(=0.001 m) + exp((2/(2-m))*LN(N) + (2/(2-m))*(LN(C*sqrt-pi) + LN(1 - 0.5*m) + LN(DS)))
    # theta[1] ~ logC, theta[2] ~ m (C is too small to control so take the natural log and process that)
    Dt_txt<-"a_0(=0.001 m) + exp((2/(2-m))*LN(N) + (2/(2-m))*(LN(C*sqrt-pi) + LN(1 - 0.5*m) + LN(DS)))"
    degradationlife <- function(theta) {
      0.001 + exp((2/(2 - theta[2]))*log(TimeDam) + (2/(2 - theta[2]))*(log(exp(theta[1])*sqrt(pi)) + log(1 - 0.5*theta[2]) + log(200)))
    }
    logdegradationlife <- function(theta) {
      log(0.001 + exp((2/(2 - theta[2]))*log(TimeDam) + (2/(2 - theta[2]))*(log(exp(theta[1])*sqrt(pi)) + log(1 - 0.5*theta[2]) + log(200))))
    }
    degfit <- function(TimeDegfit,params){
      0.001 + exp((2/(2 - params[2]))*log(TimeDegfit) + (2/(2 - params[2]))*(log(params[1]*sqrt(pi)) + log(1 - 0.5*params[2]) + log(200)))
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
    degradationlife <- function(theta) {
      0.001*exp(pi*exp(theta[1])*(200^2)*TimeDam)
    }
    logdegradationlife <- function(theta) {
      log(0.001) + (pi*exp(theta[1])*(200^2)*TimeDam)
    }
    degfit <- function(TimeDegfit,params){
      0.001*exp(pi*params[1]*(200^2)*TimeDegfit)
    }
    dl_txt<-dl
    params_txt<-c("C")
    LSQest<-c(0.1, mean(log(adtLSQ[,1])))
    sigparamno<-2
  }
  if(dl=="LinearD0"){
    # D = D0 + a*t
    # theta[1] ~ a
    # negate and log because it is relatively small

    Dt_txt<-"D = 8.5 + a*t"
    degradationlife <- function(theta,Time) {
      8.5 - Time*exp(theta[1])
    }
    logdegradationlife <- function(theta,Time) {
      log(8.5 - Time*exp(theta[1]))
    }
    degfit <- function(TimeDegfit,params){
      8.5 - TimeDegfit*exp(params[1])
    }
    dl_txt<-dl
    params_txt<-c("a")
    LSQest1 <- cbind(rep(0,length(unitnames)),log(-adtLSQ[,1]))
    LSQest<-c(0,mean(log(-adtLSQ[,1])))

    if(modelstresstype == 1){ # Replace parameter(S) = a or a(S)
      degradationlife1 <- function(theta,Time,S) {
        8.5 - Time*(paramstress(theta[2:length(theta)],S))
      }
      logdegradationlife1 <- function(theta,Time,S) {
        log(8.5 - Time*(paramstress(theta[2:length(theta)],S)))
      }
      LSQest2<-c(0,LSQest3)
      params_txt<-c(params_txt3)
    }
  }
  if(dl=="ExponentialNorm"){
    # D = 8.5*exp(a*(W/2500)*(t/5000))*exp(b*(T - 298)*(t/5000))
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"8.5*exp(a*(W/2500)*(t/5000))*exp(b*(T - 298)*(t/5000))"
    degradationlife <- function(theta,Time,Weight,Temp) {
      8.5*exp(-exp(theta[1])*(Weight/2500)*(Time/5000))*exp(-exp(theta[2])*(Temp - 298)*(Time/5000))
    }
    logdegradationlife <- function(theta,Time,Weight,Temp) {
      log(8.5) - exp(theta[1])*(Weight/2500)*(Time/5000) - exp(theta[2])*(Temp - 298)*(Time/5000)
    }
    degfit <- function(TimeDegfit,Weight,Temp,params){
      8.5*exp(-exp(params[1])*(Weight/2500)*(Time/5000))*exp(-exp(params[2])*(Temp - 298)*(TimeDegfit/5000))
    }
    dl_txt<-dl
    params_txt<-c("a","b")
    LSQest1 <- cbind(rep(0,length(unitnames)),log(-adtLSQ[,1]),log(-adtLSQ[,2]))
    LSQest<-c(0,mean(log(-adtLSQ[,1])), mean(log(-adtLSQ[,2])))
    sigparamno<-3

    if(modelstresstype == 1 || modelstresstype == 2 || modelstresstype == 5){ # Replace parameter(S) = b or b(S)
      degradationlife1 <- function(theta,Time,S) {
        8.5*exp(-exp((paramstress(theta[3:length(theta)],S)))*(S[,2]/2500)*(Time/5000))*exp(-exp(theta[2])*(S[,1] - 298)*(Time/5000))
      }
      logdegradationlife1 <- function(theta,Time,S) {
        log(8.5) - exp((paramstress(theta[3:length(theta)],S)))*(S[,2]/2500)*(Time/5000) - exp(theta[2])*(S[,1] - 298)*(Time/5000)
      }
      LSQest2<-c(0, mean(log(-adtLSQ[,2])),LSQest3)
      # params_txt<-c("b",params_txt3,"\U03BC_b","\U03C3_b")
      params_txt<-c("b",params_txt3)

    }
    if(modelstresstype == 3){ # Replace parameter(S) = a or a(S)
      degradationlife1 <- function(theta,Time,S) {
        (paramstress(theta[2:length(theta)],S)) + Time*theta[1]
      }
      logdegradationlife1 <- function(theta,Time,S) {
        log((paramstress(theta[2:length(theta)],S)) + Time*theta[1])
      }
      LSQest2<-c(0,mean(adtLSQ[,2]),LSQest3)
      # params_txt<-c("b",params_txt3,"\U03BC_a","\U03C3_a")
      params_txt<-c("b",params_txt3)
    }
  }
  if(dl=="SquareRootD0"){
    # D^(1/2) = sqrt(8.5) + a*t
    # theta[1] ~ a
    # negate and log because it is relatively small

    Dt_txt<-"D = (sqrt(8.5) + a*t)^2"
    degradationlife <- function(theta,Time) {
      (sqrt(8.5) - Time*exp(theta[1]))^2
    }
    logdegradationlife <- function(theta,Time) {
      2*log(sqrt(8.5) - Time*exp(theta[1]))
    }
    degfit <- function(TimeDegfit,params){
      (sqrt(8.5) - TimeDegfit*exp(params[1]))^2
    }
    dl_txt<-dl
    params_txt<-c("a")
    LSQest1 <- cbind(rep(0,length(unitnames)),log(-adtLSQ[,1]))
    LSQest<-c(0,mean(log(-adtLSQ[,1])))

    if(modelstresstype == 1){ # Replace parameter(S) = a or a(S)
      degradationlife1 <- function(theta,Time,S) {
        (sqrt(8.5) - Time*(paramstress(theta[2:length(theta)],S)))^2
      }
      logdegradationlife1 <- function(theta,Time,S) {
        2*log(sqrt(8.5) - Time*(paramstress(theta[2:length(theta)],S)))
      }
      LSQest2<-c(0,LSQest3)
      params_txt<-c(params_txt3)
    }
  }
  # return(list(LSQest2,adtLSQ2,params_txt))

  ## NEW (2/27/2024) Set up log-likelihoods for different cases and scenarios.  The only constant
  ## will be unit by unit MLE evaluation
  # 1. Run MLE for each unit based on dl regardless of the modelstress
  # List of uncorrected MLE and variance-covariance output
  # MLE.theta.hat.perunit<-vector(mode = "list", length = length(unitnames))
  # MLE.inv.fish.perunit<-vector(mode = "list", length = length(unitnames))
  # for(i in 1:length(unitnames)){
  #   # Since each degradationlife is dependent on a separate time stamp, a different loglikelihood will be
  #   # formed for each part of the loop
  #   # RCS 09/05/2024 - ADD WEIBULL
  #   if (dist=="Weibull") {
  #     # Remember theta[1] is logβ
  #     loglik <- function(theta){
  #       -sum(theta[1] + (exp(theta[1])-1)*log(data[which(data[,3]==unitnames[i]),2]) - exp(theta[1])*logdegradationlife(theta[2:3],data[which(data[,3]==unitnames[i]),1]) - ((data[which(data[,3]==unitnames[i]),2]/degradationlife(theta[2:3],data[which(data[,3]==unitnames[i]),1]))^exp(theta[1])))
  #     }
  #     dist_txt<-dist
  #     distparam_txt<-"\U03B2"
  #   }
  #   if (dist=="Normal" && dl!="ExponentialNorm") {
  #     # Remember theta[1] is logσ
  #     loglik <- function(theta){
  #       -sum(-theta[1] - 0.5*log(2*pi) - 0.5*(exp(theta[1])^-2)*((data[which(data[,3]==unitnames[i]),2] - degradationlife(theta[2:3],data[which(data[,3]==unitnames[i]),1]))^2))
  #     }
  #     dist_txt<-dist
  #     distparam_txt<-"\U03C3"
  #   }
  #
  #   if (dist=="Normal" && dl=="ExponentialNorm") {
  #     # Remember theta[1] is logσ
  #     loglik <- function(theta){
  #       -sum(-theta[1] - 0.5*log(2*pi) - 0.5*(exp(theta[1])^-2)*((data[which(data[,3]==unitnames[i]),2] - degradationlife(theta[2:3],data[which(data[,3]==unitnames[i]),1],data[which(data[,3]==unitnames[i]),4],data[which(data[,3]==unitnames[i]),5]))^2))
  #     }
  #     dist_txt<-dist
  #     distparam_txt<-"\U03C3"
  #   }
  #   if (dist=="Lognormal" && dl!="ExponentialNorm") {
  #     # Remember theta[1] is logσ_t
  #     loglik <- function(theta){
  #       -sum(-theta[1] - 0.5*log(2*pi) - log(data[which(data[,3]==unitnames[i]),2]) - 0.5*(exp(theta[1])^-2)*((log(data[which(data[,3]==unitnames[i]),2]) - logdegradationlife(theta[2:3],data[which(data[,3]==unitnames[i]),1]))^2))
  #     }
  #     dist_txt<-dist
  #     distparam_txt<-"\U03C3_t"
  #   }
  #   if (dist=="Lognormal" && dl=="ExponentialNorm") {
  #     # Remember theta[1] is logσ_t
  #     loglik <- function(theta){
  #       -sum(-theta[1] - 0.5*log(2*pi) - log(data[which(data[,3]==unitnames[i]),2]) - 0.5*(exp(theta[1])^-2)*((log(data[which(data[,3]==unitnames[i]),2]) - logdegradationlife(theta[2:3],data[which(data[,3]==unitnames[i]),1],data[which(data[,3]==unitnames[i]),4],data[which(data[,3]==unitnames[i]),5]))^2))
  #     }
  #     dist_txt<-dist
  #     distparam_txt<-"\U03C3_t"
  #   }
  #
  #
  #   MLEandvar <- MLE.var.covar.select(loglik,LSQest1[i,])
  #   MLE.theta.hat.perunit[[i]] <- MLEandvar[[1]]
  #   MLE.inv.fish.perunit[[i]]  <- MLEandvar[[2]]
  # }


  # 2. If correlation is selected (=1) then we do the full likelihood with parameter covariance
  # with normal distribution fit. Use mvn fit to initialize that part and numerical analysis for the rest.
  # No correlation (=0) means we compute this without correlation.  However we still need the uncorrelated
  # output as an initial LSQ for the MU and VARCOV Matrix.
  # No correlation Case
  # NOTE (2/7/2025): This is still related to a single stress condition so my answers will still have error.
  # Need to figure out the multi-stress correlation option now.
  if(is.null(modelstress) == TRUE){ # Single Stress evaluation
    if(modelstresstype == 0){
      if (dist=="Weibull") {
        loglik <- function(theta){
          -sum(theta[1] + (exp(theta[1])-1)*log(DegradationFULL) - exp(theta[1])*logdegradationlife(theta[2:3],TimeFULL) - ((DegradationFULL/degradationlife(theta[2:3],TimeFULL))^exp(theta[1])))
        }
      }
      if (dist=="Normal" && dl!="ExponentialNorm") {
        loglik <- function(theta){
          -sum(-theta[1] - 0.5*log(2*pi) - 0.5*(exp(theta[1])^-2)*((DegradationFULL - degradationlife(theta[2:3],TimeFULL))^2))
        }
        dist_txt<-dist
        distparam_txt<-"\U03C3"
      }
      if (dist=="Normal" && dl=="ExponentialNorm") {
        loglik <- function(theta){
          -sum(-theta[1] - 0.5*log(2*pi) - 0.5*(exp(theta[1])^-2)*((DegradationFULL - degradationlife(theta[2:3],TimeFULL,StresSFULL[,2],StresSFULL[,1]))^2))
        }
        dist_txt<-dist
        distparam_txt<-"\U03C3"
      }
      if (dist=="Lognormal" && dl!="ExponentialNorm") {
        loglik <- function(theta){
          -sum(-theta[1] - 0.5*log(2*pi) - log(DegradationFULL) - 0.5*(exp(theta[1])^-2)*((log(DegradationFULL) - logdegradationlife(theta[2:3],TimeFULL))^2))
        }
        dist_txt<-dist
        distparam_txt<-"\U03C3_t"
      }
      if (dist=="Lognormal" && dl=="ExponentialNorm") {
        loglik <- function(theta){
          -sum(-theta[1] - 0.5*log(2*pi) - log(DegradationFULL) - 0.5*(exp(theta[1])^-2)*((log(DegradationFULL) - logdegradationlife(theta[2:3],TimeFULL,StresSFULL[,2],StresSFULL[,1]))^2))
        }
        dist_txt<-dist
        distparam_txt<-"\U03C3_t"
      }
      MLEandvar <- MLE.var.covar.select(loglik,LSQest)
    }
  }
  # return(MLEandvar)
  if(is.null(modelstress) == FALSE){ # Multi stress evaluation
    if(modelstresstype == 1){
      if (dist=="Weibull") {
        loglik <- function(theta){
          -sum(theta[1] + (exp(theta[1])-1)*log(DegradationFULL) - exp(theta[1])*logdegradationlife1(theta[2:4],TimeFULL,StresSFULL) - ((DegradationFULL/degradationlife1(theta[2:3],TimeFULL,StresSFULL))^exp(theta[1])))
        }
      }
      if (dist=="Normal") {
        loglik <- function(theta){
          -sum(-theta[1] - 0.5*log(2*(pi^2)) + log(colSums(W_MAT*exp(-0.5*(exp(theta[1])^-2)*((Degradation_MAT - degradationlife1(theta,Time_MAT,Stress_MAT,Z_MAT))^2)))))
        }
        dist_txt<-dist
        distparam_txt<-"\U03C3"
      }
      # if (dist=="Normal" && dl!="ExponentialNorm") {
      #   loglik <- function(theta){
      #     # (RCS 02/11/2025) Will go with forming a sub function  but I may need to formulate these within the IF dist== call
      #     subfunctionNORM <- function(TIME,DEGRADATION,STRESS){
      #       exp(-0.5*(exp(theta[1])^-2)*(DEGRADATION - degradationlife1(c(1,1,1),TIME,STRESS)))*mvn_function(theta)
      #     }
      #     mapply(subfunctionNORM, TimeLIST,DegradationLIST,StresSLIST)
      #     # TimeLIST,DegradationLIST,StresSLIST
      #     # -N*theta[1] - 0.5*N*log(2*pi) + sum(1)
      #     # -sum(-theta[1] - 0.5*log(2*pi) - 0.5*(exp(theta[1])^-2)*((DegradationFULL - degradationlife1(theta,TimeFULL,StresSFULL))^2))
      #   }
      # }
      # if (dist=="Normal" && dl=="ExponentialNorm") {
      #   loglik <- function(theta){
      #     -sum(-theta[1] - 0.5*log(2*pi) - 0.5*(exp(theta[1])^-2)*((DegradationFULL - degradationlife1(theta,TimeFULL,StresSFULL))^2))
      #   }
      # }
      if (dist=="Lognormal") {
        loglik <- function(theta){
          -sum(-theta[1] - 0.5*log(2*(pi^2)) - log(DegradationFULL) + log(colSums(W_MAT*exp(-0.5*(exp(theta[1])^-2)*((log(Degradation_MAT) - logdegradationlife1(theta,Time_MAT,Stress_MAT,Z_MAT))^2)))))
        }
        dist_txt<-dist
        distparam_txt<-"\U03C3_t"
      }
      # return(LSQest2)
      # return(list(loglik,LSQest2))
      MLEandvar <- MLE.var.covar.select(loglik,LSQest2)
      # MLEandvar <- ucminf(LSQest2,loglik,hessian=1)
      # MLEandvar <- DEoptim(loglik,lower = c(-10,-10*abs(LSQest2[2:5])), upper = c(10,10*abs(LSQest2[2:5])))
    }
  }

  # return(MLEandvar)
  theta.hat <- MLEandvar[[1]]
  inv.fish  <- MLEandvar[[2]]
  loglik.hat <- -loglik(theta.hat)
  likeli.hat <- exp(loglik.hat)
  # theta.hat <- MLEandvar$par
  # inv.fish <- pinv(MLEandvar$hessian)

  # theta.hat <- MLEandvar$member$pop[1,]
  # inv.fish <- pinv(hessian(loglik,MLEandvar$member$pop[1,]))

  # return(list(theta.hat,inv.fish,loglik.hat,likeli.hat))

  # Compute the boundaries for output
  crit <- qnorm((1 + conf.level)/2)
  crit2 <- qnorm(conf.level)
  conflim<-vector(mode = "list", length = length(theta.hat))
  if(is.null(Suse)==TRUE){
    fulllimset<-vector(mode = "list", length = length(theta.hat))
  }
  if(is.null(Suse)==FALSE){
    fulllimset<-vector(mode = "list", length = (1+length(theta.hat)))
  }

  qlow <- 0.5*(1-confid)
  qhigh <- 1 - 0.5*(1-confid)

  if(is.null(modelstress) == FALSE){ # Get the confidence for each stress level (multi-stress operation)
    if(stresscount==1){ # Group unique stress levels (single)
      STRESSgroups <- unique(data[,4])
      L.mean <- rep(0,length(STRESSgroups))
      mean.degradation.STRESSgroups<-vector(mode = "list", length = length(STRESSgroups))
      mean.degradation.STRESSgroups.legend<-vector(mode = "list", length = length(STRESSgroups))
      conflim.STRESSgroups<-vector(mode = "list", length = length(STRESSgroups))
    }
    N.samples <- 100
    param_pull <- rmvnorm(N.samples,theta.hat,inv.fish)

    if(dl=="Power" || dl=="Logarithmic" || dl=="LloydLipow"){
      Time.fit <- linspace(0,max(data[,1]),1000)
    }
    if(dl=="Linear" || dl=="Exponential" || dl=="SquareRoot" || dl=="SquareRoot2"|| dl=="Mitsuom"){
      Time.fit <- linspace(0.01,max(data[,1]),1000)
    }
    for(i in 1:length(unique(stressvals))){# Need to find the end point life for the upper and lower confidence bounds of each stress level
      # if(dist=="Normal"){
      #
      # }
      L.mean[i] <- lifedegradation1(theta.hat,D0,STRESSgroups[i],Z_set[10]) # Pull mean life for each stress value
      L.samples <- rep(0,N.samples)                                         # initiate life samples
      for(j in 1:N.samples){
        # Estimate pseudo failures at D0 (these will be the stopping points for each curve plotted)
        L.samples[j] <- lifedegradation1(param_pull[j,],D0,STRESSgroups[i],Z_set[10])
      }
      conflim.out <- quantile(L.samples,c(qlow,qhigh))
      names(conflim.out) <- NULL
      conflim.STRESSgroups[[i]] <- conflim.out
      mean.degradation.STRESSgroups[[i]] <- degradationlife1(theta.hat,Time.fit,STRESSgroups[i],Z_set[10])
      if(is.null(stressunit1)==TRUE){
        mean.degradation.STRESSgroups.legend[[i]] <- rep(paste(c("Mean MLE for",STRESSgroups[i],"units"),collapse = " "),1000)
      }
      if(is.null(stressunit1)==FALSE){
        mean.degradation.STRESSgroups.legend[[i]] <- rep(paste(c("Mean MLE for",STRESSgroups[i],stressunit1),collapse = " "),1000)
      }

      # logdegradationlife1(theta.hat,Time.fit,STRESSgroups[i],Z_set[10])
    }
    if(is.null(Suse) == FALSE){ # Compute use life and use life bounds if Suse is given
      Luse <- lifedegradation1(theta.hat,D0,Suse,Z_set[10])
      L.samples <- rep(0,N.samples)
      for(j in 1:N.samples){
        # Estimate pseudo cycles-to-failure (these will be the stopping points for each curve plotted)
        L.samples[j] <- lifedegradation1(param_pull[j,],D0,Suse,Z_set[10])
      }
      conflim.out <- quantile(L.samples,c(qlow,qhigh))
      names(conflim.out) <- NULL
      conflim.use.stress <- conflim.out
      mean.degradation.use.stress <- degradationlife1(theta.hat,Time.fit,Suse,Z_set[10])
      if(is.null(stressunit1)==TRUE){
        mean.degradation.use.stress.legend <- rep(paste(c("Mean MLE for Use stress ",Suse,"units"),collapse = " "),1000)
      }
      if(is.null(stressunit1)==FALSE){
        mean.degradation.use.stress.legend <- rep(paste(c("Mean MLE for Use stress ",Suse,stressunit1),collapse = " "),1000)
      }
    }
  }
  # return(list(conflim.STRESSgroups,L.mean,mean.degradation.STRESSgroups.legend))

  # for(i in 1:length(data[,1])){  # Labeling the stresses for the plot legend
  #   if (dist=="Normal") {
  #     loglik <- function(theta){
  #       -sum(-theta[1] - 0.5*log(2*(pi^2)) + log(colSums(W_MAT*exp(-0.5*(exp(theta[1])^-2)*((Degradation_MAT - degradationlife1(theta,Time_MAT,Stress_MAT,Z_MAT))^2)))))
  #     }
  #     low_confid_curve <- function(theta){
  #       0.5*(1-confid) - (theta.hat[1]^-1)*(sqrt(2*pi^2)^-1)*sum(W_set*exp(-0.5*(exp(theta[1])^-2)*((theta[1] - degradationlife1(theta,theta[2],Stress_MAT,Z_set))^2)))
  #     }
  #   }
  # }


  for(i in 1:length(theta.hat)){
    if(sided == "twosided"){
      conflim[[i]]<-theta.hat[i] + c(-1, 1) * crit * sqrt(inv.fish[i, i])
      conflim_txt<-c(paste(c("Lower ",100*conf.level,"%"),collapse = ""),paste(c("Upper ",100*conf.level,"%"),collapse = ""))
    }
    if(sided == "onesidedlow"){
      conflim[[i]]<-theta.hat[i] - crit2 * sqrt(inv.fish[i, i])
      conflim_txt<-paste(c("One-Sided Low ",100*conf.level,"%"),collapse = "")
    }
    if(sided == "onesidedhigh"){
      conflim[[i]]<-theta.hat[i] + crit2 * sqrt(inv.fish[i, i])
      conflim_txt<-paste(c("One-Sided High ",100*conf.level,"%"),collapse = "")
    }

    # Computes back the original scale
    if((dist=="Lognormal" || dist=="Normal"  || dist=="Weibull") && i == 1){
      conflim[[i]] <- sort(exp(conflim[[i]]))
    }
    if(is.null(modelstress)==TRUE && (dl=="Exponential" || dl=="Power") && i == 3){
      conflim[[i]] <- sort(exp(conflim[[i]]))
    }
    if(is.null(modelstress)==FALSE && modelstress=="TempHumidity" && i == 2){
      conflim[[i]] <- sort(exp(conflim[[i]]))
    }
    if(is.null(modelstress)==FALSE && (modelstress=="Exponential" || modelstress=="Exponential2" || modelstress=="Arrhenius" || modelstress=="Eyring" || modelstress=="Power" || modelstress=="InversePower" || modelstress=="InversePower2") && i == 3){
      conflim[[i]] <- sort(exp(conflim[[i]]))
    }
    if (is.null(modelstress)==FALSE && modelstress=="TempNonthermal" && i == 4){
      conflim[[i]] <- sort(exp(conflim[[i]]))
    }
    if (is.null(modelstress)==FALSE && i == (1+ishift2+2)){
      conflim[[i]] <- sort(exp(conflim[[i]]))
    }
    fulllimset[[i]]<-c(theta.hat[i],conflim[[i]])
  }

  for(i in 1:length(theta.hat)){
    # Computes back the original scale
    if((dist=="Lognormal" || dist=="Normal") && i == 1){
      fulllimset[[i]] <- c(exp(theta.hat[i]),conflim[[i]])
    }
    if(is.null(modelstress)==TRUE && (dl=="Exponential" || dl=="Power") && i == 3){
      fulllimset[[i]] <- c(exp(theta.hat[i]),conflim[[i]])
    }
    if(is.null(modelstress)==FALSE && modelstress=="TempHumidity" && i == 2){
      fulllimset[[i]] <- c(exp(theta.hat[i]),conflim[[i]])
    }
    if(is.null(modelstress)==FALSE && (modelstress=="Exponential" || modelstress=="Exponential2" || modelstress=="Arrhenius" || modelstress=="Eyring" || modelstress=="Power" || modelstress=="InversePower" || modelstress=="InversePower2") && i == 3){
      fulllimset[[i]] <- c(exp(theta.hat[i]),conflim[[i]])
    }
    if (is.null(modelstress)==FALSE && modelstress=="TempNonthermal" && i == 4){
      fulllimset[[i]] <- c(exp(theta.hat[i]),conflim[[i]])
    }
    if (is.null(modelstress)==FALSE && i == (1+ishift2+2)){
      fulllimset[[i]] <- c(exp(theta.hat[i]),conflim[[i]])
    }
  }

  if(dl=="SquareRoot"){
    if (mean(adtLSQ[,2]) > 0 && min(adtLSQ[,1]) >= 0){ # Positive slope and positive intercept
      fulllimset[[2]] <- c(exp(fulllimset[[2]][1]),sort(exp(fulllimset[[2]][2:3])))
      fulllimset[[3]] <- c(exp(fulllimset[[3]][1]),sort(exp(fulllimset[[3]][2:3])))
    }
    if (mean(adtLSQ[,2]) > 0  && min(adtLSQ[,1]) < 0){ # Positive slope and positive/negative intercept
      fulllimset[[3]] <- c(exp(fulllimset[[3]][1]),sort(exp(fulllimset[[3]][2:3])))
    }
  }
  if(is.null(Suse)==TRUE){
    params_txt<-c(distparam_txt,params_txt)
  }
  if(is.null(Suse)==FALSE){
    fulllimset[[(1+length(theta.hat))]] <- c(Luse,conflim.use.stress)
    params_txt<-c(distparam_txt,params_txt,"Use Life")
  }

  AIC = 2*length(theta.hat) + 2*loglik(theta.hat)
  BIC = 2*log(length(DegradationFULL)) + 2*loglik(theta.hat)

  # Recompute necessary output
  if(dist=="Lognormal" || dist=="Normal"){
    theta.hat[1] <- exp(theta.hat[1])
  }
  if(is.null(modelstress)==TRUE && (dl=="Exponential" || dl=="Power")){
    theta.hat[3] <- exp(theta.hat[3])
  }
  if(is.null(modelstress)==FALSE && modelstress=="TempHumidity"){
    theta.hat[2] <- exp(theta.hat[2])
  }
  if(is.null(modelstress)==FALSE && (modelstress=="Exponential" || modelstress=="Exponential2" || modelstress=="Arrhenius" || modelstress=="Eyring" || modelstress=="Power" || modelstress=="InversePower" || modelstress=="InversePower2")){
    theta.hat[3] <- exp(theta.hat[3])
  }
  if (is.null(modelstress)==FALSE && modelstress=="TempNonthermal"){
    theta.hat[4] <- exp(theta.hat[4])
  }
  if (is.null(modelstress)==FALSE){
    theta.hat[1+ishift2+2] <- exp(theta.hat[1+ishift2+2])
  }
  # return(list(theta.hat,inv.fish,unlist(fulllimset),params_txt))
  # Plot the MLE Fit to the data as well as the life distribution if available
  # Set up the three cases
  if(is.null(modelstress) == TRUE){ # Single stress condition so we group by the individual units
    unitnames <- unique(data[,3]) # Pull the unit designations from column 3 of the input data
    df1 <- data.frame(timedat = data[,1], degradationdat = data[,2], group = data[,3]) # Data frame raw data

    # Plot the data
    plotout<-ggplot() +
      geom_point(data=df1, aes(timedat,degradationdat, shape = group), colour = 'black', size = 2.2) +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_shape_manual("Raw Data",values=shape_legend[1:length(unitnames)]) +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0)) +
      xlab(xlabel) +
      ylab(ylabel)
  }
  if(is.null(modelstress) == FALSE && is.null(Suse) == TRUE){ # Multi stress condition so we group by the stress (No Suse given)
    stresslabel <- rep(0,length(data[,1]))
    for(i in 1:length(data[,1])){  # Labeling the stresses for the plot legend
      stresslabel[i] <- paste(c("Data for",data[i,4],"units"),collapse = " ")
      stresslabel[i] <- paste(c("Data for",data[i,4],stressunit1),collapse = " ")
    }
    df1 <- data.frame(timedat = data[,1], degradationdat = data[,2], group = stresslabel) # Data frame raw data
    df2 <- data.frame(timedat = rep(Time.fit,length(unique(stressvals))), degradationdat = unlist(mean.degradation.STRESSgroups), group = unlist(mean.degradation.STRESSgroups.legend))

    # Plot the data
    plotout<-ggplot() +
      geom_point(data=df1, aes(timedat,degradationdat, shape = group), colour = 'black', size = 2.2) +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_shape_manual("Raw Data",values=shape_legend[1:length(unitnames)]) +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0)) +
      xlab(xlabel) +
      ylab(ylabel)

    # Plot best fit line
    plotout <- plotout + geom_line(data=df2, aes(timedat,degradationdat, colour = group), linewidth = 0.9, linetype = "dashed") +
      scale_color_manual("MLE Fit Lines",values=col_legend[1:length(unique(stressvals))])

    # return(list(MLE.point.estimate = theta.hat,var.cov.matrix = inv.fish,loglik = loglik.hat,likelihood = likeli.hat,AIC = AIC,BIC = BIC,length(unique(stresslabel)),plotoutput = plotout)) # Plot output (TEST)
  }

  if(is.null(modelstress) == FALSE && is.null(Suse) == FALSE){ # Multi stress condition so we group by the stress (Suse given)
    stresslabel <- rep(0,length(data[,1]))
    for(i in 1:length(data[,1])){  # Labeling the stresses for the plot legend
      stresslabel[i] <- paste(c("Data for",data[i,4],"units"),collapse = " ")
      stresslabel[i] <- paste(c("Data for",data[i,4],stressunit1),collapse = " ")
    }
    df1 <- data.frame(timedat = data[,1], degradationdat = data[,2], group = stresslabel) # Data frame raw data
    df2 <- data.frame(timedat = rep(Time.fit,(length(unique(stressvals))+1)), degradationdat = c(unlist(mean.degradation.STRESSgroups),mean.degradation.use.stress), group = c(unlist(mean.degradation.STRESSgroups.legend),mean.degradation.use.stress.legend))

    # Plot the data
    plotout<-ggplot() +
      geom_point(data=df1, aes(timedat,degradationdat, shape = group), colour = 'black', size = 2.2) +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_shape_manual("Raw Data",values=shape_legend[1:length(unitnames)]) +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0)) +
      xlab(xlabel) +
      ylab(ylabel)

    # Plot best fit line
    plotout <- plotout + geom_line(data=df2, aes(timedat,degradationdat, colour = group), linewidth = 0.9, linetype = "dashed") +
      scale_color_manual("MLE Fit Lines",values=col_legend[1:(length(unique(stressvals))+1)])
  }
  # return(dlmodel_txt0)



  # Produce some output text that summarizes the results
  if(is.null(modelstress) == TRUE){ # Single stress condition so we group by the individual units
    cat(c("Maximum-Likelihood estimates for the ",dl_txt,"-",dist_txt," Degradation model.\n\nD(t) = ",Dt_txt,"\n\n"),sep = "")
    print(matrix(unlist(fulllimset), nrow = length(unlist(fulllimset))/length(params_txt), ncol = length(params_txt), byrow = FALSE,dimnames = list(c("Degradation-Life Parameters Mean",conflim_txt),params_txt)))
    cat("\n")

    return(list(MLE.point.estimate = theta.hat,var.cov.matrix = inv.fish,
                loglik = loglik.hat,likelihood = likeli.hat,AIC = AIC,BIC = BIC,plotoutput = plotout)) # Plot output (TEST)
  }
  if(is.null(modelstress) == FALSE && is.null(Suse) == TRUE){ # Multi stress condition so we group by the stress (No Suse given)
    cat(c("Maximum-Likelihood estimates for the ",dl_txt,"-",ps_txt," Degradation-Stress-Life model.\n\nD(l,S) = ",dlmodel_txt0," where ",paramref_txt,param_txt2," and ",paramref_txt2,"\n\n"),sep = "")
    print(matrix(unlist(fulllimset), nrow = length(unlist(fulllimset))/length(params_txt), ncol = length(params_txt), byrow = FALSE,dimnames = list(c("Degradation-Life Parameters Mean",conflim_txt),params_txt)))
    cat("\n")
    # print(matrix(c(params_0[[1]],params_1,params_0[[2]]), nrow = length(params_0[[1]])+3, ncol = 1,byrow = FALSE,dimnames = list(c(params_txt,params_txt2,"param-stress SSE"),"EsT")))
    # cat("\n")

    return(list(MLE.point.estimate = theta.hat,var.cov.matrix = inv.fish,
                stress.level.mean.MLE.life = L.mean,
                stress.level.MLE.life.CI = conflim.STRESSgroups,
                loglik = loglik.hat,likelihood = likeli.hat,AIC = AIC,BIC = BIC,plotoutput = plotout)) # Plot output (TEST)
  }
  if(is.null(modelstress) == FALSE && is.null(Suse) == FALSE){ # Multi stress condition so we group by the stress (Suse given)
    cat(c("Maximum-Likelihood estimates for the ",dl_txt,"-",ps_txt," Degradation-Stress-Life model.\n\nD(l,S) = ",dlmodel_txt0," where ",paramref_txt,param_txt2," and ",paramref_txt2,"\n\n"),sep = "")
    print(matrix(unlist(fulllimset), nrow = length(unlist(fulllimset))/length(params_txt), ncol = length(params_txt), byrow = FALSE,dimnames = list(c("Degradation-Life Parameters Mean",conflim_txt),params_txt)))
    cat("\n")
    # print(matrix(c(params_0[[1]],params_1,params_0[[2]]), nrow = length(params_0[[1]])+3, ncol = 1,byrow = FALSE,dimnames = list(c(params_txt,params_txt2,"param-stress SSE"),"EsT")))
    # cat("\n")

    return(list(MLE.point.estimate = theta.hat,var.cov.matrix = inv.fish,
                stress.level.mean.MLE.life = L.mean,use.level.mean.MLE.life = Luse,
                stress.level.MLE.life.CI = conflim.STRESSgroups,use.level.MLE.life.CI = conflim.use.stress,
                loglik = loglik.hat,likelihood = likeli.hat,AIC = AIC,BIC = BIC,plotoutput = plotout)) # Plot output (TEST)
  }
  # cat(c("Maximum-Likelihood estimates for the ",dl_txt,"-",dist_txt," Degradation-Life-Stress model.\n\n"),sep = "")
  cat(c("Maximum-Likelihood estimates for the ",dl_txt,"-",dist_txt," Degradation model.\n\nD(t) = ",Dt_txt,"\n\n"),sep = "")
  # print(matrix(unlist(fulllimset), nrow = length(unlist(fulllimset))/length(LSQest), ncol = length(LSQest), byrow = FALSE,dimnames = list(c("Life-Stress Parameters Mean",conflim_txt),params_txt)))

  if(is.null(modelstress) == TRUE){ # Single stress condition so we group by the individual units
  }
  # cat(c("Maximum-Likelihood estimates for the ",dl_txt,"-",dist_txt," Degradation model.\n\nD(t) = ",Dt_txt,"\n\n"),sep = "")
  print(matrix(unlist(fulllimset), nrow = length(unlist(fulllimset))/length(params_txt), ncol = length(params_txt), byrow = FALSE,dimnames = list(c("Degradation-Life Parameters Mean",conflim_txt),params_txt)))
  cat(c("\n"))

  # return(list(theta.hat,inv.fish,fulllimset))
  return(list(MLE.point.estimate = theta.hat,var.cov.matrix = inv.fish,loglik = loglik.hat,likelihood = likeli.hat,AIC = AIC,BIC = BIC,plotoutput = plotout)) # Plot output (TEST)

  # return(list(theta.hat,inv.fish,conflim,loglik.hat,likeli.hat,AIC,BIC))

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
