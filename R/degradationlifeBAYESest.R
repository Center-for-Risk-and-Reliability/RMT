# Bayesian Degradation Life Estimator
# Developed by Dr. Reuel Smith, 2025

degradationlife.BAYESest <- function(pt_est,data,dl,dist="Normal",D0,modelstress=NULL,confid=0.95,SUSE=NULL,priors,nsamples=20000,burnin=10000,nchains=4,Q=20){
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

  # Check to see if dist="Exponential" so you can exclude life
  # distribution parameters.
  if (dist=="Exponential") {
    ishift<-0
  } else {
    ishift<-1
  }
  # Check to see if confidence exists
  conf.level <- confid
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

  # Pull the time and degradation data from the input
  TimeFULL <- data[,1]
  DegradationFULL <- data[,2]
  StresSFULL <- data[,4:dim(data)[2]]
  N <- length(data[,1])

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
    Ones_MAT <- matrix(data = rep(1,Q*N), nrow = Q, ncol = N,byrow = FALSE)
  }
  # ===================================================================
  # Establish the model-stress relations for parameters (usually based on the slope of linearized degradation life model)
  # ===================================================================
  # Initialize parameter-stress (modelstress) estimates for theta
  if (is.null(modelstress) == FALSE && modelstress=="Linear") {
    # theta[1] - parameter a_0, theta[2] - parameter b_0
    ishift2<-2                            # Second ishift for number of parameters
    psparams <- "real a_0; real b_0;"
    psparamsvec <- c("a_0","b_0")
    pr1<-paste(c("a_0 ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b_0 ~ ",priors[ishift+2],";"),collapse = "")
    pspriors <- paste(c(pr1,pr2),collapse = " ")

    paramstress <- "b_0 + Sf*a_0"
    logparamstress <- "log(b_0 + Sf*a_0)"
    paramstress.USE <- "b_0 + Suse*a_0"
    logparamstress.USE <- "log(b_0 + Suse*a_0)"
  }

  if (is.null(modelstress) == FALSE && modelstress=="Exponential"){
    # theta[1] - parameter a_0, theta[2] - parameter b_0
    ishift2<-2                            # Second ishift for number of parameters
    psparams <- "real a_0; real<lower=0> b_0;"
    psparamsvec <- c("a_0","b_0")
    pr1<-paste(c("a_0 ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b_0 ~ ",priors[ishift+2],";"),collapse = "")
    pspriors <- paste(c(pr1,pr2),collapse = " ")

    paramstress <- "b_0*exp(a_0*Sf)"
    logparamstress <- "log(b_0) + a_0*Sf"
    paramstress.USE <- "b_0*exp(a_0*Suse)"
    logparamstress.USE <- "log(b_0) + a_0*Suse"
  }

  if (is.null(modelstress) == FALSE && modelstress=="Exponential2"){
    # theta[1] - parameter a_0, theta[2] - parameter b_0
    ishift2<-2                            # Second ishift for number of parameters
    psparams <- "real a_0; real<lower=0> b_0;"
    psparamsvec <- c("a_0","b_0")
    pr1<-paste(c("a_0 ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b_0 ~ ",priors[ishift+2],";"),collapse = "")
    pspriors <- paste(c(pr1,pr2),collapse = " ")

    paramstress <- "b_0*exp(a_0/Sf)"
    logparamstress <- "log(b_0) + a_0/Sf"
    paramstress.USE <- "b_0*exp(a_0/Suse)"
    logparamstress.USE <- "log(b_0) + a_0/Suse"
  }

  if (is.null(modelstress) == FALSE && modelstress=="Arrhenius") {
    # lsparams[1] - parameter Ea_0, lsparams[2] - parameter b_0
    # Temperature HAS to be in Kelvin for this to work
    ishift2<-2                            # Second ishift for number of parameters
    psparams <- "real E_a_0; real<lower=0> b_0;"
    psparamsvec <- c("E_a_0","b_0")
    pr1<-paste(c("E_a_0 ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b_0 ~ ",priors[ishift+2],";"),collapse = "")
    pspriors <- paste(c(pr1,pr2),collapse = " ")

    paramstress <- "b_0*exp(E_a_0/((8.617385e-5)*Sf))"
    logparamstress <- "log(b_0) + (E_a_0/((8.617385e-5)*Sf))"
    paramstress.USE <- "b_0*exp(E_a_0/((8.617385e-5)*Suse))"
    logparamstress.USE <- "log(b_0) + (E_a_0/((8.617385e-5)*Suse))"
  }

  if (is.null(modelstress) == FALSE && modelstress=="Eyring") {
    # lsparams[1] - parameter a_0, lsparams[2] - parameter b_0
    ishift2<-2                            # Second ishift for number of parameters
    psparams <- "real a_0; real<lower=0> b_0;"
    psparamsvec <- c("a_0","b_0")
    pr1<-paste(c("a_0 ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b_0 ~ ",priors[ishift+2],";"),collapse = "")
    pspriors <- paste(c(pr1,pr2),collapse = " ")

    # NOTE (8/25/2025): Can now go back to all Eyring type functions and just set .*, ./, and .^ on vector
    # multipliers (like MATLAB)
    paramstress <- "(b_0/Sf).*exp(a_0/Sf)"
    logparamstress <- "log(b_0) - log(Sf) + (a_0/Sf)"
    paramstress.USE <- "(b_0/Suse).*exp(a_0/Suse)"
    logparamstress.USE <- "log(b_0) - log(Suse) + (a_0/Suse)"
  }

  if (is.null(modelstress) == FALSE && modelstress=="Eyring2") {
    # lsparams[1] - parameter a_0, lsparams[2] - parameter b_0
    ishift2<-2                            # Second ishift for number of parameters
    psparams <- "real a_0; real b_0;"
    psparamsvec <- c("a_0","b_0")
    pr1<-paste(c("a_0 ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b_0 ~ ",priors[ishift+2],";"),collapse = "")
    pspriors <- paste(c(pr1,pr2),collapse = " ")

    paramstress <- "(1/Sf).*exp(-(a_0 - (b_0/Sf)))"
    logparamstress <- "-log(Sf) - a_0 + (b_0/Sf)"
    paramstress.USE <- "(1/Suse).*exp(-(a_0 - (b_0/Suse)))"
    logparamstress.USE <- "-log(Suse) - a_0 + (b_0/Suse)"
  }

  if (is.null(modelstress) == FALSE && modelstress=="Power") {
    # psparams[1] - parameter a_0, psparams[2] - parameter b_0
    ishift2<-2                            # Second ishift for number of parameters
    psparams <- "real a_0; real<lower=0> b_0;"
    psparamsvec <- c("a_0","b_0")
    pr1<-paste(c("a_0 ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b_0 ~ ",priors[ishift+2],";"),collapse = "")
    pspriors <- paste(c(pr1,pr2),collapse = " ")

    paramstress <- "b_0.*(Sf.^a_0)"
    logparamstress <- "log(b_0) + a_0*log(Sf)"
    paramstress.USE <- "b_0.*(Suse.^a_0)"
    logparamstress.USE <- "log(b_0) + a_0*log(Suse)"
  }


  if (is.null(modelstress) == FALSE && modelstress=="InversePower") {
    # lsparams[1] - parameter a_0, lsparams[2] - parameter b_0
    ishift2<-2                            # Second ishift for number of parameters
    psparams <- "real a_0; real<lower=0> b_0;"
    psparamsvec <- c("a_0","b_0")
    pr1<-paste(c("a_0 ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b_0 ~ ",priors[ishift+2],";"),collapse = "")
    pspriors <- paste(c(pr1,pr2),collapse = " ")

    paramstress <- "b_0.*(Sf.^-a_0)"
    logparamstress <- "log(b_0) - a_0*log(Sf)"
    paramstress.USE <- "b_0.*(Suse^-a_0)"
    logparamstress.USE <- "log(b_0) - a_0*log(Suse)"
  }

  if (is.null(modelstress) == FALSE && modelstress=="InversePower2") {
    # lsparams[1] - parameter a_0, lsparams[2] - parameter b_0
    ishift2<-2                            # Second ishift for number of parameters
    psparams <- "real a_0; real<lower=0> b_0;"
    psparamsvec <- c("a_0","b_0")
    pr1<-paste(c("a_0 ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b_0 ~ ",priors[ishift+2],";"),collapse = "")
    pspriors <- paste(c(pr1,pr2),collapse = " ")

    paramstress <- "1./(b_0.*(Sf.^a_0))"
    logparamstress <- "-log(b_0) - a_0*log(Sf)"
    paramstress.USE <- "1./(b_0.*(Suse^a_0))"
    logparamstress.USE <- "-log(b_0) - a_0*log(Suse)"
  }

  if (is.null(modelstress) == FALSE && modelstress=="Logarithmic") {
    # lsparams[1] - parameter a_0, lsparams[2] - parameter b_0
    ishift2<-2                            # Second ishift for number of parameters
    psparams <- "real a_0; real b_0;"
    psparamsvec <- c("a_0","b_0")
    pr1<-paste(c("a_0 ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b_0 ~ ",priors[ishift+2],";"),collapse = "")
    pspriors <- paste(c(pr1,pr2),collapse = " ")

    paramstress <- "a_0*log(Sf) + b_0"
    logparamstress <- "log(a_0*log(Sf) + b_0)"
    paramstress.USE <- "a_0*log(Suse) + b_0"
    logparamstress.USE <- "log(a_0*log(Suse) + b_0)"
  }

  if (is.null(modelstress) == FALSE && modelstress=="MultiStress") {
    # CHECK THIS LAST
    # lsparams - c(ao,a1,...,an), S - c(1,S1,...,Sn)
    if((length(priors)-ishift)==2){
      ishift2<-2                            # Second ishift for number of parameters
      psparams <- "real a0_0; real a1_0; "
      psparamsvec <- c("a0_0","a1_0")
      pr1<-paste(c("a0_0 ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("a1_0 ~ ",priors[ishift+2],";"),collapse = "")
      pspriors <- paste(c(pr1,pr2),collapse = " ")

      paramstress <- "exp(a0_0 + a1_0.*Sf)"
      logparamstress <- "a0_0 + a1_0.*Sf"
      paramstress.USE <- "exp(a0_0 + a1_0.*Suse)"
      logparamstress.USE <- "a0_0 + a1_0.*Suse"
    }
    if((length(priors)-ishift)==3){
      ishift2<-3                            # Second ishift for number of parameters
      psparams <- "real a0_0; real a1_0; real a2_0; "
      psparamsvec <- c("a0_0","a1_0","a2_0")
      pr1<-paste(c("a0_0 ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("a1_0 ~ ",priors[ishift+2],";"),collapse = "")
      pr3<-paste(c("a2_0 ~ ",priors[ishift+3],";"),collapse = "")
      pspriors <- paste(c(pr1,pr2,pr3),collapse = " ")

      paramstress <- "exp(a0_0 + a1_0.*Sf[,1] + a2_0.*Sf[,2]);"
      logparamstress <- "a0_0 + a1_0.*Sf[,1] + a2_0.*Sf[,2];"
      paramstress.USE <- "exp(a0_0 + a1_0.*Suse[,1] + a2_0.*Suse[,2])"
      logparamstress.USE <- "a0_0 + a1_0.*Suse[,1] + a2_0.*Suse[,2]"
    }
    if((length(priors)-ishift)==4){
      ishift2<-4                            # Second ishift for number of parameters
      psparams <- "real a0_0; real a1_0; real a2_0; real a3_0;"
      psparamsvec <- c("a0_0","a1_0","a2_0","a3_0")
      pr1<-paste(c("a0_0 ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("a1_0 ~ ",priors[ishift+2],";"),collapse = "")
      pr3<-paste(c("a2_0 ~ ",priors[ishift+3],";"),collapse = "")
      pr4<-paste(c("a3_0 ~ ",priors[ishift+4],";"),collapse = "")
      pspriors <- paste(c(pr1,pr2,pr3,pr4),collapse = " ")

      paramstress <- "exp(a0_0 + a1_0.*Sf[,1] + a2_0.*Sf[,2] + a3_0.*Sf[,3])"
      logparamstress <- "a0_0 + a1_0.*Sf[,1] + a2_0.*Sf[,2] + a3_0.*Sf[,3]"
      paramstress.USE <- "exp(a0_0 + a1_0.*Suse[,1] + a2_0.*Suse[,2] + a3_0.*Suse[,3])"
      logparamstress.USE <- "a0_0 + a1_0.*Suse[,1] + a2_0.*Suse[,2] + a3_0.*Suse[,3]"
    }
    if((length(priors)-ishift)==5){
      ishift2<-5                            # Second ishift for number of parameters
      psparams <- "real a0_0; real a1_0; real a2_0; real a3_0; real a4_0;"
      psparamsvec <- c("a0_0","a1_0","a2_0","a3_0","a4_0")
      pr1<-paste(c("a0_0 ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("a1_0 ~ ",priors[ishift+2],";"),collapse = "")
      pr3<-paste(c("a2_0 ~ ",priors[ishift+3],";"),collapse = "")
      pr4<-paste(c("a3_0 ~ ",priors[ishift+4],";"),collapse = "")
      pr5<-paste(c("a4_0 ~ ",priors[ishift+5],";"),collapse = "")
      pspriors <- paste(c(pr1,pr2,pr3,pr4,pr5),collapse = " ")

      paramstress <- "exp(a0_0 + a1_0*Sf[,1] + a2_0*Sf[,2] + a3_0*Sf[,3] + a4_0*Sf[,4])"
      logparamstress <- "a0_0 + a1_0*Sf[,1] + a2_0*Sf[,2] + a3_0*Sf[,3] + a4_0*Sf[,4]"
      paramstress.USE <- "exp(a0_0 + a1_0.*Suse[,1] + a2_0.*Suse[,2] + a3_0.*Suse[,3] + a4_0.*Suse[,4])"
      logparamstress.USE <- "a0_0 + a1_0.*Suse[,1] + a2_0.*Suse[,2] + a3_0.*Suse[,3] + a4_0.*Suse[,4]"
    }
  }

  if (is.null(modelstress) == FALSE && modelstress=="TempHumidity") {
    # lsparams[1] - parameter A_0, lsparams[2] - parameter a_0, lsparams[3] - parameter b_0
    ishift2<-3                            # Second ishift for number of parameters
    psparams <- "real<lower=0> A_0; real a_0; real b_0;"
    psparamsvec <- c("A_0","a_0","b_0")
    pr1<-paste(c("A_0 ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("a_0 ~ ",priors[ishift+2],";"),collapse = "")
    pr3<-paste(c("b_0 ~ ",priors[ishift+3],";"),collapse = "")
    pspriors <- paste(c(pr1,pr2,pr3),collapse = " ")

    paramstress <- "A_0.*exp((a_0/Sf[,1]) + (b_0/Sf[,2]));"
    logparamstress <- "log(A_0) + (a_0/Sf[,1]) + (b_0/Sf[,2]);"
    paramstress.USE <- "A_0.*exp((a_0/Suse[,1]) + (b_0/Suse[,2]));"
    logparamstress.USE <- "log(A_0) + (a_0/Suse[,1]) + (b_0/Suse[,2]);"
  }

  if (is.null(modelstress) == FALSE && modelstress=="TempNonthermal") {
    # lsparams[1] - parameter a_0, lsparams[2] - parameter b_0, lsparams[3] - parameter c_0
    ishift2<-3                            # Second ishift for number of parameters
    psparams <- "real a_0; real b_0; real<lower=0> c_0;"
    psparamsvec <- c("a_0","b_0","c_0")
    pr1<-paste(c("a_0 ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b_0 ~ ",priors[ishift+2],";"),collapse = "")
    pr3<-paste(c("c_0 ~ ",priors[ishift+3],";"),collapse = "")
    pspriors <- paste(c(pr1,pr2,pr3),collapse = " ")

    paramstress <- "c_0./((Sf[,2].^b_0)*exp(-a./Sf[,1]))"
    logparamstress <- "log(c_0) - b_0.*log(Sf[,2]) + (a_0/Sf[,1])"
    paramstress.USE <- "c_0./((Suse[,2].^b_0)*exp(-a./Suse[,1]))"
    logparamstress.USE <- "log(c_0) - b_0.*log(Suse[,2]) + (a_0/Suse[,1])"
  }

  if (is.null(modelstress) == FALSE && modelstress=="Eyring3") {
    # psparams[1] - parameter a_0, psparams[2] - parameter b_0
    # psparams[3] - parameter c_0, psparams[4] - parameter d_0
    ishift2<-4                            # Second ishift for number of parameters
    psparams <- "real a_0; real b_0; real c_0; real d_0;"
    psparamsvec <- c("a_0","b_0","c_0","d_0")
    pr1<-paste(c("a_0 ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b_0 ~ ",priors[ishift+2],";"),collapse = "")
    pr3<-paste(c("c_0 ~ ",priors[ishift+3],";"),collapse = "")
    pr4<-paste(c("d_0 ~ ",priors[ishift+4],";"),collapse = "")
    pspriors <- paste(c(pr1,pr2,pr3,pr4),collapse = " ")

    paramstress <- "(1/Sf[,1]).*exp((a_0 + (b_0/Sf[,1])) + (c_0 + (d_0/Sf[,1])).*Sf[,2])"
    logparamstress <- "-log(Sf[,1]) + a_0 + (b_0/Sf[,1]) + (c_0 + (d_0/Sf[,1])).*Sf[,2]"
    paramstress.USE <- "(1/Suse[,1]).*exp((a_0 + (b_0/Suse[,1])) + (c_0 + (d_0/Suse[,1])).*Suse[,2])"
    logparamstress.USE <- "-log(Suse[,1]) + a_0 + (b_0/Suse[,1]) + (c_0 + (d_0/Suse[,1])).*Suse[,2]"
  }

  # ===================================================================
  # Degradation-Life (degradationlife) models
  # ===================================================================
  if(dl=="Linear"){
    if(is.null(modelstress) == TRUE){       # Degradation model text for single stress level scenario (parameters a and b)
      # D = a + b*t
      # theta[1] ~ a, theta[2] ~ b
      dlparams <- "real a; real b;"
      dlparamsvec <- c("a","b")
      pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
      dlpriors <- paste(c(pr1,pr2),collapse = " ")

      degradationlife <- "a + b*Life"
      logdegradationlife <- "log(a + b*Life)"
    } else{                                 # Degradation model text for multiple stress level scenario (parameters mu_b and sigma_b)
      dlparams <-  "real mu_b; real<lower=0> sigma_b;"
      dlparamsvec <- c("mu_b","sigma_b")
      outputdlparamset <- c("\U03BC_b","\U03C3_b")
      pr1_<-paste(c("mu_b ~ ",priors[ishift+ishift2+1],";"),collapse = "")
      pr2_<-paste(c("sigma_b ~ ",priors[ishift+ishift2+2],";"),collapse = "")
      dlpriors <- paste(c(pr1_,pr2_),collapse = " ")

      degradationlife <- paste(c("(",paramstress,") + (mu_b + (2^0.5)*sigma_b*Z)*Life"),collapse = "")
      logdegradationlife <- paste(c("log((",paramstress,") + (mu_b + (2^0.5)*sigma_b*Z)*Life)"),collapse = "")

      if(missing(SUSE)==FALSE){ # Include use life if SUSE is given
        lifeU <- paste(c("(D0./(",paramstress.USE,")).^(1./(mu_a - sigma_a.*0.2453407.*1.414214));"),collapse = "")
      }
    }
  }

  if(dl=="Exponential"){
    # D = b*exp(a*t)
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"b*exp(a*t)"
    if(is.null(modelstress) == TRUE){       # Degradation model text for single stress level scenario (parameters a and b)
      # D = b*(t^a)
      # theta[1] ~ a, theta[2] ~ b
      dlparams <- "real a; real<lower=0> b;"
      dlparamsvec <- c("a","b")
      pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
      dlpriors <- paste(c(pr1,pr2),collapse = " ")

      degradationlife <- "b*exp(Life.*a)"
      logdegradationlife <- "log(b) + a.*Life"
    } else{                                 # Degradation model text for multiple stress level scenario (parameters mu_a and sigma_a)
      dlparams <-  "real mu_a; real<lower=0> sigma_a;"
      dlparamsvec <- c("mu_a","sigma_a")
      outputdlparamset <- c("\U03BC_a","\U03C3_a")
      pr1_<-paste(c("mu_a ~ ",priors[ishift+ishift2+1],";"),collapse = "")
      pr2_<-paste(c("sigma_a ~ ",priors[ishift+ishift2+2],";"),collapse = "")
      dlpriors <- paste(c(pr1_,pr2_),collapse = " ")

      degradationlife <- paste(c("(",paramstress,").*exp(Life.*(mu_a + (2^0.5)*sigma_a*Z))"),collapse = "")
      logdegradationlife <- paste(c(logparamstress," + (mu_a + (2^0.5)*sigma_a*Z).*Life"),collapse = "")

      if(missing(SUSE)==FALSE){ # Include use life if SUSE is given
        lifeU <- paste(c("(D0./(",paramstress.USE,")).^(1./(mu_a - sigma_a.*0.2453407.*1.414214));"),collapse = "")
      }
    }
  }

  if(dl=="SquareRoot"){
    if(is.null(modelstress) == TRUE){       # Degradation model text for single stress level scenario (parameters a and b)
      # D^(1/2) = a + b*t
      # theta[1] ~ a, theta[2] ~ b
      dlparams <- "real a; real b;"
      dlparamsvec <- c("a","b")
      pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
      dlpriors <- paste(c(pr1,pr2),collapse = " ")

      degradationlife <- "(a + b*Life)^2"
      logdegradationlife <- "2*log(a + b*Life)"
    } else{                                 # Degradation model text for multiple stress level scenario (parameters mu_b and sigma_b)
      dlparams <-  "real mu_b; real<lower=0> sigma_b;"
      dlparamsvec <- c("mu_b","sigma_b")
      outputdlparamset <- c("\U03BC_b","\U03C3_b")
      pr1_<-paste(c("mu_b ~ ",priors[ishift+ishift2+1],";"),collapse = "")
      pr2_<-paste(c("sigma_b ~ ",priors[ishift+ishift2+2],";"),collapse = "")
      dlpriors <- paste(c(pr1_,pr2_),collapse = " ")

      degradationlife <- paste(c("((",paramstress,") + (mu_b + (2^0.5)*sigma_b*Z).*Life).^2"),collapse = "")
      logdegradationlife <- paste(c("2*log((",paramstress,") + (mu_b + (2^0.5)*sigma_b*Z).*Life)"),collapse = "")

      if(missing(SUSE)==FALSE){ # Include use life if SUSE is given
        lifeU <- paste(c("(D0./(",paramstress.USE,")).^(1./(mu_a - sigma_a.*0.2453407.*1.414214));"),collapse = "")
      }
    }
  }

  if(dl=="Power"){
    if(is.null(modelstress) == TRUE){       # Degradation model text for single stress level scenario (parameters a and b)
      # D = b*(t^a)
      # theta[1] ~ a, theta[2] ~ b
      dlparams <- "real a; real<lower=0> b;"
      dlparamsvec <- c("a","b")
      pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
      dlpriors <- paste(c(pr1,pr2),collapse = " ")

      degradationlife <- "b*(Life.^a)"
      logdegradationlife <- "log(b) + a*log(Life.)"
    } else{                                 # Degradation model text for multiple stress level scenario (parameters mu_a and sigma_a)
      dlparams <-  "real mu_a; real<lower=0> sigma_a;"
      dlparamsvec <- c("mu_a","sigma_a")
      outputdlparamset <- c("\U03BC_a","\U03C3_a")
      pr1_<-paste(c("mu_a ~ ",priors[ishift+ishift2+1],";"),collapse = "")
      pr2_<-paste(c("sigma_a ~ ",priors[ishift+ishift2+2],";"),collapse = "")
      dlpriors <- paste(c(pr1_,pr2_),collapse = " ")

      degradationlife <- paste(c("(",paramstress,").*(Life.^(mu_a + (2^0.5)*sigma_a*Z))"),collapse = "")
      logdegradationlife <- paste(c(logparamstress," + (mu_a + (2^0.5)*sigma_a*Z).*log(Life)"),collapse = "")

      if(missing(SUSE)==FALSE){ # Include use life if SUSE is given
        lifeU <- paste(c("(D0./(",paramstress.USE,")).^(1./(mu_a - sigma_a.*0.2453407.*1.414214));"),collapse = "")
      }
    }
  }

  if(dl=="Logarithmic"){
    if(is.null(modelstress) == TRUE){       # Degradation model text for single stress level scenario (parameters a and b)
      # D = a + b*ln(t)
      # theta[1] ~ a, theta[2] ~ b
      dlparams <- "real a; real b;"
      dlparamsvec <- c("a","b")
      pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
      dlpriors <- paste(c(pr1,pr2),collapse = " ")

      degradationlife <- "a + b*log(Life)"
      logdegradationlife <- "log(a + b*log(Life))"
    } else{                                 # Degradation model text for multiple stress level scenario (parameters mu_b and sigma_b)
      dlparams <-  "real mu_b; real<lower=0> sigma_b;"
      dlparamsvec <- c("mu_b","sigma_b")
      outputdlparamset <- c("\U03BC_b","\U03C3_b")
      pr1_<-paste(c("mu_b ~ ",priors[ishift+ishift2+1],";"),collapse = "")
      pr2_<-paste(c("sigma_b ~ ",priors[ishift+ishift2+2],";"),collapse = "")
      dlpriors <- paste(c(pr1_,pr2_),collapse = " ")

      degradationlife <- paste(c("(",paramstress,") + (mu_b + (2^0.5)*sigma_b*Z).*log(Life)"),collapse = "")
      logdegradationlife <- paste(c("log((",paramstress,") + (mu_b + (2^0.5)*sigma_b*Z)*.(Life))"),collapse = "")

      if(missing(SUSE)==FALSE){ # Include use life if SUSE is given
        lifeU <- paste(c("(D0./(",paramstress.USE,")).^(1./(mu_a - sigma_a.*0.2453407.*1.414214));"),collapse = "")
      }
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
    if(is.null(modelstress) == TRUE){       # Degradation model text for single stress level scenario (parameters a and b)
      # D = a + b/t
      # theta[1] ~ a, theta[2] ~ b
      dlparams <- "real a; real b;"
      dlparamsvec <- c("a","b")
      pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
      dlpriors <- paste(c(pr1,pr2),collapse = " ")

      degradationlife <- "a + b/Life"
      logdegradationlife <- "log(a + b/Life)"
    } else{                                 # Degradation model text for multiple stress level scenario (parameters mu_b and sigma_b)
      dlparams <-  "real mu_b; real<lower=0> sigma_b;"
      dlparamsvec <- c("mu_b","sigma_b")
      outputdlparamset <- c("\U03BC_b","\U03C3_b")
      pr1_<-paste(c("mu_b ~ ",priors[ishift+ishift2+1],";"),collapse = "")
      pr2_<-paste(c("sigma_b ~ ",priors[ishift+ishift2+2],";"),collapse = "")
      dlpriors <- paste(c(pr1_,pr2_),collapse = " ")

      degradationlife <- paste(c("(",paramstress,") + (mu_b + (2^0.5)*sigma_b*Z)/Life"),collapse = "")
      logdegradationlife <- paste(c("log((",paramstress,") + (mu_b + (2^0.5)*sigma_b*Z)/Life)"),collapse = "")

      if(missing(SUSE)==FALSE){ # Include use life if SUSE is given
        lifeU <- paste(c("(D0./(",paramstress.USE,")).^(1./(mu_a - sigma_a.*0.2453407.*1.414214));"),collapse = "")
      }
    }
  }

  if(dl=="Mitsuom"){
    # D = 1/(1 + b*(t^a))
    # theta[1] ~ a, theta[2] ~ b
    Dt_txt<-"1/(1 + b*(t^a))"

    if(is.null(modelstress) == TRUE){       # Degradation model text for single stress level scenario (parameters a and b)
      # D = 1/(1 + b*(t^a))
      # theta[1] ~ a, theta[2] ~ b
      dlparams <- "real a; real<lower=0> b;"
      dlparamsvec <- c("a","b")
      pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
      dlpriors <- paste(c(pr1,pr2),collapse = " ")

      degradationlife <- "1/(1 + b*(Life^a))"
      logdegradationlife <- "-log(1 + b*(Life^a))"
    } else{                                 # Degradation model text for multiple stress level scenario (parameters mu_a and sigma_a)
      dlparams <-  "real mu_a; real<lower=0> sigma_a;"
      dlparamsvec <- c("mu_a","sigma_a")
      outputdlparamset <- c("\U03BC_a","\U03C3_a")
      pr1_<-paste(c("mu_a ~ ",priors[ishift+ishift2+1],";"),collapse = "")
      pr2_<-paste(c("sigma_a ~ ",priors[ishift+ishift2+2],";"),collapse = "")
      dlpriors <- paste(c(pr1_,pr2_),collapse = " ")

      degradationlife <- paste(c("1./(1 + (",paramstress,").*(Life.^(mu_a + (2^0.5)*sigma_a*Z)))"),collapse = "")
      logdegradationlife <- paste(c("-log(1 + (",paramstress,").*(Life.^(mu_a + (2^0.5)*sigma_a*Z)))"),collapse = "")

      if(missing(SUSE)==FALSE){ # Include use life if SUSE is given
        lifeU <- paste(c("(D0./(",paramstress.USE,")).^(1./(mu_a - sigma_a.*0.2453407.*1.414214));"),collapse = "")
      }
    }
  }

  if(dl=="PowerVarianceFit"){
    # D = exp(b0)*(Life.^b1)
    # theta[1] ~ b0, theta[2] ~ b1
    dlparams <- "real b0; real b1;"
    dlparamsvec <- c("b0","b1")
    pr1<-paste(c("b0 ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b1 ~ ",priors[ishift+2],";"),collapse = "")
    dlpriors <- paste(c(pr1,pr2),collapse = " ")

    degradationlife <- "exp(b0)*(Life.^b1)"
    logdegradationlife <- "(b0 + b1.*log(Life))"
  }

  # return(list(degradationlife,logdegradationlife))

  # Fit to log-likelihood distributions
  if (dist=="Weibull") {
    distparam <-"real<lower=0> beta;"
    distpriors<-paste(c("beta ~ ",priors[ishift],";"),collapse = "")
    if(is.null(modelstress) == TRUE){       # Degradation model text for single stress level scenario
      loglik <- paste(c("target += weibull_lpdf(TTF | beta, ",degradationlife,");"),collapse = "")

      params <- paste(c(distparam,dlparams),collapse = " ")
      paramsvec <- c("beta",dlparamsvec)
      outputparamset <- c("\U03B2",dlparamsvec)
      priors <- paste(c(distpriors,dlpriors),collapse = " ")
    } else{                                 # Degradation model text for multiple stress level scenario
      loglik <- paste(c("target += log(beta) - 0.5*log(pi) + (beta - 1)*log(Deg_vector) + log(columns_dot_product(Wk.*(",degradationlife,"^-beta).*exp(-((Deg./",degradationlife,").^beta)),ONES_MAT));"),collapse = "")
      params <- paste(c(distparam,psparams,dlparams),collapse = " ")
      paramsvec <- c("beta",psparamsvec,dlparamsvec)
      outputparamset <- c("\U03B2",psparamsvec,outputdlparamset)
      priors <- paste(c(distpriors,pspriors,dlpriors),collapse = " ")
    }
  }
  if (dist=="Lognormal") {
    distparam <-"real<lower=0> sigma_t;"
    distpriors<-paste(c("sigma_t ~ ",priors[ishift],";"),collapse = "")
    if(is.null(modelstress) == TRUE){       # Degradation model text for single stress level scenario
      loglik <- paste(c("target += lognormal_lpdf(Deg |",logdegradationlife,", sigma_t);"),collapse = "")

      params <- paste(c(distparam,dlparams),collapse = " ")
      paramsvec <- c("sigma_t",dlparamsvec)
      outputparamset <- c("\U03C3_t",dlparamsvec)
      priors <- paste(c(distpriors,dlpriors),collapse = " ")
    } else{                                 # Degradation model text for multiple stress level scenario
      loglik <- paste(c("target += -0.5*log(2*(3.141593^2)) - log(sigma_t) - log(Deg_vector) + log(columns_dot_product(Wk.*exp(-0.5*(sigma_t^-2).*((log(Deg) - ",logdegradationlife,").^2)),ONES_MAT));"),collapse = "")

      params <- paste(c(distparam,psparams,dlparams),collapse = " ")
      paramsvec <- c("sigma_t",psparamsvec,dlparamsvec)
      outputparamset <- c("\U03C3_t",psparamsvec,outputdlparamset)
      priors <- paste(c(distpriors,pspriors,dlpriors),collapse = " ")
    }
  }
  if (dist=="Normal") {
    distparam <-"real<lower=0> sigma;"
    distpriors<-paste(c("sigma ~ ",priors[ishift],";"),collapse = "")
    if(is.null(modelstress) == TRUE){       # Degradation model text for single stress level scenario
      loglik <- paste(c("target += normal_lpdf(Deg |",degradationlife,", sigma);"),collapse = "")
      params <- paste(c(distparam,dlparams),collapse = " ")
      paramsvec <- c("sigma",dlparamsvec)
      outputparamset <- c("\U03C3",dlparamsvec)
      priors <- paste(c(distpriors,dlpriors),collapse = " ")
    } else{                                 # Degradation model text for multiple stress level scenario
      loglik <- paste(c("target += -0.5*log(2*(3.141593^2)) - log(sigma) + log(columns_dot_product(Wk.*exp(-0.5*(sigma^-2).*((Deg - ",degradationlife,").^2)),ONES_MAT));"),collapse = "")
      params <- paste(c(distparam,psparams,dlparams),collapse = " ")
      paramsvec <- c("sigma",psparamsvec,dlparamsvec)
      outputparamset <- c("\U03C3",psparamsvec,outputdlparamset)
      priors <- paste(c(distpriors,pspriors,dlpriors),collapse = " ")
    }
  }

  # if (dl=="PowerVarianceFit") {
  #   distparam <-"real alpha_0; real alpha_1;"
  #   distpriors<-paste(c("alpha_0 ~ ",priors[ishift],"; alpha_1 ~ ",priors[ishift+1],";"),collapse = "")
  #   loglik <- paste(c("target += lognormal_lpdf(Deg |",logdegradationlife,", sigma_t);"),collapse = "")
  #
  #   # loglik <- paste(c("target += -0.5*log(2*(3.141593)) - 0.5*log(exp(alpha_0)*exp(alpha_1*(Life - 3.66))) - 0.5.*((exp(alpha_0)*exp(alpha_1.*(Life - 3.66))).^-1).*((log(Deg) - ",logdegradationlife,").^2);"),collapse = "")
  #   params <- paste(c(distparam,dlparams),collapse = " ")
  #   paramsvec <- c("alpha_0","alpha_1",dlparamsvec)
  #   outputparamset <- c("\U03B1_0","\U03B1_1",dlparamsvec)
  #   priors <- paste(c(distpriors,dlpriors),collapse = " ")
  # }

  # return(lsparamsvec)
  # return(list(degradationlife,logdegradationlife,loglik,params,paramsvec,outputparamset,priors))

  # Define stancode here
  # ==========================================
  # BLOCK 1 - DATA STATEMENT AND DATA BLOCK
  # ==========================================
  if(is.null(modelstress) == TRUE && is.null(SUSE)==TRUE){ # CASE 1: One Stress Level analysis
    block1 <- "data {int<lower=0> N; vector[N] Deg; vector[N] Life; real D0;}"
    datablock <- list(N = N, Deg = DegradationFULL, Life = TimeFULL, D0 = D0)
  }
  if(is.null(modelstress) == FALSE && is.null(SUSE)==TRUE && dist=="Normal"){ # CASE 2: Multi Stress Level analysis without Use Stress given (Normal)
    block1 <- "data {int<lower=0> N; int<lower=0> Q; matrix[Q,N] Deg; matrix<lower=0>[Q,N] Life; real D0; matrix[Q,N] Sf; matrix[Q,N] Wk; matrix[Q,N] Z; matrix[Q,N] ONES_MAT;}"
    datablock <- list(N = N, Q = Q,  Deg = Degradation_MAT, Life = Time_MAT, D0 = D0, Sf = Stress_MAT, Wk = W_MAT, Z = Z_MAT, ONES_MAT = Ones_MAT)
  }
  if(is.null(modelstress) == FALSE && is.null(SUSE)==TRUE && (dist=="Lognormal" || dist=="Weibull")){ # CASE 3: Multi Stress Level analysis without Use Stress given (Lognormal or Weibull)
    block1 <- "data {int<lower=0> N; int<lower=0> Q; matrix[Q,N] Deg; vector[N] Deg_vector; matrix<lower=0>[Q,N] Life; real D0; matrix[Q,N] Sf; matrix[Q,N] Wk; matrix[Q,N] Z; matrix[Q,N] ONES_MAT;}"
    datablock <- list(N = N, Q = Q,  Deg = Degradation_MAT, Deg_vector = DegradationFULL, Life = Time_MAT, D0 = D0, Sf = Stress_MAT, Wk = W_MAT, Z = Z_MAT, ONES_MAT = Ones_MAT)
  }
  if(is.null(modelstress) == FALSE && is.null(SUSE)==FALSE && dist=="Normal"){ # CASE 4: Multi Stress Level analysis with Use Stress given (Normal)
    block1 <- "data {int<lower=0> N; int<lower=0> Q; matrix[Q,N] Deg; matrix<lower=0>[Q,N] Life; real D0; matrix[Q,N] Sf; matrix[Q,N] Wk; matrix[Q,N] Z; matrix[Q,N] ONES_MAT; real Suse;}"
    datablock <- list(N = N, Q = Q,  Deg = Degradation_MAT, Life = Time_MAT, D0 = D0, Sf = Stress_MAT, Wk = W_MAT, Z = Z_MAT, ONES_MAT = Ones_MAT, Suse = SUSE)
  }
  if(is.null(modelstress) == FALSE && is.null(SUSE)==FALSE && (dist=="Lognormal" || dist=="Weibull")){ # CASE 5: Multi Stress Level analysis with Use Stress given (Lognormal or Weibull)
    block1 <- "data {int<lower=0> N; int<lower=0> Q; matrix[Q,N] Deg; vector[N] Deg_vector; matrix<lower=0>[Q,N] Life; real D0; matrix[Q,N] Sf; matrix[Q,N] Wk; matrix[Q,N] Z; matrix[Q,N] ONES_MAT; real Suse;}"
    datablock <- list(N = N, Q = Q,  Deg = Degradation_MAT, Deg_vector = DegradationFULL, Life = Time_MAT, D0 = D0, Sf = Stress_MAT, Wk = W_MAT, Z = Z_MAT, ONES_MAT = Ones_MAT, Suse = SUSE)
  }
  # return(list(block1,datablock))

  # ==========================================
  # BLOCK 2 - PARAMETER STATEMENT
  # ==========================================
  block2 <- paste(c("parameters {",params,"}"),collapse = " ")

  # ==========================================
  # BLOCK 2b - TRANSFORMED PARAMETER STATEMENT (IGNORE FOR NOW)
  # ==========================================
  if(is.null(SUSE)==FALSE){
    block2b <- paste(c("transformed parameters { real<lower=0> Uselife; Uselife = ",lifeU,"}"),collapse = " ")
    paramsvec0 <- c(paramsvec,"Uselife")
    # pt_est <- c(pt_est,complifeU)
  }

  if(is.null(SUSE)==TRUE){
    paramsvec0 <- paramsvec
  }
  if(dl=="PowerVarianceFit"){
    block2b <- paste(c("transformed parameters { real<lower=0> time_to_failure; time_to_failure = exp((log(D0) - b0)/b1);}"),collapse = " ")
    paramsvec0 <- c(paramsvec,"time_to_failure")
    # pt_est <- c(pt_est,complifeU,compAF)
  }
  # ==========================================
  # BLOCK 3 - MODEL STATEMENT
  # ==========================================
  block3 <- paste(c("model {",priors,loglik,"}"),collapse = " ")

  stanlscode <- paste(c(block1,block2,block3),collapse=" ")
  if(is.null(SUSE)==FALSE){
    stanlscode <- paste(c(block1,block2,block2b,block3),collapse=" ")
  }
  if(dl=="PowerVarianceFit"){
    stanlscode <- paste(c(block1,block2,block2b,block3),collapse=" ")
  }
  # return(stanlscode)
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
  # return(fit)
  # Print results.  I need to get this as an output
  # ==================================================================================================
  # NOTE RCS01102026 - The following block works under cmdstanr which is not functioning at the moment
  # ==================================================================================================
  # stats <- fit$summary(variables = paramsvec0)
  # confidbounds <- mcmc_intervals_data(fit$draws(variables = paramsvec0),prob_outer = confid)
  # outputtable <- matrix(c(stats[[2]],stats[[4]],confidbounds[[5]],stats[[3]],confidbounds[[9]],stats[[8]]), nrow = length(paramsvec0), ncol = 6, byrow = FALSE,dimnames = list(paramsvec0,c("Mean","Standard Deviation",conflim_txt[1],"Median",conflim_txt[2],"R\U005E")))
  # ==================================================================================================
  # PATCH RCS01102026 - The following block works under rstan which IS functioning at the moment
  # ==================================================================================================
  stats.mean.sd <- summary(fit)$summary[,c(1,3)]
  stats.Rhat <- rhat(fit)
  confidbounds <- mcmc_intervals_data(data.frame(extract(fit, paramsvec0)),prob_outer = confid)
  outputtable <- matrix(c(unname(stats.mean.sd)[1:length(paramsvec0),1],unname(stats.mean.sd)[1:length(paramsvec0),2],confidbounds[[5]],confidbounds[[7]],confidbounds[[9]],unname(stats.Rhat)[1:length(paramsvec0)]), nrow = length(paramsvec0), ncol = 6, byrow = FALSE,dimnames = list(paramsvec0,c("Mean","Standard Deviation",conflim_txt[1],"Median",conflim_txt[2],"R\U005E")))

  # return(fit)
  # ==================================================================================================
  # NOTE RCS01102026 - The following block works under cmdstanr which is not functioning at the moment
  # ==================================================================================================
  # plot1_MCtrace <- mcmc_trace(fit$draws(paramsvec0))
  # plot2_hist <- mcmc_hist(fit$draws(paramsvec0))
  # ==================================================================================================
  # PATCH RCS01102026 - The following block works under rstan which IS functioning at the moment
  # ==================================================================================================
  # Trace the Markov Chains for each parameter
  plot1_MCtrace <- mcmc_trace(fit,paramsvec0) +
    theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4))
  # Plot histogram
  plot2_hist <- mcmc_hist(fit,paramsvec0) +
    theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4))

  # NEW Post Bayes analysis plotting of posterior (12/12/25)
  if (dist=="Normal") {
    posterior_sigma <- density(extract(fit,c("sigma"))$sigma)
    if(is.null(modelstress)==TRUE){ # CASE 1 (Disabled for now): Map the posterior degradation-life parameters only

    }
    if(is.null(modelstress)==FALSE){ # CASE 2: Map posterior for model-stress parameters and degradation life parameters
      if (modelstress=="Linear" || modelstress=="Exponential" || modelstress=="Exponential2" || modelstress=="Eyring" || modelstress=="Eyring2" ||
          modelstress=="Power" || modelstress=="InversePower" || modelstress=="InversePower2" || modelstress=="Logarithmic"){ # Parameters a and b
        posterior_a_0 <- density(extract(fit,c("a_0"))$a_0)
        posterior_b_0 <- density(extract(fit,c("b_0"))$b_0)
        X_DF_POST <- c(posterior_sigma$x,posterior_a_0$x,posterior_b_0$x)
        Y_MIN_DF_POST <- rep(0,(length(posterior_sigma$x)+length(posterior_a_0$x)+length(posterior_b_0$x)))
        Y_MAX_DF_POST <- c(posterior_sigma$y,posterior_a_0$y,posterior_b_0$y)
        DISTLABEL_DF_POST <-c(rep("σ",length(posterior_sigma$x)),rep("a_0",length(posterior_a_0$x)),rep("b_0",length(posterior_b_0$x)))
      }
    }

    if (modelstress=="Arrhenius") {
      posterior_Ea_0 <- density(extract(fit,c("E_a_0"))$E_a_0)
      posterior_b_0 <- density(extract(fit,c("b_0"))$b_0)
      X_DF_POST <- c(posterior_sigma$x,posterior_Ea_0$x,posterior_b_0$x)
      Y_MIN_DF_POST <- rep(0,(length(posterior_sigma$x)+length(posterior_Ea_0$x)+length(posterior_b_0$x)))
      Y_MAX_DF_POST <- c(posterior_sigma$y,posterior_Ea_0$y,posterior_b_0$y)
      DISTLABEL_DF_POST <- c(rep("σ",length(posterior_sigma$x)),rep("Ea_0",length(posterior_Ea_0$x)),rep("b_0",length(posterior_b_0$x)))
    }
    if (modelstress=="TempHumidity"){
      posterior_A_0 <- density(extract(fit,c("A_0"))$A_0)
      posterior_a_0 <- density(extract(fit,c("a_0"))$a_0)
      posterior_b_0 <- density(extract(fit,c("b_0"))$b_0)
      X_DF_POST <- c(posterior_sigma$x,posterior_A_0$x,posterior_a_0$x,posterior_b_0$x)
      Y_MIN_DF_POST <- rep(0,(length(posterior_sigma$x)+length(posterior_A_0$x)+length(posterior_a_0$x)+length(posterior_b_0$x)))
      Y_MAX_DF_POST <- c(posterior_sigma$y,posterior_A_0$y,posterior_a_0$y,posterior_b_0$y)
      DISTLABEL_DF_POST <- c(rep("σ",length(posterior_sigma$x)),rep("A_0",length(posterior_A_0$x)),rep("a_0",length(posterior_a_0$x)),rep("b_0",length(posterior_b_0$x)))
    }
    if (modelstress=="TempNonthermal"){
      posterior_a_0 <- density(extract(fit,c("a_0"))$a_0)
      posterior_b_0 <- density(extract(fit,c("b_0"))$b_0)
      posterior_c_0 <- density(extract(fit,c("c_0"))$c_0)
      X_DF_POST <- c(posterior_sigma$x,posterior_a_0$x,posterior_b_0$x,posterior_c_0$x)
      Y_MIN_DF_POST <- rep(0,(length(posterior_sigma$x)+length(posterior_a_0$x)+length(posterior_b_0$x)+length(posterior_c_0$x)))
      Y_MAX_DF_POST <- c(posterior_sigma$y,posterior_a_0$y,posterior_b_0$y,posterior_c_0$y)
      DISTLABEL_DF_POST <- c(rep("σ",length(posterior_sigma$x)),rep("a_0",length(posterior_a_0$x)),rep("b_0",length(posterior_b_0$x)),rep("c_0",length(posterior_c_0$x)))
    }
    if (modelstress=="Eyring3"){
      posterior_a_0 <- density(extract(fit,c("a_0"))$a_0)
      posterior_b_0 <- density(extract(fit,c("b_0"))$b_0)
      posterior_c_0 <- density(extract(fit,c("c_0"))$c_0)
      posterior_d_0 <- density(extract(fit,c("d_0"))$d_0)
      X_DF_POST <- c(posterior_sigma$x,posterior_a_0$x,posterior_b_0$x,posterior_c_0$x,posterior_d_0$x)
      Y_MIN_DF_POST <- rep(0,(length(posterior_sigma$x)+length(posterior_a_0$x)+length(posterior_b_0$x)+length(posterior_c_0$x)+length(posterior_d_0$x)))
      Y_MAX_DF_POST <- cc(posterior_sigma$y,posterior_a_0$y,posterior_b_0$y,posterior_c_0$y,posterior_d_0$y)
      DISTLABEL_DF_POST <- c(rep("σ",length(posterior_sigma$x)),rep("a_0",length(posterior_a_0$x)),rep("b_0",length(posterior_b_0$x)),rep("c_0",length(posterior_c_0$x)),rep("d_0",length(posterior_d_0$x)))
    }
    if (modelstress=="Eyring4"){
      posterior_A_0 <- density(extract(fit,c("A_0"))$A_0)
      posterior_Ea_0 <- density(extract(fit,c("E_a_0"))$E_a_0)
      posterior_b_0 <- density(extract(fit,c("b_0"))$b_0)
      X_DF_POST <- c(posterior_sigma$x,posterior_A_0$x,posterior_Ea_0$x,posterior_b_0$x)
      Y_MIN_DF_POST <- rep(0,(length(posterior_sigma$x)+length(posterior_A_0$x)+length(posterior_Ea_0$x)+length(posterior_b_0$x)))
      Y_MAX_DF_POST <- c(posterior_sigma$y,posterior_A_0$y,posterior_Ea_0$y,posterior_b_0$y)
      DISTLABEL_DF_POST <-c(rep("σ",length(posterior_sigma$x)),rep("A_0",length(posterior_A$x)),rep("Ea_0",length(posterior_Ea_0$x)),rep("b_0",length(posterior_b_0$x)))
    }
    # Finalize the posterior parameter distribution
    if(dl=="Linear" || dl=="SquareRoot" || dl=="SquareRoot2" || dl=="Logarithmic" || dl=="LloydLipow"){
      posterior_mu_b <- density(extract(fit,c("mu_b"))$mu_b)
      posterior_sigma_b <- density(extract(fit,c("sigma_b"))$sigma_b)
      X_DF_POST <- c(X_DF_POST,posterior_mu_b$x,posterior_sigma_b$x)
      Y_MIN_DF_POST <- c(Y_MIN_DF_POST,rep(0,(length(posterior_mu_b$x)+length(posterior_sigma_b$x))))
      Y_MAX_DF_POST <- c(Y_MAX_DF_POST,posterior_mu_b$y,posterior_sigma_b$y)
      DISTLABEL_DF_POST <-c(DISTLABEL_DF_POST,rep("mu_b",length(posterior_mu_b$x)),rep("sigma_b",length(posterior_sigma_b$x)))
    }
    if(dl=="Exponential" || dl=="Power" || dl=="Mitsuom"){
      posterior_mu_a <- density(extract(fit,c("mu_a"))$mu_a)
      posterior_sigma_a <- density(extract(fit,c("sigma_a"))$sigma_a)
      X_DF_POST <- c(X_DF_POST,posterior_mu_a$x,posterior_sigma_a$x)
      Y_MIN_DF_POST <- c(Y_MIN_DF_POST,rep(0,(length(posterior_mu_a$x)+length(posterior_sigma_a$x))))
      Y_MAX_DF_POST <- c(Y_MAX_DF_POST,posterior_mu_a$y,posterior_sigma_a$y)
      DISTLABEL_DF_POST <-c(DISTLABEL_DF_POST,rep("mu_a",length(posterior_mu_a$x)),rep("posterior_sigma_a",length(posterior_sigma_a$x)))
    }
    df_posterior <- data.frame(x = X_DF_POST,
                               ymin = Y_MIN_DF_POST,
                               ymax = Y_MAX_DF_POST,
                               distlabel = DISTLABEL_DF_POST)

    # Density plot for ALT parameter posterior
    plot3_density <- ggplot() + geom_ribbon(data = df_posterior, aes(x=x, ymin = ymin, ymax = ymax), fill = "red" ,alpha = 0.5) +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      facet_wrap(~distlabel, dir="v", scales = "free") +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0)) +
      xlab(" ") +
      ylab("density")
  }

  if (dist=="Lognormal") {
    posterior_sigma_t <- density(extract(fit,c("sigma_t"))$sigma_t)
    if(is.null(modelstress)==TRUE){ # CASE 1 (Disabled for now): Map the posterior degradation-life parameters only

    }
    if(is.null(modelstress)==FALSE){ # CASE 2: Map posterior for model-stress parameters and degradation life parameters
      if (modelstress=="Linear" || modelstress=="Exponential" || modelstress=="Exponential2" || modelstress=="Eyring" || modelstress=="Eyring2" ||
          modelstress=="Power" || modelstress=="InversePower" || modelstress=="InversePower2" || modelstress=="Logarithmic"){ # Parameters a and b
        posterior_a_0 <- density(extract(fit,c("a_0"))$a_0)
        posterior_b_0 <- density(extract(fit,c("b_0"))$b_0)
        X_DF_POST <- c(posterior_sigma_t$x,posterior_a_0$x,posterior_b_0$x)
        Y_MIN_DF_POST <- rep(0,(length(posterior_sigma_t$x)+length(posterior_a_0$x)+length(posterior_b_0$x)))
        Y_MAX_DF_POST <- c(posterior_sigma_t$y,posterior_a_0$y,posterior_b_0$y)
        DISTLABEL_DF_POST <-c(rep("σ_t",length(posterior_sigma_t$x)),rep("a_0",length(posterior_a_0$x)),rep("b_0",length(posterior_b_0$x)))
      }
    }
    if (modelstress=="Arrhenius") {
      posterior_Ea_0 <- density(extract(fit,c("E_a_0"))$E_a_0)
      posterior_b_0 <- density(extract(fit,c("b_0"))$b_0)
      X_DF_POST <- c(posterior_sigma_t$x,posterior_Ea_0$x,posterior_b_0$x)
      Y_MIN_DF_POST <- rep(0,(length(posterior_sigma_t$x)+length(posterior_Ea_0$x)+length(posterior_b_0$x)))
      Y_MAX_DF_POST <- c(posterior_sigma_t$y,posterior_Ea_0$y,posterior_b_0$y)
      DISTLABEL_DF_POST <- c(rep("σ_t",length(posterior_sigma_t$x)),rep("Ea_0",length(posterior_Ea_0$x)),rep("b_0",length(posterior_b_0$x)))
    }
    if (modelstress=="TempHumidity"){
      posterior_A_0 <- density(extract(fit,c("A_0"))$A_0)
      posterior_a_0 <- density(extract(fit,c("a_0"))$a_0)
      posterior_b_0 <- density(extract(fit,c("b_0"))$b_0)
      X_DF_POST <- c(posterior_sigma_t$x,posterior_A_0$x,posterior_a_0$x,posterior_b_0$x)
      Y_MIN_DF_POST <- rep(0,(length(posterior_sigma_t$x)+length(posterior_A_0$x)+length(posterior_a_0$x)+length(posterior_b_0$x)))
      Y_MAX_DF_POST <- c(posterior_sigma_t$y,posterior_A_0$y,posterior_a_0$y,posterior_b_0$y)
      DISTLABEL_DF_POST <- c(rep("σ_t",length(posterior_sigma_t$x)),rep("A_0",length(posterior_A_0$x)),rep("a_0",length(posterior_a_0$x)),rep("b_0",length(posterior_b_0$x)))
    }
    if (modelstress=="TempNonthermal"){
      posterior_a_0 <- density(extract(fit,c("a_0"))$a_0)
      posterior_b_0 <- density(extract(fit,c("b_0"))$b_0)
      posterior_c_0 <- density(extract(fit,c("c_0"))$c_0)
      X_DF_POST <- c(posterior_sigma_t$x,posterior_a_0$x,posterior_b_0$x,posterior_c_0$x)
      Y_MIN_DF_POST <- rep(0,(length(posterior_sigma_t$x)+length(posterior_a_0$x)+length(posterior_b_0$x)+length(posterior_c_0$x)))
      Y_MAX_DF_POST <- c(posterior_sigma_t$y,posterior_a_0$y,posterior_b_0$y,posterior_c_0$y)
      DISTLABEL_DF_POST <- c(rep("σ_t",length(posterior_sigma_t$x)),rep("a_0",length(posterior_a_0$x)),rep("b_0",length(posterior_b_0$x)),rep("c_0",length(posterior_c_0$x)))
    }
    if (modelstress=="Eyring3"){
      posterior_a_0 <- density(extract(fit,c("a_0"))$a_0)
      posterior_b_0 <- density(extract(fit,c("b_0"))$b_0)
      posterior_c_0 <- density(extract(fit,c("c_0"))$c_0)
      posterior_d_0 <- density(extract(fit,c("d_0"))$d_0)
      X_DF_POST <- c(posterior_sigma_t$x,posterior_a_0$x,posterior_b_0$x,posterior_c_0$x,posterior_d_0$x)
      Y_MIN_DF_POST <- rep(0,(length(posterior_sigma_t$x)+length(posterior_a_0$x)+length(posterior_b_0$x)+length(posterior_c_0$x)+length(posterior_d_0$x)))
      Y_MAX_DF_POST <- cc(posterior_sigma_t$y,posterior_a_0$y,posterior_b_0$y,posterior_c_0$y,posterior_d_0$y)
      DISTLABEL_DF_POST <- c(rep("σ_t",length(posterior_sigma_t$x)),rep("a_0",length(posterior_a_0$x)),rep("b_0",length(posterior_b_0$x)),rep("c_0",length(posterior_c_0$x)),rep("d_0",length(posterior_d_0$x)))
    }
    if (modelstress=="Eyring4"){
      posterior_A_0 <- density(extract(fit,c("A_0"))$A_0)
      posterior_Ea_0 <- density(extract(fit,c("E_a_0"))$E_a_0)
      posterior_b_0 <- density(extract(fit,c("b_0"))$b_0)
      X_DF_POST <- c(posterior_sigma_t$x,posterior_A_0$x,posterior_Ea_0$x,posterior_b_0$x)
      Y_MIN_DF_POST <- rep(0,(length(posterior_sigma_t$x)+length(posterior_A_0$x)+length(posterior_Ea_0$x)+length(posterior_b_0$x)))
      Y_MAX_DF_POST <- c(posterior_sigma_t$y,posterior_A_0$y,posterior_Ea_0$y,posterior_b_0$y)
      DISTLABEL_DF_POST <-c(rep("σ_t",length(posterior_sigma_t$x)),rep("A_0",length(posterior_A$x)),rep("Ea_0",length(posterior_Ea_0$x)),rep("b_0",length(posterior_b_0$x)))
    }
    # Finalize the posterior parameter distribution
    if(dl=="Linear" || dl=="SquareRoot" || dl=="SquareRoot2" || dl=="Logarithmic" || dl=="LloydLipow"){
      posterior_mu_b <- density(extract(fit,c("mu_b"))$mu_b)
      posterior_sigma_b <- density(extract(fit,c("sigma_b"))$sigma_b)
      X_DF_POST <- c(X_DF_POST,posterior_mu_b$x,posterior_sigma_b$x)
      Y_MIN_DF_POST <- c(Y_MIN_DF_POST,rep(0,(length(posterior_mu_b$x)+length(posterior_sigma_b$x))))
      Y_MAX_DF_POST <- c(Y_MAX_DF_POST,posterior_mu_b$y,posterior_sigma_b$y)
      DISTLABEL_DF_POST <-c(DISTLABEL_DF_POST,rep("mu_b",length(posterior_mu_b$x)),rep("sigma_b",length(posterior_sigma_b$x)))
    }
    if(dl=="Exponential" || dl=="Power" || dl=="Mitsuom"){
      posterior_mu_a <- density(extract(fit,c("mu_a"))$mu_a)
      posterior_sigma_a <- density(extract(fit,c("sigma_a"))$sigma_a)
      X_DF_POST <- c(X_DF_POST,posterior_mu_a$x,posterior_sigma_a$x)
      Y_MIN_DF_POST <- c(Y_MIN_DF_POST,rep(0,(length(posterior_mu_a$x)+length(posterior_sigma_a$x))))
      Y_MAX_DF_POST <- c(Y_MAX_DF_POST,posterior_mu_a$y,posterior_sigma_a$y)
      DISTLABEL_DF_POST <-c(DISTLABEL_DF_POST,rep("mu_a",length(posterior_mu_a$x)),rep("posterior_sigma_a",length(posterior_sigma_a$x)))
    }
    df_posterior <- data.frame(x = X_DF_POST,
                               ymin = Y_MIN_DF_POST,
                               ymax = Y_MAX_DF_POST,
                               distlabel = DISTLABEL_DF_POST)

    # Density plot for ALT parameter posterior
    plot3_density <- ggplot() + geom_ribbon(data = df_posterior, aes(x=x, ymin = ymin, ymax = ymax), fill = "red" ,alpha = 0.5) +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      facet_wrap(~distlabel, dir="v", scales = "free") +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0)) +
      xlab(" ") +
      ylab("density")
  }

  # Produce some output text that summarizes the results
  cat(c("Posterior estimates for Bayesian Analysis.\n\n"),sep = "")
  print(outputtable)
  cat(c("\n"),sep = "")

  if(is.null(modelstress)==TRUE){ # CASE 1 (Disabled for now): Map the posterior degradation-life parameters only

  }
  if(is.null(modelstress)==FALSE && is.null(SUSE)==TRUE){ # CASE 2: Map posterior for model-stress parameters and degradation life parameters
    return(list(posterior.fit=fit,post.stats=stats.mean.sd,MC.trace=plot1_MCtrace,post.histogram=plot2_hist,post.density=plot3_density,stanlscode))
  }
  if(is.null(modelstress)==FALSE && is.null(SUSE)==FALSE){ # Include plot of use life posterior
    posterior_Uselife <- density(extract(fit,c("Uselife"))$Uselife)
    df_posterior_Uselife <- data.frame(x = posterior_Uselife$x,
                                       ymin = rep(0,length(posterior_Uselife$x)),
                                       ymax = posterior_Uselife$y,
                                       distlabel = c(rep("Use Level Life",length(posterior_Uselife$x))))
    plot3A_USelife.density <- ggplot() + geom_ribbon(data = df_posterior_Uselife, aes(x=x, ymin = ymin, ymax = ymax), fill = "red" ,alpha = 0.5) +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0)) +
      xlab("Use Life") +
      ylab("density")
    return(list(posterior.fit=fit,post.stats=stats.mean.sd,MC.trace=plot1_MCtrace,post.histogram=plot2_hist,post.density=plot3_density,post.Uselife.density=plot3A_USelife.density,stanlscode))
  }
  # NOTE: Comment return above and uncomment return below if you want to view the scatterplot between posteriors.
  # This may however increase processing time.
}
