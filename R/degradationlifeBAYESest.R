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
    }
    if((length(priors)-ishift)==3){
      ishift2<-3                            # Second ishift for number of parameters
      psparams <- "real a0_0; real a1_0; real a2_0; "
      psparamsvec <- c("a0_0","a1_0","a2_0")
      pr1<-paste(c("a0_0 ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("a1_0 ~ ",priors[ishift+2],";"),collapse = "")
      pr3<-paste(c("a2_0 ~ ",priors[ishift+3],";"),collapse = "")
      pspriors <- paste(c(pr1,pr2,pr3),collapse = " ")

      paramstress <- "Lifei[i] = exp(a0_0 + a1_0*Sf[i,1] + a2_0*Sf[i,2]);"
      logparamstress <- "Lifei[i] = a0_0 + a1_0*Sf[i,1] + a2_0*Sf[i,2];"
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

      paramstress <- "Lifei[i] = exp(a0_0 + a1_0*Sf[i,1] + a2_0*Sf[i,2] + a3_0*Sf[i,3]);"
      logparamstress <- "Lifei[i] = a0_0 + a1_0*Sf[i,1] + a2_0*Sf[i,2] + a3_0*Sf[i,3];"
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

      paramstress <- "Lifei[i] = exp(a0_0 + a1_0*Sf[i,1] + a2_0*Sf[i,2] + a3_0*Sf[i,3] + a4_0*Sf[i,4]);"
      logparamstress <- "Lifei[i] = a0_0 + a1_0*Sf[i,1] + a2_0*Sf[i,2] + a3_0*Sf[i,3] + a4_0*Sf[i,4];"
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
    }
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
      paramsvec <- c("sigma_t",lsparamsvec)
      outputparamset <- c("\U03C3_t",lsparamsvec)
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

  # return(list(degradationlife,logdegradationlife,loglik,params,paramsvec,outputparamset,priors))

  # Define stancode here
  # ==========================================
  # BLOCK 1 - DATA STATEMENT AND DATA BLOCK
  # ==========================================
  if(is.null(modelstress) == TRUE && is.null(SUSE)==TRUE){ # CASE 1: One Stress Level analysis
    block1 <- "data {int<lower=0> N; vector[N] Deg; vector[N] Life;}"
    datablock <- list(N = N, Deg = DegradationFULL, Life = TimeFULL)
  }
  if(is.null(modelstress) == FALSE && is.null(SUSE)==TRUE && dist=="Normal"){ # CASE 2: Multi Stress Level analysis without Use Stress given (Normal)
    block1 <- "data {int<lower=0> N; int<lower=0> Q; matrix[Q,N] Deg; matrix<lower=0>[Q,N] Life; matrix[Q,N] Sf; matrix[Q,N] Wk; matrix[Q,N] Z; matrix[Q,N] ONES_MAT;}"
    datablock <- list(N = N, Q = Q,  Deg = Degradation_MAT, Life = Time_MAT, Sf = Stress_MAT, Wk = W_MAT, Z = Z_MAT, ONES_MAT = Ones_MAT)
  }
  if(is.null(modelstress) == FALSE && is.null(SUSE)==TRUE && (dist=="Lognormal" || dist=="Weibull")){ # CASE 3: Multi Stress Level analysis without Use Stress given (Lognormal or Weibull)
    block1 <- "data {int<lower=0> N; int<lower=0> Q; matrix[Q,N] Deg; vector[N] Deg_vector; matrix<lower=0>[Q,N] Life; matrix[Q,N] Sf; matrix[Q,N] Wk; matrix[Q,N] Z; matrix[Q,N] ONES_MAT;}"
    datablock <- list(N = N, Q = Q,  Deg = Degradation_MAT, Deg_vector = DegradationFULL, Life = Time_MAT, Sf = Stress_MAT, Wk = W_MAT, Z = Z_MAT, ONES_MAT = Ones_MAT)
  }
  if(is.null(modelstress) == FALSE && is.null(SUSE)==FALSE && dist=="Normal"){ # CASE 4: Multi Stress Level analysis with Use Stress given (Normal)
    block1 <- "data {int<lower=0> N; int<lower=0> Q; matrix[Q,N] Deg; matrix<lower=0>[Q,N] Life; matrix[Q,N] Sf; matrix[Q,N] Wk; matrix[Q,N] Z; matrix[Q,N] ONES_MAT; real Suse;}"
    datablock <- list(N = N, Q = Q,  Deg = Degradation_MAT, Life = Time_MAT, Sf = Stress_MAT, Wk = W_MAT, Z = Z_MAT, ONES_MAT = Ones_MAT, Suse = SUSE)
  }
  if(is.null(modelstress) == FALSE && is.null(SUSE)==FALSE && (dist=="Lognormal" || dist=="Weibull")){ # CASE 5: Multi Stress Level analysis with Use Stress given (Lognormal or Weibull)
    block1 <- "data {int<lower=0> N; int<lower=0> Q; matrix[Q,N] Deg; vector[N] Deg_vector; matrix<lower=0>[Q,N] Life; matrix[Q,N] Sf; matrix[Q,N] Wk; matrix[Q,N] Z; matrix[Q,N] ONES_MAT; real Suse;}"
    datablock <- list(N = N, Q = Q,  Deg = Degradation_MAT, Deg_vector = DegradationFULL, Life = Time_MAT, Sf = Stress_MAT, Wk = W_MAT, Z = Z_MAT, ONES_MAT = Ones_MAT, Suse = SUSE)
  }
  # return(list(block1,datablock))

  # ==========================================
  # BLOCK 2 - PARAMETER STATEMENT
  # ==========================================
  block2 <- paste(c("parameters {",params,"}"),collapse = " ")
  # return(list(block1,block2))

  # ==========================================
  # BLOCK 2b - TRANSFORMED PARAMETER STATEMENT (IGNORE FOR NOW)
  # ==========================================
  # if(is.null(SUSE)==FALSE){
  #   block2b <- paste(c("transformed parameters { real<lower=0> Uselife; Uselife = ",lifeU,"}"),collapse = " ")
  #   paramsvec0 <- c(paramsvec,"Uselife")
  #   # pt_est <- c(pt_est,complifeU)
  # }
  # if(is.null(SUSE)==FALSE && is.null(SACC)==FALSE){
  #   block2b <- paste(c("transformed parameters { real<lower=0> Uselife; real<lower=0> ",AFheading,"; Uselife = ",lifeU,AFheading," = ",AF,"}"),collapse = " ")
  #   paramsvec0 <- c(paramsvec,"Uselife",AFheading)
  #   # pt_est <- c(pt_est,complifeU,compAF)
  # }
  if(is.null(SUSE)==TRUE){
    paramsvec0 <- paramsvec
  }
  # ==========================================
  # BLOCK 3 - MODEL STATEMENT
  # ==========================================
  block3 <- paste(c("model {",priors,loglik,"}"),collapse = " ")
  # if (ls=="MultiStress" && dist == "Lognormal"){
  #   if(is.null(Tc)==TRUE){
  #     block3 <- paste(c("model { vector[n] Lifei; vector[m] Lifej; ",priors," for(i in 1:n){",loglifeF,"}",loglik,"}"),collapse = " ")
  #   }
  #   if(is.null(Tc)==FALSE){
  #     block3 <- paste(c("model { vector[n] Lifei; vector[m] Lifej; ",priors," for(i in 1:n){",loglifeF,"} for(j in 1:m){",loglifeC,"}",loglik,"}"),collapse = " ")
  #   }
  # }
  # if (ls=="MultiStress" && (dist == "Normal" || dist=="Weibull" || dist=="Exponential")){
  #   if(is.null(Tc)==TRUE){
  #     block3 <- paste(c("model { vector[n] Lifei; ",priors," for(i in 1:n){",lifeF,"}",loglik,"}"),collapse = " ")
  #   }
  #   if(is.null(Tc)==FALSE){
  #     block3 <- paste(c("model { vector[n] Lifei; vector[m] Lifej; ",priors," for(i in 1:n){",lifeF,"} for(j in 1:m){",lifeC,"}",loglik,"}"),collapse = " ")
  #   }
  # }
  # NOT RUN {
  stanlscode <- paste(c(block1,block2,block3),collapse=" ")
  # if(is.null(SUSE)==FALSE || is.null(SACC)==FALSE){
  #   stanlscode <- paste(c(block1,block2,block2b,block3),collapse=" ")
  # }
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
  # return(list(stanlscode,stanlsfile,init_pt_est,paramsvec0))

  # lsmod <- stan_model(model_code = stanlscode, verbose = TRUE)
  lsmod <- cmdstan_model(stanlsfile)
  # return(lsmod)
  # fit <- sampling(lsmod, data = datablock, iter = nsamples, warmup = burnin, init = pt_est)
  fit <- lsmod$sample(data = datablock, init = init_pt_est, chains = nchains, iter_warmup = burnin, iter_sampling = nsamples)
  # }
  # return(fit)
  # Print results.  I need to get this as an output
  # stats <- print(fit, pars = paramsvec, probs=c((1-confid)/2,.5,1-(1-confid)/2))
  # dataout <- fit@.MISC[["summary"]][["msd"]]
  conflim_txt<-c(paste(c("Lower ",100*conf.level,"%"),collapse = ""),paste(c("Upper ",100*conf.level,"%"),collapse = ""))
  stats <- fit$summary(variables = paramsvec0)
  # dataout <- fit$draws(format = "df")
  confidbounds <- mcmc_intervals_data(fit$draws(variables = paramsvec0),prob_outer = confid)
  outputtable <- matrix(c(stats[[2]],stats[[4]],confidbounds[[5]],stats[[3]],confidbounds[[9]],stats[[8]]), nrow = length(paramsvec0), ncol = 6, byrow = FALSE,dimnames = list(paramsvec0,c("Mean","Standard Deviation",conflim_txt[1],"Median",conflim_txt[2],"R\U005E")))


  # Trace the Markov Chains for each parameter
  # plot1_MCtrace <- traceplot(fit, pars = paramsvec, inc_warmup = TRUE, nrow = 3)
  # plot1_MCtrace <- mcmc_trace(as.matrix(fit),pars=paramsvec, facet_args = list(nrow = length(paramsvec), labeller = label_parsed))
  # plot2_hist <- stan_hist(fit)
  # plot3_density <- stan_dens(fit)
  plot1_MCtrace <- mcmc_trace(fit$draws(paramsvec0))
  plot2_hist <- mcmc_hist(fit$draws(paramsvec0))
  plot3_density <- mcmc_dens(fit$draws(paramsvec0))
  plot4_densityoverlay <- mcmc_dens_overlay(fit$draws(paramsvec0))
  plot5_scatterplot <- mcmc_pairs(fit$draws(paramsvec0))

  # NEW PLOT SET:
  # 1. Posterior curves
  # 2. Data with posterior fit and confidence

  # Produce some output text that summarizes the results
  cat(c("Posterior estimates for Bayesian Analysis.\n\n"),sep = "")
  print(outputtable)
  cat(c("\n"),sep = "")


  return(list(fit,plot1_MCtrace,plot2_hist,plot3_density))
  # NOTE: Comment return above and uncomment return below if you want to view the scatterplot between posteriors.
  # This may however increase processing time.
  # return(list(fit,plot1_MCtrace,plot2_hist,plot3_density,plot4_densityoverlay,plot5_scatterplot))
}
