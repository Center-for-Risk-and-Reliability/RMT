# Crack Propagation Calculator
# Developed by Reuel Smith, 2022

crack.propagation <- function(data,dimensions,geometry,loadconditions,units = 1,iterative = 0){
  # Computes crack propagation models and parameters and life depending on input
  # data is made up of crack length, cycle data, and propagation model.
  # (example: list(a = c(), N = c(), prop_model = "Paris").  Data will then estimate
  # crack propagation parameters of model if given (Paris is the default).  In absence of that you can state
  # crack propagation parameters if known (example: list(prop_model = "Paris", C = 10, m = -3))

  # Data: crack length (a), cycle data (N), propagation model (prop_model), model parameters, material properties (E, Sy, Kth)
  # Dimensions: W, r, initial crack a_0
  # Geometry: "edge_single_1a", "edge_double_2a", etc.
  # Loadconditions: R, P (load), DS (optional), Smax (optional)

  # Check units and set up axis labels
  if(units == 1){
    # For metric units, input must be: crack length (mm), stress (MPa), force (N), DK_th (MPa sqrt(m))
    stressunitslabel <- "MPa"
    crackunitslabel <- "mm"
    dadNunitslabel <- "mm/cycles"
    dKunitslabel <- "MPa sqrt(m)"
    a_f_corr <- 1000
  }
  if(units == 2){
    # For English units, input must be: crack length (inches), stress (ksi), force (kips), DK_th (ksi sqrt(in))
    stressunitslabel <- "ksi"
    crackunitslabel <- "inches"
    dadNunitslabel <- "inch/cycles"
    dKunitslabel <- "ksi sqrt(in)"
    a_f_corr <- 1
  }

  # Check data for propagation model and if absent, set to Paris-Edrogan
  if(length(data$prop_model) == 1){
    prop_model <- data$prop_model
  } else {
    prop_model <- "Paris"
  }

  if(length(data$prop_modelIII) == 1){
    prop_modelIII <- data$prop_modelIII
  } else {
    prop_modelIII <- "Forman"
  }

  # Establish net stress (Snet) computations based on geometry
  if(geometry == "edge_single_1a"){
    Snetcalc <- function(dimensions,loadconditions){
      Snet <- loadconditions$P/(dimensions$t*(dimensions$W - dimensions$r))
      return(Snet)
    }
    if(length(data$Sy) == 1){
      a_f1 <- (dimensions$W - loadconditions$P/(dimensions$t*data$Sy)) + dimensions$r
    }
  }
  if(geometry == "edge_double_2a" || geometry == "center_2a"){
    Snetcalc <- function(dimensions,loadconditions){
      Snet <- loadconditions$P/(dimensions$t*(dimensions$W - 2*dimensions$r))
      return(Snet)
    }
    if(length(data$Sy) == 1){
      a_f1 <- 0.5*(dimensions$W - loadconditions$P/(dimensions$t*data$Sy)) + dimensions$r
    }
  }

  # Check on whether the data is made up of crack propagation data or propagation
  # model parameters.  If the former, find the propagation da/dN and model parameters.
  # If the latter, proceed with the evaluation.
  if(length(data$a) > 0 && length(data$N) > 0){
    # CASE WITH CRACK AND CYCLE DATA
    # ==============================
    # Check to see that there is more than 5 data points each
    if(length(data$a) <= 4 || length(data$N) <= 4){
      stop('Need to have at least 5 data per axis to obtain an estimation of crack propagation model parameters.')
    }
    if(length(data$a) != length(data$N)){
      stop('Crack data (a) and cycle data (N) have to be of the same length.')
    }
    # Now compute da/dN and Delta K
    dadN <- rep(0,length(data$a) - 2)
    for(i in 1:(length(data$a) - 2)){
      params<-lm(data$a[i:(i+2)] ~ poly(data$N[i:(i+2)],2,raw=TRUE))
      C <- summary(params)$coefficients[1,1]
      B <- summary(params)$coefficients[2,1]
      A <- summary(params)$coefficients[3,1]
      dadN[i] <- sum(2*A*data$N[i+1],B)
    }
    if(length(data$fg) == 1){
      fg <- data$fg
    }
    if(length(data$fg) == 0){
      fg <- crack.correction(dimensions,geometry)
    }

    if(length(loadconditions$DS) == 0 && length(loadconditions$R) == 1){
      Smax <- Snetcalc(dimensions,loadconditions)
      DS <- Smax*(1 - loadconditions$R)
    }
    if(length(loadconditions$DS) == 1 && length(loadconditions$R) == 1){
      Smax <- Snetcalc(dimensions,loadconditions)
      DS <- loadconditions$DS
      Smax <- DS/(1 - loadconditions$R)
    }
    DK <- fg*DS*sqrt(pi*(data$a/a_f_corr))

    DKset <- DK[2:(length(data$a)-1)]
    dadNset <- dadN
    params  <- lm(log(dadNset) ~ poly(log(DKset), 1, raw=TRUE))
    R2set <- summary(params)$r.squared
    R2old <- 0
    R2new <- R2set

    # Now check to see which set of data is the most linear
    ilow <- 1
    ihigh <- length(dadN)

    for(i in 1:length(dadN)){
      # Cut lowest number and compute R2set
      R2new1 <- summary(lm(log(dadNset[(ilow+1):ihigh]) ~ poly(log(DKset[(ilow+1):ihigh]), 1, raw=TRUE)))$r.squared
      # Cut highest number and compute R2set
      R2new2 <- summary(lm(log(dadNset[ilow:(ihigh-1)]) ~ poly(log(DKset[ilow:(ihigh-1)]), 1, raw=TRUE)))$r.squared
      # Set new bounds
      if(R2new1 > R2new && R2new2 <= R2new){
        ilow <- ilow + 1
        DKset <- DKset[ilow:ihigh]
        dadNset <- dadNset[ilow:ihigh]
        R2old <- R2new
        R2new <- R2new1
      }
      if(R2new1 <= R2new && R2new2 > R2new){
        ihigh <- ihigh - 1
        DKset <- DKset[ilow:ihigh]
        dadNset <- dadNset[ilow:ihigh]
        R2old <- R2new
        R2new <- R2new2
      }
      if(R2new1 <= R2new && R2new2 <= R2new){
        break
      }
    }
    # Group dadN and DK by region (Won't need this distinction in case of exponential, NASGRO, and McClintock model)
    DK_RI <- DK[2:(length(data$a)-1)][1:ilow]
    dadN_RI <- dadN[1:ilow]
    DK_RII <- DKset
    dadN_RII <- dadNset
    DK_RIII <- DK[2:(length(data$a)-1)][ihigh:length(dadN)]
    dadN_RIII <- dadN[ihigh:length(dadN)]

    # Estimate crack propagation parameters
    if(prop_model == "Paris"){
      params  <- lm(log(dadN_RII) ~ poly(log(DK_RII), 1, raw=TRUE))
      C <- exp(summary(params)$coefficients[1,1])
      m <- summary(params)$coefficients[2,1]
      modelparams <- list(C=C,m=m)
    }
    if(prop_model == "Walker"){
      M <- matrix(c(rep(1,length(dadN_RII)),log(DK_RII),rep(-log(1 - loadconditions$R))), nrow = length(DK_RII), ncol = 3, byrow=FALSE)
      params  <- pinv(M)%*%log(dadN_RII)
      C <- exp(params[1])
      m <- params[2]
      gam <- 1 - (params[3]/params[2])
      modelparams <- list(C=C,m=m,gam=gam)
    }
  }

  if(length(data$a) == 0 && length(data$N) == 0){
    # CASE WITHOUT CRACK AND CYCLE DATA
    # =================================
    if(prop_model == "Paris" && length(data$C) == 1 && length(data$m) == 1){
      C <- data$C
      m <- data$m
      modelparams <- list(C=C,m=m)
    }
    if(prop_model == "Paris" && length(data$C) == 0 && length(data$m) == 0){
      stop('Enter Paris parameters C and m.')
    }
    if(prop_model == "Walker" && length(data$C) == 1 && length(data$m) == 1 && length(data$gam) == 1){
      C <- data$C
      m <- data$m
      gam <- data$gam
      modelparams <- list(C=C,m=m,gam=gam)
    }
    if(prop_model == "Walker" && length(data$C) == 0 && length(data$m) == 0 && length(data$gam) == 0){
      stop('Enter Walker parameters C, m, and gamma.')
    }
  }

  # Establish Crack Growth Equations
  if(prop_model == "Paris"){
    dN <- function(a1,a2,fg,delS){
      # a1 = initial, a2 = final
      a1 <- a1/a_f_corr
      a2 <- a2/a_f_corr
      DN <- ((C*(fg^m)*(delS^m)*(pi^(m/2))*(0.5*m - 1))^-1)*((a1^(1 - 0.5*m)) - (a2^(1 - 0.5*m)))
      return(DN)
    }
  }
  if(prop_model == "Walker"){
    dN <- function(a1,a2,fg,delS){
      # a1 = initial, a2 = final
      a1 <- a1/a_f_corr
      a2 <- a2/a_f_corr
      DN <- ((1 - loadconditions$R)^(m*(1-gam)))*((C*(fg^m)*(delS^m)*(pi^(m/2))*(0.5*m - 1))^-1)*((a1^(1 - 0.5*m)) - (a2^(1 - 0.5*m)))
      return(DN)
    }
  }

  # Perform Non-Iterative Crack Propagation analysis (iterative = 0)
  # TEMPORRARY SETTING OF ITERATIVE OPTION TO NON-ITERATIVE
  if(iterative == 0 || iterative == 1){
    if(length(data$fg) == 1){
      fg <- data$fg
    }
    if(length(data$fg) == 0){
      fg <- crack.correction(dimensions,geometry)
    }

    if(length(loadconditions$DS) == 0 && length(loadconditions$R) == 1){
      Smax <- Snetcalc(dimensions,loadconditions)
      DS <- Smax*(1 - loadconditions$R)
    }
    if(length(loadconditions$DS) == 1 && length(loadconditions$R) == 1){
      Smax <- Snetcalc(dimensions,loadconditions)
      DS <- loadconditions$DS
      Smax <- DS/(1 - loadconditions$R)
    }

    # Compute final crack length a_f.  Will go with Kic based evaluation to start
    # but if Sy is an input (yield stress) we need to check that so see if it is
    # less than the other evaluation or not
    a_f <- (1/pi)*((data$Kic/(fg*Smax))^2)*a_f_corr
    if(length(data$Sy) == 1){
      a_f <- min(c(a_f,a_f1))
    }
    if(length(dimensions$a_0) == 0){
      a_i <- dimensions$r
    }
    if(length(dimensions$a_0) == 1){
      a_i <- dimensions$a_0 + dimensions$r
    }
    # return(list(a_i,C))
    Nf <- dN(a_i,a_f,fg,Smax)
  }
  # Perform Iterative Crack Propagation analysis (iterative = 1)


  # Plot will be da/dN and cycles vs. crack length (against data if available)
  return(list(fg, a_i, a_f, modelparams,Nf))
}
