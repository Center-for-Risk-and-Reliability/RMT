# Crack Propagation Calculator
# Developed by Reuel Smith, (2022-2023)

crack.propagation <- function(data,dimensions,geometry,loadconditions,notcheffect = 0,units = 1,iterative = 0){
  # Computes crack propagation models and parameters and life depending on input
  # data is made up of crack length, cycle data, and propagation model.
  # (example: list(a = c(), N = c(), prop_model = "Paris").  Data will then estimate
  # crack propagation parameters of model if given (Paris is the default).  In absence of that you can state
  # crack propagation parameters if known (example: list(prop_model = "Paris", C = 10, m = -3))

  # Data: crack length (a), cycle data (N), propagation model (prop_model), model parameters, material properties (E, Sy, Kth)
  # Dimensions: W, r, initial crack a_0
  # Geometry: "edge_single_1a", "edge_double_2a", etc.
  # Loadconditions: R, P (load), DS (optional), Smax (optional)
  library(ggplot2)
  library(pracma)

  # Check units and set up axis labels
  if(units == 1){
    # For metric units, input must be: crack length (mm), stress (MPa), force (N), DK_th (MPa sqrt(m))
    stressunitslabel <- "MPa"
    crackunitslabel <- "mm"
    dadNunitslabel <- "mm/cycles"
    dKunitslabel <- "MPa sqrt(m)"
    a_f_corr <- 1000
    dadN_corr <- 1000
  }
  if(units == 2){
    # For English units, input must be: crack length (inches), stress (ksi), force (kips), DK_th (ksi sqrt(in))
    stressunitslabel <- "ksi"
    crackunitslabel <- "inches"
    dadNunitslabel <- "inch/cycles"
    dKunitslabel <- "ksi sqrt(in)"
    a_f_corr <- 1
    dadN_corr <- 1
  }

  # Establish the axes labels
  Xlab <- "Fatigue Cycles"
  Ylab <- paste(c("Fatigue Crack Length (",crackunitslabel,")"),collapse = "")
  Xlab2 <- paste(c("Stress Intensity Range (",dKunitslabel,")"),collapse = "")
  Ylab2 <- paste(c("Crack Propagation Rate (",dadNunitslabel,")"),collapse = "")

  # Check data for propagation model and if absent, set to Paris-Edrogan
  # Other options are Walker, Mechanistic, McEvily-Groeger, and NASGRO
  # Paris and Walker are RII only but need to compute RI (Threshold) and RIII
  # (Forman) as well.  The other three can be treated as all three regions
  if(length(data$prop_model) == 1){
    prop_model <- data$prop_model
  } else {
    prop_model <- "Paris"
  }

  # Establish net stress (Snet) computations based on geometry (FOR NOTCHED SPECIMENS)
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
  # Calculate or extract the Smax and/or DS
  if(length(loadconditions$DS) == 0 && length(loadconditions$R) == 1){
    Smax <- Snetcalc(dimensions,loadconditions)
    DS <- Smax*(1 - loadconditions$R)
  }
  if(length(loadconditions$DS) == 1 && length(loadconditions$R) == 1){
    # Smax <- Snetcalc(dimensions,loadconditions)
    DS <- loadconditions$DS
    Smax <- DS/(1 - loadconditions$R)
  }
  Smin <- Smax - DS

  # Predefine delta K (DK) function based on geometry
  if(geometry == "edge_single_1a" || geometry == "edge_double_2a" || geometry == "center_2a"){
    DKfunct <- function(a,DS,fg){
      DK <-fg*DS*sqrt(pi*(a/a_f_corr))
      return(DK)
    }
  }
  if(geometry == "edge_semi_circle_thick_body_1a"){
    DKfunct <- function(a,DS,fg){
      DK <-fg*2*DS*(1/pi)*sqrt(pi*(a/a_f_corr))
      return(DK)
    }
  }
  if(geometry == "corner_circle_thick_body_1a"){
    DKfunct <- function(a,DS,fg){
      DK <- (fg^2)*2*DS*(1/pi)*sqrt(pi*(a/a_f_corr))
      return(DK)
    }
  }

  # Calculate or extract a constant crack correction factor
  if(length(data$fg) == 1){
    fg <- data$fg
  }
  if(length(data$fg) == 0){
    fg <- crack.correction(dimensions,geometry)
  }
  # Compute or extract initial crack length
  if(length(dimensions$a_0) == 0 && length(data$a) == 0){
    a_0 <- 0
    a_i <- dimensions$r
  }
  if(length(dimensions$a_0) == 0 && length(data$a) > 0){
    a_0 <- min(data$a)
    a_i <- min(data$a)
  }
  if(length(dimensions$a_0) == 1){
    a_0 <- dimensions$a_0
    a_i <- dimensions$a_0 + dimensions$r
  }
  # Calculate or extract stress intensity factor threshold
  if(length(data$DKth) == 1){
    DKth <- data$DKth
  }
  if(length(data$DKth) == 0){
    # For now I will need to take DKth as based on the initial crack length if not given.  Will also need
    # to allow for different geometries here as well
    DKth <- round(DKfunct(a_i,Smax,fg),1)
  }

  # Check on whether the data is made up of crack propagation data or propagation
  # model parameters.  If the former, find the propagation da/dN and model parameters.
  # If the latter, proceed with the evaluation.
  if((length(data$a) > 0 && length(data$N) > 0 && (length(data$dadN) == 0 || length(data$dadN) > 0)) || (length(data$dadN) > 0 && length(data$DK) > 0)){
    # CASE WITH CRACK AND CYCLE DATA
    # ==============================
    # Check to see that there is more than 5 data points each
    if(length(data$a) <= 4 || length(data$N) <= 4){
      stop('Need to have at least 5 data per axis to obtain an estimation of crack propagation model parameters.')
    }
    if(length(data$a) != length(data$N)){
      stop('Crack data (a) and cycle data (N) have to be of the same length.')
    }
    if(length(data$dadN)> 0 && ((length(data$dadN) != length(data$a)-2) || (length(data$dadN) != length(data$N)-2))){
      stop('You should have two less dadN data than either the crack data or cycle data.')
    }
    # Now compute or setup da/dN and Delta K
    if(length(data$dadN) > 0){
      # When dadN is provided as data
      dadN <- data$dadN
    }
    if(length(data$dadN) == 0){
      # When dadN needs to be calculated by a and N data
      dadN <- rep(0,length(data$a) - 2)
      for(i in 1:(length(data$a) - 2)){
        params<-lm(data$a[i:(i+2)] ~ poly(data$N[i:(i+2)],2,raw=TRUE))
        C <- summary(params)$coefficients[1,1]
        B <- summary(params)$coefficients[2,1]
        A <- summary(params)$coefficients[3,1]
        dadN[i] <- sum(2*A*data$N[i+1],B)
      }
    }

    # Build Delta K vector from crack length vector
    DK <- DKfunct(data$a,Smax,fg)
    DKset <- DK[2:(length(data$a)-1)]
    dadNset <- dadN

    # Identify complete set of dadN and DK data
    DKfull <- DKset
    dadNfull <- dadNset
    if(prop_model == "Mechanistic" || prop_model == "McEvily" || prop_model == "NASGRO"){
      afull <- data$a[2:(length(data$a)-1)]
    }

    # Group dadN and DK by region (Won't need this distinction in case of exponential, NASGRO, and McClintock model)
    # Initialize by presence of RIIrange variable
    if(prop_model == "Paris" || prop_model == "Walker"){
      if(length(data$RIIrange) == 0){
        RIrange <- c(2,4)
        RIIrange <- c(2,(length(data$a)-1))
        RIIIrange <- c((length(data$a)-3),(length(data$a)-1))
      }
      if(length(data$RIIrange) == 2){
        RIrange <- c(2,min(data$RIIrange))
        RIIrange <- data$RIIrange
        RIIIrange <- c(max(data$RIIrange),(length(data$a)-1))
      }

      DK_RI <- DKset[(min(RIrange)-1):(max(RIrange)-1)]
      dadN_RI <- dadNset[(min(RIrange)-1):(max(RIrange)-1)]
      DK_RII <- DKset[(min(RIIrange)-1):(max(RIIrange)-1)]
      dadN_RII <- dadNset[(min(RIIrange)-1):(max(RIIrange)-1)]
      DK_RIII <- DKset[(min(RIIIrange)-1):(max(RIIIrange)-1)]
      dadN_RIII <- dadNset[(min(RIIIrange)-1):(max(RIIIrange)-1)]
    }

    # Estimate crack propagation parameters
    if(prop_model == "Paris"){
      params  <- lm(log(dadN_RII) ~ poly(log(DK_RII), 1, raw=TRUE))
      # Convert to m/cycles
      C <- exp(summary(params)$coefficients[1,1])/a_f_corr
      m <- summary(params)$coefficients[2,1]
      modelparams <- list(C=C,m=m)
    }
    if(prop_model == "Walker"){
      M <- matrix(c(rep(1,length(dadN_RII)),log(DK_RII),rep(-log(1 - loadconditions$R),length(dadN_RII))), nrow = length(DK_RII), ncol = 3, byrow=FALSE)
      params  <- pinv(M)%*%log(dadN_RII)
      # Convert to m/cycles
      C <- exp(params[1])/a_f_corr
      m <- params[2]
      gam <- 1 - (params[3]/params[2])
      modelparams <- list(C=C,m=m,gam=gam)
    }
    if(prop_model == "Paris" || prop_model == "Walker"){
      # Region I
      # params  <- lm(log(dadN_RI) ~ poly(log(DK_RI - DKth), 1, raw=TRUE))
      params  <- lm(log(c(dadN_RI,dadN_RII)) ~ poly(log(c(DK_RI,DK_RII) - DKth), 1, raw=TRUE))
      # Convert to m/cycles
      A <- exp(summary(params)$coefficients[1,1])/a_f_corr
      p <- summary(params)$coefficients[2,1]
      modelparamsI <- list(A=A,p=p)

      # Region III
      params  <- lm((log(c(dadN_RII,dadN_RIII)) + log(((1 - loadconditions$R)*data$Kic) - c(DK_RII,DK_RIII))) ~ poly(log(c(DK_RII,DK_RIII)), 1, raw=TRUE))
      # Convert to m/cycles
      CIII <- exp(summary(params)$coefficients[1,1])/a_f_corr
      mIII <- summary(params)$coefficients[2,1]
      modelparamsIII <- list(CIII=CIII,mIII=mIII)
    }

    # Regions I through III
    if(prop_model == "Mechanistic"){
      m <- exp(mean(log(dadNfull/a_f_corr)) - mean(log(afull/a_f_corr)))
      A_0 <- exp(mean(log(data$a/a_f_corr)) - m*mean(data$N))
      modelparams <- list(A_0=A_0,m=m)
    }
    if(prop_model == "McEvily"){
      Kmax <- DKfunct(afull,Smax,fg)
      A <- exp(mean(log(dadNfull/a_f_corr)) - mean(2*log(DKfull - DKth) + log(1 + (DKfull/(data$Kic - Kmax)))))
      modelparams <- list(A=A)
    }
    if(prop_model == "NASGRO"){
      Kmax <- DKfunct(afull,Smax,fg)
      dadNvals<-log(dadNfull/a_f_corr)
      Avals<-matrix(c(rep(1,2*length(dadNfull)),log(DKfull)-log(1 - loadconditions$R), log(1 - (DKth/DKfull)), -log(1 - (Kmax/data$Kic))),nrow = length(dadNfull), ncol = 5, byrow=FALSE)
      params  <- pinv(Avals)%*%dadNvals
      C <- exp(params[1])
      f <- 1 - exp(params[2]/params[3])
      n <- params[3]
      p <- params[4]
      q <- params[5]
      modelparams <- list(C=C,f=f,n=n,p=p,q=q)
    }
  }

  # return(modelparams)

  if(length(data$a) == 0 && length(data$N) == 0){
    # CASE WITHOUT CRACK AND CYCLE DATA
    # =================================
    if(prop_model == "Paris" && length(data$C) == 1 && length(data$m) == 1){
      C <- data$C
      m <- data$m
      modelparams <- list(C=C,m=m)
    }
    if(prop_model == "Paris" && (length(data$C) == 0  || length(data$m) == 0)){
      stop('Enter Paris parameters C and m.')
    }
    if(prop_model == "Walker" && length(data$C) == 1 && length(data$m) == 1 && length(data$gam) == 1){
      C <- data$C
      m <- data$m
      gam <- data$gam
      modelparams <- list(C=C,m=m,gam=gam)
    }
    if(prop_model == "Walker" && (length(data$C) == 0 || length(data$m) == 0 || length(data$gam) == 0)){
      stop('Enter Walker parameters C, m, and gamma.')
    }
    # Regions I through III
    if(prop_model == "Mechanistic" && length(data$m) == 1 && length(data$A_0) == 1){
      m <- data$m
      A_0 <- data$A_0
      modelparams <- list(A_0=A_0,m=m)
    }
    if(prop_model == "McEvily" && length(data$A) == 1){
      A <- data$A
      modelparams <- list(A=A)
    }
    if(prop_model == "NASGRO" && length(data$C) == 1 && length(data$f) == 1 && length(data$n) == 1 && length(data$p) == 1 && length(data$q) == 1){
      C <- data$C
      f <- data$f
      n <- data$n
      p <- data$p
      q <- data$q
      modelparams <- list(C=C,f=f,n=n,p=p,q=q)
    }

    # Establish DK set based on DKth and Kic
    DKfull <- linspace(DKth,data$Kic,100)
    if(prop_model == "Paris" || prop_model == "Walker"){
      DK_RI <- DKfull[1:20]
      DK_RII <- DKfull[20:80]
      DK_RIII <- DKfull[80:100]
    }
  }

  # Establish Crack Growth Equations
  # Region II
  if(prop_model == "Paris"){
    # Define delta N (DN) function based on geometry
    if(geometry == "edge_single_1a" || geometry == "edge_double_2a" || geometry == "center_2a"){
      dN <- function(a1,a2,fg,delS){
        # a1 = initial, a2 = final
        a1 <- a1/a_f_corr
        a2 <- a2/a_f_corr
        DN <- ((C*(fg^m)*(delS^m)*(pi^(m/2))*(0.5*m - 1))^-1)*((a1^(1 - 0.5*m)) - (a2^(1 - 0.5*m)))
        return(DN)
      }
    }
    if(geometry == "edge_semi_circle_thick_body_1a"){
      dN <- function(a1,a2,fg,delS){
        # a1 = initial, a2 = final
        a1 <- a1/a_f_corr
        a2 <- a2/a_f_corr
        DN <- ((C*(2^m)*(fg^m)*(delS^m)*(pi^(-m/2))*(0.5*m - 1))^-1)*((a1^(1 - 0.5*m)) - (a2^(1 - 0.5*m)))
        return(DN)
      }
    }
    if(geometry == "corner_circle_thick_body_1a"){
      dN <- function(a1,a2,fg,delS){
        # a1 = initial, a2 = final
        a1 <- a1/a_f_corr
        a2 <- a2/a_f_corr
        DN <- ((C*(2^m)*(fg^(2*m))*(delS^m)*(pi^(-m/2))*(0.5*m - 1))^-1)*((a1^(1 - 0.5*m)) - (a2^(1 - 0.5*m)))
        return(DN)
      }
    }
    dNalt <- function(a1,a2,delS,D,r){
      # a1 = initial, a2 = final
      a1 <- a1/a_f_corr
      a2 <- a2/a_f_corr
      DN <- ((C*(delS^m)*(pi^(m/2))*((1+7.69*sqrt(D*r))^(m/2))*(0.5*m - 1))^-1)*((a1^(1 - 0.5*m)) - (a2^(1 - 0.5*m)))
      return(DN)
    }
    dNalt2 <- function(a1,a2,fg,delS,K_t){
      # a1 = initial, a2 = final
      a1 <- a1/a_f_corr
      a2 <- a2/a_f_corr
      DN <- ((C*(fg^m)*(delS^m)*(K_t^m)*(pi^(m/2))*(0.5*m - 1))^-1)*((a1^(1 - 0.5*m)) - (a2^(1 - 0.5*m)))
      return(DN)
    }
    dadN <- function(DK){
      dadN <- dadN_corr*C*((DK)^m)
      return(dadN)
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
    dNalt <- function(a1,a2,delS,D,r){
      # a1 = initial, a2 = final
      a1 <- a1/a_f_corr
      a2 <- a2/a_f_corr
      DN <- ((1 - loadconditions$R)^(m*(1-gam)))*((C*(delS^m)*(pi^(m/2))*((1+7.69*sqrt(D*r))^(m/2))*(0.5*m - 1))^-1)*((a1^(1 - 0.5*m)) - (a2^(1 - 0.5*m)))
      return(DN)
    }
    dNalt2 <- function(a1,a2,fg,delS,K_t){
      # a1 = initial, a2 = final
      a1 <- a1/a_f_corr
      a2 <- a2/a_f_corr
      DN <- ((1 - loadconditions$R)^(m*(1-gam)))*((C*(fg^m)*(delS^m)*(K_t^m)*(pi^(m/2))*(0.5*m - 1))^-1)*((a1^(1 - 0.5*m)) - (a2^(1 - 0.5*m)))
      return(DN)
    }
    dadN <- function(DK){
      dadN <- dadN_corr*C*((DK/((1 - loadconditions$R)^(1 - gam)))^m)
      return(dadN)
    }
  }

  # Setup RI and RIII parameters if just RII parameters are given (in m/cycle)
  if(length(data$a) == 0 && length(data$N) == 0){
    dadN_RII <- dadN(DK_RII)
    # Region I
    params  <- lm(log(dadN_RII) ~ poly(log(DK_RII - DKth), 1, raw=TRUE))
    if(units==1){
      A <- exp(summary(params)$coefficients[1,1])/1000
    } else{
      A <- exp(summary(params)$coefficients[1,1])
    }
    p <- summary(params)$coefficients[2,1]
    modelparamsI <- list(A=A,p=p)

    # Region III
    params  <- lm((log(dadN_RII) + log(((1 - loadconditions$R)*data$Kic) - DK_RII)) ~ poly(log(DK_RII), 1, raw=TRUE))
    if(units==1){
      CIII <- exp(summary(params)$coefficients[1,1])/1000
    } else{
      CIII <- exp(summary(params)$coefficients[1,1])
    }
    mIII <- summary(params)$coefficients[2,1]
    modelparamsIII <- list(CIII=CIII,mIII=mIII)
  }

  if(prop_model == "Paris"|| prop_model == "Walker"){
    # Region I
    dadNI <- function(DK){
      dadN <- dadN_corr*A*((DK - DKth)^p)
      return(dadN)
    }
    # Region III Default Forman Equation
    dadNIII <- function(DK){
      dadN <- (dadN_corr*CIII*((DK)^mIII))/((1 - loadconditions$R)*data$Kic - DK)
      return(dadN)
    }
  }

  # Regions I through III
  if(prop_model == "Mechanistic"){
    dN <- function(a1,a2,fg,delS){
      # a1 = initial, a2 = final
      a1 <- a1/a_f_corr
      a2 <- a2/a_f_corr
      DN <- (log(a2) - log(a1))/m
      return(DN)
    }
    dadN <- function(a){
      dadN <- dadN_corr*m*a
      return(dadN)
    }
  }
  if(prop_model == "McEvily"){
    dN <- function(a1,a2,fg,delS){
      # a1 = initial, a2 = final
      a1 <- a1/a_f_corr
      a2 <- a2/a_f_corr
      DN <- ((1 - loadconditions$R)^(m*(1-gam)))*((C*(fg^m)*(delS^m)*(pi^(m/2))*(0.5*m - 1))^-1)*((a1^(1 - 0.5*m)) - (a2^(1 - 0.5*m)))
      return(DN)
    }
    dadN <- function(DK){
      dadN <- dadN_corr*A*((DK - DKth)^2)*(1 + (DK/(data$Kic - Kmax)))
      return(dadN)
    }
    Kmax <- DKfunct(afull,Smax,fg)
    A <- exp(mean(log(dadNfull/a_f_corr)) - mean(2*log(DKfull - DKth) + log(1 + (DKfull/(data$Kic - Kmax)))))
    modelparams <- list(A=A)
  }
  if(prop_model == "NASGRO"){
    Kmax <- DKfunct(afull,Smax,fg)
    dadNvals<-log(dadNfull/a_f_corr)
    Avals<-matrix(c(rep(1,2*length(dadNfull)),log(DKfull)-log(1 - loadconditions$R), log(1 - (DKth/DKfull)), -log(1 - (Kmax/data$Kic))),nrow = length(dadNfull), ncol = 5, byrow=FALSE)
    params  <- pinv(Avals)%*%dadNvals
    C <- exp(params[1])
    f <- 1 - exp(params[2]/params[3])
    n <- params[3]
    p <- params[4]
    q <- params[5]
    modelparams <- list(C=C,f=f,n=n,p=p,q=q)
  }
  # Compute or extract initial crack length
  if(length(dimensions$a_0) == 0){
    a_0 <- 0
    a_i <- dimensions$r
  }
  if(length(dimensions$a_0) == 1){
    a_0 <- dimensions$a_0
    a_i <- dimensions$a_0 + dimensions$r
  }

  # ================================================
  # Determine final crack length based on if this is notched (notcheffect==1) or just a crack (notcheffect==0)
  if(notcheffect == 0){
    # CASE where the sample is not-notched but simply a crack of varying shape (center, edge, circular, elliptical)
    # Compute final crack length a_f
    # Predefine delta K (DK) function based on geometry
    if(geometry == "edge_single_1a" || geometry == "edge_double_2a" || geometry == "center_2a"){
      a_f <- (1/pi)*((data$Kic/(fg*Smax))^2)*a_f_corr
    }
    if(geometry == "edge_semi_circle_thick_body_1a"){
      a_f <- (1/pi)*(((pi*data$Kic)/(2*fg*Smax))^2)*a_f_corr
    }
    if(geometry == "corner_circle_thick_body_1a"){
      a_f <- (1/pi)*(((pi*data$Kic)/(2*(fg^2)*Smax))^2)*a_f_corr
    }

    # Check to see if an a_f has already been entered and if it is less than the pre-calculated a_f
    if(length(dimensions$a_f) > 0 && (dimensions$a_f + dimensions$r) < a_f){
      a_f <- dimensions$a_f + dimensions$r
    }
    a_curve <- linspace(a_i,a_f,100)
  }
  if(notcheffect == 1){
    # CASE where the sample is notched
    # Use when transition crack length is less than initial crack
    if(geometry == "edge_single_1a"){
      Kt <- stress.concentration.factor(dimensions,"rect_1semicirc_edge")
    }
    if(geometry == "edge_double_2a" || geometry == "center_2a"){
      Kt <- stress.concentration.factor(dimensions,"rect_2semicirc_edge")
    }

    # Calculate lt as a check to see if the evaluation needs to bee adjusted
    # lt <- dimensions$r/((((1.12*Kt)/fg)^2) - 1)
    lt <- 0.13*sqrt(dimensions$r*dimensions$r)

    # Compute final crack length a_f
    if(lt < a_0){
      a_f <- (1/pi)*((data$Kic/(fg*Smax))^2)*a_f_corr
    }
    if(lt > a_0){
      a_f <- (1/pi)*((data$Kic/(fg*Smax*Kt))^2)*a_f_corr
    }

    if(length(data$Sy) == 1){
      a_f <- min(c(a_f,a_f1))
    }
    # Check to see if an a_f has already been entered and if it is less than the precalculated a_f
    if(length(dimensions$a_f) > 0 && (dimensions$a_f + dimensions$r) < a_f){
      a_f <- dimensions$a_f + dimensions$r
    }

    if(lt < a_0){
      a_curve <- linspace(a_i,a_f,100)
    }
    if(lt > a_0){
      a_curve <- linspace(a_0,a_f,100)
    }
  }

  N_curve <- rep(0,100)

  # Perform Non-Iterative Crack Propagation analysis (iterative = 0)
  # TEMPORARY SETTING OF ITERATIVE OPTION TO NON-ITERATIVE
  if(iterative == 0){
    # Non-iterative crack propagation calculation
    if(Smin < 0){
      # Just crack
      if(notcheffect == 0){
        Nf <- dN(a_i,a_f,fg,Smax)
        for(i in 1:100){
          N_curve[i] <- dN(a_i,a_curve[i],fg,Smax)
        }
        DK_curve <- linspace(DKfunct(a_i,Smax,fg),DKfunct(a_f,Smax,fg),100)
      }
      # Crack and notch
      if(notcheffect == 1){
        if(lt < a_0){
          Nf <- dN(a_i,a_f,fg,Smax)
          for(i in 1:100){
            N_curve[i] <- dN(a_i,a_curve[i],fg,Smax)
          }
          DK_curve <- linspace(DKfunct(a_i,Smax,fg),DKfunct(a_f,Smax,fg),100)
        }
        if(lt > a_0){
          Nf <- dNalt2(a_0,a_f,fg,Smax,Kt)
          for(i in 1:100){
            N_curve[i] <- dNalt2(a_0,a_curve[i],fg,Smax,Kt)
          }
          DK_curve <- linspace(DKfunct(a_0,Smax,fg),data$Kic,100)
        }
      }
    }
    if(Smin >= 0){
      # Just crack
      if(notcheffect == 0){
        Nf <- dN(a_i,a_f,fg,DS)
        for(i in 1:100){
          N_curve[i] <- dN(a_i,a_curve[i],fg,DS)
        }
        DK_curve <- linspace(DKfunct(a_i,DS,fg),DKfunct(a_f,DS,fg),100)
      }
      # Crack and notch
      if(notcheffect == 1){
        if(lt < a_0){
          Nf <- dN(a_i,a_f,fg,DS)
          for(i in 1:100){
            N_curve[i] <- dN(a_i,a_curve[i],fg,DS)
          }
          DK_curve <- linspace(DKfunct(a_i,DS,fg),DKfunct(a_f,DS,fg),100)
        }
        if(lt > a_0){
          Nf <- dNalt2(a_0,a_f,fg,DS,Kt)
          for(i in 1:100){
            N_curve[i] <- dNalt2(a_0,a_curve[i],fg,DS,Kt)
          }
          DK_curve <- linspace(DKfunct(a_0,DS,fg),data$Kic,100)
        }
      }
    }

    if(prop_model == "Paris" || prop_model == "Walker"){
      # Region I
      DKI_curve <- linspace(DKth,min(DK_RII),51)
      DKI_curve <- DKI_curve[2:51]
      dadNI_curve <- dadNI(DKI_curve)
      # Region II
      DKII_curve <- linspace(min(DK_RII),max(DK_RII),100)
      dadNII_curve <- dadN(DKII_curve)
      # Region III
      DKIII_curve <- linspace(max(DK_RII),data$Kic,51)
      DKIII_curve <- DKIII_curve[1:50]
      dadNIII_curve <- dadNIII(DKIII_curve)
    }
  }
  # Perform Iterative Crack Propagation analysis (iterative == 1)
  # Base this from the previous calculation if you want to iterate
  fg_v <- rep(0,100)
  dadN_avg <- rep(0,99)
  if(iterative == 1){
    # Take a_curve and iterate further
    if(length(data$fg) == 1){
      fg_v <- rep(data$fg,100)
    }
    if(lt < a_0){
      dimensions_alt <- dimensions
      for(i in 1:100){
        dimensions_alt$a_0 <- a_curve[i] - dimensions$r
        if(length(data$fg) == 0){
          fg_v[i] <- crack.correction(dimensions_alt,geometry)
        }
      }
      DK_v <- fg_v*Smax*sqrt(pi*(a_curve/a_f_corr))
      dadN_v <- dadN(DK_v)
      for(i in 1:99){
        dadN_avg[i] <- 0.5*sum(dadN_v[i:(i+1)])
      }
      dN <- (((a_f - a_i)/99)/a_f_corr)/dadN_avg
      N_curve <- cumsum(c(0,dN))
    }
    DK_curve <- DKfunct(min(a_curve),Smax,fg_v)
    dadN_curve <- dadN(DK_curve)
    dadNIII_curve <- dadNIII(DK_curve)
    Nf <- N_curve[100]
    # maxDK <-max(DK_v)
  }


  # Plot #1: cycles vs. crack length
  df <- data.frame(a = a_curve, N = N_curve)
  plotout1<-ggplot() +
    geom_path(data=df, aes(N, a), colour = 'blue', size = 0.9, linetype = "dashed") +
    labs(x=Xlab,y=Ylab)
  if(length(data$a) > 0 && length(data$N) > 0){
    df1 <- data.frame(adata = data$a, Ndata = data$N, dataset = rep("data",length(data$N)))
    plotout1 <- plotout1 + geom_point(data=df1, aes(Ndata, adata, colour=dataset), size = 1.9)
  }
  # Plot #2: da/dN and cycles vs. crack length (for now only for input with crack, cycle, or growth data)
  df2 <- data.frame(DK = c(DKI_curve,DKII_curve,DKIII_curve), dadN = c(dadNI_curve,dadNII_curve,dadNIII_curve), Region = c(rep("Region I",50),rep("Region II",100),rep("Region III",50)))

  plotout2<-ggplot() +
    geom_path(data=df2, aes(DK, dadN, colour = Region), size = 0.9, linetype = "dashed") +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    annotation_logticks() +
    labs(x=Xlab2,y=Ylab2)

  if((length(data$a) > 0 && length(data$N) > 0) || (length(data$dadN) > 0 && length(data$DK) > 0)){
    df3 <- data.frame(DKdata = DKfull, dadNdata = dadNfull, dataset = rep("data",length(dadNfull)))
    plotout2 <- plotout2 + geom_point(data=df3, aes(DKdata, dadNdata, colour=dataset), size = 1.9)
  }


  # return(list(crackpropplot = plotout1,growthratecurve = plotout2,fg, a_i, a_f, modelparams,Nf,Smax,DK_RI,DK_RII,DK_RIII,dadN_RI,dadN_RII,dadN_RIII))
  # return(list(crackpropplot = plotout1,fg, a_i, a_f, modelparams,Smax,Nf,DKth,DKset,dadNset))
  return(list(crackpropplot = plotout1,growthratecurve = plotout2,fg, Smax, Smin, a_i, a_f, modelparams,Nf,DKth,DK_RII,dadN_RII))
}
