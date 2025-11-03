# Bayesian Life-Stress Estimator
# Developed by Dr. Reuel Smith, 2021-2023

lifestress.BAYESest <- function(pt_est,ls,dist,TTF,SF,Tc=NULL,Sc=NULL,SUSE=NULL,SACC=NULL,confid=0.95,priors,nsamples=20000,burnin=10000,nchains=4,param2=NULL){
  #Load pracma library for erf
  library(pracma)
  library(StanHeaders)
  library(rstan)
  library(ggplot2)
  library(shinystan)
  library(cmdstanr)
  library(bayesplot)

  # (UPDATE 11/14/2023): Adding an IPL form of l = 1/(b x S^a) and Exponential form b x exp(a/S)

  # Add input to this to include prior estimates for LS parameters.
  # Example: priors<-c("normal(3,4)","normal(1,4)", "lognormal(-2,3)")
  # I will have to cite the Rstan text for distributions in the code.  Use lookup("") for the translation.
  # The code takes these and separates them so that they are written into the stan file.
  # Then the code will run the program and compute the Bayes estimation

  # Check to see if dist="Exponential" so you can exclude life
  # distribution parameters.
  if (dist=="Exponential"  || (dist=="Weibull" && is.null(param2) == FALSE)) {
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

  # Initialize life-stress parameter estimates for theta
  if (ls=="Linear") {
    # theta[1] - parameter a, theta[2] - parameter b
    lsparams <- "real a; real b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "b + Sf*a"
    loglifeF <- "log(b + Sf*a)"
    if(is.null(Tc)==FALSE){
      lifeC <- "b + Sc*a"
      loglifeC <- "log(b + Sc*a)"
    }
    if(missing(SUSE)==FALSE){
      lifeU <- "b + Suse*a;"
      complifeU <- pt_est[ishift+2] + SUSE*pt_est[ishift+1]
    }
    if(missing(SUSE)==FALSE && missing(SACC)==FALSE){
      # AFheading <- paste(c("AFat",SACC),collapse = "")
      AFheading <- paste(c("AFatSACC"),collapse = "")
      AF <- "(b + Suse*a)/(b + Sacc*a);"
      complifeU <- pt_est[ishift+2] + SUSE*pt_est[ishift+1]
      compAF <- (pt_est[ishift+2] + SUSE*pt_est[ishift+1])/(pt_est[ishift+2] + SACC*pt_est[ishift+1])
    }
  }

  if (ls=="Exponential"){
    # theta[1] - parameter a, theta[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "b*exp(a*Sf)"
    loglifeF <- "log(b) + a*Sf"
    if(is.null(Tc)==FALSE){
      lifeC <- "b*exp(a*Sc)"
      loglifeC <- "log(b) + a*Sc"
    }
    if(missing(SUSE)==FALSE){
      lifeU <- "b*exp(a*Suse);"
      complifeU <- pt_est[ishift+2]*exp(SUSE*pt_est[ishift+1])
    }
    if(missing(SUSE)==FALSE && missing(SACC)==FALSE){
      # AFheading <- paste(c("AFat",SACC),collapse = "")
      AFheading <- paste(c("AFatSACC"),collapse = "")
      AF <- "exp(a*(Suse - Sacc));"
      complifeU <- pt_est[ishift+2]*exp(SUSE*pt_est[ishift+1])
      compAF <- exp((SUSE - SACC)*pt_est[ishift+1])
    }
  }

  if (ls=="Exponential2"){
    # theta[1] - parameter a, theta[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "b*exp(a/Sf)"
    loglifeF <- "log(b) + a/Sf"
    if(is.null(Tc)==FALSE){
      lifeC <- "b*exp(a/Sc)"
      loglifeC <- "log(b) + a/Sc"
    }
    if(missing(SUSE)==FALSE){
      lifeU <- "b*exp(a/Suse);"
      complifeU <- pt_est[ishift+2]*exp(pt_est[ishift+1]/SUSE)
    }
    if(missing(SUSE)==FALSE && missing(SACC)==FALSE){
      # AFheading <- paste(c("AFat",SACC),collapse = "")
      AFheading <- paste(c("AFatSACC"),collapse = "")
      AF <- "exp(a*((1/Suse) - (1/Sacc)));"
      complifeU <- pt_est[ishift+2]*exp(pt_est[ishift+1]/SUSE)
      compAF <- exp(((1/SUSE) - (1/SACC))*pt_est[ishift+1])
    }
  }

  if (ls=="Arrhenius") {
    # lsparams[1] - parameter Ea, lsparams[2] - parameter b
    # Temperature HAS to be in Kelvin for this to work
    lsparams <- "real Ea; real<lower=0> b;"
    lsparamsvec <- c("Ea","b")
    pr1<-paste(c("Ea ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "b*exp(Ea/((8.617385e-5)*Sf))"
    loglifeF <- "log(b) + (Ea/((8.617385e-5)*Sf))"
    if(is.null(Tc)==FALSE){
      lifeC <- "b*exp(Ea/((8.617385e-5)*Sc))"
      loglifeC <- "log(b) + (Ea/((8.617385e-5)*Sc))"
    }
    if(missing(SUSE)==FALSE){
      lifeU <- "b*exp(Ea/((8.617385e-5)*Suse));"
      complifeU <- pt_est[ishift+2]*exp(pt_est[ishift+1]/(8.617385e-5*SUSE))
    }
    if(missing(SUSE)==FALSE && missing(SACC)==FALSE){
      # AFheading <- paste(c("AFat",SACC),collapse = "")
      AFheading <- paste(c("AFatSACC"),collapse = "")
      AF <- "exp((Ea/8.617385e-5)*((1/Suse) - (1/Sacc)));"
      complifeU <- pt_est[ishift+2]*exp(pt_est[ishift+1]/(8.617385e-5*SUSE))
      compAF <- exp(((1/SUSE) - (1/SACC))*(pt_est[ishift+1]/8.617385e-5))
    }
  }

  if (ls=="Eyring") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    # NOTE (8/25/2025): Can now go back to all Eyring type functions and just set .*, ./, and .^ on vector
    # multipliers (like MATLAB)
    lifeF <- "(b/Sf).*exp(a/Sf)"
    loglifeF <- "log(b) - log(Sf) + (a/Sf)"
    if(is.null(Tc)==FALSE){
      lifeC <- "(b/Sc).*exp(a/Sc)"
      loglifeC <- "log(b) - log(Sc) + (a/Sc)"
    }
    if(missing(SUSE)==FALSE){
      lifeU <- "(b/Suse)*exp(a/Suse);"
      complifeU <- (pt_est[ishift+2]/SUSE)*exp(pt_est[ishift+1]/SUSE)
    }
    if(missing(SUSE)==FALSE && missing(SACC)==FALSE){
      # AFheading <- paste(c("AFat",SACC),collapse = "")
      AFheading <- paste(c("AFatSACC"),collapse = "")
      AF <- "(Sacc/Suse)*exp(a*((1/Suse) - (1/Sacc)));"
      complifeU <- (pt_est[ishift+2]/SUSE)*exp(pt_est[ishift+1]/SUSE)
      compAF <- (SACC/SUSE)*exp(((1/SUSE) - (1/SACC))*pt_est[ishift+1])
    }
  }

  if (ls=="Eyring2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "(1/Sf).*exp(-(a - (b/Sf)))"
    loglifeF <- "-log(Sf) - a + (b/Sf)"

    if(is.null(Tc)==FALSE){
      lifeC <- "(1/Sc).*exp(-(a - (b/Sc)))"
      loglifeC <- "-log(Sc) - a + (b/Sc)"
    }
    if(missing(SUSE)==FALSE){
      lifeU <- "(1/Suse)*exp(-(a - (b/Suse)));"
      complifeU <- (1/SUSE)*exp(-(pt_est[ishift+1] - (pt_est[ishift+2]/SUSE)))
    }
    if(missing(SUSE)==FALSE && missing(SACC)==FALSE){
      # AFheading <- paste(c("AFat",SACC),collapse = "")
      AFheading <- paste(c("AFatSACC"),collapse = "")
      AF <- "(Sacc/Suse)*exp(b*((1/Suse) - (1/Sacc)));"
      complifeU <- (1/SUSE)*exp(-(pt_est[ishift+1] - (pt_est[ishift+2]/SUSE)))
      compAF <- (SACC/SUSE)*exp(((1/SUSE) - (1/SACC))*pt_est[ishift+2])
    }
  }

  if (ls=="Power") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "b*(Sf^a)"
    loglifeF <- "log(b) + a*log(Sf)"
    if(is.null(Tc)==FALSE){
      lifeC <- "b*(Sc^a)"
      loglifeC <- "log(b) + a*log(Sc)"
    }
    if(missing(SUSE)==FALSE){
      lifeU <- "b*(Suse^a);"
      complifeU <- pt_est[ishift+2]*(SUSE^pt_est[ishift+1])
    }
    if(missing(SUSE)==FALSE && missing(SACC)==FALSE){
      # AFheading <- paste(c("AFat",SACC),collapse = "")
      AFheading <- paste(c("AFatSACC"),collapse = "")
      AF <- "(Suse/Sacc)^a;"
      complifeU <- pt_est[ishift+2]*(SUSE^pt_est[ishift+1])
      compAF <- (SUSE/SACC)^pt_est[ishift+1]
    }
  }


  if (ls=="InversePower") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "b*(Sf^-a)"
    loglifeF <- "log(b) - a*log(Sf)"
    if(is.null(Tc)==FALSE){
      lifeC <- "b*(Sc^-a)"
      loglifeC <- "log(b) - a*log(Sc)"
    }
    if(missing(SUSE)==FALSE){
      lifeU <- "b*(Suse^-a);"
      complifeU <- pt_est[ishift+2]*(SUSE^-pt_est[ishift+1])
    }
    if(missing(SUSE)==FALSE && missing(SACC)==FALSE){
      # AFheading <- paste(c("AFat",SACC),collapse = "")
      AFheading <- paste(c("AFatSACC"),collapse = "")
      AF <- "(Sacc/Suse)^a;"
      complifeU <- pt_est[ishift+2]*(SUSE^-pt_est[ishift+1])
      compAF <- (SACC/SUSE)^pt_est[ishift+1]
    }
  }

  if (ls=="InversePower2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "1/(b*(Sf^a))"
    loglifeF <- "-log(b) - a*log(Sf)"
    if(is.null(Tc)==FALSE){
      lifeC <- "1/(b*(Sc^a))"
      loglifeC <- "-log(b) - a*log(Sc)"
    }
    if(missing(SUSE)==FALSE){
      lifeU <- "1/(b*(Suse^a));"
      complifeU <- 1/(pt_est[ishift+2]*(SUSE^pt_est[ishift+1]))
    }
    if(missing(SUSE)==FALSE && missing(SACC)==FALSE){
      # AFheading <- paste(c("AFat",SACC),collapse = "")
      AFheading <- paste(c("AFatSACC"),collapse = "")
      AF <- "(Sacc/Suse)^a;"
      complifeU <- 1/(pt_est[ishift+2]*(SUSE^pt_est[ishift+1]))
      compAF <- (SACC/SUSE)^pt_est[ishift+1]
    }
  }

  if (ls=="Logarithmic") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- "a*log(Sf) + b"
    loglifeF <- "log(a*log(Sf) + b)"
    if(is.null(Tc)==FALSE){
      lifeC <- "a*log(Sc) + b"
      loglifeC <- "log(a*log(Sc) + b)"
    }
    if(missing(SUSE)==FALSE){
      lifeU <- "a*log(Suse) + b;"
      complifeU <- pt_est[ishift+2] + log(SUSE)*pt_est[ishift+1]
    }
    if(missing(SUSE)==FALSE && missing(SACC)==FALSE){
      # AFheading <- paste(c("AFat",SACC),collapse = "")
      AFheading <- paste(c("AFatSACC"),collapse = "")
      AF <- "(a*log(Suse) + b)/(a*log(Sacc) + b);"
      complifeU <- pt_est[ishift+2] + log(SUSE)*pt_est[ishift+1]
      compAF <- (pt_est[ishift+2] + log(SUSE)*pt_est[ishift+1])/(pt_est[ishift+2] + log(SACC)*pt_est[ishift+1])
    }
  }

  if (ls=="MultiStress") {
    # CHECK THIS LAST
    # lsparams - c(ao,a1,...,an), S - c(1,S1,...,Sn)
    if((length(priors)-ishift)==2){
      lsparams <- "real a0; real a1; "
      lsparamsvec <- c("a0","a1")
      pr1<-paste(c("a0 ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("a1 ~ ",priors[ishift+2],";"),collapse = "")
      lspriors <- paste(c(pr1,pr2),collapse = " ")

      lifeF <- "exp(a0 + a1.*Sf)"
      loglifeF <- "a0 + a1.*Sf"
      if(is.null(Tc)==FALSE){
        lifeC <- "exp(a0 + a1.*Sc)"
        loglifeC <- "a0 + a1.*Sc"
      }
      if(missing(SUSE)==FALSE){
        lifeU <- "exp(a0 + a1*Suse);"
        complifeU <- exp(pt_est[ishift+1] + pt_est[ishift+2]*SUSE)
      }
      if(missing(SUSE)==FALSE && missing(SACC)==FALSE){
        # AFheading <- paste(c("AFat",SACC),collapse = "")
        AFheading <- paste(c("AFatSACC"),collapse = "")
        AF <- "exp(a1*(Suse - Sacc));"
        complifeU <- exp(pt_est[ishift+1] + pt_est[ishift+2]*SUSE)
        compAF <- exp((SUSE - SACC)*pt_est[ishift+2])
      }
    }
    if((length(priors)-ishift)==3){
      lsparams <- "real a0; real a1; real a2; "
      lsparamsvec <- c("a0","a1","a2")
      pr1<-paste(c("a0 ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("a1 ~ ",priors[ishift+2],";"),collapse = "")
      pr3<-paste(c("a2 ~ ",priors[ishift+3],";"),collapse = "")
      lspriors <- paste(c(pr1,pr2,pr3),collapse = " ")

      lifeF <- "Lifei[i] = exp(a0 + a1*Sf[i,1] + a2*Sf[i,2]);"
      loglifeF <- "Lifei[i] = a0 + a1*Sf[i,1] + a2*Sf[i,2];"
      if(is.null(Tc)==FALSE){
        lifeC <- "Lifej[j] = exp(a0 + a1*Sc[j,1] + a2*Sc[j,2]);"
        loglifeC <- "Lifej[j] = a0 + a1*Sc[j,1] + a2*Sc[j,2];"
      }
      if(missing(SUSE)==FALSE){
        lifeU <- "exp(a0 + a1*Suse[1] + a2*Suse[2]);"
        complifeU <- exp(pt_est[ishift+1] + pt_est[ishift+2]*SUSE[1] + pt_est[ishift+3]*SUSE[2])
      }
      if(missing(SUSE)==FALSE && missing(SACC)==FALSE){
        # AFheading <- paste(c("AFat",SACC),collapse = "")
        AFheading <- paste(c("AFatSACC"),collapse = "")
        AF <- "exp(a1*(Suse[1] - Sacc[1]) + a2*(Suse[2] - Sacc[2]));"
        complifeU <- exp(pt_est[ishift+1] + pt_est[ishift+2]*SUSE[1] + pt_est[ishift+3]*SUSE[2])
        compAF <- exp((SUSE[1] - SACC[1])*pt_est[ishift+2] + (SUSE[2] - SACC[2])*pt_est[ishift+3])
      }
    }
    if((length(priors)-ishift)==4){
      lsparams <- "real a0; real a1; real a2; real a3;"
      lsparamsvec <- c("a0","a1","a2","a3")
      pr1<-paste(c("a0 ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("a1 ~ ",priors[ishift+2],";"),collapse = "")
      pr3<-paste(c("a2 ~ ",priors[ishift+3],";"),collapse = "")
      pr4<-paste(c("a3 ~ ",priors[ishift+4],";"),collapse = "")
      lspriors <- paste(c(pr1,pr2,pr3,pr4),collapse = " ")

      lifeF <- "Lifei[i] = exp(a0 + a1*Sf[i,1] + a2*Sf[i,2] + a3*Sf[i,3]);"
      loglifeF <- "Lifei[i] = a0 + a1*Sf[i,1] + a2*Sf[i,2] + a3*Sf[i,3];"
      if(is.null(Tc)==FALSE){
        lifeC <- "Lifej[j] = exp(a0 + a1*Sc[j,1] + a2*Sc[j,2] + a3*Sc[j,3]);"
        loglifeC <- "Lifej[j] = a0 + a1*Sc[j,1] + a2*Sc[j,2] + a3*Sc[j,3];"
      }
      if(missing(SUSE)==FALSE){
        lifeU <- "exp(a0 + a1*Suse[1] + a2*Suse[2] + a3*Suse[3]);"
        complifeU <- exp(pt_est[ishift+1] + pt_est[ishift+2]*SUSE[1] + pt_est[ishift+3]*SUSE[2] + pt_est[ishift+4]*SUSE[3])
      }
      if(missing(SUSE)==FALSE && missing(SACC)==FALSE){
        AFheading <- paste(c("AFat"),collapse = "")
        AF <- "exp(a1*(Suse[1] - Sacc[1]) + a2*(Suse[2] - Sacc[2]) + a3*(Suse[3] - Sacc[3]));"
        complifeU <- exp(pt_est[ishift+1] + pt_est[ishift+2]*SUSE[1] + pt_est[ishift+3]*SUSE[2] + pt_est[ishift+4]*SUSE[3])
        compAF <- exp((SUSE[1] - SACC[1])*pt_est[ishift+2] + (SUSE[2] - SACC[2])*pt_est[ishift+3] + (SUSE[3] - SACC[3])*pt_est[ishift+4])
      }
    }
    if((length(priors)-ishift)==5){
      lsparams <- "real a0; real a1; real a2; real a3; real a4;"
      lsparamsvec <- c("a0","a1","a2","a3")
      pr1<-paste(c("a0 ~ ",priors[ishift+1],";"),collapse = "")
      pr2<-paste(c("a1 ~ ",priors[ishift+2],";"),collapse = "")
      pr3<-paste(c("a2 ~ ",priors[ishift+3],";"),collapse = "")
      pr4<-paste(c("a3 ~ ",priors[ishift+4],";"),collapse = "")
      pr5<-paste(c("a4 ~ ",priors[ishift+5],";"),collapse = "")
      lspriors <- paste(c(pr1,pr2,pr3,pr4,pr5),collapse = " ")

      lifeF <- "Lifei[i] = exp(a0 + a1*Sf[i,1] + a2*Sf[i,2] + a3*Sf[i,3] + a4*Sf[i,4]);"
      loglifeF <- "Lifei[i] = a0 + a1*Sf[i,1] + a2*Sf[i,2] + a3*Sf[i,3] + a4*Sf[i,4];"
      if(is.null(Tc)==FALSE){
        lifeC <- "Lifej[j] = exp(a0 + a1*Sc[j,1] + a2*Sc[j,2] + a3*Sc[j,3] + a4*Sc[j,4]);"
        loglifeC <- "Lifej[j] = a0 + a1*Sc[j,1] + a2*Sc[j,2] + a3*Sc[j,3] + a4*Sc[j,4];"
      }
      if(missing(SUSE)==FALSE){
        lifeU <- "exp(a0 + a1*Suse[1] + a2*Suse[2] + a3*Suse[3] + a4*Suse[4]);"
        complifeU <- exp(pt_est[ishift+1] + pt_est[ishift+2]*SUSE[1] + pt_est[ishift+3]*SUSE[2] + pt_est[ishift+4]*SUSE[3] + pt_est[ishift+5]*SUSE[4])
      }
      if(missing(SUSE)==FALSE && missing(SACC)==FALSE){
        # AFheading <- paste(c("AFat",SACC),collapse = "")
        AFheading <- paste(c("AFatSACC"),collapse = "")
        AF <- "exp(a1*(Suse[1] - Sacc[1]) + a2*(Suse[2] - Sacc[2]) + a3*(Suse[3] - Sacc[3]) + a4*(Suse[4] - Sacc[4]));"
        complifeU <- exp(pt_est[ishift+1] + pt_est[ishift+2]*SUSE[1] + pt_est[ishift+3]*SUSE[2] + pt_est[ishift+4]*SUSE[3] + pt_est[ishift+5]*SUSE[4])
        compAF <- exp((SUSE[1] - SACC[1])*pt_est[ishift+2] + (SUSE[2] - SACC[2])*pt_est[ishift+3] + (SUSE[3] - SACC[3])*pt_est[ishift+4] + (SUSE[4] - SACC[4])*pt_est[ishift+5])
      }
    }
  }

  if (ls=="TempHumidity") {
    # lsparams[1] - parameter A, lsparams[2] - parameter a, lsparams[3] - parameter b
    lsparams <- "real<lower=0> A; real a; real b;"
    lsparamsvec <- c("A","a","b")
    pr1<-paste(c("A ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("a ~ ",priors[ishift+2],";"),collapse = "")
    pr3<-paste(c("b ~ ",priors[ishift+3],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2,pr3),collapse = " ")

    lifeF <- "A.*exp((a/Sf[,1]) + (b/Sf[,2]));"
    loglifeF <- "log(A) + (a/Sf[,1]) + (b/Sf[,2]);"
    if(is.null(Tc)==FALSE){
      lifeC <- "A.*exp((a/Sc[1]) + (b/Sc[2]))"
      loglifeC <- "log(A) + (a/Sc[1]) + (b/Sc[2])"
    }
    if(missing(SUSE)==FALSE){
      lifeU <- "A*exp((a/Suse[1]) + (b/Suse[2]));"
      complifeU <- pt_est[ishift+1]*exp((pt_est[ishift+2]/SUSE[1]) + (pt_est[ishift+3]/SUSE[2]))
    }
    if(missing(SUSE)==FALSE && missing(SACC)==FALSE){
      # AFheading <- paste(c("AFat",SACC),collapse = "")
      AFheading <- paste(c("AFatSACC"),collapse = "")
      # AFheading <- paste(c("AFat",SACC[1],"_",SACC[2]),collapse = "")
      AF <- "exp(a*((1/Suse[1]) - (1/Sacc[1])) + b*((1/Suse[2]) - (1/Sacc[2])));"
      complifeU <- pt_est[ishift+1]*exp((pt_est[ishift+2]/SUSE[1]) + (pt_est[ishift+3]/SUSE[2]))
      compAF <- exp(pt_est[ishift+2]*(1/SUSE[1] - 1/SACC[1]) + pt_est[ishift+3]*(1/SUSE[2] - 1/SACC[2]))
    }
  }

  if (ls=="TempNonthermal") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
    lsparams <- "real a; real b; real<lower=0> c;"
    lsparamsvec <- c("a","b","c")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    pr3<-paste(c("c ~ ",priors[ishift+3],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2,pr3),collapse = " ")

    lifeF <- "c./((Sf[,2].^b)*exp(-a./Sf[,1]))"
    loglifeF <- "log(c) - b.*log(Sf[,2]) + (a/Sf[,1])"
    if(is.null(Tc)==FALSE){
      lifeC <- "c/((Sc[,2]^b)*exp(-a/Sc[,1]))"
      loglifeC <- "log(c) - b*log(Sc[,2]) + (a/Sc[,1])"
    }
    if(missing(SUSE)==FALSE){
      lifeU <- "c/((Suse[2]^b)*exp(-a/Suse[1]));"
      complifeU <- pt_est[ishift+3]/((SUSE[2]^pt_est[ishift+2])*exp(-pt_est[ishift+1]/SUSE[1]))
    }
    if(missing(SUSE)==FALSE && missing(SACC)==FALSE){
      AFheading <- paste(c("AFatSACC"),collapse = "")
      AF <- "((Sacc[2]/Suse[2])^b)*exp(-a*((1/Suse[1]) - (1/Sacc[1])));"
      complifeU <- pt_est[ishift+3]/((SUSE[2]^pt_est[ishift+2])*exp(-pt_est[ishift+1]/SUSE[1]))
      compAF <- ((SACC[2]/SUSE[2])^pt_est[ishift+2])*exp(-((1/SUSE[1]) - (1/SACC[1]))*pt_est[ishift+1])
    }
  }

  if (ls=="Eyring3") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # lsparams[3] - parameter c, lsparams[4] - parameter d
    lsparams <- "real a; real b; real c; real d;"
    lsparamsvec <- c("a","b","c","d")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    pr3<-paste(c("c ~ ",priors[ishift+3],";"),collapse = "")
    pr4<-paste(c("d ~ ",priors[ishift+4],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2,pr3,pr4),collapse = " ")

    lifeF <- "(1/Sf[,1]).*exp((a + (b/Sf[,1])) + (c + (d/Sf[,1])).*Sf[,2])"
    loglifeF <- "-log(Sf[,1]) + a + (b/Sf[,1]) + (c + (d/Sf[,1])).*Sf[,2]"
    if(is.null(Tc)==FALSE){
      lifeC <- "(1/Sc[,1]).*exp((a + (b/Sc[,1])) + (c + (d/Sc[,1])).*Sc[,2])"
      loglifeC <- "-log(Sc[,1]) + a + (b/Sc[,1]) + (c + (d/Sc[,1])).*Sc[,2]"
    }
    if(missing(SUSE)==FALSE){
      lifeU <- "(1/Suse[1])*exp((a + (b/Suse[1])) + (c + (d/Suse[1]))*Suse[2]);"
      complifeU <- (1/SUSE[1])*exp((pt_est[ishift+1] + (pt_est[ishift+2]/SUSE[1])) + (pt_est[ishift+3] + (pt_est[ishift+4]/SUSE[1]))*SUSE[2])
    }
    if(missing(SUSE)==FALSE && missing(SACC)==FALSE){
      # AFheading <- paste(c("AFat",SACC),collapse = "")
      AFheading <- paste(c("AFatSACC"),collapse = "")
      AF <- "(Sacc[1]/Suse[1])*exp(b*((1/Suse[1]) - (1/Sacc[1])) + c*(Suse[2] - Sacc[2]) + d*((Suse[2]/Suse[1]) - (Sacc[2]/Sacc[1])));"
      complifeU <- (1/SUSE[1])*exp((pt_est[ishift+1] + (pt_est[ishift+2]/SUSE[1])) + (pt_est[ishift+3] + (pt_est[ishift+4]/SUSE[1]))*SUSE[2])
      compAF <- (SACC[1]/SUSE[1])*exp(pt_est[ishift+2]*((1/SUSE[1]) - (1/SACC[1])) + pt_est[ishift+3]*(SUSE[2] - SACC[2]) + pt_est[ishift+4]*((SUSE[2]/SUSE[1]) - (SACC[2]/SACC[1])))
    }
  }



  # Fit to log-likelihood distributions
  if (dist=="Weibull") {
    if(is.null(param2) == TRUE){
      distparam <-"real<lower=0> beta;"
      distpriors<-paste(c("beta ~ ",priors[ishift],";"),collapse = "")

      if(ls=="MultiStress"){
        if(is.null(Tc)==TRUE){
          loglik <- paste(c("target += weibull_lpdf(TTF | beta, Lifei);"),collapse = "")
        } else{
          loglik <- paste(c("target += weibull_lpdf(TTF | beta, Lifei) + weibull_lccdf(TTS | beta, Lifej);"),collapse = "")
        }
      } else{
        if(is.null(Tc)==TRUE){
          loglik <- paste(c("target += weibull_lpdf(TTF | beta,",lifeF,");"),collapse = "")
        } else{
          loglik <- paste(c("target += weibull_lpdf(TTF | beta,",lifeF,") + weibull_lccdf(TTS | beta,",lifeC,");"),collapse = "")
        }
      }
      params <- paste(c(distparam,lsparams),collapse = " ")
      paramsvec <- c("beta",lsparamsvec)
      outputparamset <- c("\U03B2",lsparamsvec)
      priors <- paste(c(distpriors,lspriors),collapse = " ")
    }
    if(is.null(param2) == FALSE){
      if(ls=="MultiStress"){
        if(is.null(Tc)==TRUE){
          loglik <- paste(c("target += weibull_lpdf(TTF | ",param2,", Lifei);"),collapse = "")
        } else{
          loglik <- paste(c("target += weibull_lpdf(TTF | ",param2,", Lifei) + weibull_lccdf(TTS | ",param2,", Lifej);"),collapse = "")
        }
      } else{
        if(is.null(Tc)==TRUE){
          loglik <- paste(c("target += weibull_lpdf(TTF | ",param2,",",lifeF,");"),collapse = "")
        } else{
          loglik <- paste(c("target += weibull_lpdf(TTF | ",param2,",",lifeF,") + weibull_lccdf(TTS | ",param2,",",lifeC,");"),collapse = "")
        }
      }
      params <- paste(c(lsparams),collapse = " ")
      paramsvec <- lsparamsvec
      outputparamset <- lsparamsvec
      priors <- paste(c(lspriors),collapse = " ")
    }

  }
  if (dist=="3PWeibull") {
    distparam <-"real<lower=0> beta; real gamma;"
    distpriors<-paste(c("beta ~ ",priors[ishift],";","gamma ~ ",priors[ishift+1],";"),collapse = "")

    if(ls=="MultiStress"){
      if(is.null(Tc)==TRUE){
        loglik <- paste(c("target += weibull_lpdf(TTF | beta, Lifei);"),collapse = "")
      } else{
        loglik <- paste(c("target += weibull_lpdf(TTF | beta, Lifei) + weibull_lccdf(TTS | beta, Lifej);"),collapse = "")
      }
    } else{
      if(is.null(Tc)==TRUE){
        loglik <- paste(c("target += weibull_lpdf(TTF | beta,",lifeF,");"),collapse = "")
      } else{
        loglik <- paste(c("target += weibull_lpdf(TTF | beta,",lifeF,") + weibull_lccdf(TTS | beta,",lifeC,");"),collapse = "")
      }
    }
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("beta","gamma",lsparamsvec)
    outputparamset <- c("\U03B2","\U03B3",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")
  }
  if (dist=="Lognormal") {
    distparam <-"real<lower=0> sigma_t;"
    distpriors<-paste(c("sigma_t ~ ",priors[ishift],";"),collapse = "")

    if(ls=="MultiStress"){
      if(is.null(Tc)==TRUE){
        loglik <- paste(c("target += lognormal_lpdf(TTF |Lifei, sigma_t);"),collapse = "")
      } else{
        loglik <- paste(c("target += lognormal_lpdf(TTF |Lifei, sigma_t) + lognormal_lccdf(TTS |Lifej, sigma_t);"),collapse = "")
      }
    } else{
      if(is.null(Tc)==TRUE){
        loglik <- paste(c("target += lognormal_lpdf(TTF |",loglifeF,", sigma_t);"),collapse = "")
      } else{
        loglik <- paste(c("target += lognormal_lpdf(TTF |",loglifeF,", sigma_t) + lognormal_lccdf(TTS |",loglifeC,", sigma_t);"),collapse = "")
      }
    }

    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("sigma_t",lsparamsvec)
    outputparamset <- c("\U03C3_t",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")
  }
  if (dist=="Normal") {
    distparam <-"real<lower=0> sigma;"
    distpriors<-paste(c("sigma ~ ",priors[ishift],";"),collapse = "")
    if(ls=="MultiStress"){
      if(is.null(Tc)==TRUE){
        loglik <- paste(c("target += normal_lpdf(TTF | Lifei, sigma);"),collapse = "")
      } else{
        loglik <- paste(c("target += normal_lpdf(TTF | Lifei, sigma) + normal_lccdf(TTS | Lifej, sigma);"),collapse = "")
      }
    } else{
      if(is.null(Tc)==TRUE){
        loglik <- paste(c("target += normal_lpdf(TTF |",lifeF,", sigma);"),collapse = "")
      } else{
        loglik <- paste(c("target += normal_lpdf(TTF |",lifeF,", sigma) + normal_lccdf(TTS |",lifeC,", sigma);"),collapse = "")
      }
    }
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("sigma",lsparamsvec)
    outputparamset <- c("\U03C3",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")
  }
  if (dist=="Exponential") {
    if(ls=="MultiStress"){
      if(is.null(Tc)==TRUE){
        loglik <- paste(c("target += exponential_lpdf(TTF | 1/(Lifei));"),collapse = "")
      } else{
        loglik <- paste(c("target += exponential_lpdf(TTF | 1/(Lifei)) + exponential_lccdf(TTS 1/(Lifej));"),collapse = "")
      }
    } else{
      if(is.null(Tc)==TRUE){
        loglik <- paste(c("target += exponential_lpdf(TTF | 1/(",lifeF,"));"),collapse = "")
      } else{
        loglik <- paste(c("target += exponential_lpdf(TTF | 1/(",lifeF,")) + exponential_lccdf(TTS 1/(",lifeC,"));"),collapse = "")
      }
    }
    params <- lsparams
    paramsvec <- lsparamsvec
    outputparamset <- lsparamsvec
    priors <- lspriors
  }
  if (dist=="2PExponential") {
    distparam <-"real<lower=0> sigma;"
    distpriors<-paste(c("sigma ~ ",priors[ishift],";"),collapse = "")
    if(ls=="MultiStress"){
      if(is.null(Tc)==TRUE){
        loglik <- paste(c("target += double_exponential_lpdf(TTF | Lifei, sigma);"),collapse = "")
      } else{
        loglik <- paste(c("target += double_exponential_lpdf(TTF | Lifei, sigma) + double_exponential_lccdf(TTS | Lifej, sigma);"),collapse = "")
      }
    } else{
      if(is.null(Tc)==TRUE){
        loglik <- paste(c("target += double_exponential_lpdf(TTF |",lifeF,", sigma);"),collapse = "")
      } else{
        loglik <- paste(c("target += double_exponential_lpdf(TTF |",lifeF,", sigma) + double_exponential_lccdf(TTS |",lifeC,", sigma);"),collapse = "")
      }
    }
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("sigma",lsparamsvec)
    outputparamset <- c("\U03C3",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")
  }

  # Define stancode here
  if(is.null(Tc) == TRUE){
    # Identify if the stress Sf is single or not
    if(is.null(dim(SF))==TRUE){ # Single stress case Sf treated as vector
      block1 <- "data {int<lower=0> n; vector[n] TTF; vector[n] Sf;}"
      datablock <- list(n = length(TTF), TTF = TTF, Sf = SF)
      if(is.null(SUSE)==FALSE){
        block1 <- "data {int<lower=0> n; vector[n] TTF; vector[n] Sf; real Suse;}"
        datablock <- list(n = length(TTF), TTF = TTF, Sf = SF, Suse = SUSE)
      }
      if(is.null(SUSE)==FALSE && is.null(SACC)==FALSE){
        block1 <- "data {int<lower=0> n; vector[n] TTF; vector[n] Sf; real Suse; real Sacc;}"
        datablock <- list(n = length(TTF), TTF = TTF, Sf = SF, Suse = SUSE, Sacc = SACC)
      }
    }
    if(is.null(dim(SF))==FALSE){ # multi-stress case, Sf treated as matrix
      block1 <- "data {int<lower=0> n; int<lower=0> n2; vector[n] TTF; matrix[n,n2] Sf;}"
      datablock <- list(n = length(TTF),  n2 = dim(SF)[2], TTF = TTF, Sf = SF)
      if(is.null(SUSE)==FALSE){
        block1 <- "data {int<lower=0> n; int<lower=0> n2; vector[n] TTF; matrix[n,n2] Sf; vector[n2] Suse;}"
        datablock <- list(n = length(TTF), n2 = dim(SF)[2], TTF = TTF, Sf = SF, Suse = SUSE)
      }
      if(is.null(SUSE)==FALSE && is.null(SACC)==FALSE){
        block1 <- "data {int<lower=0> n; int<lower=0> n2; vector[n] TTF; matrix[n,n2] Sf; vector[n2] Suse; vector[n2] Sacc;}"
        datablock <- list(n = length(TTF), n2 = dim(SF)[2], TTF = TTF, Sf = SF, Suse = SUSE, Sacc = SACC)
      }
    }

  } else{
    if(is.null(dim(SF))==TRUE){
      block1 <- "data {int<lower=0> n; int<lower=0> m; vector[n] TTF; vector[m] TTS; vector[n] Sf; vector[m] Sc;}"
      datablock <- list(n = length(TTF), m = length(Tc), TTF = TTF, Sf = SF, TTS = Tc, Sc = Sc)
      if(is.null(SUSE)==FALSE){
        block1 <- "data {int<lower=0> n; int<lower=0> m; vector[n] TTF; vector[m] TTS; vector[n] Sf; vector[m] Sc; real Suse;}"
        datablock <- list(n = length(TTF), m = length(Tc), TTF = TTF, Sf = SF, TTS = Tc, Sc = Sc, Suse = SUSE)
      }
      if(is.null(SUSE)==FALSE && is.null(SACC)==FALSE){
        block1 <- "data {int<lower=0> n; int<lower=0> m; vector[n] TTF; vector[m] TTS; vector[n] Sf; vector[m] Sc; real Suse; real Sacc;}"
        datablock <- list(n = length(TTF), m = length(Tc), TTF = TTF, Sf = SF, TTS = Tc, Sc = Sc, Suse = SUSE, Sacc = SACC)
      }
    }
    if(is.null(dim(SF))==FALSE){
      block1 <- "data {int<lower=0> n; int<lower=0> m; vector[n] TTF; vector[m] TTS; vector[n] Sf; vector[m] Sc;}"
      datablock <- list(n = length(TTF), m = length(Tc), TTF = TTF, Sf = SF, TTS = Tc, Sc = Sc)
      if(is.null(SUSE)==FALSE){
        block1 <- "data {int<lower=0> n; int<lower=0> m; vector[n] TTF; vector[m] TTS; vector[n] Sf; vector[m] Sc; real Suse;}"
        datablock <- list(n = length(TTF), m = length(Tc), TTF = TTF, Sf = SF, TTS = Tc, Sc = Sc, Suse = SUSE)
      }
      if(is.null(SUSE)==FALSE && is.null(SACC)==FALSE){
        block1 <- "data {int<lower=0> n; int<lower=0> m; vector[n] TTF; vector[m] TTS; vector[n] Sf; vector[m] Sc; real Suse; real Sacc;}"
        datablock <- list(n = length(TTF), m = length(Tc), TTF = TTF, Sf = SF, TTS = Tc, Sc = Sc, Suse = SUSE, Sacc = SACC)
      }
    }

  }
  block2 <- paste(c("parameters {",params,"}"),collapse = " ")
  if(is.null(SUSE)==FALSE){
    block2b <- paste(c("transformed parameters { real<lower=0> Uselife; Uselife = ",lifeU,"}"),collapse = " ")
    paramsvec0 <- c(paramsvec,"Uselife")
    # pt_est <- c(pt_est,complifeU)
  }
  if(is.null(SUSE)==FALSE && is.null(SACC)==FALSE){
    block2b <- paste(c("transformed parameters { real<lower=0> Uselife; real<lower=0> ",AFheading,"; Uselife = ",lifeU,AFheading," = ",AF,"}"),collapse = " ")
    paramsvec0 <- c(paramsvec,"Uselife",AFheading)
    # pt_est <- c(pt_est,complifeU,compAF)
  }
  if(is.null(SUSE)==TRUE && is.null(SACC)==TRUE){
    paramsvec0 <- paramsvec
  }
  block3 <- paste(c("model {",priors,loglik,"}"),collapse = " ")
  if (ls=="MultiStress" && dist == "Lognormal"){
    if(is.null(Tc)==TRUE){
      block3 <- paste(c("model { vector[n] Lifei; vector[m] Lifej; ",priors," for(i in 1:n){",loglifeF,"}",loglik,"}"),collapse = " ")
    }
    if(is.null(Tc)==FALSE){
      block3 <- paste(c("model { vector[n] Lifei; vector[m] Lifej; ",priors," for(i in 1:n){",loglifeF,"} for(j in 1:m){",loglifeC,"}",loglik,"}"),collapse = " ")
    }
  }
  if (ls=="MultiStress" && (dist == "Normal" || dist=="Weibull" || dist=="Exponential")){
    if(is.null(Tc)==TRUE){
      block3 <- paste(c("model { vector[n] Lifei; ",priors," for(i in 1:n){",lifeF,"}",loglik,"}"),collapse = " ")
    }
    if(is.null(Tc)==FALSE){
      block3 <- paste(c("model { vector[n] Lifei; vector[m] Lifej; ",priors," for(i in 1:n){",lifeF,"} for(j in 1:m){",lifeC,"}",loglik,"}"),collapse = " ")
    }
  }
  # NOT RUN {
  stanlscode <- paste(c(block1,block2,block3),collapse=" ")
  if(is.null(SUSE)==FALSE || is.null(SACC)==FALSE){
    stanlscode <- paste(c(block1,block2,block2b,block3),collapse=" ")
  }
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
  # return(list(stanlscode,stanlsfile))

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

  # Produce some output text that summarizes the results
  cat(c("Posterior estimates for Bayesian Analysis.\n\n"),sep = "")
  print(outputtable)
  cat(c("\n"),sep = "")


  return(list(fit,plot1_MCtrace,plot2_hist,plot3_density))
  # NOTE: Comment return above and uncomment return below if you want to view the scatterplot between posteriors.
  # This may however increase processing time.
  # return(list(fit,plot1_MCtrace,plot2_hist,plot3_density,plot4_densityoverlay,plot5_scatterplot))
}
