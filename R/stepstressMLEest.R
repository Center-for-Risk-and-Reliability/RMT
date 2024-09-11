# MLE Step-Stress Estimator
# Developed by Dr. Reuel Smith, 2021-2022

stepstress.MLEest <- function(LSQest,data,stepstresstable,ls,dist,confid=0.95,sided="twosided"){
  #Load pracma library for erf
  library(pracma)
  library(dplyr)
  library(matrixcalc)
  library(stringr)
  # UPDATE(11/14/2023): Going to set the more troublesome parameters (the ones that need to be positive)
  # to log form like the others

  # Check to see if dist="Exponential" so you can exclude life
  # distribution parameters.
  if (dist=="Exponential") {
    ishift<-0
  } else {
    ishift<-1
  }
  # Check to see if confidence exists
  conf.level <- confid

  # Setup positivity check vector for parameters
  positivity_v<-rep(0,length(LSQest))

  if(is.list(data[[1]])==FALSE){
    # FOR ONE STEP-STRESS TEST
    # Re-sort the input data and then separate censored from non-censored data
    stpstrdatsort<-stepstress.data.cum(data,stepstresstable)
    xirctot<-sort.xircstressdata(stpstrdatsort[[1]])
    # The failure and survival stress levels are assigned a vector of the same
    # sizes as the adjusted time vectors for the MLE operation.
    tF <- xirctot[[1]]
    tc <- xirctot[[2]]
    SF <- xirctot[[3]]
    Sc <- xirctot[[4]]
    tF_end <- stpstrdatsort[[3]]
    tc_end <- stpstrdatsort[[4]]
    # Obtain the cumulative sum of the end times
    t_end <- cumsum(stepstresstable[,dim(stepstresstable)[[2]]])
    # Compute number of stresses
    Nstress<-dim(stepstresstable)[2]-1
  } else{
    # FOR SET OF STEP-STRESS TESTS TO BE COMPUTED TOGETHER
    stpstrdatsort<-vector(mode = "list", length(data))
    xirctot<-vector(mode = "list", length(data))
    tF<-vector(mode = "list", length(data))
    tc<-vector(mode = "list", length(data))
    SF<-vector(mode = "list", length(data))
    Sc<-vector(mode = "list", length(data))
    tF_end <- vector(mode = "list", length(data))
    tc_end <- vector(mode = "list", length(data))
    t_end<-vector(mode = "list", length(data))
    S_end<-vector(mode = "list", length(data))
    for(i in 1:length(data)){
      # Re-sort the input data and then separate censored from non-censored data
      stpstrdatsort[[i]]<-stepstress.data.cum(data[[i]],stepstresstable[[i]])
      xirctot[[i]]<-sort.xircstressdata(stpstrdatsort[[i]][[1]])
      # The failure and survival stress levels are assigned a vector of the same
      # sizes as the adjusted time vectors for the MLE operation.
      tF[[i]] <- xirctot[[i]][[1]]
      tc[[i]] <- xirctot[[i]][[2]]
      SF[[i]] <- xirctot[[i]][[3]]
      Sc[[i]] <- xirctot[[i]][[4]]
      tF_end[[i]] <- stpstrdatsort[[i]][[3]]
      tc_end[[i]] <- stpstrdatsort[[i]][[4]]
      # Obtain the cumulative sum of the end times
      t_end[[i]] <- cumsum(stepstresstable[[i]][,dim(stepstresstable[[i]])[[2]]])
      S_end[[i]] <- stepstresstable[[i]][,1:(dim(stepstresstable[[i]])[[2]]-1)]
    }
    t_S_end <- mapply(c,t_end,S_end, SIMPLIFY=FALSE)
    tF<-unlist(tF)
    tc<-unlist(tc)
    SF<-unlist(SF)
    Sc<-unlist(Sc)
    tF_end <- unlist(tF_end)
    tc_end <- unlist(tc_end)
    # t_end<-unlist(t_end)
    Nstress<-dim(stepstresstable[[1]])[2]-1
  }

  # return(list(t_end,S_end,tF_end,tc_end))

  # Initialize life-stress parameter estimates for theta
  if (ls=="Linear") {
    # theta[1] - parameter a, theta[2] - parameter b

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      (theta[ishift+2] + S1*theta[ishift+1])/(theta[ishift+2] + S2*theta[ishift+1])
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      theta[ishift+1]*(S2 - S1)
    }

    # Life functions
    lifeF <- function(theta) {
      theta[ishift+2] + SF*theta[ishift+1]
    }
    loglifeF <- function(theta) {
      log(theta[ishift+2] + SF*theta[ishift+1])
    }
    lifeC <- function(theta) {
      theta[ishift+2] + Sc*theta[ishift+1]
    }
    loglifeC <- function(theta) {
      log(theta[ishift+2] + Sc*theta[ishift+1])
    }
    ls_txt<-ls
    params_txt<-c("a","b")
  }

  if (ls=="Exponential"){
    # theta[1] - parameter a, theta[2] - parameter b
    # shift b to log b
    # positivity_v[ishift+2]<-1
    LSQest[ishift+2] <- log(LSQest[ishift+2])

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      exp(theta[ishift+1]*(S1 - S2))
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      exp(theta[ishift+2])*((exp(theta[ishift+1]*S2)) - (exp(theta[ishift+1]*S1)))
    }

    # Life functions
    lifeF <- function(theta) {
      exp(theta[ishift+2])*exp(SF*theta[ishift+1])
    }
    loglifeF <- function(theta) {
      theta[ishift+2] + SF*theta[ishift+1]
    }
    lifeC <- function(theta) {
      exp(theta[ishift+2])*exp(Sc*theta[ishift+1])
    }
    loglifeC <- function(theta) {
      theta[ishift+2] + Sc*theta[ishift+1]
    }
    ls_txt<-ls
    params_txt<-c("a","b")
  }

  if (ls=="Exponential2"){
    # theta[1] - parameter a, theta[2] - parameter b
    # shift b to log b
    # positivity_v[ishift+2]<-1
    LSQest[ishift+2] <- log(LSQest[ishift+2])

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      exp(theta[ishift+1]*((1/S1) - (1/S2)))
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      exp(theta[ishift+2])*((exp(theta[ishift+1]/S2)) - (exp(theta[ishift+1]/S1)))
    }

    # Life functions
    lifeF <- function(theta) {
      exp(theta[ishift+2])*exp(theta[ishift+1]/SF)
    }
    loglifeF <- function(theta) {
      theta[ishift+2] + (theta[ishift+1]/SF)
    }
    lifeC <- function(theta) {
      exp(theta[ishift+2])*exp(theta[ishift+1]/Sc)
    }
    loglifeC <- function(theta) {
      theta[ishift+2] + (theta[ishift+1]/Sc)
    }
    ls_txt<-ls
    params_txt<-c("a","b")
  }

  if (ls=="Arrhenius") {
    # lsparams[1] - parameter Ea, lsparams[2] - parameter b
    # Temperature HAS to be in Kelvin for this to work
    # shift b to log b
    # positivity_v[ishift+2]<-1
    LSQest[ishift+2] <- log(LSQest[ishift+2])

    K<-8.617385e-5

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      exp((theta[ishift+1]/K)*((1/S1) - (1/S2)))
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      exp(theta[ishift+2])*((exp(theta[ishift+1]/(K*S2))) - (exp(theta[ishift+1]/(K*S1))))
    }

    # Life functions
    lifeF <- function(theta) {
      exp(theta[ishift+2])*exp(theta[ishift+1]/(K*SF))
    }
    loglifeF <- function(theta) {
      theta[ishift+2] + theta[ishift+1]/(K*SF)
    }
    lifeC <- function(theta) {
      exp(theta[ishift+2])*exp(theta[ishift+1]/(K*Sc))
    }
    loglifeC <- function(theta) {
      theta[ishift+2] + theta[ishift+1]/(K*Sc)
    }
    ls_txt<-ls
    params_txt<-c("E_a","b")
  }

  if (ls=="Eyring") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      (S2/S1)*exp(theta[ishift+1]*((1/S1) - (1/S2)))
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      theta[ishift+2]*(((1/S2)*exp(theta[ishift+1]/S2)) - (1/S1)*exp(theta[ishift+1]/S1))
    }

    # Life functions
    lifeF <- function(theta) {
      (theta[ishift+2]/SF)*exp(theta[ishift+1]/SF)
    }
    loglifeF <- function(theta) {
      log(theta[ishift+2]) - log(SF) + theta[ishift+1]/SF
    }
    lifeC <- function(theta) {
      (theta[ishift+2]/Sc)*exp(theta[ishift+1]/Sc)
    }
    loglifeC <- function(theta) {
      log(theta[ishift+2]) - log(Sc) + theta[ishift+1]/Sc
    }
    ls_txt<-ls
    params_txt<-c("a","b")
  }

  if (ls=="Eyring2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      (S2/S1)*exp(theta[ishift+2]*((1/S1) - (1/S2)))
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      ((1/S2)*exp(-(theta[ishift+1] - (theta[ishift+2]/S2)))) - ((1/S1)*exp(-(theta[ishift+1] - (theta[ishift+2]/S1))))
    }

    # Life functions
    lifeF <- function(theta) {
      (1/SF)*exp(-(theta[ishift+1] - (theta[ishift+2]/SF)))
    }
    loglifeF <- function(theta) {
      -log(SF) - theta[ishift+1] + theta[ishift+2]/SF
    }
    lifeC <- function(theta) {
      (1/Sc)*exp(-(theta[ishift+1] - (theta[ishift+2]/Sc)))
    }
    loglifeC <- function(theta) {
      - log(Sc) - theta[ishift+1] + theta[ishift+2]/Sc
    }
    ls_txt<-ls
    params_txt<-c("a","b")
  }

  if (ls=="Power") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # First redefine b parameter as logb
    LSQest[ishift+2] <- log(LSQest[ishift+2])
    # positivity_v[ishift+2]<-1

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      (S1/S2)^theta[ishift+1]
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      exp(theta[ishift+2])*((S2^theta[ishift+1])-(S1^theta[ishift+1]))
    }

    # Life functions
    lifeF <- function(theta) {
      exp(theta[ishift+2])*(SF^theta[ishift+1])
    }
    loglifeF <- function(theta) {
      theta[ishift+2] + theta[ishift+1]*log(SF)
    }
    lifeC <- function(theta) {
      exp(theta[ishift+2])*(Sc^theta[ishift+1])
    }
    loglifeC <- function(theta) {
      theta[ishift+2] + theta[ishift+1]*log(Sc)
    }
    ls_txt<-ls
    params_txt<-c("a","b")
  }

  if (ls=="InversePower") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # First redefine b parameter as logb
    LSQest[ishift+2] <- log(LSQest[ishift+2])
    # positivity_v[ishift+2]<-1

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      (S2/S1)^theta[ishift+1]
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      exp(theta[ishift+2])*((S2^-theta[ishift+1])-(S1^-theta[ishift+1]))
    }

    # Life functions
    lifeF <- function(theta) {
      exp(theta[ishift+2])*(SF^-theta[ishift+1])
    }
    loglifeF <- function(theta) {
      theta[ishift+2] - theta[ishift+1]*log(SF)
    }
    lifeC <- function(theta) {
      exp(theta[ishift+2])*(Sc^-theta[ishift+1])
    }
    loglifeC <- function(theta) {
      theta[ishift+2] - theta[ishift+1]*log(Sc)
    }
    ls_txt<-"Inverse-Power"
    params_txt<-c("a","b")
  }

  if (ls=="InversePower2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # First redefine b parameter as logb
    LSQest[ishift+2] <- log(LSQest[ishift+2])
    # positivity_v[ishift+2]<-1

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      (S2/S1)^theta[ishift+1]
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      (1/exp(theta[ishift+2]))*((S2^-theta[ishift+1])-(S1^-theta[ishift+1]))
    }

    # Life functions
    lifeF <- function(theta) {
      1/(exp(theta[ishift+2])*(SF^theta[ishift+1]))
    }
    loglifeF <- function(theta) {
      -(theta[ishift+2]) - theta[ishift+1]*log(SF)
    }
    lifeC <- function(theta) {
      1/(exp(theta[ishift+2])*(Sc^theta[ishift+1]))
    }
    loglifeC <- function(theta) {
      -theta[ishift+2] - theta[ishift+1]*log(Sc)
    }
    ls_txt<-"Inverse-Power"
    params_txt<-c("a","b")
  }

  if (ls=="InversePower3") {
    # NOTE 12/19/23: This particular example sets parameter a at 1/0.09 and solves for b
    # I may want to have an option where you can hold input a known parameter and just solve
    # for the other one.  ALTA/Weibull++ apparently has that option.
    # parameter a given as 1/0.09, lsparams[1] - parameter b
    # First redefine b parameter as logb
    LSQest[ishift+1] <- log(LSQest[ishift+1])
    # positivity_v[ishift+2]<-1

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      (S2/S1)^(1/0.09)
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      exp(theta[ishift+1])*((S2^-(1/0.09))-(S1^-(1/0.09)))
    }

    # Life functions
    lifeF <- function(theta) {
      exp(theta[ishift+1])*(SF^-(1/0.09))
    }
    loglifeF <- function(theta) {
      theta[ishift+1] - (1/0.09)*log(SF)
    }
    lifeC <- function(theta) {
      exp(theta[ishift+1])*(Sc^-(1/0.09))
    }
    loglifeC <- function(theta) {
      theta[ishift+1] - (1/0.09)*log(Sc)
    }
    ls_txt<-"Inverse-Power"
    params_txt<-c("b")
  }

  if (ls=="Logarithmic") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      (theta[ishift+2] + log(S1)*theta[ishift+1])/(theta[ishift+2] + log(S2)*theta[ishift+1])
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      theta[ishift+1]*(log(S2) - log(S1))
    }

    # Life functions
    lifeF <- function(theta) {
      theta[ishift+1]*log(SF) + theta[ishift+2]
    }
    loglifeF <- function(theta) {
      log(theta[ishift+1]*log(SF) + theta[ishift+2])
    }
    lifeC <- function(theta) {
      theta[ishift+1]*log(Sc) + theta[ishift+2]
    }
    loglifeC <- function(theta) {
      log(theta[ishift+1]*log(Sc) + theta[ishift+2])
    }
    ls_txt<-ls
    params_txt<-c("a","b")
  }

  if (ls=="MultiStress") {
    # lsparams - c(ao,a1,...,an), S - c(1,S1,...,Sn)

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      exp((S1 - S2)%*%theta[ishift+2:length(LSQest)])
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      exp(c(1,S2)%*%theta[ishift+1:length(LSQest)]) - exp(c(1,S1)%*%theta[ishift+1:length(LSQest)])
    }

    # Life functions
    lifeF <- function(theta) {
      if(dim(SF)[2]==2){
        eqn1<-exp(theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2])
      }
      if(dim(SF)[2]==3){
        eqn1<-exp(theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3])
      }
      if(dim(SF)[2]==4){
        eqn1<-exp(theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4])
      }
      if(dim(SF)[2]==5){
        eqn1<-exp(theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+4]*SF[,5])
      }
      if(dim(SF)[2]==6){
        eqn1<-exp(theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+4]*SF[,5]+theta[ishift+5]*SF[,6])
      }
      if(dim(SF)[2]==7){
        eqn1<-exp(theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+4]*SF[,5]+theta[ishift+5]*SF[,6]+theta[ishift+6]*SF[,7])
      }
      if(dim(SF)[2]==8){
        eqn1<-exp(theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+4]*SF[,5]+theta[ishift+5]*SF[,6]+theta[ishift+6]*SF[,7]+theta[ishift+7]*SF[,8])
      }
      if(dim(SF)[2]==9){
        eqn1<-exp(theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+4]*SF[,5]+theta[ishift+5]*SF[,6]+theta[ishift+6]*SF[,7]+theta[ishift+7]*SF[,8]+theta[ishift+8]*SF[,9])
      }
      if(dim(SF)[2]==10){
        eqn1<-exp(theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+4]*SF[,5]+theta[ishift+5]*SF[,6]+theta[ishift+6]*SF[,7]+theta[ishift+7]*SF[,8]+theta[ishift+8]*SF[,9]+theta[ishift+9]*SF[,10])
      }
      return(eqn1)
    }
    loglifeF <- function(theta) {
      if(dim(SF)[2]==2){
        eqn2<-theta[ishift+1]*rep(1,length(SF[,1]))+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]
      }
      if(dim(SF)[2]==3){
        eqn2<-theta[ishift+1]*rep(1,length(SF[,1]))+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]
      }
      if(dim(SF)[2]==4){
        eqn2<-theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]
      }
      if(dim(SF)[2]==5){
        eqn2<-theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+4]*SF[,5]
      }
      if(dim(SF)[2]==6){
        eqn2<-theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+4]*SF[,5]+theta[ishift+5]*SF[,6]
      }
      if(dim(SF)[2]==7){
        eqn2<-theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+4]*SF[,5]+theta[ishift+5]*SF[,6]+theta[ishift+6]*SF[,7]
      }
      if(dim(SF)[2]==8){
        eqn2<-theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+4]*SF[,5]+theta[ishift+5]*SF[,6]+theta[ishift+6]*SF[,7]+theta[ishift+7]*SF[,8]
      }
      if(dim(SF)[2]==9){
        eqn2<-theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+4]*SF[,5]+theta[ishift+5]*SF[,6]+theta[ishift+6]*SF[,7]+theta[ishift+7]*SF[,8]+theta[ishift+8]*SF[,9]
      }
      if(dim(SF)[2]==10){
        eqn2<-theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+4]*SF[,5]+theta[ishift+5]*SF[,6]+theta[ishift+6]*SF[,7]+theta[ishift+7]*SF[,8]+theta[ishift+8]*SF[,9]+theta[ishift+9]*SF[,10]
      }
      return(eqn2)
    }
    if(missing(Tc)==FALSE){
      lifeC <- function(theta) {
        if(dim(Sc)[2]==2){
          eqn3<-exp(theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2])
        }
        if(dim(Sc)[2]==3){
          eqn3<-exp(theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3])
        }
        if(dim(Sc)[2]==4){
          eqn3<-exp(theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4])
        }
        if(dim(Sc)[2]==5){
          eqn3<-exp(theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+4]*Sc[,5])
        }
        if(dim(Sc)[2]==6){
          eqn3<-exp(theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+4]*Sc[,5]+theta[ishift+5]*Sc[,6])
        }
        if(dim(Sc)[2]==7){
          eqn3<-exp(theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+4]*Sc[,5]+theta[ishift+5]*Sc[,6]+theta[ishift+6]*Sc[,7])
        }
        if(dim(Sc)[2]==8){
          eqn3<-exp(theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+4]*Sc[,5]+theta[ishift+5]*Sc[,6]+theta[ishift+6]*Sc[,7]+theta[ishift+7]*Sc[,8])
        }
        if(dim(Sc)[2]==9){
          eqn3<-exp(theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+4]*Sc[,5]+theta[ishift+5]*Sc[,6]+theta[ishift+6]*Sc[,7]+theta[ishift+7]*Sc[,8]+theta[ishift+8]*Sc[,9])
        }
        if(dim(Sc)[2]==10){
          eqn3<-exp(theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+4]*Sc[,5]+theta[ishift+5]*Sc[,6]+theta[ishift+6]*Sc[,7]+theta[ishift+7]*Sc[,8]+theta[ishift+8]*Sc[,9]+theta[ishift+9]*Sc[,10])
        }
        return(eqn3)
      }
      loglifeC <- function(theta) {
        if(dim(Sc)[2]==2){
          eqn4<-theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]
        }
        if(dim(Sc)[2]==3){
          eqn4<-theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]
        }
        if(dim(Sc)[2]==4){
          eqn4<-theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]
        }
        if(dim(Sc)[2]==5){
          eqn4<-theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+4]*Sc[,5]
        }
        if(dim(Sc)[2]==6){
          eqn4<-theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+4]*Sc[,5]+theta[ishift+5]*Sc[,6]
        }
        if(dim(Sc)[2]==7){
          eqn4<-theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+4]*Sc[,5]+theta[ishift+5]*Sc[,6]+theta[ishift+6]*Sc[,7]
        }
        if(dim(Sc)[2]==8){
          eqn4<-theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+4]*Sc[,5]+theta[ishift+5]*Sc[,6]+theta[ishift+6]*Sc[,7]+theta[ishift+7]*Sc[,8]
        }
        if(dim(Sc)[2]==9){
          eqn4<-theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+4]*Sc[,5]+theta[ishift+5]*Sc[,6]+theta[ishift+6]*Sc[,7]+theta[ishift+7]*Sc[,8]+theta[ishift+8]*Sc[,9]
        }
        if(dim(Sc)[2]==10){
          eqn4<-theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+4]*Sc[,5]+theta[ishift+5]*Sc[,6]+theta[ishift+6]*Sc[,7]+theta[ishift+7]*Sc[,8]+theta[ishift+8]*Sc[,9]+theta[ishift+9]*Sc[,10]
        }
        return(eqn4)
      }
    }
    ls_txt<-ls
    params_txt<-c("a","b")
  }

  if (ls=="TempHumidity") {
    # lsparams[1] - parameter A, lsparams[2] - parameter a, lsparams[3] - parameter b
    # First redefine A parameter as logA
    LSQest[ishift+1] <- log(LSQest[ishift+1])
    # positivity_v[ishift+1]<-1

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      exp((theta[ishift+2]*((1/S1[1]) - (1/S2[1]))) + (theta[ishift+3]*((1/S1[2]) - (1/S2[2]))))
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      exp(theta[ishift+1])*((exp((theta[ishift+2]/S2[1])- (theta[ishift+3]/S2[2]))) - (exp((theta[ishift+2]/S1[1])- (theta[ishift+3]/S1[2]))))
    }

    # Life functions
    lifeF <- function(theta) {
      exp(theta[ishift+1])*exp((theta[ishift+2]/SF[,1]) + (theta[ishift+3]/SF[,2]))
    }
    loglifeF <- function(theta) {
      theta[ishift+1] + (theta[ishift+2]/SF[,1]) + (theta[ishift+3]/SF[,2])
    }
    lifeC <- function(theta) {
      exp(theta[ishift+1])*exp((theta[ishift+2]/Sc[,1]) + (theta[ishift+3]/Sc[,2]))
    }
    loglifeC <- function(theta) {
      theta[ishift+1] + (theta[ishift+2]/Sc[,1]) + (theta[ishift+3]/Sc[,2])
    }
    ls_txt<-"Temperature-Humidity"
    params_txt<-c("A","a","b")
  }

  if (ls=="TempNonthermal") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
    # First redefine c parameter as logc
    LSQest[ishift+3] <- log(LSQest[ishift+3])

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      ((S2[1]/S1[1])^theta[ishift+2])*exp(-theta[ishift+1]*((1/S1[2]) - (1/S2[2])))
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      exp(theta[ishift+3])*((1/((S2[1]^theta[ishift+2])*exp(-theta[ishift+1]/S2[2]))) - (1/((S1[1]^theta[ishift+2])*exp(-theta[ishift+1]/S1[2]))))
    }

    # Life functions
    lifeF <- function(theta) {
      exp(theta[ishift+3])/((SF[,2]^theta[ishift+2])*exp(-theta[ishift+1]/SF[,1]))
    }
    loglifeF <- function(theta) {
      theta[ishift+3] - theta[ishift+2]*log(SF[,2]) + (theta[ishift+1]/SF[,1])
    }
    lifeC <- function(theta) {
      exp(theta[ishift+3])/((Sc[,2]^theta[ishift+2])*exp(-theta[ishift+1]/Sc[,1]))
    }
    loglifeC <- function(theta) {
      theta[ishift+3] - theta[ishift+2]*log(Sc[,2]) + (theta[ishift+1]/Sc[,1])
    }

    ls_txt<-"Temperature-Non-thermal"
    params_txt<-c("a","b","c")
  }

  if (ls=="Eyring3") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # lsparams[3] - parameter c, lsparams[4] - parameter d

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      (S2[1]/S1[1])*exp(theta[ishift+2]*((1/S1[1])-(1/S2[1])) + theta[ishift+3]*(S1[2] - S2[2]) + theta[ishift+4]*((S1[2]/S1[1]) - (S2[2]/S2[1])))
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      theta[ishift+1]*((exp((theta[ishift+2]/S2[1])- (theta[ishift+3]/S2[2]))) - (exp((theta[ishift+2]/S1[1])- (theta[ishift+3]/S1[2]))))
    }

    # Life functions
    lifeF <- function(theta) {
      (1/SF[,1])*exp((theta[ishift+1] + (theta[ishift+2]/SF[,1])) + (theta[ishift+3] + (theta[ishift+4]/SF[,1]))*SF[,2])
    }
    loglifeF <- function(theta) {
      -log(SF[,1]) + theta[ishift+1] + (theta[ishift+2]/SF[,1]) + (theta[ishift+3] + (theta[ishift+4]/SF[,1]))*SF[,2]
    }
    lifeC <- function(theta) {
      (1/Sc[,1])*exp((theta[ishift+1] + (theta[ishift+2]/Sc[,1])) + (theta[ishift+3] + (theta[ishift+4]/Sc[,1]))*Sc[,2])
    }
    loglifeC <- function(theta) {
      -log(Sc[,1]) + theta[ishift+1] + (theta[ishift+2]/Sc[,1]) + (theta[ishift+3] + (theta[ishift+4]/Sc[,1]))*Sc[,2]
    }
    ls_txt<-ls
    params_txt<-c("a","b","c","d")
  }

  if (ls=="PH1") {
    # lsparams[1] - parameter beta0, lsparams[2] - parameter beta1, lsparams[3] - parameter beta2
    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      exp((theta[ishift+2]*((1/S2[1]) - (1/S1[1]))) + (theta[ishift+3]*((1/S2[2]) - (1/S1[2]))))
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      exp(-theta[ishift+1])*((exp(-(theta[ishift+2]/S2[1]) - (theta[ishift+3]/S2[2]))) - (exp(-(theta[ishift+2]/S1[1]) - (theta[ishift+3]/S1[2]))))
    }

    # Life functions
    lifeF <- function(theta) {
      exp(-theta[ishift+1])*exp(-(theta[ishift+2]/SF[,1]) - (theta[ishift+3]/SF[,2]))
    }
    loglifeF <- function(theta) {
      -theta[ishift+1] - (theta[ishift+2]/SF[,1]) - (theta[ishift+3]/SF[,2])
    }
    lifeC <- function(theta) {
      exp(-theta[ishift+1])*exp(-(theta[ishift+2]/Sc[,1]) - (theta[ishift+3]/Sc[,2]))
    }
    loglifeC <- function(theta) {
      -theta[ishift+1] - (theta[ishift+2]/Sc[,1]) - (theta[ishift+3]/Sc[,2])
    }
    ls_txt<-ls
    params_txt<-c("\U03B2\U2080","\U03B2\U2081","\U03B2\U2082")
  }
  # Initialize tau_i and tau_j functions and the step numbers (taui_steps and tauj_steps) to look for
  if(is.list(data[[1]])==FALSE){
    # FOR ONE STEP-STRESS TEST
    taui_steps<-unique(stpstrdatsort[[6]])
    tauj_steps<-unique(stpstrdatsort[[7]])
  } else{
    # FOR SET OF STEP-STRESS TESTS TO BE COMPUTED TOGETHER
    taui_steps<-vector(mode = "list", length(data))
    tauj_steps<-vector(mode = "list", length(data))
    for(i in 1:length(data)){
      taui_steps[[i]]<-unique(stpstrdatsort[[i]][[6]])
      tauj_steps[[i]]<-unique(stpstrdatsort[[i]][[7]])
    }
  }

  # return(list(xirctot,stpstrdatsort,tF,tF_end,SF,tc,tc_end,Sc))


  if(is.list(data[[1]])==FALSE){
    # FOR ONE STEP-STRESS TEST
    # ===================================
    # Fit to log-likelihood distributions
    if (dist=="Weibull") {
      positivity_v[1]<-1

      # Set up the tau functions for steps 2 through 10
      tau1 <- function(theta){
        t_end[1]*(1/L1dL2(theta,as.numeric(stepstresstable[1,1:Nstress]),as.numeric(stepstresstable[2,1:Nstress])))
      }
      tau2 <- function(theta){
        (t_end[2]-t_end[1]+tau1(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[2,1:Nstress]),as.numeric(stepstresstable[3,1:Nstress])))
      }
      tau3 <- function(theta){
        (t_end[3]-t_end[2]+tau2(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[3,1:Nstress]),as.numeric(stepstresstable[4,1:Nstress])))
      }
      tau4 <- function(theta){
        (t_end[4]-t_end[3]+tau3(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[4,1:Nstress]),as.numeric(stepstresstable[5,1:Nstress])))
      }
      tau5 <- function(theta){
        (t_end[5]-t_end[4]+tau4(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[5,1:Nstress]),as.numeric(stepstresstable[6,1:Nstress])))
      }
      tau6 <- function(theta){
        (t_end[6]-t_end[5]+tau5(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[6,1:Nstress]),as.numeric(stepstresstable[7,1:Nstress])))
      }
      tau7 <- function(theta){
        (t_end[7]-t_end[6]+tau6(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[7,1:Nstress]),as.numeric(stepstresstable[8,1:Nstress])))
      }
      tau8 <- function(theta){
        (t_end[8]-t_end[7]+tau7(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[8,1:Nstress]),as.numeric(stepstresstable[9,1:Nstress])))
      }
      tau9 <- function(theta){
        (t_end[9]-t_end[8]+tau8(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[9,1:Nstress]),as.numeric(stepstresstable[10,1:Nstress])))
      }
      tau10 <- function(theta){
        (t_end[10]-t_end[9]+tau9(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[10,1:Nstress]),as.numeric(stepstresstable[11,1:Nstress])))
      }
      tau11 <- function(theta){
        (t_end[11]-t_end[10]+tau10(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[11,1:Nstress]),as.numeric(stepstresstable[12,1:Nstress])))
      }
      tau12 <- function(theta){
        (t_end[12]-t_end[11]+tau11(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[12,1:Nstress]),as.numeric(stepstresstable[13,1:Nstress])))
      }
      tau13 <- function(theta){
        (t_end[13]-t_end[12]+tau12(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[13,1:Nstress]),as.numeric(stepstresstable[14,1:Nstress])))
      }
      tau14 <- function(theta){
        (t_end[14]-t_end[13]+tau13(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[14,1:Nstress]),as.numeric(stepstresstable[15,1:Nstress])))
      }
      tau15 <- function(theta){
        (t_end[15]-t_end[14]+tau14(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[15,1:Nstress]),as.numeric(stepstresstable[16,1:Nstress])))
      }
      tau16 <- function(theta){
        (t_end[16]-t_end[15]+tau15(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[16,1:Nstress]),as.numeric(stepstresstable[17,1:Nstress])))
      }
      tau17 <- function(theta){
        (t_end[17]-t_end[16]+tau16(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[17,1:Nstress]),as.numeric(stepstresstable[18,1:Nstress])))
      }
      tau18 <- function(theta){
        (t_end[18]-t_end[17]+tau17(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[18,1:Nstress]),as.numeric(stepstresstable[19,1:Nstress])))
      }
      tau19 <- function(theta){
        (t_end[19]-t_end[18]+tau18(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[19,1:Nstress]),as.numeric(stepstresstable[20,1:Nstress])))
      }
      tau20 <- function(theta){
        (t_end[20]-t_end[19]+tau19(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[20,1:Nstress]),as.numeric(stepstresstable[21,1:Nstress])))
      }
    }
    if (dist=="Lognormal") {
      # Set up the tau functions for steps 2 through 10
      tau1 <- function(theta){
        t_end[1]*(1/L1dL2(theta,as.numeric(stepstresstable[1,1:Nstress]),as.numeric(stepstresstable[2,1:Nstress])))
      }
      tau2 <- function(theta){
        (t_end[2]-t_end[1]+tau1(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[2,1:Nstress]),as.numeric(stepstresstable[3,1:Nstress])))
      }
      tau3 <- function(theta){
        (t_end[3]-t_end[2]+tau2(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[3,1:Nstress]),as.numeric(stepstresstable[4,1:Nstress])))
      }
      tau4 <- function(theta){
        (t_end[4]-t_end[3]+tau3(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[4,1:Nstress]),as.numeric(stepstresstable[5,1:Nstress])))
      }
      tau5 <- function(theta){
        (t_end[5]-t_end[4]+tau4(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[5,1:Nstress]),as.numeric(stepstresstable[6,1:Nstress])))
      }
      tau6 <- function(theta){
        (t_end[6]-t_end[5]+tau5(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[6,1:Nstress]),as.numeric(stepstresstable[7,1:Nstress])))
      }
      tau7 <- function(theta){
        (t_end[7]-t_end[6]+tau6(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[7,1:Nstress]),as.numeric(stepstresstable[8,1:Nstress])))
      }
      tau8 <- function(theta){
        (t_end[8]-t_end[7]+tau7(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[8,1:Nstress]),as.numeric(stepstresstable[9,1:Nstress])))
      }
      tau9 <- function(theta){
        (t_end[9]-t_end[8]+tau8(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[9,1:Nstress]),as.numeric(stepstresstable[10,1:Nstress])))
      }
      tau10 <- function(theta){
        (t_end[10]-t_end[9]+tau9(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[10,1:Nstress]),as.numeric(stepstresstable[11,1:Nstress])))
      }
      tau11 <- function(theta){
        (t_end[11]-t_end[10]+tau10(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[11,1:Nstress]),as.numeric(stepstresstable[12,1:Nstress])))
      }
      tau12 <- function(theta){
        (t_end[12]-t_end[11]+tau11(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[12,1:Nstress]),as.numeric(stepstresstable[13,1:Nstress])))
      }
      tau13 <- function(theta){
        (t_end[13]-t_end[12]+tau12(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[13,1:Nstress]),as.numeric(stepstresstable[14,1:Nstress])))
      }
      tau14 <- function(theta){
        (t_end[14]-t_end[13]+tau13(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[14,1:Nstress]),as.numeric(stepstresstable[15,1:Nstress])))
      }
      tau15 <- function(theta){
        (t_end[15]-t_end[14]+tau14(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[15,1:Nstress]),as.numeric(stepstresstable[16,1:Nstress])))
      }
      tau16 <- function(theta){
        (t_end[16]-t_end[15]+tau15(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[16,1:Nstress]),as.numeric(stepstresstable[17,1:Nstress])))
      }
      tau17 <- function(theta){
        (t_end[17]-t_end[16]+tau16(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[17,1:Nstress]),as.numeric(stepstresstable[18,1:Nstress])))
      }
      tau18 <- function(theta){
        (t_end[18]-t_end[17]+tau17(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[18,1:Nstress]),as.numeric(stepstresstable[19,1:Nstress])))
      }
      tau19 <- function(theta){
        (t_end[19]-t_end[18]+tau18(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[19,1:Nstress]),as.numeric(stepstresstable[20,1:Nstress])))
      }
      tau20 <- function(theta){
        (t_end[20]-t_end[19]+tau19(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[20,1:Nstress]),as.numeric(stepstresstable[21,1:Nstress])))
      }
    }
    if (dist=="Normal") {
      # Set up the tau functions for steps 2 through 10
      tau1 <- function(theta){
        t_end[1] + L2mL1(theta,as.numeric(stepstresstable[1,1:Nstress]),as.numeric(stepstresstable[2,1:Nstress]))
      }
      tau2 <- function(theta){
        t_end[2]-t_end[1]+tau1(theta) + L2mL1(theta,as.numeric(stepstresstable[2,1:Nstress]),as.numeric(stepstresstable[3,1:Nstress]))
      }
      tau3 <- function(theta){
        t_end[3]-t_end[2]+tau2(theta) + L2mL1(theta,as.numeric(stepstresstable[3,1:Nstress]),as.numeric(stepstresstable[4,1:Nstress]))
      }
      tau4 <- function(theta){
        t_end[4]-t_end[3]+tau3(theta) + L2mL1(theta,as.numeric(stepstresstable[4,1:Nstress]),as.numeric(stepstresstable[5,1:Nstress]))
      }
      tau5 <- function(theta){
        t_end[5]-t_end[4]+tau4(theta) + L2mL1(theta,as.numeric(stepstresstable[5,1:Nstress]),as.numeric(stepstresstable[6,1:Nstress]))
      }
      tau6 <- function(theta){
        t_end[6]-t_end[5]+tau5(theta) + L2mL1(theta,as.numeric(stepstresstable[6,1:Nstress]),as.numeric(stepstresstable[7,1:Nstress]))
      }
      tau7 <- function(theta){
        t_end[7]-t_end[6]+tau6(theta) + L2mL1(theta,as.numeric(stepstresstable[7,1:Nstress]),as.numeric(stepstresstable[8,1:Nstress]))
      }
      tau8 <- function(theta){
        t_end[8]-t_end[7]+tau7(theta) + L2mL1(theta,as.numeric(stepstresstable[8,1:Nstress]),as.numeric(stepstresstable[9,1:Nstress]))
      }
      tau9 <- function(theta){
        t_end[9]-t_end[8]+tau8(theta) + L2mL1(theta,as.numeric(stepstresstable[9,1:Nstress]),as.numeric(stepstresstable[10,1:Nstress]))
      }
      tau10 <- function(theta){
        t_end[10]-t_end[9]+tau9(theta) + L2mL1(theta,as.numeric(stepstresstable[10,1:Nstress]),as.numeric(stepstresstable[11,1:Nstress]))
      }
      tau11 <- function(theta){
        t_end[11]-t_end[10]+tau10(theta) + L2mL1(theta,as.numeric(stepstresstable[11,1:Nstress]),as.numeric(stepstresstable[12,1:Nstress]))
      }
      tau12 <- function(theta){
        t_end[12]-t_end[11]+tau11(theta) + L2mL1(theta,as.numeric(stepstresstable[12,1:Nstress]),as.numeric(stepstresstable[13,1:Nstress]))
      }
      tau13 <- function(theta){
        t_end[13]-t_end[12]+tau12(theta) + L2mL1(theta,as.numeric(stepstresstable[13,1:Nstress]),as.numeric(stepstresstable[14,1:Nstress]))
      }
      tau14 <- function(theta){
        t_end[14]-t_end[13]+tau13(theta) + L2mL1(theta,as.numeric(stepstresstable[14,1:Nstress]),as.numeric(stepstresstable[15,1:Nstress]))
      }
      tau15 <- function(theta){
        t_end[15]-t_end[14]+tau14(theta) + L2mL1(theta,as.numeric(stepstresstable[15,1:Nstress]),as.numeric(stepstresstable[16,1:Nstress]))
      }
      tau16 <- function(theta){
        t_end[16]-t_end[15]+tau15(theta) + L2mL1(theta,as.numeric(stepstresstable[16,1:Nstress]),as.numeric(stepstresstable[17,1:Nstress]))
      }
      tau17 <- function(theta){
        t_end[17]-t_end[16]+tau16(theta) + L2mL1(theta,as.numeric(stepstresstable[17,1:Nstress]),as.numeric(stepstresstable[18,1:Nstress]))
      }
      tau18 <- function(theta){
        t_end[18]-t_end[17]+tau17(theta) + L2mL1(theta,as.numeric(stepstresstable[18,1:Nstress]),as.numeric(stepstresstable[19,1:Nstress]))
      }
      tau19 <- function(theta){
        t_end[19]-t_end[18]+tau18(theta) + L2mL1(theta,as.numeric(stepstresstable[19,1:Nstress]),as.numeric(stepstresstable[20,1:Nstress]))
      }
      tau20 <- function(theta){
        (t_end[20]-t_end[19]+tau19(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[20,1:Nstress]),as.numeric(stepstresstable[21,1:Nstress])))
      }
    }
    if (dist=="Exponential") {
      # Set up the tau functions for steps 2 through 10
      tau1 <- function(theta){
        t_end[1]*(1/L1dL2(theta,as.numeric(stepstresstable[1,1:Nstress]),as.numeric(stepstresstable[2,1:Nstress])))
      }
      tau2 <- function(theta){
        (t_end[2]-t_end[1]+tau1(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[2,1:Nstress]),as.numeric(stepstresstable[3,1:Nstress])))
      }
      tau3 <- function(theta){
        (t_end[3]-t_end[2]+tau2(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[3,1:Nstress]),as.numeric(stepstresstable[4,1:Nstress])))
      }
      tau4 <- function(theta){
        (t_end[4]-t_end[3]+tau3(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[4,1:Nstress]),as.numeric(stepstresstable[5,1:Nstress])))
      }
      tau5 <- function(theta){
        (t_end[5]-t_end[4]+tau4(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[5,1:Nstress]),as.numeric(stepstresstable[6,1:Nstress])))
      }
      tau6 <- function(theta){
        (t_end[6]-t_end[5]+tau5(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[6,1:Nstress]),as.numeric(stepstresstable[7,1:Nstress])))
      }
      tau7 <- function(theta){
        (t_end[7]-t_end[6]+tau6(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[7,1:Nstress]),as.numeric(stepstresstable[8,1:Nstress])))
      }
      tau8 <- function(theta){
        (t_end[8]-t_end[7]+tau7(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[8,1:Nstress]),as.numeric(stepstresstable[9,1:Nstress])))
      }
      tau9 <- function(theta){
        (t_end[9]-t_end[8]+tau8(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[9,1:Nstress]),as.numeric(stepstresstable[10,1:Nstress])))
      }
      tau10 <- function(theta){
        (t_end[10]-t_end[9]+tau9(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[10,1:Nstress]),as.numeric(stepstresstable[11,1:Nstress])))
      }
      tau11 <- function(theta){
        (t_end[11]-t_end[10]+tau10(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[11,1:Nstress]),as.numeric(stepstresstable[12,1:Nstress])))
      }
      tau12 <- function(theta){
        (t_end[12]-t_end[11]+tau11(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[12,1:Nstress]),as.numeric(stepstresstable[13,1:Nstress])))
      }
      tau13 <- function(theta){
        (t_end[13]-t_end[12]+tau12(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[13,1:Nstress]),as.numeric(stepstresstable[14,1:Nstress])))
      }
      tau14 <- function(theta){
        (t_end[14]-t_end[13]+tau13(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[14,1:Nstress]),as.numeric(stepstresstable[15,1:Nstress])))
      }
      tau15 <- function(theta){
        (t_end[15]-t_end[14]+tau14(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[15,1:Nstress]),as.numeric(stepstresstable[16,1:Nstress])))
      }
      tau16 <- function(theta){
        (t_end[16]-t_end[15]+tau15(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[16,1:Nstress]),as.numeric(stepstresstable[17,1:Nstress])))
      }
      tau17 <- function(theta){
        (t_end[17]-t_end[16]+tau16(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[17,1:Nstress]),as.numeric(stepstresstable[18,1:Nstress])))
      }
      tau18 <- function(theta){
        (t_end[18]-t_end[17]+tau17(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[18,1:Nstress]),as.numeric(stepstresstable[19,1:Nstress])))
      }
      tau19 <- function(theta){
        (t_end[19]-t_end[18]+tau18(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[19,1:Nstress]),as.numeric(stepstresstable[20,1:Nstress])))
      }
      tau20 <- function(theta){
        (t_end[20]-t_end[19]+tau19(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[20,1:Nstress]),as.numeric(stepstresstable[21,1:Nstress])))
      }
    }
    if (dist=="2PExponential") {
      # Set up the tau functions for steps 2 through 10
      tau1 <- function(theta){
        t_end[1] + L2mL1(theta,as.numeric(stepstresstable[1,1:Nstress]),as.numeric(stepstresstable[2,1:Nstress]))
      }
      tau2 <- function(theta){
        t_end[2]-t_end[1]+tau1(theta) + L2mL1(theta,as.numeric(stepstresstable[2,1:Nstress]),as.numeric(stepstresstable[3,1:Nstress]))
      }
      tau3 <- function(theta){
        t_end[3]-t_end[2]+tau2(theta) + L2mL1(theta,as.numeric(stepstresstable[3,1:Nstress]),as.numeric(stepstresstable[4,1:Nstress]))
      }
      tau4 <- function(theta){
        t_end[4]-t_end[3]+tau3(theta) + L2mL1(theta,as.numeric(stepstresstable[4,1:Nstress]),as.numeric(stepstresstable[5,1:Nstress]))
      }
      tau5 <- function(theta){
        t_end[5]-t_end[4]+tau4(theta) + L2mL1(theta,as.numeric(stepstresstable[5,1:Nstress]),as.numeric(stepstresstable[6,1:Nstress]))
      }
      tau6 <- function(theta){
        t_end[6]-t_end[5]+tau5(theta) + L2mL1(theta,as.numeric(stepstresstable[6,1:Nstress]),as.numeric(stepstresstable[7,1:Nstress]))
      }
      tau7 <- function(theta){
        t_end[7]-t_end[6]+tau6(theta) + L2mL1(theta,as.numeric(stepstresstable[7,1:Nstress]),as.numeric(stepstresstable[8,1:Nstress]))
      }
      tau8 <- function(theta){
        t_end[8]-t_end[7]+tau7(theta) + L2mL1(theta,as.numeric(stepstresstable[8,1:Nstress]),as.numeric(stepstresstable[9,1:Nstress]))
      }
      tau9 <- function(theta){
        t_end[9]-t_end[8]+tau8(theta) + L2mL1(theta,as.numeric(stepstresstable[9,1:Nstress]),as.numeric(stepstresstable[10,1:Nstress]))
      }
      tau10 <- function(theta){
        t_end[10]-t_end[9]+tau9(theta) + L2mL1(theta,as.numeric(stepstresstable[10,1:Nstress]),as.numeric(stepstresstable[11,1:Nstress]))
      }
      tau11 <- function(theta){
        t_end[11]-t_end[10]+tau10(theta) + L2mL1(theta,as.numeric(stepstresstable[11,1:Nstress]),as.numeric(stepstresstable[12,1:Nstress]))
      }
      tau12 <- function(theta){
        t_end[12]-t_end[11]+tau11(theta) + L2mL1(theta,as.numeric(stepstresstable[12,1:Nstress]),as.numeric(stepstresstable[13,1:Nstress]))
      }
      tau13 <- function(theta){
        t_end[13]-t_end[12]+tau12(theta) + L2mL1(theta,as.numeric(stepstresstable[13,1:Nstress]),as.numeric(stepstresstable[14,1:Nstress]))
      }
      tau14 <- function(theta){
        t_end[14]-t_end[13]+tau13(theta) + L2mL1(theta,as.numeric(stepstresstable[14,1:Nstress]),as.numeric(stepstresstable[15,1:Nstress]))
      }
      tau15 <- function(theta){
        t_end[15]-t_end[14]+tau14(theta) + L2mL1(theta,as.numeric(stepstresstable[15,1:Nstress]),as.numeric(stepstresstable[16,1:Nstress]))
      }
      tau16 <- function(theta){
        t_end[16]-t_end[15]+tau15(theta) + L2mL1(theta,as.numeric(stepstresstable[16,1:Nstress]),as.numeric(stepstresstable[17,1:Nstress]))
      }
      tau17 <- function(theta){
        t_end[17]-t_end[16]+tau16(theta) + L2mL1(theta,as.numeric(stepstresstable[17,1:Nstress]),as.numeric(stepstresstable[18,1:Nstress]))
      }
      tau18 <- function(theta){
        t_end[18]-t_end[17]+tau17(theta) + L2mL1(theta,as.numeric(stepstresstable[18,1:Nstress]),as.numeric(stepstresstable[19,1:Nstress]))
      }
      tau19 <- function(theta){
        t_end[19]-t_end[18]+tau18(theta) + L2mL1(theta,as.numeric(stepstresstable[19,1:Nstress]),as.numeric(stepstresstable[20,1:Nstress]))
      }
      tau20 <- function(theta){
        (t_end[20]-t_end[19]+tau19(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[20,1:Nstress]),as.numeric(stepstresstable[21,1:Nstress])))
      }
    }

    # ===========================================================
    # Initialize main tau_i function for failed step-stress units
    tauifunct<-integer(0)
    taujfunct<-integer(0)

    taui <- function(theta){
      for(i in 1:length(taui_steps)){
        if(length(tauifunct)==0){
          if(taui_steps[i]==1){
            tauifunct<-rep(0,length(which(stpstrdatsort[[6]]==1)))
          }
          if(taui_steps[i]==2){
            tauifunct<-rep(tau1(theta),length(which(stpstrdatsort[[6]]==2)))
          }
          if(taui_steps[i]==3){
            tauifunct<-rep(tau2(theta),length(which(stpstrdatsort[[6]]==3)))
          }
          if(taui_steps[i]==4){
            tauifunct<-rep(tau3(theta),length(which(stpstrdatsort[[6]]==4)))
          }
          if(taui_steps[i]==5){
            tauifunct<-rep(tau4(theta),length(which(stpstrdatsort[[6]]==5)))
          }
          if(taui_steps[i]==6){
            tauifunct<-rep(tau5(theta),length(which(stpstrdatsort[[6]]==6)))
          }
          if(taui_steps[i]==7){
            tauifunct<-rep(tau6(theta),length(which(stpstrdatsort[[6]]==7)))
          }
          if(taui_steps[i]==8){
            tauifunct<-rep(tau7(theta),length(which(stpstrdatsort[[6]]==8)))
          }
          if(taui_steps[i]==9){
            tauifunct<-rep(tau8(theta),length(which(stpstrdatsort[[6]]==9)))
          }
          if(taui_steps[i]==10){
            tauifunct<-rep(tau9(theta),length(which(stpstrdatsort[[6]]==10)))
          }
          if(taui_steps[i]==11){
            tauifunct<-rep(tau10(theta),length(which(stpstrdatsort[[6]]==11)))
          }
          if(taui_steps[i]==12){
            tauifunct<-rep(tau11(theta),length(which(stpstrdatsort[[6]]==12)))
          }
          if(taui_steps[i]==13){
            tauifunct<-rep(tau12(theta),length(which(stpstrdatsort[[6]]==13)))
          }
          if(taui_steps[i]==14){
            tauifunct<-rep(tau13(theta),length(which(stpstrdatsort[[6]]==14)))
          }
          if(taui_steps[i]==15){
            tauifunct<-rep(tau14(theta),length(which(stpstrdatsort[[6]]==15)))
          }
          if(taui_steps[i]==16){
            tauifunct<-rep(tau15(theta),length(which(stpstrdatsort[[6]]==16)))
          }
          if(taui_steps[i]==17){
            tauifunct<-rep(tau16(theta),length(which(stpstrdatsort[[6]]==17)))
          }
          if(taui_steps[i]==18){
            tauifunct<-rep(tau17(theta),length(which(stpstrdatsort[[6]]==18)))
          }
          if(taui_steps[i]==19){
            tauifunct<-rep(tau18(theta),length(which(stpstrdatsort[[6]]==19)))
          }
          if(taui_steps[i]==20){
            tauifunct<-rep(tau19(theta),length(which(stpstrdatsort[[6]]==20)))
          }
          if(taui_steps[i]==21){
            tauifunct<-rep(tau20(theta),length(which(stpstrdatsort[[6]]==21)))
          }
        } else {
          if(taui_steps[i]==1){
            tauifunct<-c(tauifunct,rep(0,length(which(stpstrdatsort[[6]]==1))))
          }
          if(taui_steps[i]==2){
            tauifunct<-c(tauifunct,rep(tau1(theta),length(which(stpstrdatsort[[6]]==2))))
          }
          if(taui_steps[i]==3){
            tauifunct<-c(tauifunct,rep(tau2(theta),length(which(stpstrdatsort[[6]]==3))))
          }
          if(taui_steps[i]==4){
            tauifunct<-c(tauifunct,rep(tau3(theta),length(which(stpstrdatsort[[6]]==4))))
          }
          if(taui_steps[i]==5){
            tauifunct<-c(tauifunct,rep(tau4(theta),length(which(stpstrdatsort[[6]]==5))))
          }
          if(taui_steps[i]==6){
            tauifunct<-c(tauifunct,rep(tau5(theta),length(which(stpstrdatsort[[6]]==6))))
          }
          if(taui_steps[i]==7){
            tauifunct<-c(tauifunct,rep(tau6(theta),length(which(stpstrdatsort[[6]]==7))))
          }
          if(taui_steps[i]==8){
            tauifunct<-c(tauifunct,rep(tau7(theta),length(which(stpstrdatsort[[6]]==8))))
          }
          if(taui_steps[i]==9){
            tauifunct<-c(tauifunct,rep(tau8(theta),length(which(stpstrdatsort[[6]]==9))))
          }
          if(taui_steps[i]==10){
            tauifunct<-c(tauifunct,rep(tau9(theta),length(which(stpstrdatsort[[6]]==10))))
          }
          if(taui_steps[i]==11){
            tauifunct<-c(tauifunct,rep(tau10(theta),length(which(stpstrdatsort[[6]]==11))))
          }
          if(taui_steps[i]==12){
            tauifunct<-c(tauifunct,rep(tau11(theta),length(which(stpstrdatsort[[6]]==12))))
          }
          if(taui_steps[i]==13){
            tauifunct<-c(tauifunct,rep(tau12(theta),length(which(stpstrdatsort[[6]]==13))))
          }
          if(taui_steps[i]==14){
            tauifunct<-c(tauifunct,rep(tau13(theta),length(which(stpstrdatsort[[6]]==14))))
          }
          if(taui_steps[i]==15){
            tauifunct<-c(tauifunct,rep(tau14(theta),length(which(stpstrdatsort[[6]]==15))))
          }
          if(taui_steps[i]==16){
            tauifunct<-c(tauifunct,rep(tau15(theta),length(which(stpstrdatsort[[6]]==16))))
          }
          if(taui_steps[i]==17){
            tauifunct<-c(tauifunct,rep(tau16(theta),length(which(stpstrdatsort[[6]]==17))))
          }
          if(taui_steps[i]==18){
            tauifunct<-c(tauifunct,rep(tau17(theta),length(which(stpstrdatsort[[6]]==18))))
          }
          if(taui_steps[i]==19){
            tauifunct<-c(tauifunct,rep(tau18(theta),length(which(stpstrdatsort[[6]]==19))))
          }
          if(taui_steps[i]==20){
            tauifunct<-c(tauifunct,rep(tau19(theta),length(which(stpstrdatsort[[6]]==20))))
          }
          if(taui_steps[i]==21){
            tauifunct<-c(tauifunct,rep(tau20(theta),length(which(stpstrdatsort[[6]]==21))))
          }
        }
      }
      return(tauifunct)
    }

    # Initialize main tau_j function for censored step-stress units
    tauj <- function(theta){
      for(i in 1:length(tauj_steps)){
        if(length(taujfunct)==0){
          if(tauj_steps[i]==1){
            taujfunct<-rep(0,length(which(stpstrdatsort[[7]]==1)))
          }
          if(tauj_steps[i]==2){
            taujfunct<-rep(tau1(theta),length(which(stpstrdatsort[[7]]==2)))
          }
          if(tauj_steps[i]==3){
            taujfunct<-rep(tau2(theta),length(which(stpstrdatsort[[7]]==3)))
          }
          if(tauj_steps[i]==4){
            taujfunct<-rep(tau3(theta),length(which(stpstrdatsort[[7]]==4)))
          }
          if(tauj_steps[i]==5){
            taujfunct<-rep(tau4(theta),length(which(stpstrdatsort[[7]]==5)))
          }
          if(tauj_steps[i]==6){
            taujfunct<-rep(tau5(theta),length(which(stpstrdatsort[[7]]==6)))
          }
          if(tauj_steps[i]==7){
            taujfunct<-rep(tau6(theta),length(which(stpstrdatsort[[7]]==7)))
          }
          if(tauj_steps[i]==8){
            taujfunct<-rep(tau7(theta),length(which(stpstrdatsort[[7]]==8)))
          }
          if(tauj_steps[i]==9){
            taujfunct<-rep(tau8(theta),length(which(stpstrdatsort[[7]]==9)))
          }
          if(tauj_steps[i]==10){
            taujfunct<-rep(tau9(theta),length(which(stpstrdatsort[[7]]==10)))
          }
          if(tauj_steps[i]==11){
            taujfunct<-rep(tau10(theta),length(which(stpstrdatsort[[7]]==11)))
          }
          if(tauj_steps[i]==12){
            taujfunct<-rep(tau11(theta),length(which(stpstrdatsort[[7]]==12)))
          }
          if(tauj_steps[i]==13){
            taujfunct<-rep(tau12(theta),length(which(stpstrdatsort[[7]]==13)))
          }
          if(tauj_steps[i]==14){
            taujfunct<-rep(tau13(theta),length(which(stpstrdatsort[[7]]==14)))
          }
          if(tauj_steps[i]==15){
            taujfunct<-rep(tau14(theta),length(which(stpstrdatsort[[7]]==15)))
          }
          if(tauj_steps[i]==16){
            taujfunct<-rep(tau15(theta),length(which(stpstrdatsort[[7]]==16)))
          }
          if(tauj_steps[i]==17){
            taujfunct<-rep(tau16(theta),length(which(stpstrdatsort[[7]]==17)))
          }
          if(tauj_steps[i]==18){
            taujfunct<-rep(tau17(theta),length(which(stpstrdatsort[[7]]==18)))
          }
          if(tauj_steps[i]==19){
            taujfunct<-rep(tau18(theta),length(which(stpstrdatsort[[7]]==19)))
          }
          if(tauj_steps[i]==20){
            taujfunct<-rep(tau19(theta),length(which(stpstrdatsort[[7]]==20)))
          }
          if(tauj_steps[i]==21){
            taujfunct<-rep(tau20(theta),length(which(stpstrdatsort[[7]]==21)))
          }
        } else {
          if(tauj_steps[i]==1){
            taujfunct<-c(taujfunct,rep(0,length(which(stpstrdatsort[[7]]==1))))
          }
          if(tauj_steps[i]==2){
            taujfunct<-c(taujfunct,rep(tau1(theta),length(which(stpstrdatsort[[7]]==2))))
          }
          if(tauj_steps[i]==3){
            taujfunct<-c(taujfunct,rep(tau2(theta),length(which(stpstrdatsort[[7]]==3))))
          }
          if(tauj_steps[i]==4){
            taujfunct<-c(taujfunct,rep(tau3(theta),length(which(stpstrdatsort[[7]]==4))))
          }
          if(tauj_steps[i]==5){
            taujfunct<-c(taujfunct,rep(tau4(theta),length(which(stpstrdatsort[[7]]==5))))
          }
          if(tauj_steps[i]==6){
            taujfunct<-c(taujfunct,rep(tau5(theta),length(which(stpstrdatsort[[7]]==6))))
          }
          if(tauj_steps[i]==7){
            taujfunct<-c(taujfunct,rep(tau6(theta),length(which(stpstrdatsort[[7]]==7))))
          }
          if(tauj_steps[i]==8){
            taujfunct<-c(taujfunct,rep(tau7(theta),length(which(stpstrdatsort[[7]]==8))))
          }
          if(tauj_steps[i]==9){
            taujfunct<-c(taujfunct,rep(tau8(theta),length(which(stpstrdatsort[[7]]==9))))
          }
          if(tauj_steps[i]==10){
            taujfunct<-c(taujfunct,rep(tau9(theta),length(which(stpstrdatsort[[7]]==10))))
          }
          if(tauj_steps[i]==11){
            taujfunct<-c(taujfunct,rep(tau10(theta),length(which(stpstrdatsort[[7]]==11))))
          }
          if(tauj_steps[i]==12){
            taujfunct<-c(taujfunct,rep(tau11(theta),length(which(stpstrdatsort[[7]]==12))))
          }
          if(tauj_steps[i]==13){
            taujfunct<-c(taujfunct,rep(tau12(theta),length(which(stpstrdatsort[[7]]==13))))
          }
          if(tauj_steps[i]==14){
            taujfunct<-c(taujfunct,rep(tau13(theta),length(which(stpstrdatsort[[7]]==14))))
          }
          if(tauj_steps[i]==15){
            taujfunct<-c(taujfunct,rep(tau14(theta),length(which(stpstrdatsort[[7]]==15))))
          }
          if(tauj_steps[i]==16){
            taujfunct<-c(taujfunct,rep(tau15(theta),length(which(stpstrdatsort[[7]]==16))))
          }
          if(tauj_steps[i]==17){
            taujfunct<-c(taujfunct,rep(tau16(theta),length(which(stpstrdatsort[[7]]==17))))
          }
          if(tauj_steps[i]==18){
            taujfunct<-c(taujfunct,rep(tau17(theta),length(which(stpstrdatsort[[7]]==18))))
          }
          if(tauj_steps[i]==19){
            taujfunct<-c(taujfunct,rep(tau18(theta),length(which(stpstrdatsort[[7]]==19))))
          }
          if(tauj_steps[i]==20){
            taujfunct<-c(taujfunct,rep(tau9(theta),length(which(stpstrdatsort[[7]]==20))))
          }
          if(tauj_steps[i]==21){
            taujfunct<-c(taujfunct,rep(tau20(theta),length(which(stpstrdatsort[[7]]==21))))
          }
        }
      }
      return(taujfunct)
    }
  } else{
    # FOR SET OF STEP-STRESS TESTS TO BE COMPUTED TOGETHER
    taui <- function(theta){
      # To start with, theta has to be set into each of the lapply equations so we must define EVERYTHING
      # within the taui and tauj functions including the life-stress functions as a function of theta and their
      # respective step-stress groups.  The first test will be strictly with the Inverse Power-Weibull life-stress
      # relation and when it works I will add the other groups accordingly.

      # Re-Initialize life-stress parameter estimates for theta
      if (ls=="InversePower") {
        # lsparams[1] - parameter a, lsparams[2] - parameter b
        # Life ratio L1/L2
        L1dL2 <- function(theta,S1,S2) {
          (S2/S1)^theta[ishift+1]
        }
        # Life difference L2-L1
        L2mL1 <- function(theta,S1,S2) {
          theta[ishift+2]*((S2^-theta[ishift+1])-(S1^-theta[ishift+1]))
        }

        # Life functions
        lifeF <- function(theta) {
          theta[ishift+2]*(SF^-theta[ishift+1])
        }
        loglifeF <- function(theta) {
          log(theta[ishift+2]) - theta[ishift+1]*log(SF)
        }
        lifeC <- function(theta) {
          theta[ishift+2]*(Sc^-theta[ishift+1])
        }
        loglifeC <- function(theta) {
          log(theta[ishift+2]) - theta[ishift+1]*log(Sc)
        }
      }
      # Fit to log-likelihood distributions
      if (dist=="Weibull") {
        positivity_v[1]<-1

        # Set up the tau functions for steps 2 through 10
        tau1 <- function(theta){
          t_end[1]*(1/L1dL2(theta,as.numeric(stepstresstable[1,1:Nstress]),as.numeric(stepstresstable[2,1:Nstress])))
        }
        tau2 <- function(theta){
          (t_end[2]-t_end[1]+tau1(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[2,1:Nstress]),as.numeric(stepstresstable[3,1:Nstress])))
        }
        tau3 <- function(theta){
          (t_end[3]-t_end[2]+tau2(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[3,1:Nstress]),as.numeric(stepstresstable[4,1:Nstress])))
        }
        tau4 <- function(theta){
          (t_end[4]-t_end[3]+tau3(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[4,1:Nstress]),as.numeric(stepstresstable[5,1:Nstress])))
        }
        tau5 <- function(theta){
          (t_end[5]-t_end[4]+tau4(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[5,1:Nstress]),as.numeric(stepstresstable[6,1:Nstress])))
        }
        tau6 <- function(theta){
          (t_end[6]-t_end[5]+tau5(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[6,1:Nstress]),as.numeric(stepstresstable[7,1:Nstress])))
        }
        tau7 <- function(theta){
          (t_end[7]-t_end[6]+tau6(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[7,1:Nstress]),as.numeric(stepstresstable[8,1:Nstress])))
        }
        tau8 <- function(theta){
          (t_end[8]-t_end[7]+tau7(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[8,1:Nstress]),as.numeric(stepstresstable[9,1:Nstress])))
        }
        tau9 <- function(theta){
          (t_end[9]-t_end[8]+tau8(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[9,1:Nstress]),as.numeric(stepstresstable[10,1:Nstress])))
        }
        tau10 <- function(theta){
          (t_end[10]-t_end[9]+tau9(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[10,1:Nstress]),as.numeric(stepstresstable[11,1:Nstress])))
        }
        tau11 <- function(theta){
          (t_end[11]-t_end[10]+tau10(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[11,1:Nstress]),as.numeric(stepstresstable[12,1:Nstress])))
        }
        tau12 <- function(theta){
          (t_end[12]-t_end[11]+tau11(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[12,1:Nstress]),as.numeric(stepstresstable[13,1:Nstress])))
        }
        tau13 <- function(theta){
          (t_end[13]-t_end[12]+tau12(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[13,1:Nstress]),as.numeric(stepstresstable[14,1:Nstress])))
        }
        tau14 <- function(theta){
          (t_end[14]-t_end[13]+tau13(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[14,1:Nstress]),as.numeric(stepstresstable[15,1:Nstress])))
        }
        tau15 <- function(theta){
          (t_end[15]-t_end[14]+tau14(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[15,1:Nstress]),as.numeric(stepstresstable[16,1:Nstress])))
        }
        tau16 <- function(theta){
          (t_end[16]-t_end[15]+tau15(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[16,1:Nstress]),as.numeric(stepstresstable[17,1:Nstress])))
        }
        tau17 <- function(theta){
          (t_end[17]-t_end[16]+tau16(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[17,1:Nstress]),as.numeric(stepstresstable[18,1:Nstress])))
        }
        tau18 <- function(theta){
          (t_end[18]-t_end[17]+tau17(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[18,1:Nstress]),as.numeric(stepstresstable[19,1:Nstress])))
        }
        tau19 <- function(theta){
          (t_end[19]-t_end[18]+tau18(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[19,1:Nstress]),as.numeric(stepstresstable[20,1:Nstress])))
        }
        tau20 <- function(theta){
          (t_end[20]-t_end[19]+tau19(theta))*(1/L1dL2(theta,as.numeric(stepstresstable[20,1:Nstress]),as.numeric(stepstresstable[21,1:Nstress])))
        }
      }
      # Final tau function
      taui_function <- unlist(lapply(stepstressGROUP,taui_vectors))
      return(taui_function)
    }
    tauj <- function(theta){
      # To start with, theta has to be set into each of the lapply equations so we must define EVERYTHING
      # within the taui and tauj functions including the life-stress functions as a function of theta and their
      # respective step-stress groups
    }

  }
  # return(list(xirctot,stpstrdatsort,tF,tF_end,SF,taui,tau1.1,tau1.2,tau1.3,tau1.4,tau1.5,tc,tc_end,Sc,tauj))


  taui_full_funct <- integer(0)
  tauj_full_funct <- integer(0)
  taui_full <- function(theta){
    for(i in 1:length(data)){
      if(length(taui_full_funct)==0){
        taui_full_funct <- get(taui[i])(theta)
      } else{
        taui_full_funct <- c(taui_full_funct,get(taui[i])(theta))
      }
    }
    return(taui_full_funct)
  }
  tauj_full <- function(theta){
    for(i in 1:length(data)){
      if(length(tauj_full_funct)==0){
        tauj_full_funct <- get(tauj[i])(theta)
      } else{
        tauj_full_funct <- c(tauj_full_funct,get(tauj[i])(theta))
      }
    }
    return(tauj_full_funct)
  }

  # Form the failure and censored time functions
  ti <- function(theta) {
    if(is.list(data[[1]])==FALSE){
      # FOR ONE STEP-STRESS TEST
      tF - tF_end + taui(theta)
    } else{
      tF - tF_end + taui_full(theta)
    }
  }
  tj <- function(theta) {
    if(is.list(data[[1]])==FALSE){
      # FOR ONE STEP-STRESS TEST
      tc - tc_end + tauj(theta)
    } else{
      tc - tc_end + tauj_full(theta)
    }
  }

  # return(list(xirctot,stpstrdatsort,tF,tF_end,SF,taui,taui_full,taui_full(LSQest),tc,tc_end,Sc,tauj,tauj_full,tauj_full(LSQest),t_end,taui_full_funct))
  # Fit to log-likelihood distributions
  if (dist=="Weibull") {
    positivity_v[1]<-1

    loglik <- function(theta){
      -sum(log(theta[1]) + (theta[1]-1)*log(ti(theta)) - theta[1]*loglifeF(theta) - ((ti(theta)/lifeF(theta))^theta[1])) - sum(- ((tj(theta)/lifeC(theta))^theta[1]))
    }

    dist_txt<-dist
    distparam_txt<-"\U03B2"
  }
  if (dist=="Lognormal") {
    positivity_v[1]<-1

    loglik <- function(theta){
      -sum(-log(theta[1]) - 0.5*log(2*pi) - log(ti(theta)) - 0.5*(theta[1]^-2)*((log(ti(theta)) - loglifeF(theta))^2)) - sum(log(0.5 - 0.5*erf((2^-0.5)*(theta[1]^-1)*(log(tj(theta)) - loglifeC(theta)))))
    }
    dist_txt<-dist
    distparam_txt<-"\U03C3_t"
  }
  if (dist=="Normal") {
    positivity_v[1]<-1

    loglik <- function(theta){
      -sum(-log(theta[1]) - 0.5*log(2*pi) - 0.5*(theta[1]^-2)*((ti(theta) - lifeF(theta))^2)) - sum(log(0.5 - 0.5*erf((2^-0.5)*(theta[1]^-1)*(tj(theta) - lifeC(theta)))))
    }
    dist_txt<-dist
    distparam_txt<-"\U03C3"
  }
  if (dist=="Exponential") {
    loglik <- function(theta){
      -sum(-loglifeF(theta) - ti(theta)/lifeF(theta)) - sum(-tj(theta)/lifeC(theta))
    }
    dist_txt<-dist
  }
  if (dist=="2PExponential") {
    positivity_v[1]<-1

    loglik <- function(theta){
      -sum(-log(theta[1]) - (theta[1]^-1)*(ti(theta) - lifeF(theta)) - 1) - sum(-(theta[1])*(tj(theta) - lifeC(theta)) - 1)
    }
    dist_txt<-"Two-Parameter Exponential"
    distparam_txt<-"\U03C3"
  }

  if(dist=="Exponential") {
    params_txt<-params_txt
  }
  else {
    params_txt<-c(distparam_txt,params_txt)
  }

  # return(list(loglik,LSQest))

  MLEandvar <- MLE.var.covar.select(loglik,LSQest)
  theta.hat <- MLEandvar[[1]]
  inv.fish  <- MLEandvar[[2]]

  # return(list(theta.hat,inv.fish))

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
      if(ls=="Exponential" && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(ls=="Exponential2" && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(ls=="Arrhenius" && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if((ls=="InversePower") && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if((ls=="InversePower2") && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if((ls=="InversePower3") && i == (ishift+1)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(ls=="TempHumidity" && i == (ishift+1)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if (ls=="TempNonthermal" && i == (ishift+3)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      conflim_txt<-c(paste(c("Lower ",100*conf.level,"%"),collapse = ""),paste(c("Upper ",100*conf.level,"%"),collapse = ""))
    }
    if(sided == "onesidedlow"){
      if(positivity_v[i]==0){
        conflim[[i]]<-theta.hat[i] - crit2 * sqrt(inv.fish[i, i])
      } else if(positivity_v[i]==1){
        conflim[[i]]<-theta.hat[i]*exp(-crit2 * (sqrt(inv.fish[i, i])/theta.hat[i]))
      }
      if(ls=="Exponential" && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(ls=="Exponential2" && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(ls=="Arrhenius" && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if((ls=="Power") && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if((ls=="InversePower") && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if((ls=="InversePower2") && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if((ls=="InversePower3") && i == (ishift+1)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(ls=="TempHumidity" && i == (ishift+1)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if (ls=="TempNonthermal" && i == (ishift+3)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      conflim_txt<-paste(c("One-Sided Low ",100*conf.level,"%"),collapse = "")
    }
    if(sided == "onesidedhigh"){
      if(positivity_v[i]==0){
        conflim[[i]]<-theta.hat[i] + crit2 * sqrt(inv.fish[i, i])
      } else if(positivity_v[i]==1){
        conflim[[i]]<-theta.hat[i]*exp(crit2 * (sqrt(inv.fish[i, i])/theta.hat[i]))
      }
      if(ls=="Exponential" && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(ls=="Exponential2" && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(ls=="Arrhenius" && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if((ls=="Power") && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if((ls=="InversePower") && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if((ls=="InversePower2") && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if((ls=="InversePower3") && i == (ishift+1)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(ls=="TempHumidity" && i == (ishift+1)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if (ls=="TempNonthermal" && i == (ishift+3)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      conflim_txt<-paste(c("One-Sided High ",100*conf.level,"%"),collapse = "")
    }
    fulllimset[[i]]<-c(theta.hat[i],conflim[[i]])

    if((ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius"|| ls=="Power" || ls=="InversePower" || ls=="InversePower2" || ls=="InversePower3") && i == (ishift+2)){
      fulllimset[[i]] <- c(exp(theta.hat[i]),conflim[[i]])
    }
    if(ls=="TempHumidity" && i == (ishift+1)){
      fulllimset[[i]] <- c(exp(theta.hat[i]),conflim[[i]])
    }
    if (ls=="TempNonthermal" && i == (ishift+3)){
      fulllimset[[i]] <- c(exp(theta.hat[i]),conflim[[i]])
    }
  }

  AIC = 2*length(theta.hat) + 2*loglik(theta.hat)
  BIC = 2*log(length(tF)+length(tc)) + 2*loglik(theta.hat)

  # if(is.null(Suse) == FALSE){
  #   # Generate some data
  # }

  # Produce some output text that summarizes the results
  cat(c("Maximum-Likelihood estimates for the ",ls_txt,"-",dist_txt," Step-Stress Life-Stress model.\n\n"),sep = "")
  print(matrix(unlist(fulllimset), nrow = length(unlist(fulllimset))/length(LSQest), ncol = length(LSQest), byrow = FALSE,dimnames = list(c("Step-Stress Life-Stress Parameters Mean",conflim_txt),params_txt)))

  #return(list(ti(LSQest),tj(LSQest),theta.hat,conflim))
  # Add a plot for the LSQ and MLE here as well (CDF)

  # Recompute necessary output
  if(ls=="Exponential"){
    theta.hat[ishift+2] <- exp(theta.hat[ishift+2])
  }
  if(ls=="Exponential2"){
    theta.hat[ishift+2] <- exp(theta.hat[ishift+2])
  }
  if(ls=="Arrhenius"){
    theta.hat[ishift+2] <- exp(theta.hat[ishift+2])
  }
  if(ls=="Power"){
    theta.hat[ishift+2] <- exp(theta.hat[ishift+2])
  }
  if(ls=="InversePower"){
    theta.hat[ishift+2] <- exp(theta.hat[ishift+2])
  }
  if(ls=="InversePower2"){
    theta.hat[ishift+2] <- exp(theta.hat[ishift+2])
  }
  if(ls=="InversePower3"){
    theta.hat[ishift+1] <- exp(theta.hat[ishift+1])
  }
  if(ls=="TempHumidity"){
    theta.hat[ishift+1] <- exp(theta.hat[ishift+1])
  }
  if (ls=="TempNonthermal"){
    theta.hat[ishift+3] <- exp(theta.hat[ishift+3])
  }

  return(list(theta.hat,inv.fish,conflim,AIC=AIC,BIC=BIC))
}
