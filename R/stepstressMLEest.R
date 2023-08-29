# MLE Step-Stress Estimator
# Developed by Dr. Reuel Smith, 2021-2022

stepstress.MLEest <- function(LSQest,data,stepstresstable,ls,dist,confid,sided){
  #Load pracma library for erf
  library(pracma)
  library(dplyr)
  library(matrixcalc)

  # Check to see if dist="Exponential" so you can exclude life
  # distribution parameters.
  if (dist=="Exponential") {
    ishift<-0
  } else {
    ishift<-1
  }
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

  # Setup positivity check vector for parameters
  positivity_v<-rep(0,length(LSQest))

  # Re-sort the input data and then separate censored from non-censored data
  stpstrdatsort<-stepstress.data.cum(data,stepstresstable)
  xirctot<-sort.xircstressdata(stpstrdatsort[[1]])

  # The failure and survival stress levels are assigned a vector of the same
  # sizes as the adjusted time vectors for the MLE operation.
  SF <- xirctot[[3]]
  Sc <- xirctot[[4]]


  # Obtain the cumulative sum of the end times
  t_end <- cumsum(stepstresstable[,dim(stepstresstable)[[2]]])

  # Compute number of stresses
  Nstress<-dim(stepstresstable)[2]-1

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
  }

  if (ls=="Exponential"){
    # theta[1] - parameter a, theta[2] - parameter b
    positivity_v[ishift+2]<-1

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      exp(theta[ishift+1]*(S1 - S2))
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      theta[ishift+2]*((exp(theta[ishift+1]*S2)) - (exp(theta[ishift+1]*S1)))
    }

    # Life functions
    lifeF <- function(theta) {
      theta[ishift+2]*exp(SF*theta[ishift+1])
    }
    loglifeF <- function(theta) {
      log(theta[ishift+2]) + SF*theta[ishift+1]
    }
    lifeC <- function(theta) {
      theta[ishift+2]*exp(Sc*theta[ishift+1])
    }
    loglifeC <- function(theta) {
      log(theta[ishift+2]) + Sc*theta[ishift+1]
    }
  }

  if (ls=="Arrhenius") {
    # lsparams[1] - parameter Ea, lsparams[2] - parameter b
    # Temperature HAS to be in Kelvin for this to work
    positivity_v[ishift+1]<-1
    positivity_v[ishift+2]<-1

    K<-8.617385e-5

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      exp((theta[ishift+1]/K)*((1/S1) - (1/S2)))
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      theta[ishift+2]*((exp(theta[ishift+1]/(K*S2))) - (exp(theta[ishift+1]/(K*S1))))
    }

    # Life functions
    lifeF <- function(theta) {
      theta[ishift+2]*exp(theta[ishift+1]/(K*SF))
    }
    loglifeF <- function(theta) {
      log(theta[ishift+2]) + theta[ishift+1]/(K*SF)
    }
    lifeC <- function(theta) {
      theta[ishift+2]*exp(theta[ishift+1]/(K*Sc))
    }
    loglifeC <- function(theta) {
      log(theta[ishift+2]) + theta[ishift+1]/(K*Sc)
    }
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
  }

  if (ls=="Power") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    positivity_v[ishift+2]<-1

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      (S1/S2)^theta[ishift+1]
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      theta[ishift+2]*((S2^theta[ishift+1])-(S1^theta[ishift+1]))
    }

    # Life functions
    lifeF <- function(theta) {
      theta[ishift+2]*(SF^theta[ishift+1])
    }
    loglifeF <- function(theta) {
      log(theta[ishift+2]) + theta[ishift+1]*log(SF)
    }
    lifeC <- function(theta) {
      theta[ishift+2]*(Sc^theta[ishift+1])
    }
    loglifeC <- function(theta) {
      log(theta[ishift+2]) + theta[ishift+1]*log(Sc)
    }
  }

  if (ls=="InversePower") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    positivity_v[ishift+2]<-1

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
  }

  if (ls=="TempHumidity") {
    # lsparams[1] - parameter A, lsparams[2] - parameter a, lsparams[3] - parameter b
    positivity_v[ishift+1]<-1

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      exp((theta[ishift+2]*((1/S1[1]) - (1/S2[1]))) + (theta[ishift+3]*((1/S1[2]) - (1/S2[2]))))
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      theta[ishift+1]*((exp((theta[ishift+2]/S2[1])- (theta[ishift+3]/S2[2]))) - (exp((theta[ishift+2]/S1[1])- (theta[ishift+3]/S1[2]))))
    }

    # Life functions
    lifeF <- function(theta) {
      theta[ishift+1]*exp((theta[ishift+2]/SF[,1]) + (theta[ishift+3]/SF[,2]))
    }
    loglifeF <- function(theta) {
      log(theta[ishift+1]) + (theta[ishift+2]/SF[,1]) + (theta[ishift+3]/SF[,2])
    }
    lifeC <- function(theta) {
      theta[ishift+1]*exp((theta[ishift+2]/Sc[,1]) + (theta[ishift+3]/Sc[,2]))
    }
    loglifeC <- function(theta) {
      log(theta[ishift+1]) + (theta[ishift+2]/Sc[,1]) + (theta[ishift+3]/Sc[,2])
    }
  }

  if (ls=="TempNonthermal") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
    positivity_v[ishift+3]<-1

    # Life ratio L1/L2
    L1dL2 <- function(theta,S1,S2) {
      ((S2[1]/S1[1])^theta[ishift+2])*exp(-theta[ishift+1]*((1/S1[2]) - (1/S2[2])))
    }
    # Life difference L2-L1
    L2mL1 <- function(theta,S1,S2) {
      theta[ishift+3]*((1/((S2[1]^theta[ishift+2])*exp(-theta[ishift+1]/S2[2]))) - (1/((S1[1]^theta[ishift+2])*exp(-theta[ishift+1]/S1[2]))))
    }

    # Life functions
    lifeF <- function(theta) {
      theta[ishift+3]/((SF[,2]^theta[ishift+2])*exp(-theta[ishift+1]/SF[,1]))
    }
    loglifeF <- function(theta) {
      log(theta[ishift+3]) - theta[ishift+2]*log(SF[,2]) + (theta[ishift+1]/SF[,1])
    }
    lifeC <- function(theta) {
      theta[ishift+3]/((Sc[,2]^theta[ishift+2])*exp(-theta[ishift+1]/Sc[,1]))
    }
    loglifeC <- function(theta) {
      log(theta[ishift+3]) - theta[ishift+2]*log(Sc[,2]) + (theta[ishift+1]/Sc[,1])
    }
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
  }
  tauifunct<-integer(0)
  taujfunct<-integer(0)
  taui_steps<-unique(stpstrdatsort[[6]])
  tauj_steps<-unique(stpstrdatsort[[7]])

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
      }
    }
    return(tauifunct)
  }

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
      }
    }
    return(taujfunct)
  }



  # Form the failure and censored time functions
  ti <- function(theta) {
    xirctot[[1]] - stpstrdatsort[[3]] + taui(theta)
  }
  tj <- function(theta) {
    xirctot[[2]] - stpstrdatsort[[4]] + tauj(theta)
  }

  # Fit to log-likelihood distributions
  if (dist=="Weibull") {
    positivity_v[1]<-1

    loglik <- function(theta){
      -sum(log(theta[1]) + (theta[1]-1)*log(ti(theta)) - theta[1]*loglifeF(theta) - ((ti(theta)/lifeF(theta))^theta[1])) - sum(- ((tj(theta)/lifeC(theta))^theta[1]))
    }
  }
  if (dist=="Lognormal") {
    positivity_v[1]<-1

    loglik <- function(theta){
      -sum(-log(theta[1]) - 0.5*log(2*pi) - log(ti(theta)) - 0.5*(theta[1]^-2)*((log(ti(theta)) - loglifeF(theta))^2)) - sum(log(0.5 - 0.5*erf((2^-0.5)*(theta[1]^-1)*(log(tj(theta)) - loglifeC(theta)))))
    }
  }
  if (dist=="Normal") {
    positivity_v[1]<-1

    loglik <- function(theta){
      -sum(-log(theta[1]) - 0.5*log(2*pi) - 0.5*(theta[1]^-2)*((ti(theta) - lifeF(theta))^2)) - sum(log(0.5 - 0.5*erf((2^-0.5)*(theta[1]^-1)*(tj(theta) - lifeC(theta)))))
    }
  }
  if (dist=="Exponential") {
    loglik <- function(theta){
      -sum(-loglifeF(theta) - ti(theta)/lifeF(theta)) - sum(-tj(theta)/lifeC(theta))
    }
  }
  if (dist=="2PExponential") {
    positivity_v[1]<-1

    loglik <- function(theta){
      -sum(-log(theta[1]) - (theta[1]^-1)*(ti(theta) - lifeF(theta)) - 1) - sum(-(theta[1])*(tj(theta) - lifeC(theta)) - 1)
    }
  }

  MLEandvar <- MLE.var.covar.select(loglik,LSQest)
  theta.hat <- MLEandvar[[1]]
  inv.fish  <- MLEandvar[[2]]

  crit <- qnorm((1 + conf.level)/2)
  crit2 <- qnorm(conf.level)
  conflim<-vector(mode = "list", length = length(LSQest))

  for(i in 1:length(LSQest)){
    if(sided == "twosided"){
      if(positivity_v[i]==0){
        conflim[[i]]<-theta.hat[i] + c(-1, 1) * crit * sqrt(inv.fish[i, i])
      } else if(positivity_v[i]==1){
        conflim[[i]]<-theta.hat[i]*exp(c(-1, 1) * crit * (sqrt(inv.fish[i, i])/theta.hat[i]))
      }
    }
    if(sided == "onesidedlow"){
      if(positivity_v[i]==0){
        conflim[[i]]<-theta.hat[i] - crit2 * sqrt(inv.fish[i, i])
      } else if(positivity_v[i]==1){
        conflim[[i]]<-theta.hat[i]*exp(-crit2 * (sqrt(inv.fish[i, i])/theta.hat[i]))
      }
    }
    if(sided == "onesidedhigh"){
      if(positivity_v[i]==0){
        conflim[[i]]<-theta.hat[i] + crit2 * sqrt(inv.fish[i, i])
      } else if(positivity_v[i]==1){
        conflim[[i]]<-theta.hat[i]*exp(crit2 * (sqrt(inv.fish[i, i])/theta.hat[i]))
      }
    }
  }
  #return(list(ti(LSQest),tj(LSQest),theta.hat,conflim))
  return(list(theta.hat,inv.fish,conflim))
}
