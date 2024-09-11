# Maximum Likelihood Life-Stress Estimator
# Developed by Dr. Reuel Smith, 2021-2022

lifestress.MLEest <- function(LSQest,ls,dist,TTF,SF,Tc=NULL,Sc=NULL,Suse=NULL,confid=0.95,sided="twosided"){
  # (NOTE -11/13/2023 RS) Add Exponential2 model that uses inverse stress as input

  #Load pracma library for erf
  library(pracma)
  library(matrixcalc)
  library(ucminf)
  library(MASS)

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

  # Initialize life-stress parameter estimates for theta
  if (ls=="Linear") {
    # theta[1] - parameter a, theta[2] - parameter b
    lifeF <- function(theta) {
      theta[ishift+2] + SF*theta[ishift+1]
    }
    loglifeF <- function(theta) {
      log(theta[ishift+2] + SF*theta[ishift+1])
    }
    if(is.null(Tc)==FALSE){
      lifeC <- function(theta) {
        theta[ishift+2] + Sc*theta[ishift+1]
      }
      loglifeC <- function(theta) {
        log(theta[ishift+2] + Sc*theta[ishift+1])
      }
    }
    ls_txt<-ls
    params_txt<-c("a","b")
  }

  if (ls=="Exponential"){
    # theta[1] - parameter a, theta[2] - parameter b
    # shift b to log b
    # positivity_v[ishift+2]<-1
    LSQest[ishift+2] <- log(LSQest[ishift+2])

    lifeF <- function(theta) {
      exp(theta[ishift+2])*exp(SF*theta[ishift+1])
    }
    loglifeF <- function(theta) {
      theta[ishift+2] + SF*theta[ishift+1]
    }
    if(is.null(Tc)==FALSE){
      lifeC <- function(theta) {
        exp(theta[ishift+2])*exp(Sc*theta[ishift+1])
      }
      loglifeC <- function(theta) {
        theta[ishift+2] + Sc*theta[ishift+1]
      }
    }
    ls_txt<-ls
    params_txt<-c("a","b")
  }

  # (NOTE -11/13/2023 RS) Add Exponential2 model that uses inverse stress as input
  if (ls=="Exponential2"){
    # theta[1] - parameter a, theta[2] - parameter b
    # shift b to log b
    # positivity_v[ishift+2]<-1
    LSQest[ishift+2] <- log(LSQest[ishift+2])

    lifeF <- function(theta) {
      exp(theta[ishift+2])*exp(theta[ishift+1]/SF)
    }
    loglifeF <- function(theta) {
      theta[ishift+2] + theta[ishift+1]/SF
    }
    if(is.null(Tc)==FALSE){
      lifeC <- function(theta) {
        exp(theta[ishift+2])*exp(theta[ishift+1]/Sc)
      }
      loglifeC <- function(theta) {
        theta[ishift+2] + theta[ishift+1]/Sc
      }
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
    lifeF <- function(theta) {
      exp(theta[ishift+2])*exp(theta[ishift+1]/(K*SF))
    }
    loglifeF <- function(theta) {
      theta[ishift+2] + theta[ishift+1]/(K*SF)
    }
    if(is.null(Tc)==FALSE){
      lifeC <- function(theta) {
        exp(theta[ishift+2])*exp(theta[ishift+1]/(K*Sc))
      }
      loglifeC <- function(theta) {
        theta[ishift+2] + theta[ishift+1]/(K*Sc)
      }
    }
    ls_txt<-ls
    params_txt<-c("E_a","b")
  }

  if (ls=="Eyring") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # shift b to log b
    # positivity_v[ishift+2]<-1
    LSQest[ishift+2] <- log(LSQest[ishift+2])

    lifeF <- function(theta) {
      (exp(theta[ishift+2])/SF)*exp(theta[ishift+1]/SF)
    }
    loglifeF <- function(theta) {
      theta[ishift+2] - log(SF) + theta[ishift+1]/SF
    }
    if(is.null(Tc)==FALSE){
      lifeC <- function(theta) {
        (exp(theta[ishift+2])/Sc)*exp(theta[ishift+1]/Sc)
      }
      loglifeC <- function(theta) {
        theta[ishift+2] - log(Sc) + theta[ishift+1]/Sc
      }
    }
    ls_txt<-ls
    params_txt<-c("a","b")
  }

  if (ls=="Eyring2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lifeF <- function(theta) {
      (1/SF)*exp(-(theta[ishift+1] - (theta[ishift+2]/SF)))
    }
    loglifeF <- function(theta) {
      -log(SF) - theta[ishift+1] + theta[ishift+2]/SF
    }
    if(is.null(Tc)==FALSE){
      lifeC <- function(theta) {
        (1/Sc)*exp(-(theta[ishift+1] - (theta[ishift+2]/Sc)))
      }
      loglifeC <- function(theta) {
        - log(Sc) - theta[ishift+1] + theta[ishift+2]/Sc
      }
    }
    ls_txt<-"Eyring (Type-2)"
    params_txt<-c("a","b")
  }

  if (ls=="Power") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # First redefine b parameter as logb
    LSQest[ishift+2] <- log(LSQest[ishift+2])
    # positivity_v[ishift+2]<-1

    lifeF <- function(theta) {
      exp(theta[ishift+2])*(SF^theta[ishift+1])
    }
    loglifeF <- function(theta) {
      theta[ishift+2] + theta[ishift+1]*log(SF)
    }
    if(is.null(Tc)==FALSE){
      lifeC <- function(theta) {
        exp(theta[ishift+2])*(Sc^theta[ishift+1])
      }
      loglifeC <- function(theta) {
        theta[ishift+2] + theta[ishift+1]*log(Sc)
      }
    }
    ls_txt<-ls
    params_txt<-c("a","b")
  }

  if (ls=="InversePower") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # First redefine b parameter as logb
    LSQest[ishift+2] <- log(LSQest[ishift+2])
    # positivity_v[ishift+2]<-1

    lifeF <- function(theta) {
      exp(theta[ishift+2])*(SF^-theta[ishift+1])
    }
    loglifeF <- function(theta) {
      theta[ishift+2] - theta[ishift+1]*log(SF)
    }
    if(is.null(Tc)==FALSE){
      lifeC <- function(theta) {
        exp(theta[ishift+2])*(Sc^-theta[ishift+1])
      }
      loglifeC <- function(theta) {
        theta[ishift+2] - theta[ishift+1]*log(Sc)
      }
    }
    ls_txt<-"Inverse Power"
    params_txt<-c("a","b")
  }

  if (ls=="InversePower2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # First redefine b parameter as logb
    LSQest[ishift+2] <- log(LSQest[ishift+2])
    # positivity_v[ishift+2]<-1

    lifeF <- function(theta) {
      1/(exp(theta[ishift+2])*(SF^theta[ishift+1]))
    }
    loglifeF <- function(theta) {
      -theta[ishift+2] - theta[ishift+1]*log(SF)
    }
    if(is.null(Tc)==FALSE){
      lifeC <- function(theta) {
        1/(exp(theta[ishift+2])*(Sc^theta[ishift+1]))
      }
      loglifeC <- function(theta) {
        -theta[ishift+2] - theta[ishift+1]*log(Sc)
      }
    }
    ls_txt<-"Inverse Power"
    params_txt<-c("a","b")
  }

  if (ls=="Logarithmic") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b

    lifeF <- function(theta) {
      theta[ishift+1]*log(SF) + theta[ishift+2]
    }
    loglifeF <- function(theta) {
      log(theta[ishift+1]*log(SF) + theta[ishift+2])
    }
    if(is.null(Tc)==FALSE){
      lifeC <- function(theta) {
        theta[ishift+1]*log(Sc) + theta[ishift+2]
      }
      loglifeC <- function(theta) {
        log(theta[ishift+1]*log(Sc) + theta[ishift+2])
      }
    }
    ls_txt<-ls
    params_txt<-c("a","b")
  }

  if (ls=="MultiStress") {
    # lsparams - c(ao,a1,...,an), S - c(1,S1,...,Sn)
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
    if(is.null(Tc)==FALSE){
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
    ls_txt<-"Multi-Stress"
    params_txt<-c("a","b")
  }

  if (ls=="TempHumidity") {
    # lsparams[1] - parameter A, lsparams[2] - parameter a, lsparams[3] - parameter b
    # (NOTE: 11/13/2023) 'A' seems to be overly sensitive to error messages so I may use 'log A' here instead
    # along with 'exp(log A)' in place of A
    # First redefine A parameter as logA
    LSQest[ishift+1] <- log(LSQest[ishift+1])
    # positivity_v[ishift+1]<-1

    lifeF <- function(theta) {
      # theta[ishift+1]*exp((theta[ishift+2]/SF[,1]) + (theta[ishift+3]/SF[,2]))
      exp(theta[ishift+1])*exp((theta[ishift+2]/SF[,1]) + (theta[ishift+3]/SF[,2]))
    }
    loglifeF <- function(theta) {
      # log(theta[ishift+1]) + (theta[ishift+2]/SF[,1]) + (theta[ishift+3]/SF[,2])
      theta[ishift+1] + (theta[ishift+2]/SF[,1]) + (theta[ishift+3]/SF[,2])
    }
    if(is.null(Tc)==FALSE){
      lifeC <- function(theta) {
        # theta[ishift+1]*exp((theta[ishift+2]/Sc[,1]) + (theta[ishift+3]/Sc[,2]))
        exp(theta[ishift+1])*exp((theta[ishift+2]/Sc[,1]) + (theta[ishift+3]/Sc[,2]))
      }
      loglifeC <- function(theta) {
        # log(theta[ishift+1]) + (theta[ishift+2]/Sc[,1]) + (theta[ishift+3]/Sc[,2])
        theta[ishift+1] + (theta[ishift+2]/Sc[,1]) + (theta[ishift+3]/Sc[,2])
      }
    }
    ls_txt<-"Temperature-Humidity"
    params_txt<-c("A","a","b")
  }

  if (ls=="TempNonthermal") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
    # First redefine c parameter as logc
    LSQest[ishift+3] <- log(LSQest[ishift+3])
    # positivity_v[ishift+3]<-1

    lifeF <- function(theta) {
      exp(theta[ishift+3])/((SF[,2]^theta[ishift+2])*exp(-theta[ishift+1]/SF[,1]))
    }
    loglifeF <- function(theta) {
      theta[ishift+3] - theta[ishift+2]*log(SF[,2]) + (theta[ishift+1]/SF[,1])
    }
    if(is.null(Tc)==FALSE){
      lifeC <- function(theta) {
        exp(theta[ishift+3])/((Sc[,2]^theta[ishift+2])*exp(-theta[ishift+1]/Sc[,1]))
      }
      loglifeC <- function(theta) {
        theta[ishift+3] - theta[ishift+2]*log(Sc[,2]) + (theta[ishift+1]/Sc[,1])
      }
    }
    ls_txt<-"Temperature-Non-thermal"
    params_txt<-c("a","b","c")
  }

  if (ls=="Eyring3") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # lsparams[3] - parameter c, lsparams[4] - parameter d
    lifeF <- function(theta) {
      (1/SF[,1])*exp((theta[ishift+1] + (theta[ishift+2]/SF[,1])) + (theta[ishift+3] + (theta[ishift+4]/SF[,1]))*SF[,2])
    }
    loglifeF <- function(theta) {
      -log(SF[,1]) + theta[ishift+1] + (theta[ishift+2]/SF[,1]) + (theta[ishift+3] + (theta[ishift+4]/SF[,1]))*SF[,2]
    }
    if(is.null(Tc)==FALSE){
      lifeC <- function(theta) {
        (1/Sc[,1])*exp((theta[ishift+1] + (theta[ishift+2]/Sc[,1])) + (theta[ishift+3] + (theta[ishift+4]/Sc[,1]))*Sc[,2])
      }
      loglifeC <- function(theta) {
        -log(Sc[,1]) + theta[ishift+1] + (theta[ishift+2]/Sc[,1]) + (theta[ishift+3] + (theta[ishift+4]/Sc[,1]))*Sc[,2]
      }
    }
    ls_txt<-"Eyring (Type 3)"
    params_txt<-c("a","b","c","d")
  }

  if (ls=="Eyring4") {
    # lsparams[1] - parameter A, lsparams[2] - parameter b
    # lsparams[3] - parameter Ea
    # Temperature HAS to be in Kelvin for this to work
    K<-8.617385e-5
    lifeF <- function(theta) {
      theta[ishift+1]*exp(theta[ishift+3]/(K*SF[,1]))*(SF[,2]^-theta[ishift+2])
    }
    loglifeF <- function(theta) {
      log(theta[ishift+1]) + theta[ishift+3]/(K*SF[,1]) - theta[ishift+2]*log(SF[,2])
    }
    if(is.null(Tc)==FALSE){
      lifeC <- function(theta) {
        theta[ishift+1]*exp(theta[ishift+3]/(K*Sc[,1]))*(Sc[,2]^-theta[ishift+2])
      }
      loglifeC <- function(theta) {
        log(theta[ishift+1]) + theta[ishift+3]/(K*Sc[,1]) - theta[ishift+2]*log(Sc[,2])
      }
    }
    ls_txt<-"Eyring (Type 4)"
    params_txt<-c("A","b","E_a")
  }

  # Fit to log-likelihood distributions
  if (dist=="Weibull") {
    # positivity_v[1]<-1
    # First redefine beta parameter as logbeta
    LSQest[1] <- log(LSQest[1])

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(theta[1] + (exp(theta[1])-1)*log(TTF) - exp(theta[1])*loglifeF(theta) - ((TTF/lifeF(theta))^exp(theta[1])))
      }
    } else{
      loglik <- function(theta){
        -sum(theta[1] + (exp(theta[1])-1)*log(TTF) - exp(theta[1])*loglifeF(theta) - ((TTF/lifeF(theta))^exp(theta[1]))) - sum(- ((Tc/lifeC(theta))^exp(theta[1])))
      }
    }
    dist_txt<-dist
    distparam_txt<-"\U03B2"
  }

  if (dist=="3PWeibull") {
    # positivity_v[1]<-1
    # First redefine beta parameter as logbeta
    LSQest[1] <- log(LSQest[1])

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(theta[1] + (exp(theta[1])-1)*log(TTF-theta[2]) - exp(theta[1])*loglifeF(theta) - (((TTF-theta[2])/lifeF(theta))^exp(theta[1])))
      }
    } else{
      loglik <- function(theta){
        -sum(theta[1] + (exp(theta[1])-1)*log(TTF-theta[2]) - exp(theta[1])*loglifeF(theta) - (((TTF-theta[2])/lifeF(theta))^exp(theta[1]))) - sum(- (((Tc-theta[1])/lifeC(theta))^exp(theta[1])))
      }
    }
    dist_txt<-"Three-Parameter Weibull"
    distparam_txt<-c("\U03B2","\U03B3")
  }

  if (dist=="Lognormal") {
    # positivity_v[1]<-1
    # First redefine sigmat parameter as logsigmat
    LSQest[1] <- log(LSQest[1])

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(-theta[1] - 0.5*log(2*pi) - log(TTF) - 0.5*(exp(theta[1])^-2)*((log(TTF) - loglifeF(theta))^2))
      }
    } else{
      loglik <- function(theta){
        -sum(-theta[1] - 0.5*log(2*pi) - log(TTF) - 0.5*(exp(theta[1])^-2)*((log(TTF) - loglifeF(theta))^2)) - sum(log(0.5 - 0.5*erf((2^-0.5)*(exp(theta[1])^-1)*(log(Tc) - loglifeC(theta)))))
      }
    }
    dist_txt<-dist
    distparam_txt<-"\U03C3_t"
  }
  if (dist=="Normal") {
    # positivity_v[1]<-1
    # First redefine sigma parameter as logsigma
    LSQest[1] <- log(LSQest[1])

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(-theta[1] - 0.5*log(2*pi) - 0.5*(exp(theta[1])^-2)*((TTF - lifeF(theta))^2))
      }
    } else{
      loglik <- function(theta){
        -sum(-theta[1] - 0.5*log(2*pi) - 0.5*(exp(theta[1])^-2)*((TTF - lifeF(theta))^2)) - sum(log(0.5 - 0.5*erf((2^-0.5)*(exp(theta[1])^-1)*(Tc - lifeC(theta)))))
      }
    }
    dist_txt<-dist
    distparam_txt<-"\U03C3"
  }
  if (dist=="Exponential") {
    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(-loglifeF(theta) - TTF/lifeF(theta))
      }
    } else{
      loglik <- function(theta){
        -sum(-loglifeF(theta) - TTF/lifeF(theta)) - sum(-Tc/lifeC(theta))
      }
    }
    dist_txt<-dist
  }
  if (dist=="2PExponential") {
    # positivity_v[1]<-1
    # First redefine sigmat parameter as logsigmat
    LSQest[1] <- log(LSQest[1])

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(-theta[1] - (exp(theta[1])^-1)*(TTF - lifeF(theta)) - 1)
      }
    } else{
      loglik <- function(theta){
        -sum(-theta[1] - (exp(theta[1])^-1)*(TTF - lifeF(theta)) - 1) - sum(-(theta[1])*(Tc - lifeC(theta)) - 1)
      }
    }
    dist_txt<-"Two-Parameter Exponential"
    distparam_txt<-"\U03C3"
  }

  # Group all parameters
  # Check to see if any distribution parameters were tabulated
  if(dist=="Exponential") {
    params_txt<-params_txt
  }
  else {
    params_txt<-c(distparam_txt,params_txt)
  }

  return(list(loglik,LSQest))

  MLEandvar <- MLE.var.covar.select(loglik,LSQest)
  # return(MLEandvar)
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
      # Computes back the original scale
      if((dist=="Weibull" || dist=="3PWeibull" || dist=="Lognormal" || dist=="Normal" || dist=="2PExponential") && i == 1){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(ls=="TempHumidity" && i == (ishift+1)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if((ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" || ls=="Power" || ls=="InversePower" || ls=="InversePower2") && i == (ishift+2)){
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
      if((dist=="Weibull" || dist=="3PWeibull" || dist=="Lognormal" || dist=="Normal" || dist=="2PExponential") && i == 1){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(ls=="TempHumidity" && i == (ishift+1)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if((ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" || ls=="Power" || ls=="InversePower" || ls=="InversePower2") && i == (ishift+2)){
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
      if((dist=="Weibull" || dist=="3PWeibull" || dist=="Lognormal" || dist=="Normal" || dist=="2PExponential") && i == 1){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if(ls=="TempHumidity" && i == (ishift+1)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if((ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" ||  ls=="Power" || ls=="InversePower" || ls=="InversePower2") && i == (ishift+2)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      if (ls=="TempNonthermal" && i == (ishift+3)){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      conflim_txt<-paste(c("One-Sided High ",100*conf.level,"%"),collapse = "")
    }
    fulllimset[[i]]<-c(theta.hat[i],conflim[[i]])

    if((dist=="Weibull" || dist=="3PWeibull" || dist=="Lognormal" || dist=="Normal" || dist=="2PExponential") && i == 1){
      fulllimset[[i]] <- c(exp(theta.hat[i]),conflim[[i]])
    }
    if(ls=="TempHumidity" && i == (ishift+1)){
      fulllimset[[i]] <- c(exp(theta.hat[i]),conflim[[i]])
    }
    if((ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" || ls=="Power" || ls=="InversePower" || ls=="InversePower2") && i == (ishift+2)){
      fulllimset[[i]] <- c(exp(theta.hat[i]),conflim[[i]])
    }
    if (ls=="TempNonthermal" && i == (ishift+3)){
      fulllimset[[i]] <- c(exp(theta.hat[i]),conflim[[i]])
    }
  }

  AIC = 2*length(theta.hat) + 2*loglik(theta.hat)
  BIC = 2*log(length(TTF)+length(Tc)) + 2*loglik(theta.hat)

  if(is.null(Suse) == FALSE){
    # Generate some data
  }

  # Produce some output text that summarizes the results
  cat(c("Maximum-Likelihood estimates for the ",ls_txt,"-",dist_txt," Life-Stress model.\n\n"),sep = "")
  print(matrix(unlist(fulllimset), nrow = length(unlist(fulllimset))/length(LSQest), ncol = length(LSQest), byrow = FALSE,dimnames = list(c("Life-Stress Parameters Mean",conflim_txt),params_txt)))

  # Recompute necessary output
  if(dist=="Weibull" || dist=="3PWeibull" || dist=="Lognormal" || dist=="Normal" || dist=="2PExponential"){
    theta.hat[1] <- exp(theta.hat[1])
  }
  if(ls=="TempHumidity"){
    theta.hat[ishift+1] <- exp(theta.hat[ishift+1])
  }
  if(ls=="Exponential"){
    theta.hat[ishift+2] <- exp(theta.hat[ishift+2])
  }
  if(ls=="Exponential2"){
    theta.hat[ishift+2] <- exp(theta.hat[ishift+2])
  }
  if(ls=="Arrhenius"){
    theta.hat[ishift+2] <- exp(theta.hat[ishift+2])
  }
  if(ls=="Power"  || ls=="Eyring" || ls=="InversePower" || ls=="InversePower2"){
    theta.hat[ishift+2] <- exp(theta.hat[ishift+2])
  }
  if (ls=="TempNonthermal"){
    theta.hat[ishift+3] <- exp(theta.hat[ishift+3])
  }
  return(list(theta.hat,inv.fish,conflim,AIC,BIC))
}
