# Maximum Likelihood Life-Stress Estimator
# Developed by Dr. Reuel Smith, 2021-2026

lifestress.MLEest <- function(data,ls,dist,Suse=NULL,confid=0.95,
                              sided="twosided",pp="Blom",xlabel1="X",Llab=NULL,Slab=NULL,SLab2=NULL,param2=NULL,
                              stressunit1 = NULL, stressunit2 = NULL){
  # (NOTE -11/13/2023 RS) Add Exponential2 model that uses inverse stress as input
  # UPDATE (12/17/2025 RCS) Simplified lifestress.MLEest so that user only enters data table and the tool divides data into
  # TTF, SF, Tc, and Sc itself along with obtaining the LSQest

  #Load pracma library for erf
  library(pracma)
  library(matrixcalc)
  library(ucminf)
  library(MASS)

  # Check to see if dist="Exponential" so you can exclude life
  # distribution parameters.
  if (dist=="Exponential" || (dist=="Weibull" && is.null(param2) == FALSE)) {
    ishift<-0
  } else {
    ishift<-1
  }
  # Check to see if confidence exists
  conf.level <- confid

  # Break up data into TTF, SF, Tc, and Sc (if Tc and Sc exist)
  data.segments <- sort.xircstressdata(data)
  TTF <- data.segments[[1]]   # Time to failure data
  Tc <- data.segments[[2]]    # censored time data (if exists)
  SF <- data.segments[[3]]    # Stress at failure data
  Sc <- data.segments[[4]]    # Stress at censored time

  output.1 <- lifestress.LSQest(data,ls,dist = dist,pp=pp,xlabel1=xlabel1,Llab=Llab,Suse=Suse,Slab=Slab,stressunit1 = stressunit1)[c(1,3)] # Compute least squares estimate LSQest and pull stress levels
  LSQest <- output.1[[2]]     # Least squares estimate to compute MLE estimate
  S <- output.1[[1]]          # Unit stresses to compute life L based on MLE

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
    if(is.null(Suse) == FALSE){
      lifeUSE <- function(theta) {
        theta[ishift+2] + Suse*theta[ishift+1]
      }
      lifeUSEVAR <- function(theta,VARCOV) {
        VARCOV[ishift+2,ishift+2] +
          (Suse^2)*VARCOV[ishift+1,ishift+1] +
          2*Suse*VARCOV[ishift+1,ishift+2]
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
    if(is.null(Suse) == FALSE){
      lifeUSE <- function(theta) {
        exp(theta[ishift+2])*exp(Suse*theta[ishift+1])
      }
      lifeUSEVAR <- function(theta,VARCOV) {
        (lifeUSE(theta)^2)*VARCOV[ishift+2,ishift+2] +
          (Suse^2)*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+1] +
          2*Suse*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+2]
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
    if(is.null(Suse) == FALSE){
      lifeUSE <- function(theta) {
        exp(theta[ishift+2])*exp(theta[ishift+1]/Suse)
      }
      lifeUSEVAR <- function(theta,VARCOV) {
        (lifeUSE(theta)^2)*VARCOV[ishift+2,ishift+2] +
          (1/Suse^2)*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+1] +
          2*(1/Suse)*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+2]
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
    if(is.null(Suse) == FALSE){
      lifeUSE <- function(theta) {
        exp(theta[ishift+2])*exp(theta[ishift+1]/(K*Suse))
      }
      lifeUSEVAR <- function(theta,VARCOV) {
        (lifeUSE(theta)^2)*VARCOV[ishift+2,ishift+2] +
          (1/(K*Suse)^2)*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+1] +
          2*(1/(K*Suse))*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+2]
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
    if(is.null(Suse) == FALSE){
      lifeUSE <- function(theta) {
        (exp(theta[ishift+2])/Suse)*exp(theta[ishift+1]/Suse)
      }
      lifeUSEVAR <- function(theta,VARCOV) {
        (lifeUSE(theta)^2)*VARCOV[ishift+2,ishift+2] +
          (1/Suse^2)*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+1] +
          2*(1/Suse)*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+2]
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
    if(is.null(Suse) == FALSE){
      lifeUSE <- function(theta) {
        (1/Suse)*exp(-(theta[ishift+1] - (theta[ishift+2]/Suse)))
      }
      lifeUSEVAR <- function(theta,VARCOV) {
        (1/Suse^2)*(lifeUSE(theta)^2)*VARCOV[ishift+2,ishift+2] +
          (lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+1] -
          2*(1/Suse)*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+2]
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
    if(is.null(Suse) == FALSE){
      lifeUSE <- function(theta) {
        exp(theta[ishift+2])*(Suse^theta[ishift+1])
      }
      lifeUSEVAR <- function(theta,VARCOV) {
        (lifeUSE(theta)^2)*VARCOV[ishift+2,ishift+2] +
          (log(Suse)^2)*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+1] +
          2*log(Suse)*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+2]
      }
    }
    ls_txt<-ls
    params_txt<-c("a","b")
  }

  if (ls=="PowerwithBias") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
    # First redefine b parameter as logb
    # and c parameter as logc
    LSQest[ishift+2] <- log(LSQest[ishift+2])
    LSQest[ishift+3] <- log(LSQest[ishift+3])

    lifeF <- function(theta) {
      exp(theta[ishift+3]) + exp(theta[ishift+2])*(SF^theta[ishift+1])
    }
    loglifeF <- function(theta) {
      log(exp(theta[ishift+3]) + exp(theta[ishift+2])*(SF^theta[ishift+1]))
    }
    if(is.null(Tc)==FALSE){
      lifeC <- function(theta) {
        exp(theta[ishift+3]) + exp(theta[ishift+2])*(Sc^theta[ishift+1])
      }
      loglifeC <- function(theta) {
        log(exp(theta[ishift+3]) + exp(theta[ishift+2])*(Sc^theta[ishift+1]))
      }
    }
    if(is.null(Suse) == FALSE){
      lifeUSE <- function(theta) {
        exp(theta[ishift+3]) + exp(theta[ishift+2])*(Suse^theta[ishift+1])
      }
      lifeUSEVAR <- function(theta,VARCOV) {
        ((exp(theta[ishift+3]))^2)*VARCOV[ishift+3,ishift+3] +
          ((exp(theta[ishift+2])*(Suse^theta[ishift+1]))^2)*VARCOV[ishift+2,ishift+2] +
          ((log(Suse)*exp(theta[ishift+2])*(Suse^theta[ishift+1]))^2)*VARCOV[ishift+1,ishift+1] -
          2*(log(Suse)*exp(theta[ishift+2])*(Suse^theta[ishift+1]))*(exp(theta[ishift+2])*(Suse^theta[ishift+1]))*VARCOV[ishift+1,ishift+2] -
          2*(log(Suse)*exp(theta[ishift+2])*(Suse^theta[ishift+1]))*(exp(theta[ishift+3]))*VARCOV[ishift+1,ishift+3] -
          2*(exp(theta[ishift+2])*(Suse^theta[ishift+1]))*(exp(theta[ishift+3]))*VARCOV[ishift+2,ishift+3]
      }
    }
    ls_txt<-ls
    params_txt<-c("a","b","c")
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
    if(is.null(Suse) == FALSE){
      lifeUSE <- function(theta) {
        exp(theta[ishift+2])*(Suse^-theta[ishift+1])
      }
      lifeUSEVAR <- function(theta,VARCOV) {
        (lifeUSE(theta)^2)*VARCOV[ishift+2,ishift+2] +
          (log(Suse)^2)*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+1] -
          2*log(Suse)*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+2]
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
    if(is.null(Suse) == FALSE){
      lifeUSE <- function(theta) {
        1/(exp(theta[ishift+2])*(Suse^theta[ishift+1]))
      }
      lifeUSEVAR <- function(theta,VARCOV) {
        (lifeUSE(theta)^2)*VARCOV[ishift+2,ishift+2] +
          (log(Suse)^2)*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+1] +
          2*log(Suse)*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+2]
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
    if(is.null(Suse) == FALSE){
      lifeUSE <- function(theta) {
        theta[ishift+1]*log(Suse) + theta[ishift+2]
      }
      lifeUSEVAR <- function(theta,VARCOV) {
        VARCOV[ishift+2,ishift+2] +
          (log(Suse)^2)*VARCOV[ishift+1,ishift+1] +
          2*log(Suse)*VARCOV[ishift+1,ishift+2]
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
        eqn1<-exp(theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+6]*SF[,5])
      }
      if(dim(SF)[2]==6){
        eqn1<-exp(theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+6]*SF[,5]+theta[ishift+7]*SF[,6])
      }
      if(dim(SF)[2]==7){
        eqn1<-exp(theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+6]*SF[,5]+theta[ishift+7]*SF[,6]+theta[ishift+8]*SF[,7])
      }
      if(dim(SF)[2]==8){
        eqn1<-exp(theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+6]*SF[,5]+theta[ishift+7]*SF[,6]+theta[ishift+8]*SF[,7]+theta[ishift+9]*SF[,8])
      }
      if(dim(SF)[2]==9){
        eqn1<-exp(theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+6]*SF[,5]+theta[ishift+7]*SF[,6]+theta[ishift+8]*SF[,7]+theta[ishift+9]*SF[,8]+theta[ishift+10]*SF[,9])
      }
      if(dim(SF)[2]==10){
        eqn1<-exp(theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+6]*SF[,5]+theta[ishift+7]*SF[,6]+theta[ishift+8]*SF[,7]+theta[ishift+9]*SF[,8]+theta[ishift+10]*SF[,9]+theta[ishift+11]*SF[,10])
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
        eqn2<-theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+6]*SF[,5]
      }
      if(dim(SF)[2]==6){
        eqn2<-theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+6]*SF[,5]+theta[ishift+7]*SF[,6]
      }
      if(dim(SF)[2]==7){
        eqn2<-theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+6]*SF[,5]+theta[ishift+7]*SF[,6]+theta[ishift+8]*SF[,7]
      }
      if(dim(SF)[2]==8){
        eqn2<-theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+6]*SF[,5]+theta[ishift+7]*SF[,6]+theta[ishift+8]*SF[,7]+theta[ishift+9]*SF[,8]
      }
      if(dim(SF)[2]==9){
        eqn2<-theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+6]*SF[,5]+theta[ishift+7]*SF[,6]+theta[ishift+8]*SF[,7]+theta[ishift+9]*SF[,8]+theta[ishift+10]*SF[,9]
      }
      if(dim(SF)[2]==10){
        eqn2<-theta[ishift+1]+theta[ishift+2]*SF[,1]+theta[ishift+3]*SF[,2]+theta[ishift+4]*SF[,3]+theta[ishift+5]*SF[,4]+theta[ishift+6]*SF[,5]+theta[ishift+7]*SF[,6]+theta[ishift+8]*SF[,7]+theta[ishift+9]*SF[,8]+theta[ishift+10]*SF[,9]+theta[ishift+11]*SF[,10]
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
          eqn3<-exp(theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+6]*Sc[,5])
        }
        if(dim(Sc)[2]==6){
          eqn3<-exp(theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+6]*Sc[,5]+theta[ishift+7]*Sc[,6])
        }
        if(dim(Sc)[2]==7){
          eqn3<-exp(theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+6]*Sc[,5]+theta[ishift+7]*Sc[,6]+theta[ishift+8]*Sc[,7])
        }
        if(dim(Sc)[2]==8){
          eqn3<-exp(theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+6]*Sc[,5]+theta[ishift+7]*Sc[,6]+theta[ishift+8]*Sc[,7]+theta[ishift+9]*Sc[,8])
        }
        if(dim(Sc)[2]==9){
          eqn3<-exp(theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+6]*Sc[,5]+theta[ishift+7]*Sc[,6]+theta[ishift+8]*Sc[,7]+theta[ishift+9]*Sc[,8]+theta[ishift+10]*Sc[,9])
        }
        if(dim(Sc)[2]==10){
          eqn3<-exp(theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+6]*Sc[,5]+theta[ishift+7]*Sc[,6]+theta[ishift+8]*Sc[,7]+theta[ishift+9]*Sc[,8]+theta[ishift+10]*Sc[,9]+theta[ishift+11]*Sc[,10])
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
          eqn4<-theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+6]*Sc[,5]
        }
        if(dim(Sc)[2]==6){
          eqn4<-theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+6]*Sc[,5]+theta[ishift+7]*Sc[,6]
        }
        if(dim(Sc)[2]==7){
          eqn4<-theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+6]*Sc[,5]+theta[ishift+7]*Sc[,6]+theta[ishift+8]*Sc[,7]
        }
        if(dim(Sc)[2]==8){
          eqn4<-theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+6]*Sc[,5]+theta[ishift+7]*Sc[,6]+theta[ishift+8]*Sc[,7]+theta[ishift+9]*Sc[,8]
        }
        if(dim(Sc)[2]==9){
          eqn4<-theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+6]*Sc[,5]+theta[ishift+7]*Sc[,6]+theta[ishift+8]*Sc[,7]+theta[ishift+9]*Sc[,8]+theta[ishift+10]*Sc[,9]
        }
        if(dim(Sc)[2]==10){
          eqn4<-theta[ishift+1]+theta[ishift+2]*Sc[,1]+theta[ishift+3]*Sc[,2]+theta[ishift+4]*Sc[,3]+theta[ishift+5]*Sc[,4]+theta[ishift+6]*Sc[,5]+theta[ishift+7]*Sc[,6]+theta[ishift+8]*Sc[,7]+theta[ishift+9]*Sc[,8]+theta[ishift+10]*Sc[,9]+theta[ishift+11]*Sc[,10]
        }
        return(eqn4)
      }
    }
    if(is.null(Suse) == FALSE){
      lifeUSE <- function(theta) {
        if(length(Suse)==2){
          eqn5<-exp(theta[ishift+1]+theta[ishift+2]*Suse[1]+theta[ishift+3]*Suse[2])
        }
        if(length(Suse)==3){
          eqn5<-exp(theta[ishift+1]+theta[ishift+2]*Suse[1]+theta[ishift+3]*Suse[2]+theta[ishift+4]*Suse[3])
        }
        if(length(Suse)==4){
          eqn5<-exp(theta[ishift+1]+theta[ishift+2]*Suse[1]+theta[ishift+3]*Suse[2]+theta[ishift+4]*Suse[3]+theta[ishift+5]*Suse[4])
        }
        if(length(Suse)==5){
          eqn5<-exp(theta[ishift+1]+theta[ishift+2]*Suse[1]+theta[ishift+3]*Suse[2]+theta[ishift+4]*Suse[3]+theta[ishift+5]*Suse[4]+theta[ishift+6]*Suse[5])
        }
        if(length(Suse)==6){
          eqn5<-exp(theta[ishift+1]+theta[ishift+2]*Suse[1]+theta[ishift+3]*Suse[2]+theta[ishift+4]*Suse[3]+theta[ishift+5]*Suse[4]+theta[ishift+6]*Suse[5]+theta[ishift+7]*Suse[6])
        }
        if(length(Suse)==7){
          eqn5<-exp(theta[ishift+1]+theta[ishift+2]*Suse[1]+theta[ishift+3]*Suse[2]+theta[ishift+4]*Suse[3]+theta[ishift+5]*Suse[4]+theta[ishift+6]*Suse[5]+theta[ishift+7]*Suse[6]+theta[ishift+8]*Suse[7])
        }
        if(length(Suse)==8){
          eqn5<-exp(theta[ishift+1]+theta[ishift+2]*Suse[1]+theta[ishift+3]*Suse[2]+theta[ishift+4]*Suse[3]+theta[ishift+5]*Suse[4]+theta[ishift+6]*Suse[5]+theta[ishift+7]*Suse[6]+theta[ishift+8]*Suse[7]+theta[ishift+9]*Suse[8])
        }
        if(length(Suse)==9){
          eqn5<-exp(theta[ishift+1]+theta[ishift+2]*Suse[1]+theta[ishift+3]*Suse[2]+theta[ishift+4]*Suse[3]+theta[ishift+5]*Suse[4]+theta[ishift+6]*Suse[5]+theta[ishift+7]*Suse[6]+theta[ishift+8]*Suse[7]+theta[ishift+9]*Suse[8]+theta[ishift+10]*Suse[9])
        }
        if(length(Suse)==10){
          eqn5<-exp(theta[ishift+1]+theta[ishift+2]*Suse[1]+theta[ishift+3]*Suse[2]+theta[ishift+4]*Suse[3]+theta[ishift+5]*Suse[4]+theta[ishift+6]*Suse[5]+theta[ishift+7]*Suse[6]+theta[ishift+8]*Suse[7]+theta[ishift+9]*Suse[8]+theta[ishift+10]*Suse[9]+theta[ishift+11]*Suse[10])
        }
        return(eqn5)
      }
      lifeUSEVAR <- function(theta,VARCOV) {
        VARCOV[ishift+2,ishift+2] +
          (log(Suse)^2)*VARCOV[ishift+1,ishift+1] +
          2*log(Suse)*VARCOV[ishift+1,ishift+2]
        if(length(Suse)==2){
          eqn6<-(Suse[2]^2)*VARCOV[ishift+3,ishift+3] +
            (Suse[1]^2)*VARCOV[ishift+2,ishift+2] +
            VARCOV[ishift+1,ishift+1] +
            2*Suse[1]*VARCOV[ishift+1,ishift+2] + 2*Suse[2]*VARCOV[ishift+1,ishift+3] +
            2*Suse[1]*Suse[2]*VARCOV[ishift+2,ishift+3]
        }
        if(length(Suse)==3){
          eqn6<-(Suse[3]^2)*VARCOV[ishift+4,ishift+4] +
            (Suse[2]^2)*VARCOV[ishift+3,ishift+3] +
            (Suse[1]^2)*VARCOV[ishift+2,ishift+2] +
            VARCOV[ishift+1,ishift+1] +
            2*Suse[1]*VARCOV[ishift+1,ishift+2] + 2*Suse[2]*VARCOV[ishift+1,ishift+3] + 2*Suse[3]*VARCOV[ishift+1,ishift+4] +
            2*Suse[1]*Suse[2]*VARCOV[ishift+2,ishift+3] + 2*Suse[1]*Suse[3]*VARCOV[ishift+2,ishift+4] +
            2*Suse[2]*Suse[3]*VARCOV[ishift+3,ishift+4]
        }
        if(length(Suse)==4){
          eqn6<-(Suse[4]^2)*VARCOV[ishift+5,ishift+5] +
            (Suse[3]^2)*VARCOV[ishift+4,ishift+4] +
            (Suse[2]^2)*VARCOV[ishift+3,ishift+3] +
            (Suse[1]^2)*VARCOV[ishift+2,ishift+2] +
            VARCOV[ishift+1,ishift+1] +
            2*Suse[1]*VARCOV[ishift+1,ishift+2] +
            2*Suse[2]*VARCOV[ishift+1,ishift+3] +
            2*Suse[3]*VARCOV[ishift+1,ishift+4] +
            2*Suse[4]*VARCOV[ishift+1,ishift+5] +
            2*Suse[1]*Suse[2]*VARCOV[ishift+2,ishift+3] +
            2*Suse[1]*Suse[3]*VARCOV[ishift+2,ishift+4] +
            2*Suse[1]*Suse[4]*VARCOV[ishift+2,ishift+5] +
            2*Suse[2]*Suse[3]*VARCOV[ishift+3,ishift+4] +
            2*Suse[2]*Suse[4]*VARCOV[ishift+3,ishift+5] +
            2*Suse[3]*Suse[4]*VARCOV[ishift+4,ishift+5]
        }
        if(length(Suse)==5){
          eqn6<-(Suse[5]^2)*VARCOV[ishift+6,ishift+6] +
            (Suse[4]^2)*VARCOV[ishift+5,ishift+5] +
            (Suse[3]^2)*VARCOV[ishift+4,ishift+4] +
            (Suse[2]^2)*VARCOV[ishift+3,ishift+3] +
            (Suse[1]^2)*VARCOV[ishift+2,ishift+2] +
            VARCOV[ishift+1,ishift+1] +
            2*Suse[1]*VARCOV[ishift+1,ishift+2] +
            2*Suse[2]*VARCOV[ishift+1,ishift+3] +
            2*Suse[3]*VARCOV[ishift+1,ishift+4] +
            2*Suse[4]*VARCOV[ishift+1,ishift+5] +
            2*Suse[5]*VARCOV[ishift+1,ishift+6] +
            2*Suse[1]*Suse[2]*VARCOV[ishift+2,ishift+3] +
            2*Suse[1]*Suse[3]*VARCOV[ishift+2,ishift+4] +
            2*Suse[1]*Suse[4]*VARCOV[ishift+2,ishift+5] +
            2*Suse[1]*Suse[5]*VARCOV[ishift+2,ishift+6] +
            2*Suse[2]*Suse[3]*VARCOV[ishift+3,ishift+4] +
            2*Suse[2]*Suse[4]*VARCOV[ishift+3,ishift+5] +
            2*Suse[2]*Suse[5]*VARCOV[ishift+3,ishift+6] +
            2*Suse[3]*Suse[4]*VARCOV[ishift+4,ishift+5] +
            2*Suse[3]*Suse[5]*VARCOV[ishift+4,ishift+6] +
            2*Suse[4]*Suse[5]*VARCOV[ishift+5,ishift+6]
        }
        if(length(Suse)==6){
          eqn6<-(Suse[6]^2)*VARCOV[ishift+7,ishift+7] +
            (Suse[5]^2)*VARCOV[ishift+6,ishift+6] +
            (Suse[4]^2)*VARCOV[ishift+5,ishift+5] +
            (Suse[3]^2)*VARCOV[ishift+4,ishift+4] +
            (Suse[2]^2)*VARCOV[ishift+3,ishift+3] +
            (Suse[1]^2)*VARCOV[ishift+2,ishift+2] +
            VARCOV[ishift+1,ishift+1] +
            2*Suse[1]*VARCOV[ishift+1,ishift+2] +
            2*Suse[2]*VARCOV[ishift+1,ishift+3] +
            2*Suse[3]*VARCOV[ishift+1,ishift+4] +
            2*Suse[4]*VARCOV[ishift+1,ishift+5] +
            2*Suse[5]*VARCOV[ishift+1,ishift+6] +
            2*Suse[6]*VARCOV[ishift+1,ishift+7] +
            2*Suse[1]*Suse[2]*VARCOV[ishift+2,ishift+3] +
            2*Suse[1]*Suse[3]*VARCOV[ishift+2,ishift+4] +
            2*Suse[1]*Suse[4]*VARCOV[ishift+2,ishift+5] +
            2*Suse[1]*Suse[5]*VARCOV[ishift+2,ishift+6] +
            2*Suse[1]*Suse[6]*VARCOV[ishift+2,ishift+7] +
            2*Suse[2]*Suse[3]*VARCOV[ishift+3,ishift+4] +
            2*Suse[2]*Suse[4]*VARCOV[ishift+3,ishift+5] +
            2*Suse[2]*Suse[5]*VARCOV[ishift+3,ishift+6] +
            2*Suse[2]*Suse[6]*VARCOV[ishift+3,ishift+7] +
            2*Suse[3]*Suse[4]*VARCOV[ishift+4,ishift+5] +
            2*Suse[3]*Suse[5]*VARCOV[ishift+4,ishift+6] +
            2*Suse[3]*Suse[6]*VARCOV[ishift+4,ishift+7] +
            2*Suse[4]*Suse[5]*VARCOV[ishift+5,ishift+6] +
            2*Suse[4]*Suse[6]*VARCOV[ishift+5,ishift+7] +
            2*Suse[5]*Suse[6]*VARCOV[ishift+6,ishift+7]
        }
        if(length(Suse)==7){
          eqn6<-(Suse[7]^2)*VARCOV[ishift+8,ishift+8] +
            (Suse[6]^2)*VARCOV[ishift+7,ishift+7] +
            (Suse[5]^2)*VARCOV[ishift+6,ishift+6] +
            (Suse[4]^2)*VARCOV[ishift+5,ishift+5] +
            (Suse[3]^2)*VARCOV[ishift+4,ishift+4] +
            (Suse[2]^2)*VARCOV[ishift+3,ishift+3] +
            (Suse[1]^2)*VARCOV[ishift+2,ishift+2] +
            VARCOV[ishift+1,ishift+1] +
            2*Suse[1]*VARCOV[ishift+1,ishift+2] +
            2*Suse[2]*VARCOV[ishift+1,ishift+3] +
            2*Suse[3]*VARCOV[ishift+1,ishift+4] +
            2*Suse[4]*VARCOV[ishift+1,ishift+5] +
            2*Suse[5]*VARCOV[ishift+1,ishift+6] +
            2*Suse[6]*VARCOV[ishift+1,ishift+7] +
            2*Suse[7]*VARCOV[ishift+1,ishift+8] +
            2*Suse[1]*Suse[2]*VARCOV[ishift+2,ishift+3] +
            2*Suse[1]*Suse[3]*VARCOV[ishift+2,ishift+4] +
            2*Suse[1]*Suse[4]*VARCOV[ishift+2,ishift+5] +
            2*Suse[1]*Suse[5]*VARCOV[ishift+2,ishift+6] +
            2*Suse[1]*Suse[6]*VARCOV[ishift+2,ishift+7] +
            2*Suse[1]*Suse[7]*VARCOV[ishift+2,ishift+8] +
            2*Suse[2]*Suse[3]*VARCOV[ishift+3,ishift+4] +
            2*Suse[2]*Suse[4]*VARCOV[ishift+3,ishift+5] +
            2*Suse[2]*Suse[5]*VARCOV[ishift+3,ishift+6] +
            2*Suse[2]*Suse[6]*VARCOV[ishift+3,ishift+7] +
            2*Suse[2]*Suse[7]*VARCOV[ishift+3,ishift+8] +
            2*Suse[3]*Suse[4]*VARCOV[ishift+4,ishift+5] +
            2*Suse[3]*Suse[5]*VARCOV[ishift+4,ishift+6] +
            2*Suse[3]*Suse[6]*VARCOV[ishift+4,ishift+7] +
            2*Suse[3]*Suse[7]*VARCOV[ishift+4,ishift+8] +
            2*Suse[4]*Suse[5]*VARCOV[ishift+5,ishift+6] +
            2*Suse[4]*Suse[6]*VARCOV[ishift+5,ishift+7] +
            2*Suse[4]*Suse[7]*VARCOV[ishift+5,ishift+8] +
            2*Suse[5]*Suse[6]*VARCOV[ishift+6,ishift+7] +
            2*Suse[5]*Suse[7]*VARCOV[ishift+6,ishift+8] +
            2*Suse[6]*Suse[7]*VARCOV[ishift+7,ishift+8]
        }
        if(length(Suse)==8){
          eqn6<-(Suse[8]^2)*VARCOV[ishift+9,ishift+9] +
            (Suse[7]^2)*VARCOV[ishift+8,ishift+8] +
            (Suse[6]^2)*VARCOV[ishift+7,ishift+7] +
            (Suse[5]^2)*VARCOV[ishift+6,ishift+6] +
            (Suse[4]^2)*VARCOV[ishift+5,ishift+5] +
            (Suse[3]^2)*VARCOV[ishift+4,ishift+4] +
            (Suse[2]^2)*VARCOV[ishift+3,ishift+3] +
            (Suse[1]^2)*VARCOV[ishift+2,ishift+2] +
            VARCOV[ishift+1,ishift+1] +
            2*Suse[1]*VARCOV[ishift+1,ishift+2] +
            2*Suse[2]*VARCOV[ishift+1,ishift+3] +
            2*Suse[3]*VARCOV[ishift+1,ishift+4] +
            2*Suse[4]*VARCOV[ishift+1,ishift+5] +
            2*Suse[5]*VARCOV[ishift+1,ishift+6] +
            2*Suse[6]*VARCOV[ishift+1,ishift+7] +
            2*Suse[7]*VARCOV[ishift+1,ishift+8] +
            2*Suse[8]*VARCOV[ishift+1,ishift+9] +
            2*Suse[1]*Suse[2]*VARCOV[ishift+2,ishift+3] +
            2*Suse[1]*Suse[3]*VARCOV[ishift+2,ishift+4] +
            2*Suse[1]*Suse[4]*VARCOV[ishift+2,ishift+5] +
            2*Suse[1]*Suse[5]*VARCOV[ishift+2,ishift+6] +
            2*Suse[1]*Suse[6]*VARCOV[ishift+2,ishift+7] +
            2*Suse[1]*Suse[7]*VARCOV[ishift+2,ishift+8] +
            2*Suse[1]*Suse[8]*VARCOV[ishift+2,ishift+9] +
            2*Suse[2]*Suse[3]*VARCOV[ishift+3,ishift+4] +
            2*Suse[2]*Suse[4]*VARCOV[ishift+3,ishift+5] +
            2*Suse[2]*Suse[5]*VARCOV[ishift+3,ishift+6] +
            2*Suse[2]*Suse[6]*VARCOV[ishift+3,ishift+7] +
            2*Suse[2]*Suse[7]*VARCOV[ishift+3,ishift+8] +
            2*Suse[2]*Suse[8]*VARCOV[ishift+3,ishift+9] +
            2*Suse[3]*Suse[4]*VARCOV[ishift+4,ishift+5] +
            2*Suse[3]*Suse[5]*VARCOV[ishift+4,ishift+6] +
            2*Suse[3]*Suse[6]*VARCOV[ishift+4,ishift+7] +
            2*Suse[3]*Suse[7]*VARCOV[ishift+4,ishift+8] +
            2*Suse[3]*Suse[8]*VARCOV[ishift+4,ishift+9] +
            2*Suse[4]*Suse[5]*VARCOV[ishift+5,ishift+6] +
            2*Suse[4]*Suse[6]*VARCOV[ishift+5,ishift+7] +
            2*Suse[4]*Suse[7]*VARCOV[ishift+5,ishift+8] +
            2*Suse[4]*Suse[8]*VARCOV[ishift+5,ishift+9] +
            2*Suse[5]*Suse[6]*VARCOV[ishift+6,ishift+7] +
            2*Suse[5]*Suse[7]*VARCOV[ishift+6,ishift+8] +
            2*Suse[5]*Suse[8]*VARCOV[ishift+6,ishift+9] +
            2*Suse[6]*Suse[7]*VARCOV[ishift+7,ishift+8] +
            2*Suse[6]*Suse[8]*VARCOV[ishift+7,ishift+9] +
            2*Suse[7]*Suse[8]*VARCOV[ishift+8,ishift+9]
        }
        if(length(Suse)==9){
          eqn6<-(Suse[9]^2)*VARCOV[ishift+10,ishift+10] +
            (Suse[8]^2)*VARCOV[ishift+9,ishift+9] +
            (Suse[7]^2)*VARCOV[ishift+8,ishift+8] +
            (Suse[6]^2)*VARCOV[ishift+7,ishift+7] +
            (Suse[5]^2)*VARCOV[ishift+6,ishift+6] +
            (Suse[4]^2)*VARCOV[ishift+5,ishift+5] +
            (Suse[3]^2)*VARCOV[ishift+4,ishift+4] +
            (Suse[2]^2)*VARCOV[ishift+3,ishift+3] +
            (Suse[1]^2)*VARCOV[ishift+2,ishift+2] +
            VARCOV[ishift+1,ishift+1] +
            2*Suse[1]*VARCOV[ishift+1,ishift+2] +
            2*Suse[2]*VARCOV[ishift+1,ishift+3] +
            2*Suse[3]*VARCOV[ishift+1,ishift+4] +
            2*Suse[4]*VARCOV[ishift+1,ishift+5] +
            2*Suse[5]*VARCOV[ishift+1,ishift+6] +
            2*Suse[6]*VARCOV[ishift+1,ishift+7] +
            2*Suse[7]*VARCOV[ishift+1,ishift+8] +
            2*Suse[8]*VARCOV[ishift+1,ishift+9] +
            2*Suse[9]*VARCOV[ishift+1,ishift+10] +
            2*Suse[1]*Suse[2]*VARCOV[ishift+2,ishift+3] +
            2*Suse[1]*Suse[3]*VARCOV[ishift+2,ishift+4] +
            2*Suse[1]*Suse[4]*VARCOV[ishift+2,ishift+5] +
            2*Suse[1]*Suse[5]*VARCOV[ishift+2,ishift+6] +
            2*Suse[1]*Suse[6]*VARCOV[ishift+2,ishift+7] +
            2*Suse[1]*Suse[7]*VARCOV[ishift+2,ishift+8] +
            2*Suse[1]*Suse[8]*VARCOV[ishift+2,ishift+9] +
            2*Suse[1]*Suse[9]*VARCOV[ishift+2,ishift+10] +
            2*Suse[2]*Suse[3]*VARCOV[ishift+3,ishift+4] +
            2*Suse[2]*Suse[4]*VARCOV[ishift+3,ishift+5] +
            2*Suse[2]*Suse[5]*VARCOV[ishift+3,ishift+6] +
            2*Suse[2]*Suse[6]*VARCOV[ishift+3,ishift+7] +
            2*Suse[2]*Suse[7]*VARCOV[ishift+3,ishift+8] +
            2*Suse[2]*Suse[8]*VARCOV[ishift+3,ishift+9] +
            2*Suse[2]*Suse[9]*VARCOV[ishift+3,ishift+10] +
            2*Suse[3]*Suse[4]*VARCOV[ishift+4,ishift+5] +
            2*Suse[3]*Suse[5]*VARCOV[ishift+4,ishift+6] +
            2*Suse[3]*Suse[6]*VARCOV[ishift+4,ishift+7] +
            2*Suse[3]*Suse[7]*VARCOV[ishift+4,ishift+8] +
            2*Suse[3]*Suse[8]*VARCOV[ishift+4,ishift+9] +
            2*Suse[3]*Suse[9]*VARCOV[ishift+4,ishift+10] +
            2*Suse[4]*Suse[5]*VARCOV[ishift+5,ishift+6] +
            2*Suse[4]*Suse[6]*VARCOV[ishift+5,ishift+7] +
            2*Suse[4]*Suse[7]*VARCOV[ishift+5,ishift+8] +
            2*Suse[4]*Suse[8]*VARCOV[ishift+5,ishift+9] +
            2*Suse[4]*Suse[9]*VARCOV[ishift+5,ishift+10] +
            2*Suse[5]*Suse[6]*VARCOV[ishift+6,ishift+7] +
            2*Suse[5]*Suse[7]*VARCOV[ishift+6,ishift+8] +
            2*Suse[5]*Suse[8]*VARCOV[ishift+6,ishift+9] +
            2*Suse[6]*Suse[7]*VARCOV[ishift+7,ishift+8] +
            2*Suse[6]*Suse[8]*VARCOV[ishift+7,ishift+9] +
            2*Suse[6]*Suse[9]*VARCOV[ishift+7,ishift+10] +
            2*Suse[7]*Suse[8]*VARCOV[ishift+8,ishift+9] +
            2*Suse[7]*Suse[9]*VARCOV[ishift+8,ishift+10] +
            2*Suse[8]*Suse[9]*VARCOV[ishift+9,ishift+10]
        }
        if(length(Suse)==10){
          eqn6<-(Suse[10]^2)*VARCOV[ishift+11,ishift+11] +
            (Suse[9]^2)*VARCOV[ishift+10,ishift+10] +
            (Suse[8]^2)*VARCOV[ishift+9,ishift+9] +
            (Suse[7]^2)*VARCOV[ishift+8,ishift+8] +
            (Suse[6]^2)*VARCOV[ishift+7,ishift+7] +
            (Suse[5]^2)*VARCOV[ishift+6,ishift+6] +
            (Suse[4]^2)*VARCOV[ishift+5,ishift+5] +
            (Suse[3]^2)*VARCOV[ishift+4,ishift+4] +
            (Suse[2]^2)*VARCOV[ishift+3,ishift+3] +
            (Suse[1]^2)*VARCOV[ishift+2,ishift+2] +
            VARCOV[ishift+1,ishift+1] +
            2*Suse[1]*VARCOV[ishift+1,ishift+2] +
            2*Suse[2]*VARCOV[ishift+1,ishift+3] +
            2*Suse[3]*VARCOV[ishift+1,ishift+4] +
            2*Suse[4]*VARCOV[ishift+1,ishift+5] +
            2*Suse[5]*VARCOV[ishift+1,ishift+6] +
            2*Suse[6]*VARCOV[ishift+1,ishift+7] +
            2*Suse[7]*VARCOV[ishift+1,ishift+8] +
            2*Suse[8]*VARCOV[ishift+1,ishift+9] +
            2*Suse[9]*VARCOV[ishift+1,ishift+10] +
            2*Suse[10]*VARCOV[ishift+1,ishift+11] +
            2*Suse[1]*Suse[2]*VARCOV[ishift+2,ishift+3] +
            2*Suse[1]*Suse[3]*VARCOV[ishift+2,ishift+4] +
            2*Suse[1]*Suse[4]*VARCOV[ishift+2,ishift+5] +
            2*Suse[1]*Suse[5]*VARCOV[ishift+2,ishift+6] +
            2*Suse[1]*Suse[6]*VARCOV[ishift+2,ishift+7] +
            2*Suse[1]*Suse[7]*VARCOV[ishift+2,ishift+8] +
            2*Suse[1]*Suse[8]*VARCOV[ishift+2,ishift+9] +
            2*Suse[1]*Suse[9]*VARCOV[ishift+2,ishift+10] +
            2*Suse[1]*Suse[10]*VARCOV[ishift+2,ishift+11] +
            2*Suse[2]*Suse[3]*VARCOV[ishift+3,ishift+4] +
            2*Suse[2]*Suse[4]*VARCOV[ishift+3,ishift+5] +
            2*Suse[2]*Suse[5]*VARCOV[ishift+3,ishift+6] +
            2*Suse[2]*Suse[6]*VARCOV[ishift+3,ishift+7] +
            2*Suse[2]*Suse[7]*VARCOV[ishift+3,ishift+8] +
            2*Suse[2]*Suse[8]*VARCOV[ishift+3,ishift+9] +
            2*Suse[2]*Suse[9]*VARCOV[ishift+3,ishift+10] +
            2*Suse[2]*Suse[10]*VARCOV[ishift+3,ishift+11] +
            2*Suse[3]*Suse[4]*VARCOV[ishift+4,ishift+5] +
            2*Suse[3]*Suse[5]*VARCOV[ishift+4,ishift+6] +
            2*Suse[3]*Suse[6]*VARCOV[ishift+4,ishift+7] +
            2*Suse[3]*Suse[7]*VARCOV[ishift+4,ishift+8] +
            2*Suse[3]*Suse[8]*VARCOV[ishift+4,ishift+9] +
            2*Suse[3]*Suse[9]*VARCOV[ishift+4,ishift+10] +
            2*Suse[3]*Suse[10]*VARCOV[ishift+4,ishift+11] +
            2*Suse[4]*Suse[5]*VARCOV[ishift+5,ishift+6] +
            2*Suse[4]*Suse[6]*VARCOV[ishift+5,ishift+7] +
            2*Suse[4]*Suse[7]*VARCOV[ishift+5,ishift+8] +
            2*Suse[4]*Suse[8]*VARCOV[ishift+5,ishift+9] +
            2*Suse[4]*Suse[9]*VARCOV[ishift+5,ishift+10] +
            2*Suse[4]*Suse[10]*VARCOV[ishift+5,ishift+11] +
            2*Suse[5]*Suse[6]*VARCOV[ishift+6,ishift+7] +
            2*Suse[5]*Suse[7]*VARCOV[ishift+6,ishift+8] +
            2*Suse[5]*Suse[8]*VARCOV[ishift+6,ishift+9] +
            2*Suse[6]*Suse[7]*VARCOV[ishift+7,ishift+8] +
            2*Suse[6]*Suse[8]*VARCOV[ishift+7,ishift+9] +
            2*Suse[6]*Suse[9]*VARCOV[ishift+7,ishift+10] +
            2*Suse[6]*Suse[10]*VARCOV[ishift+7,ishift+11] +
            2*Suse[7]*Suse[8]*VARCOV[ishift+8,ishift+9] +
            2*Suse[7]*Suse[9]*VARCOV[ishift+8,ishift+10] +
            2*Suse[8]*Suse[9]*VARCOV[ishift+9,ishift+10] +
            2*Suse[8]*Suse[10]*VARCOV[ishift+9,ishift+11] +
            2*Suse[9]*Suse[10]*VARCOV[ishift+10,ishift+11]
        }
        return(eqn6)
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
    if(is.null(Suse) == FALSE){
      lifeUSE <- function(theta) {
        exp(theta[ishift+1])*exp((theta[ishift+2]/Suse[1]) + (theta[ishift+3]/Suse[2]))
      }
      lifeUSEVAR <- function(theta,VARCOV) {
        (1/Suse[2]^2)*(lifeUSE(theta)^2)*VARCOV[ishift+3,ishift+3] +
          (1/Suse[1]^2)*(lifeUSE(theta)^2)*VARCOV[ishift+2,ishift+2] +
          (lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+1] +
          2*(1/Suse[1])*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+2] +
          2*(1/Suse[2])*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+3] +
          2*(1/Suse[1])*(1/Suse[2])*(lifeUSE(theta)^2)*VARCOV[ishift+2,ishift+3]
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
    if(is.null(Suse) == FALSE){
      lifeUSE <- function(theta) {
        exp(theta[ishift+3])/((Suse[2]^theta[ishift+2])*exp(-theta[ishift+1]/Suse[1]))
      }
      lifeUSEVAR <- function(theta,VARCOV) {
        c(c(1/Suse[1]*lifeUSE(theta),-log(Suse[2])*lifeUSE(theta),lifeUSE(theta))%*%
          VARCOV[((ishift+1):(ishift+3)),((ishift+1):(ishift+3))]%*%
          c(1/Suse[1]*lifeUSE(theta),-log(Suse[2])*lifeUSE(theta),lifeUSE(theta)))
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
    if(is.null(Suse) == FALSE){
      lifeUSE <- function(theta) {
        (1/Suse[1])*exp((theta[ishift+1] + (theta[ishift+2]/Suse[1])) + (theta[ishift+3] + (theta[ishift+4]/Suse[1]))*Suse[2])
      }
      lifeUSEVAR <- function(theta,VARCOV) {
        ((Suse[2]/Suse[1])^2)*(lifeUSE(theta)^2)*VARCOV[ishift+4,ishift+4] +
          (Suse[2]^2)*VARCOV[ishift+3,ishift+3] +
          (1/Suse[1]^2)*(lifeUSE(theta)^2)*VARCOV[ishift+2,ishift+2] +
          (lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+1] +
          2*(1/Suse[1])*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+2] +
          2*Suse[2]*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+3] +
          2*(Suse[2]/Suse[1])*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+3] +
          2*(Suse[2]/Suse[1])*(lifeUSE(theta)^2)*VARCOV[ishift+2,ishift+3] +
          2*(Suse[2]/(Suse[1]^2))*(lifeUSE(theta)^2)*VARCOV[ishift+2,ishift+4] +
          2*((Suse[2]^2)/Suse[1])*(lifeUSE(theta)^2)*VARCOV[ishift+3,ishift+4]
      }
    }
    ls_txt<-"Eyring (Type 3)"
    params_txt<-c("a","b","c","d")
  }

  if (ls=="Eyring4") {
    # lsparams[1] - parameter A, lsparams[2] - parameter b
    # lsparams[3] - parameter Ea
    # shift A to log A
    # positivity_v[ishift+1]<-1
    LSQest[ishift+1] <- log(LSQest[ishift+1])
    # Temperature HAS to be in Kelvin for this to work
    K<-8.617385e-5
    lifeF <- function(theta) {
      exp(theta[ishift+1])*exp(theta[ishift+3]/(K*SF[,1]))*(SF[,2]^-theta[ishift+2])
    }
    loglifeF <- function(theta) {
      theta[ishift+1] + theta[ishift+3]/(K*SF[,1]) - theta[ishift+2]*log(SF[,2])
    }
    if(is.null(Tc)==FALSE){
      lifeC <- function(theta) {
        exp(theta[ishift+1])*exp(theta[ishift+3]/(K*Sc[,1]))*(Sc[,2]^-theta[ishift+2])
      }
      loglifeC <- function(theta) {
        theta[ishift+1] + theta[ishift+3]/(K*Sc[,1]) - theta[ishift+2]*log(Sc[,2])
      }
    }
    if(is.null(Suse) == FALSE){
      lifeUSE <- function(theta) {
        exp(theta[ishift+1])*exp(theta[ishift+3]/(K*Suse[1]))*(Suse[2]^-theta[ishift+2])
      }
      lifeUSEVAR <- function(theta,VARCOV) {
        (1/(K*Suse[1])^2)*(lifeUSE(theta)^2)*VARCOV[ishift+3,ishift+3] +
          (log(Suse[2])^2)*(lifeUSE(theta)^2)*VARCOV[ishift+2,ishift+2] +
          (lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+1] -
          2*log(Suse[2])*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+2] +
          2*(1/(K*Suse[2]))*(lifeUSE(theta)^2)*VARCOV[ishift+1,ishift+3] -
          2*(1/(K*Suse[1]))*log(Suse[2])*(lifeUSE(theta)^2)*VARCOV[ishift+2,ishift+3]
      }
    }
    ls_txt<-"Eyring (Type 4)"
    params_txt<-c("A","b","E_a")
  }


  # UPDATE (10/22/2024) - Form data matrix for MLE probability plot
  if(is.null(Tc) == TRUE){
    data <- cbind(TTF,rep(1,length(TTF)),SF)
  }
  if(is.null(Tc) == FALSE && is.null(dim(Sc)) == TRUE){
    data <- cbind(c(TTF,Tc),c(rep(1,length(TTF)),rep(0,length(Tc))),c(SF,Sc))
  }
  if(is.null(Tc) == FALSE && is.null(dim(Sc)) == FALSE){
    data <- cbind(c(TTF,Tc),c(rep(1,length(TTF)),rep(0,length(Tc))),rbind(SF,Sc))
  }

  # Fit to log-likelihood distributions
  if (dist=="Weibull") {
    # For standard beta estimation
    if(is.null(param2) == TRUE){
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
    }
    # For when beta is known
    if(is.null(param2) == FALSE){
      # Store full LSQest for later
      LSQest0 <- LSQest
      # Set new LSQest for everything but the known beta
      LSQest <- LSQest0[2:length(LSQest0)]
      # First redefine beta parameter as logbeta
      if(is.null(Tc)){
        loglik <- function(theta){
          -sum(log(param2) + (param2-1)*log(TTF) - param2*loglifeF(theta) - ((TTF/lifeF(theta))^param2))
        }
      } else{
        loglik <- function(theta){
          -sum(log(param2) + (param2-1)*log(TTF) - param2*loglifeF(theta) - ((TTF/lifeF(theta))^param2)) - sum(- ((Tc/lifeC(theta))^param2))
        }
      }
    }

    plotoutput <- probplot.wbl(data,pp,xlabel1,MLE_i = 1,stressunit1 = stressunit1,stressunit2 = stressunit2)$prob_plot
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
    plotoutput <- probplot.wbl3P(data,pp,xlabel1,MLE_i = 1,stressunit1 = stressunit1,stressunit2 = stressunit2)$prob_plot
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
    plotoutput <- probplot.logn(data,pp,xlabel1,MLE_i = 1,stressunit1 = stressunit1,stressunit2 = stressunit2)$prob_plot
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
    plotoutput <- probplot.nor(data,pp,xlabel1,MLE_i = 1,stressunit1 = stressunit1,stressunit2 = stressunit2)$prob_plot
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
    plotoutput <- probplot.exp(data,pp,xlabel1,MLE_i = 1,stressunit1 = stressunit1,stressunit2 = stressunit2)$prob_plot
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
    plotoutput <- probplot.exp2P(data,pp,xlabel1,MLE_i = 1,stressunit1 = stressunit1,stressunit2 = stressunit2)$prob_plot
    dist_txt<-"Two-Parameter Exponential"
    distparam_txt<-"\U03C3"
  }
  if (dist=="Gamma") {
    # Re-parameterize alpha and beta to mu and lambda
    LSQest0 <- LSQest
    # mu = ln beta + ln alpha (-Inf to Inf)
    # LSQest0[1] <- log(LSQest[2]) + log(LSQest[1])
    # lambda = 1/sqrt(alpha) (positive so make loglambda)
    if(is.null(Tc)){
      LSQest0[1] <- mean(log(1/sqrt(exp(loglifeF(LSQest) - log(LSQest[1])))))
    } else{
      LSQest0[1] <- mean(log(1/sqrt(exp((loglifeF(LSQest) + loglifeC(LSQest)) - log(LSQest[1])))))
    }
    LSQest <- LSQest0

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(-log(gamma(1/(exp(theta[1]))^2)) - (1/(exp(theta[1]))^2)*(loglifeF(theta) + log((exp(theta[1]))^2)) + ((1/(exp(theta[1]))^2) - 1)*log(TTF) - (1/(exp(theta[1]))^2)*exp(log(TTF) - loglifeF(theta)))
      }
    } else{
      loglik <- function(theta){
        -sum(-log(gamma(1/(exp(theta[1]))^2)) - (1/(exp(theta[1]))^2)*(loglifeF(theta) + log((exp(theta[1]))^2)) + ((1/(exp(theta[1]))^2) - 1)*log(TTF) - (1/(exp(theta[1]))^2)*exp(log(TTF) - loglifeF(theta))) - sum(log(1 - Rgamma(1/(exp(theta[1]))^2,exp(log(Tc) - loglifeC(theta))/(exp(theta[1]))^2),lower=TRUE))
      }
    }
  }
  if (dist=="3PGamma") {
    # Re-parameterize alpha, beta, and gamma to mu, lambda, and sigma
    LSQest0 <- c(0,0,0)
    # mu = ln beta + (1/gamma)*ln alpha (-Inf to Inf)
    LSQest0[1] <- log(LSQest[2]) + (1/LSQest[3])*log(LSQest[1])
    # lambda = 1/sqrt(alpha) (positive so make loglambda)
    LSQest0[2] <- log(1/sqrt(LSQest[1]))
    # sigma = 1/gamma*sqrt(alpha) (positive so make logsigma)
    LSQest0[3] <- log(1/(LSQest[3]*sqrt(LSQest[1])))
    LSQest <- LSQest0

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(log(abs(exp(theta[2]))) - log(exp(theta[3])) - log(gamma(1/(exp(theta[2]))^2)) - (1/(exp(theta[2]))^2)*(((theta[1]*exp(theta[2]))/exp(theta[3])) + log((exp(theta[2]))^2)) + ((1/(exp(theta[2])*exp(theta[3]))) - 1)*log(TTF) - (1/(exp(theta[2]))^2)*exp((exp(theta[2])/exp(theta[3]))*(log(TTF) - theta[1])))
      }
    } else{
      loglik <- function(theta){
        -sum(log(abs(exp(theta[2]))) - log(exp(theta[3])) - log(gamma(1/(exp(theta[2]))^2)) - (1/(exp(theta[2]))^2)*(((theta[1]*exp(theta[2]))/exp(theta[3])) + log((exp(theta[2]))^2)) + ((1/(exp(theta[2])*exp(theta[3]))) - 1)*log(TTF) - (1/(exp(theta[2]))^2)*exp((exp(theta[2])/exp(theta[3]))*(log(TTF) - theta[1]))) - sum(log(1 - Rgamma(1/(exp(theta[2]))^2,exp((exp(theta[2])/exp(theta[3]))*(log(Tc) - theta[1]))/(exp(theta[2]))^2),lower=TRUE))
      }
    }
  }
  if (dist=="Logistic") {
    # positivity_v[2]<-1
    LSQest[1] <- log(LSQest[1])

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(((TTF - lifeF(theta))/exp(theta[1])) - theta[1] - 2*log(1 + exp((TTF - lifeF(theta))/exp(theta[1]))))
      }
    } else{
      loglik <- function(theta){
        -sum(((TTF - lifeF(theta))/exp(theta[1])) - theta[1] - 2*log(1 + exp((TTF - lifeF(theta))/exp(theta[1])))) + sum(log(1 + exp((Tc - lifeC(theta))/exp(theta[1]))))
      }
    }
  }
  if (dist=="Loglogistic") {
    # positivity_v[2]<-1
    LSQest[1] <- log(LSQest[1])

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(((log(TTF) - loglifeF(theta))/exp(theta[1])) - theta[1] - log(TTF*((1 + exp((log(TTF) - loglifeF(theta))/exp(theta[1])))^2)))
      }
    } else{
      loglik <- function(theta){
        -sum(((log(TTF) - loglifeF(theta))/exp(theta[1])) - theta[1] - log(TTF*((1 + exp((log(TTF) - loglifeF(theta))/exp(theta[1])))^2))) + sum(log(1 + exp((log(Tc) - loglifeC(theta))/exp(theta[1]))))
      }
    }
  }
  if (dist=="Gumbel") {
    # positivity_v[2]<-1
    # shift sigma to log sigma
    LSQest[1] <- log(LSQest[1])

    if(is.null(Tc)){
      loglik <- function(theta){
        -sum(-theta[1] + ((TTF - lifeF(theta))/exp(theta[1])) - exp((TTF - lifeF(theta))/exp(theta[1])))
      }
    } else{
      loglik <- function(theta){
        -sum(-theta[1] + ((TTF - lifeF(theta))/exp(theta[1])) - exp((TTF - lifeF(theta))/exp(theta[1]))) + sum(exp((Tc - lifeC(theta))/exp(theta[1])))
      }
    }
  }

  # Group all parameters
  # Check to see if any distribution parameters were tabulated
  if(dist=="Exponential") {
    params_txt<-params_txt
  }
  else {
    params_txt<-c(distparam_txt,params_txt)
  }

  # return(list(loglik,LSQest))

  MLEandvar <- MLE.var.covar.select(loglik,LSQest)
  # return(MLEandvar)
  theta.hat <- MLEandvar[[1]]
  inv.fish  <- MLEandvar[[2]]
  loglik.hat <- -loglik(theta.hat)
  likeli.hat <- exp(loglik.hat)

  if(is.null(Suse) == FALSE){
    # Compute use life and variance with untransformed parameters
    uselife <- lifeUSE(theta.hat)
    uselife_VAR <- lifeUSEVAR(theta.hat,inv.fish)
    uselifelim <- vector(mode = "list", length = 1)
  }
  # Reimplement beta in estimate if Weibull beta was known
  if(is.null(param2) == FALSE && dist == "Weibull"){
    # Reset LSQest to original
    LSQest <- LSQest0
    ishift <- 1
    # Add setbeta to theta.hat
    theta.hat <- c(log(param2),theta.hat)
    # Add zero variance for beta to VARCOV
    i2 <- size(inv.fish)[1]
    inv.fish <- cbind(rep(0,i2),inv.fish)
    inv.fish <- rbind(rep(0,(i2+1)),inv.fish)
  }
  # return(list(uselife,uselife_VAR))


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
      if (ls=="PowerwithBias" && (i == (ishift+2) || i == (ishift+3))){
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
      if (ls=="PowerwithBias" && (i == (ishift+2) || i == (ishift+3))){
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
      if (ls=="PowerwithBias" && (i == (ishift+2) || i == (ishift+3))){
        conflim[[i]] <- sort(exp(conflim[[i]]))
      }
      conflim_txt<-paste(c("One-Sided High ",100*conf.level,"%"),collapse = "")
    }
    fulllimset[[i]]<-c(theta.hat[i],conflim[[i]])

    if((dist=="Weibull"  || dist=="3PWeibull" || dist=="Lognormal" || dist=="Normal" || dist=="2PExponential") && i == 1){
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
    if (ls=="PowerwithBias" && (i == (ishift+2) || i == (ishift+3))){
      fulllimset[[i]] <- c(exp(theta.hat[i]),conflim[[i]])
    }
  }

  AIC = 2*length(theta.hat) + 2*loglik(theta.hat)
  BIC = 2*log(length(TTF)+length(Tc)) + 2*loglik(theta.hat)

  # return(ishift)

  if(is.null(Suse) == FALSE){
    if(sided == "twosided"){
      uselifelim <- uselife + c(-1, 1) * crit * sqrt(uselife_VAR)
      if(min(uselifelim) < 0){
        uselifelim <- uselife*exp(c(-1, 1) * crit * (sqrt(uselife_VAR)/uselife))
      }
    }
    if(sided == "onesidedlow"){
      uselifelim <- uselife - crit2 * sqrt(uselife_VAR)
      if(uselifelim < 0){
        uselifelim <- uselife*exp(-crit2 * (sqrt(uselife_VAR)/uselife))
      }
    }
    if(sided == "onesidedhigh"){
      uselifelim <- uselife + crit2 * sqrt(uselife_VAR)
    }
  }
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
  if (ls=="PowerwithBias"){
    theta.hat[ishift+2] <- exp(theta.hat[ishift+2])
    theta.hat[ishift+3] <- exp(theta.hat[ishift+3])
  }

  # Calculate
  # Set all distribution parameters for each stress level as equal
  if(ls!="Exponential"){
    L <- lifestress.select(ls)[[1]](theta.hat[2:length(theta.hat)],S)
    distparams <- rep(theta.hat[1],length(S))
  } else{
    L <- lifestress.select(ls)[[1]](theta.hat,S)
    distparams <- NULL
  }
  # return(list(S,L,distparams))

  if(ls=="Linear" || ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" || ls=="Eyring2" || ls=="Power" || ls=="PowerwithBias" || ls=="InversePower" || ls=="InversePower2" || ls=="InversePower2" || ls=="Logarithmic"){
    if(is.null(Suse)==TRUE){
      relplotoutput <- lifestress.relationplot.LSQ.2(data,ls,dist,theta.hat,S,L,distparams=distparams,stressunit1 = stressunit1)
    }
    if(is.null(Suse)==FALSE){
      relplotoutput <- lifestress.relationplot.LSQ.2(data,ls,dist,theta.hat,S,L,distparams=distparams,Suse = Suse,stressunit1 = stressunit1)
    }
    if(is.null(Llab)==FALSE && is.null(Slab)==FALSE){
      relplotoutput <- lifestress.relationplot.LSQ.2(data,ls,dist,theta.hat,S,L,distparams=distparams,Suse = Suse,Llab = Llab,Slab = Slab,stressunit1 = stressunit1)
    }
  }

  # Produce some output text that summarizes the results
  cat(c("Maximum-Likelihood estimates for the ",ls_txt,"-",dist_txt," Life-Stress model.\n\n"),sep = "")
  if(is.null(Suse) == TRUE){
    print(matrix(unlist(fulllimset), nrow = length(unlist(fulllimset))/length(LSQest), ncol = length(LSQest), byrow = FALSE,dimnames = list(c("Life-Stress Parameters Mean",conflim_txt),params_txt)))
  }
  if(is.null(Suse) == FALSE){
    # Add column for use life mean and confidence bounds
    fulllimset2<-fulllimset
    fulllimset2[[length(LSQest)+1]] <- c(uselife,uselifelim)
    print(matrix(unlist(fulllimset2), nrow = length(unlist(fulllimset))/length(LSQest), ncol = (length(LSQest)+1), byrow = FALSE,dimnames = list(c("Life-Stress Parameters Mean",conflim_txt),c(params_txt,"Use Life"))))
  }


  if(is.null(Suse) == TRUE){
    return(list(MLE.point.estimate = theta.hat,var.cov.matrix = inv.fish,param.CI = conflim,loglik = loglik.hat,ikelihood = likeli.hat,AIC = AIC,BIC = BIC,plotoutput=plotoutput,relplotoutput=relplotoutput$relationplot))
  }
  if(is.null(Suse) == FALSE){
    return(list(MLE.point.estimate = theta.hat,var.cov.matrix = inv.fish,use.life = uselife,param.CI = conflim,use.life.CI = uselifelim,loglik = loglik.hat,likelihood = likeli.hat,AIC = AIC,BIC = BIC,plotoutput=plotoutput,relplotoutput=relplotoutput$relationplot))
  }
}
