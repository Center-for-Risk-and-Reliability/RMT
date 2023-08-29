# MLE Variance-Covariance Selector by Life-Stress Model
# Developed by Dr. Reuel Smith, 2021-2022

MLEvarcov.select.by.lifestress <- function(ls,dist) {
  if (ls=="Linear") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsmodel <- function(lsparams,S) {
      lsparams[2] + S*lsparams[1]
    }
    lin_lsmodel <- function(lsparams,S) {
      log(lsparams[2] + S*lsparams[1])
    }
  }
  if (ls=="Exponential") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsmodel <- function(lsparams,S) {
      lsparams[2]*exp(S*lsparams[1])
    }
    lin_lsmodel <- function(lsparams,S) {
      log(lsparams[2]) + S*lsparams[1]
    }
  }
  if (ls=="Arrhenius") {
    # lsparams[1] - parameter Ea, lsparams[2] - parameter b
    # Temperature HAS to be in Kelvin for this to work
    K<-8.617385e-5
    lsmodel <- function(lsparams,S) {
      lsparams[2]*exp(lsparams[1]/(K*S))
    }
    lin_lsmodel <- function(lsparams,S) {
      log(lsparams[2]) + lsparams[1]/(K*S)
    }
  }
  if (ls=="Eyring") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsmodel <- function(lsparams,S) {
      (lsparams[2]/S)*exp(lsparams[1]/S)
    }
    lin_lsmodel <- function(lsparams,S) {
      log(lsparams[2]) - log(S) + lsparams[1]/S
    }
  }
  if (ls=="Eyring2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsmodel <- function(lsparams,S) {
      (1/S)*exp(-(lsparams[1] - (lsparams[2]/S)))
    }
    lin_lsmodel <- function(lsparams,S) {
      -log(S) - log(lsparams[1]) + lsparams[2]/S
    }
  }
  if (ls=="Power") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsmodel <- function(lsparams,S) {
      lsparams[2]*(S^lsparams[1])
    }
    lin_lsmodel <- function(lsparams,S) {
      log(lsparams[2]) + lsparams[1]*log(S)
    }
  }
  if (ls=="InversePower") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsmodel <- function(lsparams,S) {
      lsparams[2]*(S^-lsparams[1])
    }
    lin_lsmodel <- function(lsparams,S) {
      log(lsparams[2]) - lsparams[1]*log(S)
    }
  }
  if (ls=="Logarithmic") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsmodel <- function(lsparams,S) {
      lsparams[1]*log(S) + lsparams[2]
    }
    lin_lsmodel <- function(lsparams,S) {
      log(lsparams[1]*log(S) + lsparams[2])
    }
  }
  if (ls=="MultiStress") {
    # lsparams - c(ao,a1,...,an), S - c(1,S1,...,Sn)
    lsmodel <- function(lsparams,S) {
      exp(lsparams%*%c(1,S))
    }
    lin_lsmodel <- function(lsparams,S) {
      lsparams%*%c(1,S)
    }
  }
  if (ls=="TempHumidity") {
    # lsparams[1] - parameter A, lsparams[2] - parameter a, lsparams[3] - parameter b
    lsmodel <- function(lsparams,S,H) {
      lsparams[1]*exp((lsparams[2]/S) + (lsparams[3]/H))
    }
    lin_lsmodel <- function(lsparams,S,H) {
      log(lsparams[1]) + (lsparams[2]/S) + (lsparams[3]/H)
    }
  }
  if (ls=="TempNonthermal") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
    lsmodel <- function(lsparams,S,U) {
      lsparams[3]/((U^lsparams[2])*exp(-lsparams[1]/S))
    }
    lin_lsmodel <- function(lsparams,S,U) {
      log(lsparams[3]) - lsparams[2]*log(U) + (lsparams[1]/S)
    }
  }
  if (ls=="Eyring3") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # lsparams[3] - parameter c, lsparams[4] - parameter d
    lsmodel <- function(lsparams,S,U) {
      (1/S)*exp((lsparams[1] + (lsparams[2]/S)) + (lsparams[3] + (lsparams[4]/S))*U)
    }
    lin_lsmodel <- function(lsparams,S,U) {
      -log(S) + lsparams[1] + (lsparams[2]/S) + (lsparams[3] + (lsparams[4]/S))*U
    }
  }

  return(list(lsmodel,lin_lsmodel))
}
