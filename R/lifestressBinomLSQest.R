# Least-Squares Binomial Accelerated Life-Stress Estimator
# Developed by Mohammad Modarres and Reuel Smith, 2022

lifestress.Binom.LSQest <- function(data,interact_stress,weight0){
  #Load pracma library for pseudoinverse
  library(pracma)
  # Load port library to cap LSQ estimation bounds such that 0 < p_0 < 1 and
  # theta > 0
  # library(port)

  # Sets the default for weight use to zero if input is not used
  if(missing(weight0)){
    weight0<-0
  }

  # Check and see that there are multiple stress levels
  if(dim(data)[2]<3 && weight0==0) {
    stop('Need more than one stress level to generate estimates')
  }
  if(dim(data)[2]<4 && weight0==1) {
    stop('Need more than one stress level to generate estimates')
  }

  # Check to see if there are more than one set of data
  if(dim(data)[1]==1) {
    stop('Need more than one data entry to generate estimates')
  }

  # Isolate binomial data n and k
  ni <- data[,1]
  ki <- data[,2]

  # Isolate the stress differentials
  Dstressi <- data[,3:(dim(data)[2]-weight0)]
  names(Dstressi) <- NULL

  # Sets the default computation for stresses to simply single evaluation
  # and not check for interactions
  if(missing(interact_stress)){
    interact_stress<-0
  }


  # Writeup for the output text
  ls_txt<-"Binomial Acceleration"
  pdf_txt<-"(n|k) \U2219 [1 - (1 - p_0) \U2219 exp(-\U03A3_(j=1)^m \U03B8_j \U0394S_j)]^k \U2219 [(1 - p_0) \U2219 exp(-\U03A3_(j=1)^m \U03B8_j \U0394S_j)]^(n - k)"
  distparam_txt<-"p_0"

  # Set y in non-linear regression analysis
  y<-log(ki/ni)-(lfactorial(ni) - lfactorial(ki) - lfactorial(ni-ki))

  # Set up Matrix and Vector to perform pseudo-inverse operation to get LSQ parameters
  # theta_0, ... theta_n.
  V_main <- log(1 - ki/ni)
  # V_main <- log(ki/ni)-(lfactorial(ni) - lfactorial(ki) - lfactorial(ni-ki))

  Dstress_txt<-paste("\U03B8", 1:(dim(data)[2]-2-weight0), sep = "_")

  if(is.null(dim(Dstressi))==TRUE){
    # One Stress Differential case
    # Dstressi_m <- matrix(c(rep(1,dim(data)[1]),Dstressi), nrow = dim(data)[1], ncol = 2, byrow = FALSE)
    Dstressi_m<-Dstressi
    LSQout <- nls(y ~ ki*log(1 - (1 - p0)*exp(-c(theta%*%Dstressi))) + (ni-ki)*log(1 - p0) - (ni-ki)*c(theta%*%Dstressi), start = list(p0 = 0.5, theta = 0.1), algorithm = "port", lower = list(p0 = 0, theta = 0), upper = list(p0 = 1, theta = 100),  control=list(maxit=100000,tol=1e-10,warnOnly=T,minFactor=1e-10))
     } else {
    # Multi Stress Differential
    for(i in 1:(dim(data)[2]-2-weight0)){
      if(i==1){
        Dstressi_v<-Dstressi[,i]
      } else {
        Dstressi_v<-c(Dstressi_v,Dstressi[,i])
      }
    }
    if(interact_stress==0){
      # No Stress dependency check
      Dstressi_m <- matrix(Dstressi_v, nrow = dim(data)[1], ncol = length(Dstressi_v)/dim(data)[1], byrow = FALSE)
      LSQout <- nls(y ~ ki*log(1 - (1 - p0)*exp(-c(Dstressi_m%*%c(theta1,theta2)))) + (ni-ki)*log(1 - p0) - (ni-ki)*c(Dstressi_m%*%c(theta1,theta2)), start = list(p0 = 0.5, theta1 = 0.1, theta2 = 0.1), algorithm = "port", lower = list(p0 = 0, theta1 = 0, theta2 = 0), upper = list(p0 = 1, theta1 = 100, theta2 = 100),  control=list(maxit=100000,tol=1e-10,warnOnly=T,minFactor=1e-10))
    }
    if(interact_stress==1){
      # Stress dependency check
      # Permutations of stress pairs
      Dstresspairs <- combn(dim(Dstressi)[2],2)
      for(i in 1:dim(Dstresspairs)[2]){
        if(i==1){
          Dstresspair_v <- Dstressi[, Dstresspairs[1,i]]*Dstressi[, Dstresspairs[2,i]]
          Dstressjoint_txt<-paste("\U03B8", paste(Dstresspairs[,i],collapse=""), sep = "_")
        } else{
          Dstresspair_v <- c(Dstresspair_v,Dstressi[, Dstresspairs[1,i]]*Dstressi[, Dstresspairs[2,i]])
          Dstressjoint_txt<-c(Dstressjoint_txt,paste("\U03B8", paste(Dstresspairs[,i],collapse=""), sep = "_"))
        }
      }
      Dstressi_m <- matrix(c(Dstressi_v,Dstresspair_v), nrow = dim(data)[1], ncol = (length(Dstressi_v)/dim(data)[1])+dim(Dstresspairs)[2], byrow = FALSE)
      LSQout <- nls(y ~ ki*log(1 - (1 - p0)*exp(-c(Dstressi_m%*%c(theta1,theta2,theta12)))) + (ni-ki)*log(1 - p0) - (ni-ki)*c(Dstressi_m%*%c(theta1,theta2,theta12)), start = list(p0 = 0.5, theta1 = 0.1, theta2 = 0.1, theta12 = 0.1), algorithm = "port", lower = list(p0 = 0, theta1 = 0, theta2 = 0, theta12 = 0), upper = list(p0 = 1, theta1 = 100, theta2 = 100, theta12 = 100),  control=list(maxit=100000,tol=1e-10,warnOnly=T,minFactor=1e-10))
      Dstress_txt<-c(Dstress_txt,Dstressjoint_txt)
    }
    # fx <- ki*log(1 - (1 - p0)*exp(-c(Dstressi_m%*%theta))) + (ni-ki)*log(1 - p0) - (ni-ki)*c(Dstressi_m%*%theta)
    # x0 <- list(p0 = 0.5, theta = rep(1,dim(Dstressi_m)[2]))
    # LSQout <- nls(y ~ ki*log(1 - (1 - p0)*exp(-c(Dstressi_m%*%theta))) + (ni-ki)*log(1 - p0) - (ni-ki)*c(Dstressi_m%*%theta), start = list(p0 = 0.5, theta = rep(0.1,dim(Dstressi_m)[2])), algorithm = "port", lower = list(p0 = 0, theta = rep(0,dim(Dstressi_m)[2])), upper = list(p0 = 1, theta = rep(1,dim(Dstressi_m)[2])),  control=list(maxit=100000,tol=1e-10,warnOnly=T,minFactor=1e-10))
    # LSQout <- nls(y ~ ki*log(1 - (1 - p0)*exp(-c(Dstressi_m%*%exp(theta)))) + (ni-ki)*log(1 - p0) - (ni-ki)*c(Dstressi_m%*%exp(theta)), start = list(p0 = 0.5, theta = log(rep(1,dim(Dstressi_m)[2]))))
  }


  # Compute the theta LSQ parameters
  # LSQout <- nls(y ~ fx, start = x0)
  LSQ  <- LSQout$m$getAllPars()
  names(LSQ) <- NULL

  params_txt<-c(distparam_txt,Dstress_txt)

  # Check to see if p_0 is a valid estimate (between 0 and 1)
  if(LSQ[1]<0 || LSQ[1]>1) {
    stop('Estimate for p_0 is out of bounds.  Check data for outlier input or correctness.')
  }

  # Produce some output text that summariZes the results
  cat(c("Least-Squares estimates for the ",ls_txt," Life-Stress model.\n\nPr(k|p_0,\U03B8_j) = ",pdf_txt,"\n\n"),sep = "")
  print(matrix(c(LSQ), nrow = 1, ncol = length(LSQ), byrow = TRUE,dimnames = list(c("Life-Stress Parameters"),params_txt)))
  cat("\n\n")

  # Return parameter list
  return(list(LSQ,Dstressi_m,params_txt))
  # return(list(S,L,distparams,V_main,M_main))
}
