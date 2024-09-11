# Maximum Likelihood Variance-Covariance Matrix Select
# Developed by Dr. Reuel Smith, 2022-2024

MLE.var.covar.select <- function(loglik,LSQest,bounded=NULL){
  library(pracma)
  library(matrixcalc)
  library(ucminf)
  library(numDeriv)
  library(dfoptim)
  library(DEoptim)

  # There is still some question as to which R script and library to use for MLE estimation.
  # nlm provides the mean estimate at all times, however there are cases that produce a negative
  # variance in some variables.  Further study will need to be done to find the best option.  We may
  # need to have multiple options on standby just in case.  There is the 'hessian' function from
  # the numDeriv library as a first alternative.  A second alternative is the 'ucminf' function
  # from the library of the same name
  # fish <- hessian(loglik,out$estimate)

  # Optimization Option Number as of 7/19/2024
  N <- 3

  theta.hat.check <- vector(mode = "list", length = N)
  fisher.mat.check <- vector(mode = "list", length = N)
  var.covar.mat.check <- vector(mode = "list", length = N)
  real.theta.hat.check <- rep(0,N)
  real.fisher.mat.check <- rep(0,N)
  pos.variance.check <- rep(0,N)

  if(is.null(bounded)==TRUE){
    # =======================================================
    # OPTION 1: Based on nlm (default)
    # =======================================================
    out1 <- nlm(loglik, theta <- LSQest, hessian=TRUE)
    theta.hat1 <- out1$estimate
    fish1 <- out1$hessian

    theta.hat.check[[1]] <- theta.hat1
    fisher.mat.check[[1]] <- fish1
    if(is.null(theta.hat1) == FALSE){
      real.theta.hat.check[1] <- 1
    }
    # return(fish1)
    if(is.null(fish1) == FALSE  && is.nan(sum(fish1)) == FALSE && is.infinite(sum(fish1)) == FALSE){
      real.fisher.mat.check[1] <- 1
      var.covar.mat.check[[1]] <- pinv(fish1)
      if(min(diag(pinv(fish1))) >= 0){
        pos.variance.check[1] <- 1
      }
    }
    # =======================================================
    # OPTION 2: Based on ucminf
    # =======================================================
    out2 <- suppressWarnings(ucminf(LSQest,loglik,hessian=1))
    theta.hat2 <- out2$par
    fish2 <- out2$hessian
    theta.hat.check[[2]] <- theta.hat2
    fisher.mat.check[[2]] <- fish2
    if(is.null(theta.hat2) == FALSE){
      real.theta.hat.check[2] <- 1
    }
    if(is.null(fish2) == FALSE && is.nan(sum(fish2)) == FALSE && is.infinite(sum(fish2)) == FALSE){
      real.fisher.mat.check[2] <- 1
      var.covar.mat.check[[2]] <- pinv(fish2)
      if(min(diag(pinv(fish2))) >= 0){
        pos.variance.check[2] <- 1
      }
    }
    # =======================================================
    # OPTION 3: Based on optim
    # =======================================================
    out3 <- suppressWarnings(optim(LSQest,loglik,hessian=1))
    theta.hat3 <- out3$par
    fish3 <- out3$hessian

    theta.hat.check[[3]] <- theta.hat3
    fisher.mat.check[[3]] <- fish3
    if(is.null(theta.hat3) == FALSE){
      real.theta.hat.check[3] <- 1
    }
    if(is.null(fish3) == FALSE  && is.nan(sum(fish3)) == FALSE && is.infinite(sum(fish3)) == FALSE){
      real.fisher.mat.check[3] <- 1
      var.covar.mat.check[[3]] <- pinv(fish3)
      if(min(diag(pinv(fish3))) >= 0){
        pos.variance.check[3] <- 1
      }
    }
  }

  if(is.null(bounded)==FALSE){
    # =======================================================
    # OPTION 4: Based on nlminb
    # =======================================================
    if(length(LSQest)==2){
      out4 <- suppressWarnings(nlminb(LSQest,loglik,lower = c(-1000,-100), upper = c(bounded,100)))
    }
    if(length(LSQest)==3){
      out4 <- suppressWarnings(nlminb(LSQest,loglik,lower = c(-100,-100,-1000), upper = c(100,100,bounded)))
    }
    theta.hat4 <- out4$par
    fish4 <- hessian(loglik,out4$par)
    # Findings at the edge will result in NaNs for fisher matrices.  Need to set those at zero to prevent error.
    if(length(LSQest)==2 && is.nan(sum(fish4)) == TRUE){
      fish4[1,]<-rep(0,2)
      fish4[,1]<-rep(0,2)
    }
    if(length(LSQest)==3 && is.nan(sum(fish4)) == TRUE){
      fish4[3,]<-rep(0,3)
      fish4[,3]<-rep(0,3)
    }

    theta.hat.check[[1]] <- theta.hat4
    fisher.mat.check[[1]] <- fish4
    if(is.null(theta.hat4) == FALSE){
      real.theta.hat.check[1] <- 1
    }
    if(is.null(fish4) == FALSE && is.infinite(sum(fish4)) == FALSE){
      real.fisher.mat.check[1] <- 1
      var.covar.mat.check[[1]] <- pinv(fish4)
      if(min(diag(pinv(fish4))) >= 0){
        pos.variance.check[1] <- 1
      }
    }
    # =======================================================
    # OPTION 5: Based on DEoptim
    # =======================================================
    if(length(LSQest)==2){
      out5 <- suppressWarnings(DEoptim(loglik,lower = c(-1000,-100), upper = c(bounded,100)))
    }
    if(length(LSQest)==3){
      out5 <- suppressWarnings(DEoptim(loglik,lower = c(-100,-100,-1000), upper = c(100,100,bounded)))
    }

    theta.hat5 <- out5$member$pop[1,]
    fish5 <- hessian(loglik,out5$member$pop[1,])
    # Findings at the edge will result in NaNs for fisher matrices.  Need to set those at zero to prevent error.
    if(length(LSQest)==2 && is.nan(sum(fish5)) == TRUE){
      fish5[1,]<-rep(0,2)
      fish5[,1]<-rep(0,2)
    }
    if(length(LSQest)==3 && is.nan(sum(fish5)) == TRUE){
      fish5[3,]<-rep(0,3)
      fish5[,3]<-rep(0,3)
    }

    theta.hat.check[[2]] <- theta.hat5
    fisher.mat.check[[2]] <- fish5
    if(is.null(theta.hat5) == FALSE){
      real.theta.hat.check[2] <- 1
    }
    if(is.null(fish5) == FALSE && is.infinite(sum(fish5)) == FALSE){
      real.fisher.mat.check[2] <- 1
      var.covar.mat.check[[2]] <- pinv(fish5)
      if(min(diag(pinv(fish5))) >= 0){
        pos.variance.check[2] <- 1
      }
    }
    # =======================================================
    # OPTION 6: Based on nmkb
    # =======================================================
    if(length(LSQest)==2){
      out6 <- suppressWarnings(nmkb(LSQest,loglik,lower = c(-1000,-100), upper = c(bounded,100)))
    }
    if(length(LSQest)==3){
      out6 <- suppressWarnings(nmkb(LSQest,loglik,lower = c(-100,-100,-1000), upper = c(100,100,bounded)))
    }
    theta.hat6 <- out6$par
    fish6 <- hessian(loglik,out6$par)
    # Findings at the edge will result in NaNs for fisher matrices.  Need to set those at zero to prevent error.
    if(length(LSQest)==2 && is.nan(sum(fish6)) == TRUE){
      fish6[1,]<-rep(0,2)
      fish6[,1]<-rep(0,2)
    }
    if(length(LSQest)==3 && is.nan(sum(fish6)) == TRUE){
      fish6[3,]<-rep(0,3)
      fish6[,3]<-rep(0,3)
    }

    theta.hat.check[[3]] <- theta.hat6
    fisher.mat.check[[3]] <- fish6
    if(is.null(theta.hat6) == FALSE){
      real.theta.hat.check[3] <- 1
    }
    if(is.null(fish6) == FALSE && is.infinite(sum(fish6)) == FALSE){
      real.fisher.mat.check[3] <- 1
      var.covar.mat.check[[3]] <- pinv(fish6)
      if(min(diag(pinv(fish6))) >= 0){
        pos.variance.check[3] <- 1
      }
    }
  }

  # return(pos.variance.check)
  # Return response for the solution based on checks
  for(i in 1:N){
    if((real.theta.hat.check[i] + pos.variance.check[i] + real.fisher.mat.check[i]) == 3){
      theta.hat <- theta.hat.check[[i]]
      inv.fish <- var.covar.mat.check[[i]]
      return(list(theta.hat,inv.fish,fulllog = list(thetaest = theta.hat.check,fishermatest=fisher.mat.check,varcovmat = var.covar.mat.check,check_realtheta = real.theta.hat.check,check_realfishermat = real.fisher.mat.check,check_positivevariance = pos.variance.check)))
    }
  }
  for(i in 1:N){
    if((real.theta.hat.check[i] + pos.variance.check[i]) == 2){
      theta.hat <- theta.hat.check[[i]]
      inv.fish <- var.covar.mat.check[[i]]
      return(list(theta.hat,inv.fish,fulllog = list(thetaest = theta.hat.check,fishermatest=fisher.mat.check,varcovmat = var.covar.mat.check,check_realtheta = real.theta.hat.check,check_realfishermat = real.fisher.mat.check,check_positivevariance = pos.variance.check)))
    }
  }
  for(i in 1:N){
    if(real.theta.hat.check[i] == 1 && pos.variance.check[i] == 0){
      theta.hat <- theta.hat.check[[i]]
      inv.fish <- matrix(rep(0,length(LSQest)^2),length(LSQest),length(LSQest))
      return(list(theta.hat,inv.fish,fulllog = list(thetaest = theta.hat.check,fishermatest=fisher.mat.check,varcovmat = var.covar.mat.check,check_realtheta = real.theta.hat.check,check_realfishermat = real.fisher.mat.check,check_positivevariance = pos.variance.check)))
    }
  }
  # If no results work, then return a NULL estimate for theta.hat and fisher matrix
  theta.hat <- NULL
  inv.fish <- NULL

  return(list(theta.hat,inv.fish,fulllog = list(thetaest = theta.hat.check,fishermatest=fisher.mat.check,varcovmat = var.covar.mat.check,check_realtheta = real.theta.hat.check,check_realfishermat = real.fisher.mat.check,check_positivevariance = pos.variance.check)))
}
