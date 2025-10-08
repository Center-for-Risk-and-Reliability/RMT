# Gamma Probability Plot Parameters
# Developed by Dr. Reuel Smith, 2022-2025

probplotparam.gam <- function(xi,F,CDFrangesetting = 1) {
  library(pracma)
  library(nls.multstart)
  library(zipfR)

  # Set up the x and y data for fitting
  if(F[length(F)]==1){
    yfit <- (F[1:length(F)-1])
    xfit <- log(xi[1:length(F)-1])
  } else {
    yfit <- (F)
    xfit <- log(xi)
  }

  # Calculate least-squares estimate for gamma parameters
  if(length(xi) >= 2){
    # RCS06252024 - Estimate by establishing Sum of square error (SSE)
    # theta[1] - mu
    # theta[2] - loglambda
    gammaSSE <- function(theta){
      SSE <- sum((yfit - Rgamma(1/(exp(theta[2])^2),exp(xfit - theta[1])/(exp(theta[2])^2), lower = TRUE)))^2
      return(SSE)
    }
    # RCS06252024 - Get initial estimate based on lognormal fit (estimate when lambda = 0 or close to 0)
    x0 <- c(probplotparam.logn(xi,F)[[3]])
    x0[2] <- log(x0[2])

    # RCS06262024 - Refine estimate of mu and lambda
    # https://help.reliasoft.com/reference/life_data_analysis/lda/the_generalized_gamma_distribution.html
    LSQ <- nlminb(x0,gammaSSE,hessian=TRUE,lower = c(-Inf,-Inf),upper = c(Inf,Inf))$par
    # Output parameter estimates
    alpest <- 1/exp(LSQ[2])^2
    betest <- exp(LSQ[1] - log(alpest))
    SSE <- gammaSSE(LSQ)
    SST <- sum((yfit - mean(yfit))^2)
    R2 <- 1 - (SSE/SST)

    # Calculate upper and lower bound of best fit line
    ttfc <- exp(LSQ[1] + log((1/alpest)*Rgamma.inv(alpest,c(0.001,0.999),lower = TRUE)))

    # ttfc <- betest*Rgamma.inv(alpest,c(0.001,0.999),lower = TRUE)
  }
  if(length(xi) == 1){
    alpest <- 1
    betest <- xi
    R2 <- 0
    ttfc <- rep(xi,2)
  }
  # Upper and lower bounds of the Percent Failure axis in percent
  fcB <- Rgamma.inv(alpest,c(0.001,0.999),lower = TRUE)

  gamresults <- matrix(c(alpest,betest), nrow = 1, ncol = 2, byrow = TRUE,dimnames = list(c("Gamma Parameters"),c("alpha", "beta")))
  return(list(ttfc,fcB,gamresults,R2,SSE=SSE))
}
