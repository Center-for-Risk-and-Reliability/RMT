# Generalized Three-Parameter Gamma Probability Plot Parameters
# Developed by Dr. Reuel Smith, 2024

probplotparam.gam3P <- function(xi,F,CDFrangesetting = 1) {
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
    # theta[2] - logsigma
    # theta[3] - loglambda
    gammaSSE <- function(theta){
      SSE <- sum((yfit - Rgamma(1/(exp(theta[3])^2),exp((exp(theta[3])/exp(theta[2]))*(xfit - theta[1]))/(exp(theta[3])^2), lower = TRUE)))^2
      return(SSE)
    }
    # RCS06252024 - Get initial estimate based on lognormal fit (estimate when lambda = 0 or close to 0)
    x0 <- c(probplotparam.logn(xi,F)[[3]])
    x0[2] <- log(x0[2])
    x0 <- c(x0,0)

    # RCS06262024 - Refine estimate of mu and lambda
    # https://help.reliasoft.com/reference/life_data_analysis/lda/the_generalized_gamma_distribution.html
    LSQ <- nlminb(x0,gammaSSE,hessian=TRUE,lower = c(-Inf,-Inf,-Inf),upper = c(Inf,Inf,Inf))$par
    # Output parameter estimates
    alpest <- 1/exp(LSQ[3])^2
    gamest <- exp(LSQ[3])/exp(LSQ[2])
    betest <- exp(LSQ[1] - (1/gamest)*log(alpest))
    SSE <- gammaSSE(LSQ)
    SST <- sum((yfit - mean(yfit))^2)
    R2 <- 1 - (SSE/SST)

    # Calculate upper and lower bound of best fit line
    ttfc <- exp(LSQ[1] + (1/gamest)*log((1/alpest)*Rgamma.inv(alpest,c(0.001,0.999),lower = TRUE)))
    if(CDFrangesetting == 1){ # Minitab range 1% to 99%
      ttfc <- exp(LSQ[1] + (1/gamest)*log((1/alpest)*Rgamma.inv(alpest,c(0.01,0.99),lower = TRUE)))
    }
    if(CDFrangesetting == 2){ # Weibull++ range 0.1% to 99.9%
      ttfc <- exp(LSQ[1] + (1/gamest)*log((1/alpest)*Rgamma.inv(alpest,c(0.001,0.999),lower = TRUE)))
    }
  }
  if(length(xi) == 1){
    alpest <- 1
    betest <- xi
    gamest <- 0
    R2 <- 0
    ttfc <- rep(xi,2)
  }
  # Upper and lower bounds of the Percent Failure axis in percent
  if(CDFrangesetting == 1){ # Minitab range 1% to 99%
    fcB <- Rgamma.inv(alpest,c(0.01,0.99),lower = TRUE)
  }
  if(CDFrangesetting == 2){ # Weibull++ range 0.1% to 99.9%
    fcB <- Rgamma.inv(alpest,c(0.001,0.999),lower = TRUE)
  }

  gamresults <- matrix(c(alpest,betest,gamest), nrow = 1, ncol = 3, byrow = TRUE,dimnames = list(c("Generalized Gamma Parameters"),c("alpha", "beta","gamma")))
  return(list(ttfc,fcB,gamresults,R2,SSE=SSE))
}
