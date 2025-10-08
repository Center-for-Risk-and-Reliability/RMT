# Three Parameter Weibull Probability Plot Parameters
# Developed by Dr. Reuel Smith, 2021-2022

probplotparam.wbl3P <- function(xi,R,CDFrangesetting = 1) {
  library(nls.multstart)
  # Upper and lower bounds of the Percent Failure axis in percent
  if(CDFrangesetting == 1){ # Minitab range 1% to 99%
    fcB <- log(log(1/(1-c(.01,0.99))))
  }
  if(CDFrangesetting == 2){ # Weibull++ range 0.1% to 99.9%
    fcB <- log(log(1/(1-c(.001,0.999))))
  }
  if(R[length(R)]==0){
    yfit <- log(log(1/R[1:length(R)-1]))
    xfit <- xi[1:length(R)-1]
  } else {
    yfit <- log(log(1/R))
    xfit <- xi
  }

  # Gamma estimates
  if(length(xi) >= 2){
    # Define linearized equation
    # wbl3Plineq <- function(alp0,bet0,gam0){
    #   F_x <- bet0*log(xfit - gam0) - bet0*log(alp0)
    #   return(F_x)
    # }
    # RCS06102024 - Estimate by establishing Sum of square error (SSE)
    wbl3PSSE <- function(theta){
      SSE <- sum((yfit - theta[2]*log(xfit - theta[3]) + theta[2]*theta[1])^2)
      return(SSE)
    }
    # Define dataset and run analysis

    # alp0bet0 <- c(probplotparam.wbl(xi,R)[[3]])
    # df <- data.frame(y = yfit, xfit = xfit)
    # LSQout <- nls_multstart(y ~ wbl3Plineq(alp0,bet0,gam0), data = df, iter = 5000, start_lower = c(alp0 = 0.1, bet0 = 0.1, gam0 = -0.99*min(xi)),start_upper = c(alp0 = 10*alp0bet0[1], bet0 = 10*alp0bet0[2], gam0 = 0.99*min(xi)), lower = c(alp0 = 0.01, bet0 = 0.01, gam0 = -1000))
    # # return(LSQout)

    # RCS06102024 - Define upper bound of gamma as 0.99*x_min and obtain estimate of alpha and beta
    gamest0 <- 0.5*min(xi)
    alpbetest0 <- c(probplotparam.wbl(xi-gamest0,R)[[3]])
    # LSQ <- nlm(wbl3PSSE, c(log(alpbetest0[1]),alpbetest0[2],gamest0), hessian=TRUE)$estimate
    LSQ <- nlminb(c(log(alpbetest0[1]),alpbetest0[2],gamest0),wbl3PSSE,hessian=TRUE,lower = c(-Inf,0,-Inf),upper = c(Inf,Inf,0.99*min(xi)))$par
    # Output parameter estimates
    alpha <- exp(LSQ[1])
    beta <- LSQ[2]
    gammaparam <- LSQ[3]
    SSE <- wbl3PSSE(LSQ)

    # LSQ <- LSQout$m$getPars()
    # names(LSQ) <- NULL
    # alpha <- LSQ[1]
    # beta <- LSQ[2]
    # gammaparam <- LSQ[3]
    # Calculate upper and lower bound of best fit line
    ttfc <- gammaparam + exp(log(alpha) + (fcB/beta))
    fcB <- fcB
    # Coefficient of Determination
    R2 <- probplotparam.wbl(xi-gammaparam,R)[[4]]
  }
  if(length(xi) == 1){
    alpha <- xi
    beta <- NA
    gammaparam <- 0
    R2 <- 0
    ttfc <- rep(xi,2)
  }

  wbl3Presults <- matrix(c(alpha,beta,gammaparam), nrow = 1, ncol = 3, byrow = TRUE,dimnames = list(c("3P Wbl Parameters"),c("alpha", "beta", "Gamma")))
  return(list(ttfc,fcB,wbl3Presults,R2,SSE=SSE))
}
