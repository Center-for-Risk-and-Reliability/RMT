# Weibull Probability Plot Parameters
# Developed by Dr. Reuel Smith, 2021-2025

probplotparam.wbl <- function(xi,R,CDFrangesetting = 1) {
  # Upper and lower bounds of the Percent Failure axis in percent
  if(CDFrangesetting == 1){ # Minitab range 1% to 99%
    fcB <- log(log(1/(1-c(.01,0.99))))
  }
  if(CDFrangesetting == 2){ # Weibull++ range 0.1% to 99.9%
    fcB <- log(log(1/(1-c(.001,0.999))))
  }
  # Set up the x and y data for fitting
  if(R[length(R)]==0){
    yfit <- log(log(1/R[1:length(R)-1]))
    xfit <- log(xi[1:length(R)-1])
  } else {
    yfit <- log(log(1/R))
    xfit <- log(xi)
  }
  # Calculate least-squares estimate for Weibull parameters
  if(length(xi) >= 2){
    wblparm  <- lm(yfit ~ poly(xfit, 1, raw=TRUE))
    beta <- summary(wblparm)$coefficients[2,1]
    alpha <- exp(-summary(wblparm)$coefficients[1,1]/beta)
    SSE <- sum((fitted(wblparm) - yfit)^2)
    R2 <- summary(wblparm)$r.squared
    # Calculate upper and lower bound of best fit line
    intercept <- summary(wblparm)$coefficients[1,1]
    if(CDFrangesetting == 1){ # Minitab range 1% to 99%
      ttfc <- c(exp((log(log(1./(1-0.01)))-intercept)/beta),exp((log(log(1./(1-0.99)))-intercept)/beta))
    }
    if(CDFrangesetting == 2){ # Weibull++ range 0.1% to 99.9%
      ttfc <- c(exp((log(log(1./(1-0.001)))-intercept)/beta),exp((log(log(1./(1-0.999)))-intercept)/beta))
    }

  }
  if(length(xi) == 1){
    alpha <- xi
    beta <- NA
    SSE <- 0
    R2 <- 0
    ttfc <- rep(xi,2)
  }
  wblresults <- matrix(c(alpha,beta), nrow = 1, ncol = 2, byrow = TRUE,dimnames = list(c("Wbl Parameters"),c("alpha", "beta")))

  return(list(ttfc,fcB,wblresults,R2,SSE=SSE))
}
