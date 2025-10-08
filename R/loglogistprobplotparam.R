# Log-logistic Probability Plot Parameters
# Developed by Dr. Reuel Smith, 2023

probplotparam.loglogist <- function(xi,F,CDFrangesetting = 1) {
  # Upper and lower bounds of the Percent Failure axis in percent
  if(CDFrangesetting == 1){ # Minitab range 1% to 99%
    fcB <- log((1/c(.01,0.99))-1)
  }
  if(CDFrangesetting == 2){ # Weibull++ range 0.1% to 99.9%
    fcB <- log((1/c(.001,0.999))-1)
  }
  # Set up the x and y data for fitting
  if(F[length(F)]==0){
    yfit <- log((1/F[1:length(F)-1])-1)
    xfit <- log(xi[1:length(F)-1])
  } else {
    yfit <- log((1/F)-1)
    xfit <- log(xi)
  }
  # Calculate least-squares estimate for Weibull parameters
  if(length(xi) >= 2){
    loglogistparm  <- lm(yfit ~ poly(xfit, 1, raw=TRUE))
    sigma <- 1/summary(loglogistparm)$coefficients[2,1]
    mu <- -summary(loglogistparm)$coefficients[1,1]*sigma
    SSE <- sum((fitted(loglogistparm) - yfit)^2)
    R2 <- summary(loglogistparm)$r.squared
    # Calculate upper and lower bound of best fit line
    intercept <- summary(loglogistparm)$coefficients[1,1]
    if(CDFrangesetting == 1){ # Minitab range 1% to 99%
      ttfc <- c(exp((log((1/0.01)-1)-intercept)*sigma),exp((log((1/0.99)-1)-intercept)*sigma))
    }
    if(CDFrangesetting == 2){ # Weibull++ range 0.1% to 99.9%
      ttfc <- c(exp((log((1/0.001)-1)-intercept)*sigma),exp((log((1/0.999)-1)-intercept)*sigma))
    }
  }
  if(length(xi) == 1){
    mu <- xi
    sigma <- NA
    SSE <- 0
    R2 <- 0
    ttfc <- rep(xi,2)
  }
  loglogistresults <- matrix(c(mu,-sigma), nrow = 1, ncol = 2, byrow = TRUE,dimnames = list(c("Log-logistic Parameters"),c("mu", "sigma")))

  return(list(ttfc,fcB,loglogistresults,R2,SSE=SSE))
}
