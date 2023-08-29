# Logistic Probability Plot Parameters
# Developed by Dr. Reuel Smith, 2023

probplotparam.logist <- function(xi,F) {
  # Upper and lower bounds of the Percent Failure axis in percent
  fcB <- log((1/c(.001,0.999))-1)
  # Set up the x and y data for fitting
  if(F[length(F)]==0){
    yfit <- log((1/F[1:length(F)-1])-1)
    xfit <- xi[1:length(F)-1]
  } else {
    yfit <- log((1/F)-1)
    xfit <- xi
  }
  # Calculate least-squares estimate for Weibull parameters
  if(length(xi) >= 2){
    logistparm  <- lm(yfit ~ poly(xfit, 1, raw=TRUE))
    sigma <- 1/summary(logistparm)$coefficients[2,1]
    mu <- -summary(logistparm)$coefficients[1,1]*sigma
    R2 <- summary(logistparm)$r.squared
    # Calculate upper and lower bound of best fit line
    intercept <- summary(logistparm)$coefficients[1,1]
    ttfc <- c((log((1/0.001)-1)-intercept)*sigma,(log((1/0.999)-1)-intercept)*sigma)
  }
  if(length(xi) == 1){
    mu <- xi
    sigma <- NA
    R2 <- 0
    ttfc <- rep(xi,2)
  }
  logistresults <- matrix(c(mu,-sigma), nrow = 1, ncol = 2, byrow = TRUE,dimnames = list(c("Logistic Parameters"),c("mu", "sigma")))

  return(list(ttfc,fcB,logistresults,R2))
}
