# Gumbel Probability Plot Parameters
# Developed by Dr. Reuel Smith, 2023

probplotparam.gumb <- function(xi,R) {
  # Upper and lower bounds of the Percent Failure axis in percent
  fcB <- log(-log(1-c(.001,0.999)))
  # Set up the x and y data for fitting
  if(R[length(R)]==0){
    yfit <- log(-log(R[1:length(R)-1]))
    xfit <- xi[1:length(R)-1]
  } else {
    yfit <- log(-log(R))
    xfit <- xi
  }
  # Calculate least-squares estimate for Weibull parameters
  if(length(xi) >= 2){
    gumbparm  <- lm(yfit ~ poly(xfit, 1, raw=TRUE))
    sigma <- 1/summary(gumbparm)$coefficients[2,1]
    mu <- -summary(gumbparm)$coefficients[1,1]*sigma
    R2 <- summary(gumbparm)$r.squared
    # Calculate upper and lower bound of best fit line
    intercept <- summary(gumbparm)$coefficients[1,1]
    ttfc <- c((log(-log(1-0.001))-intercept)*sigma,(log(-log(1-0.999))-intercept)*sigma)
  }
  if(length(xi) == 1){
    mu <- xi
    sigma <- NA
    R2 <- 0
    ttfc <- rep(xi,2)
  }
  gumbresults <- matrix(c(mu,sigma), nrow = 1, ncol = 2, byrow = TRUE,dimnames = list(c("Gumbel Parameters"),c("mu", "sigma")))

  return(list(ttfc,fcB,gumbresults,R2))
}
