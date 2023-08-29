# Two Parameter Exponential Probability Plot Parameters
# Developed by Dr. Reuel Smith, 2021-2022

probplotparam.exp2P <- function(xi,F) {
  # Upper and lower bounds of the Percent Failure axis in percent
  fcB <- -log(1-c(.001,0.999))
  # Set up the x and y data for fitting
  if(F[length(F)]==1){
    yfit <- -log(1-F[1:length(F)-1])
    xfit <- xi[1:length(F)-1]
  } else {
    yfit <- -log(1-F)
    xfit <- xi
  }
  # Calculate least-squares estimate for 2P exponential parameter
  if(length(xi) >= 2){
    pb2 <- lm(yfit ~ poly(xfit, 1, raw=TRUE))
    sigma <- 1/summary(pb2)$coefficients[2,1]
    theta <- -sigma*summary(pb2)$coefficients[1,1]
    R2 <- summary(pb2)$r.squared
    # Calculate upper and lower bound of best fit line
    ttfc <- c((theta-sigma*log(1-0.001)),(theta-sigma*log(1-0.999)))
  }
  if(length(xi) == 1){
    sigma <- 1/xi
    theta <- 0
    R2 <- 0
    ttfc <- rep(xi,2)
  }

  exp2Presults <- matrix(c(theta,sigma), nrow = 1, ncol = 2, byrow = TRUE,dimnames = list(c("2P Exponential Parameter"),c("theta","sigma")))
  return(list(ttfc,fcB,exp2Presults,R2))
}
