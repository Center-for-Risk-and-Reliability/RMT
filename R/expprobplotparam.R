# Exponential Probability Plot Parameters
# Developed by Dr. Reuel Smith, 2021-2022

probplotparam.exp <- function(xi,F) {
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
  # Calculate least-squares estimate for exponential parameter
  pb2 <- lm(yfit ~ xfit + 0)
  lambda <- summary(pb2)$coefficients[,1]
  SSE <- sum((fitted(pb2) - yfit)^2)
  R2 <- summary(pb2)$r.squared
  expresults <- matrix(c(lambda), nrow = 1, ncol = 1, byrow = TRUE,dimnames = list(c("Exponential Parameter"),c("lambda")))
  # Calculate upper and lower bound of best fit line
  ttfc <- c(-log(1-0.001)/lambda,-log(1-0.999)/lambda)
  return(list(ttfc,fcB,expresults,R2,SSE=SSE))
}
