# Gamma Probability Plot Parameters
# Developed by Dr. Reuel Smith, 2022-2023

probplotparam.gam <- function(xi,F) {
  library(pracma)
  library(nls.multstart)
  library(zipfR)

  # Set up the x and y data for fitting
  if(F[length(F)]==1){
    yprefit <- F[1:length(F)-1]
    xfit <- xi[1:length(F)-1]
  } else {
    yprefit <- F
    xfit <- xi
  }
  gammaeqn1 <- function(alp,bet){
    F_x <- bet*Rgamma.inv(alp,yprefit)
    return(F_x)
  }

  # Calculate least-squares estimate for gamma parameters
  if(length(xi) >= 2){
    df2 <- data.frame(xi = xfit, Fi = yprefit)
    LSQout <- nls_multstart(xi ~ gammaeqn1(alp,bet), data = df2, iter = 5000, start_lower = c(alp = 0.1,bet = 0.1),start_upper = c(alp = 10,bet = 10), lower = c(alp = 0,bet = 0))
    alpest <- LSQout$m$getPars()[[1]]
    betest <- LSQout$m$getPars()[[2]]

    yfit <- Igamma.inv(alpest,yprefit)

    pb2 <- lm(yfit ~ poly(xfit, 1, raw=TRUE))
    intercept <- summary(pb2)$coefficients[1,1]
    slope <- summary(pb2)$coefficients[2,1]
    R2 <- summary(pb2)$r.squared

    gamresults <- matrix(c(alpest,betest), nrow = 1, ncol = 2, byrow = TRUE,dimnames = list(c("Gamma Parameters"),c("alpha", "beta")))

    # Calculate upper and lower bound of best fit line
    ttfc <- c(betest*Rgamma.inv(alpest,0.001),betest*Rgamma.inv(alpest,1-0.999,lower = FALSE))
  }
  if(length(xi) == 1){
    alpest <- 1
    betest <- xi
    R2 <- 0
    gamresults <- matrix(c(alpest,betest), nrow = 1, ncol = 2, byrow = TRUE,dimnames = list(c("Gamma Parameters"),c("alpha", "beta")))
    ttfc <- rep(xi,2)
  }
  # Upper and lower bounds of the Percent Failure axis in percent
  fcB <- c(Rgamma.inv(alpest,0.001),Rgamma.inv(alpest,1-0.999,lower = FALSE))

  return(list(ttfc,fcB,gamresults,R2))
}
