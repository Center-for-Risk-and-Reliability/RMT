# Three Parameter Weibull Probability Plot Parameters
# Developed by Dr. Reuel Smith, 2021-2022

probplotparam.wbl3P <- function(xi,R) {
  # Upper and lower bounds of the Percent Failure axis in percent
  fcB <- log(log(1/(1-c(.001,0.999))))
  #fcB <- log(log(1/(1-seq(from = 0.001, to = 0.999, by = 0.001))))

  # Gamma estimates
  if(length(xi) >= 2){
    gamest<-seq(from = 0.001, to = 0.999, by = 0.01)
    R2up<-0
    for(i in 1:(length(gamest)-1)){
      # Set up the x and y data for fitting
      if(R[length(R)]==0){
        yfit <- log(log(1/R[1:length(R)-1]))
        xfit <- log(xi[1:length(R)-1] - gamest[i])
      } else {
        yfit <- log(log(1/R))
        xfit <- log(xi - gamest[i])
      }
      # Calculate least-squares estimate for Weibull parameters
      wblparm  <- lm(yfit ~ poly(xfit, 1, raw=TRUE))
      R2 <- summary(wblparm)$r.squared
      if(R2>=R2up){
        R2up <- R2
        wblparmup <- wblparm
        gamparam <- gamest[i]
      }
    }
    beta <- summary(wblparmup)$coefficients[2,1]
    alpha <- exp(-summary(wblparmup)$coefficients[1,1]/beta)
    gammaparam <- gamparam
    R2 <- summary(wblparmup)$r.squared
    # Calculate upper and lower bound of best fit line
    intercept <- summary(wblparmup)$coefficients[1,1]
    ttfc <- c(exp((log(log(1./(1-0.001)))-intercept)/beta)+gammaparam,exp((log(log(1./(1-0.999)))-intercept)/beta)+gammaparam)
    #ttfc <- exp((log(log(1./(1-fcB)))-intercept)/beta)+gammaparam
  }
  if(length(xi) == 1){
    alpha <- xi
    beta <- NA
    gammaparam <- 0
    R2 <- 0
    ttfc <- rep(xi,2)
  }

  wblresults <- matrix(c(alpha,beta,gammaparam), nrow = 1, ncol = 3, byrow = TRUE,dimnames = list(c("3P Wbl Parameters"),c("alpha", "beta", "Gamma")))
  return(list(ttfc,fcB,wblresults,R2))
}
