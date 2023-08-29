# Lognormal Probability Plot Parameters
# Developed by Dr. Reuel Smith, 2021-2022

probplotparam.logn <- function(xi,F) {
  # Upper and lower bounds of the Percent Failure axis in percent
  fcB <- qnorm(c(.001,0.999),mean=0,sd=1)
  # Set up the x and y data for fitting
  if(F[length(F)]==1){
    yfit <- qnorm(F[1:length(F)-1],mean=0,sd=1)
    #yfit <- F*100
    xfit <- log(xi[1:length(F)-1])
  } else {
    yfit <- qnorm(F,mean=0,sd=1)
    #yfit <- F*100
    xfit <- log(xi)
  }

  # Calculate least-squares estimate for lognormal parameters
  if(length(xi) >= 2){
    pb2 <- lm(yfit ~ poly(xfit, 1, raw=TRUE))
    intercept <- summary(pb2)$coefficients[1,1]
    slope <- summary(pb2)$coefficients[2,1]
    sigmat <- 1/slope
    meant <- -intercept*sigmat
    R2 <- summary(pb2)$r.squared
    #meant <- (50-slope)/intercept
    #t84 <- (84-slope)/intercept
    #sigmat <- t84-meant
    lognresults <- matrix(c(meant,sigmat), nrow = 1, ncol = 2, byrow = TRUE,dimnames = list(c("Lognormal Parameters"),c("mean_t", "sigma_t")))
    # Calculate upper and lower bound of best fit line
    xfit2 <- c(meant,meant+sigmat)
    yfit2 <- c(qnorm(0.5,mean=0,sd=1),qnorm(0.84,mean=0,sd=1))
    pb2_B <- lm(yfit2 ~ poly(xfit2, 1, raw=TRUE))
    min_mu_sig <- summary(pb2_B)$coefficients[2,1]
    siginv <- summary(pb2_B)$coefficients[1,1]
    ttfc <- c(exp((qnorm(0.001,mean=0,sd=1)-siginv)/min_mu_sig),exp((qnorm(0.999,mean=0,sd=1)-siginv)/min_mu_sig))
  }
  if(length(xi) == 1){
    sigmat <- 0
    meant <- log(xi)
    R2 <- 0
    lognresults <- matrix(c(meant,sigmat), nrow = 1, ncol = 2, byrow = TRUE,dimnames = list(c("Lognormal Parameters"),c("mean_t", "sigma_t")))
    ttfc <- rep(xi,2)
  }
  return(list(ttfc,fcB,lognresults,R2))
}
